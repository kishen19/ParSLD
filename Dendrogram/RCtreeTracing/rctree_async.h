#pragma once

// A cross-linked implementation of RC-tree construction.

#include "concurrent_table.h"
#include "gbbs/gbbs.h"
#include "rctree_utils.h"

namespace gbbs {

template <class Graph>
auto build_rctree_async(Graph& GA) {
  using W = typename Graph::weight_type;

  timer t;
  t.start();
  auto n = GA.n;
  auto m = GA.m / 2;

  // Stores the degrees as nodes are peeled.
  auto deg = sequence<uintE>::from_function(
     n, [&](uintE u) { return GA.get_vertex(u).out_degree(); });
  auto edges = sequence<std::tuple<uintE, uintE, uintE>>::uninitialized(2 * m);

  auto offsets = parlay::scan(deg).first;
  auto offset_f = [&](const uintE& src, const uintE& dst, const W& wgh,
                      const uintE& ind) {
    edges[offsets[src] + ind] = {std::min(src, dst), std::max(src, dst), offsets[src] + ind};
    std::cout << "Edge: " << src << " " << dst << " " << wgh << std::endl;
    // ngh's offset in src is ind.
  };
  parallel_for(0, n, [&](uintE u) {
    GA.get_vertex(u).out_neighbors().map_with_index(offset_f);
  }, 100000000);
  parlay::sort_inplace(edges);
  t.next("Sort edges");

  using edge_info = std::tuple<uintE, uintE, uintE, W>;
  // Neighbors of vertex i at offsets[i] -- offsets[i+1]:
  // We store:
  // - neighbor_id
  // - index_in_neighbor
  // - edge_index of this edge (the original identity of this edge)
  // - weight of the edge (redundant, but saves a cache miss)
  auto neighbors = parlay::sequence<edge_info>::uninitialized(2*m);

  parlay::parallel_for(0, m, [&](uintE i) {
       auto [u, v, ind1] = edges[2*i];
       auto [_u, _v, ind2] = edges[2*i+1];

       auto [_, wgh] = GA.get_vertex(u).out_neighbors().get_ith_neighbor(ind1-offsets[u]);
       // smaller of the two inds corresponds to the smaller vertex
       // id.
       // G.get_vertex(u).get_ith_neighbor(ind1)
       neighbors[ind1] = {v, ind2-offsets[v], i, wgh};
       neighbors[ind2] = {u, ind1-offsets[u], i, wgh};
  });
  t.next("Initialize Neighbors 4-tuples");

  edges.clear();

  // Step 1: Compute RC Tree
  auto rctree = parlay::sequence<RCtree_node<W>>::uninitialized(n);
  // The rc-tree node where this edge was contracted / merged. This
  // is used in the query-path when we traverse up the RC tree for an
  // edge.
  parallel_for(0, n, [&](uintE i) {
    rctree[i].parent = i;
    // to ensure unique buckets, even if the input tree is disconnected
    rctree[i].edge_index = m + i;
  });

  parlay::random r(42);
  auto priorities = sequence<uint64_t>::from_function(
     n, [&](uintE i) { return r.ith_rand(i); });

  t.next("Preprocess (initialize RCTree Nodes and Priorities");

  auto get_neighbor = [&](uintE src) -> edge_info {
    uintE offset = offsets[src];
    uintE deg = ((src == n - 1) ? (2 * m) : offsets[src + 1]) - offset;

    for (uintE i = 0; i < deg; ++i) {
      if (std::get<0>(neighbors[offset + i]) != UINT_E_MAX) {
        return neighbors[offset + i];
      }
    }
    std::cerr << "Calling get single neighbor on a bad source: " << src
              << std::endl;
    assert(false);
    exit(-1);
  };

  auto get_both_neighbors = [&](uintE src) -> std::vector<edge_info> {
    uintE offset = offsets[src];
    uintE deg = ((src == n - 1) ? (2 * m) : offsets[src + 1]) - offset;
    std::vector<edge_info> ret;

    for (uintE i = 0; i < deg; ++i) {
      if (std::get<0>(neighbors[offset + i]) != UINT_E_MAX) {
        ret.push_back(neighbors[offset + i]);
      }
    }
    assert(ret.size() == 2);
    return ret;
  };

  // Implements the rake operation.
  auto rake_f = [&](uintE src) -> bool {
    // We know that src has degree 1. We just need to get its neighbor
    // now by scanning neighbors[offsets[src]].
    auto [dst, index_in_dst, edge_index, wgh] = get_neighbor(src);
    assert(deg[src] == 1);

    // Need to symmetry break two degree-1 nodes here.
    if (deg[dst] > 1 || (deg[dst] == 1 && src < dst)) {
      deg[src]--;
      rctree[src].parent = dst;
      rctree[src].edge_index = edge_index;
      rctree[src].wgh = wgh;
      if (deg[dst] != 1) {
        // Mark this edge as deleted.
        std::get<0>(neighbors[offsets[dst] + index_in_dst]) = UINT_E_MAX;
      }
      return true;
    } else {
      return false;
    }
  };

  auto pri_greater = [&] (uintE src, uintE dst) {
    uintE deg_dst = deg[dst];
    return deg_dst != 2 || (priorities[src] > priorities[dst]);
  };

  // // Implements the compress operation.
  //    Merge into the cluster with min-weight edge
  auto compress = [&](const uintE& src) {
    // Invariant: src has degree 2
    uintE deg_src = deg[src];
    auto our_neighbors = get_both_neighbors(src);
    if (our_neighbors.size() != 2 || deg_src != 2) {
      std::cerr << "Bad compress: num neighbors = " << our_neighbors.size() << std::endl;
      assert(false);
      exit(-1);
    }
    auto [dst1, index_in_dst1, edge_index1, wgh1] = our_neighbors[0];
    auto [dst2, index_in_dst2, edge_index2, wgh2] = our_neighbors[1];

    // Check whether we want to compress this node (check neighbor
    // priorities).
    rctree[src].parent = UINT_E_MAX - 3 + static_cast<uint16_t>(pri_greater(src, dst1)) + static_cast<uint16_t>(pri_greater(src, dst2));
  };

  // We separate compress and finish_compress here to avoid a
  // race-condition (since src and dst1 could both have compress
  // called on them, e.g., in a path).
  auto finish_compress = [&](uintE src) {
    auto cur = src;
    while (gbbs::atomic_compare_and_swap(&rctree[cur].parent, UINT_E_MAX - 1, UINT_E_MAX)) {
      auto our_neighbors = get_both_neighbors(cur);
      if (our_neighbors.size() != 2) {
        assert(false);
        exit(-1);
      }
      deg[cur] -= 2;
      auto [dst1, index_in_dst1, edge_index1, wgh1] = our_neighbors[0];
      auto [dst2, index_in_dst2, edge_index2, wgh2] = our_neighbors[1];

      neighbors[offsets[cur]] = our_neighbors[0];
      neighbors[offsets[cur] + 1] = our_neighbors[1];

      // Check which edge is smaller; break ties using the indices.
      // Set the parent to be the neighbor along the smaller-weight
      // edge.
      if (std::make_pair(wgh1, edge_index1) <
          std::make_pair(wgh2, edge_index2)) {
        // The edge to 1 is smaller.
        // TODO(): combine the cases into one by std::swap?
        rctree[cur].edge_index = edge_index1;
        rctree[cur].wgh = wgh1;
        // Now, update neighbors for dst1 and dst2 to cross-link with
        // each other.
        neighbors[offsets[dst1] + index_in_dst1] = {dst2, index_in_dst2,
                                                    edge_index2, wgh2};
        neighbors[offsets[dst2] + index_in_dst2] = {dst1, index_in_dst1,
                                                    edge_index2, wgh2};
        // Set the parent.
        rctree[cur].parent = dst1;
      } else {
        rctree[cur].edge_index = edge_index2;
        rctree[cur].wgh = wgh2;
        // Now, update neighbors for dst1 and dst2 to cross-link with
        // each other.
        neighbors[offsets[dst1] + index_in_dst1] = {dst2, index_in_dst2,
                                                    edge_index1, wgh1};
        neighbors[offsets[dst2] + index_in_dst2] = {dst1, index_in_dst1,
                                                    edge_index1, wgh1};
        // Set the parent.
        rctree[cur].parent = dst2;
      }

      if (pri_greater(dst1, dst2) && deg[dst1]==2) {
        gbbs::write_add(&rctree[dst1].parent, 1);
        // We know that dst1 will be compressed in this super-round,
        // and it has higher priority than dst2, so dst1 will be our
        // parent.
        cur = dst1;
      } else if (pri_greater(dst2, dst1) && deg[dst2]==2) {
        gbbs::write_add(&rctree[dst2].parent, 1);
        // Ditto, in the other direction.
        cur = dst2;
      }
    }
  };

  auto delayed_n = parlay::delayed_seq<uintE>(n, [&] (uintE i) { return i; });
  auto degree_one =
     parlay::filter(delayed_n, [&](uintE i) { return deg[i] == 1; });
  auto degree_two =
     parlay::filter(delayed_n, [&](uintE i) { return deg[i] == 2; });
  t.next("Before While loop (filters)");
  size_t round = 0;
  while (degree_one.size() > 0) {
    // timer rt;
    // rt.start();
    // (1) filtering out the initial frontier (degree 1 vertices)
    // Peel:
    // (a) emit the id of the neighbor / peeled (degree1) vertex
    // (b) histogram the ids (do this using semisort / sort)
    // (c) update the degrees, and filter out those that become degree 1 / 2
    std::cout << "Degree_one.size = " << degree_one.size() << std::endl;
    std::cout << "Degree_two.size = " << degree_two.size() << std::endl;

   // Compress
   if (degree_two.size() > 0) {
     // Identify an independent set of degree 2 nodes to compress
     parallel_for(0, degree_two.size(),
                 [&](uintE i) { compress(degree_two[i]); });

     // Async'ly compress the rest of the degree 2 nodes.
     parallel_for(0, degree_two.size(),
                 [&](uintE i) { finish_compress(degree_two[i]); });
   }
   // rt.next("Compress time");

    // Rake
    parallel_for(0, degree_one.size(), [&](uintE i) {
      auto outp = rake_f(degree_one[i]);
      degree_one[i] = outp ? rctree[degree_one[i]].parent : n;
    });
    // rt.next("Rake time");

    parlay::sort_inplace(parlay::make_slice(degree_one));
    auto flags = parlay::delayed_seq<bool>(degree_one.size(), [&](uintE i) {
      return (i == 0) || (degree_one[i] != degree_one[i - 1]);
    });
    auto indices = parlay::pack_index(flags);

    parlay::parallel_for(0, indices.size(), [&](uintE i) {
      auto start = indices[i];
      auto end = (i == indices.size() - 1) ? degree_one.size() : indices[i + 1];
      uintE dst = degree_one[start];
      uintE degree_lost = end - start;
      if (dst < n) {
        deg[dst] -= degree_lost;
      }
    });
    auto unique_keys = parlay::delayed_seq<uintE>(
       indices.size(),
       [&](uintE i) { return degree_one[indices[i]]; });

    // Update degree one and two.
    auto degree_one_new = parlay::filter(
       unique_keys, [&](uintE id) { return (id < n) ? deg[id] == 1 : false; });
    auto degree_two_new = parlay::filter(unique_keys, [&](uintE id) {
      // Check if this node was symmetry broken to live.
      return (id < n) ? (deg[id] == 2) : false;
    });
    degree_one = std::move(degree_one_new);
    degree_two = std::move(degree_two_new);

    round++;
    // rt.next("Prepare next round time");
    // std::cout << rem << std::endl;
  }
  std::cout << "# Rounds = " << round << std::endl;
  t.next("RC Tree Time");

  return std::make_tuple(std::move(rctree), std::move(offsets), std::move(neighbors));
}

}   // namespace gbbs
