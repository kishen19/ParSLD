#pragma once

// A cross-linked implementation of RC-tree construction.

#include "concurrent_table.h"
#include "gbbs/gbbs.h"
#include "rctree_utils.h"

namespace gbbs {

template <class Graph>
auto build_rctree_crosslink(Graph& GA) {
  using W = typename Graph::weight_type;

  // Key in the hash table (lex-ordered edge).
  using K = std::pair<uintE, uintE>;
  // Value in the hash table.
  using V = std::pair<uintE, W>;

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
    // ngh's offset in src is ind.
  };
  parallel_for(0, n, [&](uintE u) {
    GA.get_vertex(u).out_neighbors().map_with_index(offset_f);
  });
  parlay::sort_inplace(edges);

  // Sort the triples, then drop the first endpoint (project to pairs).
  // parlay::sort_inplace(triples);
  t.next("Sort edges");

  auto get_key = [](uintE x, uintE y) {
    return std::make_pair(std::min(x, y), std::max(x, y));
  };

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

  // Step 1: Compute RC Tree
  auto rctree = parlay::sequence<RCtree_node<W>>::uninitialized(n);

  // The rc-tree node where this edge was contracted / merged. This
  // is used in the query-path when we traverse up the RC tree for an
  // edge.
  parallel_for(0, n, [&](uintE i) {
    rctree[i].parent = i;
    // to ensure unique buckets, even if the input tree is disconnected
    rctree[i].edge_index = m + i;
    rctree[i].alt = UINT_E_MAX;
    rctree[i].round = UINT_E_MAX;
  });

  parlay::random r(42);
  auto priorities = sequence<uint64_t>::from_function(
     n, [&](uintE i) { return r.ith_rand(i); });
  uintE round = 0;

//  auto priorities = sequence<uint32_t>::from_function(
//     n, [&](uintE i) { return parlay::hash32_3(i); });
//  uintE rem = m, round = 0;

  t.next("Preprocess (initialize RCTree Nodes and Priorities");

  // Implements the rake operation.
  // If deg 1 cluster merges into a deg>1 cluster:
  //    simply merge and unite the clusters
  //

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

  auto rake_f = [&](uintE src) -> bool {
    // auto entries = neighbors[src].entries(UINT_E_MAX - 1);

    // We know that src has degree 1. We just need to get its neighbor
    // now by scanning neighbors[offsets[src]].

    auto [dst, index_in_dst, edge_index, wgh] = get_neighbor(src);

    // Need to symmetry break two degree-1 nodes here.
    if (deg[dst] > 1 || (deg[dst] == 1 && src < dst)) {

      deg[src]--;
      rctree[src].parent = dst;
      rctree[src].round = round;
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

  // // Implements the compress operation.
  //    Merge into the cluster with min-weight edge
  auto compress = [&](const uintE& src) {
    auto our_neighbors = get_both_neighbors(src);
    if (our_neighbors.size() != 2) {
      std::cerr << "Bad compress" << std::endl;
      assert(false);
      exit(-1);
    }
    auto [dst1, index_in_dst1, edge_index1, wgh1] = our_neighbors[0];
    auto [dst2, index_in_dst2, edge_index2, wgh2] = our_neighbors[1];

    // Check whether we want to compress this node.
    if (priorities[src] > priorities[dst1] &&
        priorities[src] > priorities[dst2]) {
      deg[src] -= 2;
      rctree[src].round = round;
    }
  };

  // We separate compress and finish_compress here to avoid a
  // race-condition (since src and dst1 could both have compress
  // called on them, e.g., in a path).
  auto finish_compress = [&](uintE src) {
    if (rctree[src].round != UINT_E_MAX) {

      auto our_neighbors = get_both_neighbors(src);
      if (our_neighbors.size() != 2) {
        std::cerr << "Bad compress" << std::endl;
        assert(false);
        exit(-1);
      }
      auto [dst1, index_in_dst1, edge_index1, wgh1] = our_neighbors[0];
      auto [dst2, index_in_dst2, edge_index2, wgh2] = our_neighbors[1];

      // Any order here is fine, since we post-process and check both.
      rctree[src].parent = dst1;
      rctree[src].alt = dst2;

      // Check which edge is smaller; break ties using the indices.
      if (std::make_pair(wgh1, edge_index1) <
          std::make_pair(wgh2, edge_index2)) {
        // The edge to 1 is smaller.
        // TODO(): combine the cases into one by std::swap?
        rctree[src].edge_index = edge_index1;
        rctree[src].wgh = wgh1;

        // Now, update neighbors for dst1 and dst2 to cross-link with
        // each other.
        neighbors[offsets[dst1] + index_in_dst1] = {dst2, index_in_dst2,
                                                    edge_index2, wgh2};
        neighbors[offsets[dst2] + index_in_dst2] = {dst1, index_in_dst1,
                                                    edge_index2, wgh2};
      } else {
        rctree[src].edge_index = edge_index2;
        rctree[src].wgh = wgh2;

        // Now, update neighbors for dst1 and dst2 to cross-link with
        // each other.
        neighbors[offsets[dst1] + index_in_dst1] = {dst2, index_in_dst2,
                                                    edge_index1, wgh1};
        neighbors[offsets[dst2] + index_in_dst2] = {dst1, index_in_dst1,
                                                    edge_index1, wgh1};
      }
    }
  };

  auto delayed_n = parlay::delayed_seq<uintE>(n, [&] (uintE i) { return i; });
  auto degree_one =
     parlay::filter(delayed_n, [&](uintE i) { return deg[i] == 1; });
  auto degree_two =
     parlay::filter(delayed_n, [&](uintE i) { return deg[i] == 2; });
  t.next("Before While loop (filters)");
  while (degree_one.size() > 0) {
    timer rt;
    rt.start();
    // Note that we need both the degree 1 and degree 2 buckets.
    // Only two buckets we care about: {1, 2, everything_else}
    // Dense iterations do help us here...

    // (1) filtering out the initial frontier (degree 1 vertices)
    // Peel:
    // (a) emit the id of the neighbor / peeled (degree1) vertex
    // (b) histogram the ids (do this using semisort / sort)
    // (c) update the degrees, and filter out those that become degree 1 / 2

    std::cout << "Degree_one.size = " << degree_one.size() << std::endl;
    std::cout << "Degree_two.size = " << degree_two.size() << std::endl;

    // Compress
    if (degree_two.size() > 0) {
      // Adds new edges
      parallel_for(0, degree_two.size(),
                  [&](uintE i) { compress(degree_two[i]); });
      // sync, then finish_compress.
      parallel_for(0, degree_two.size(),
                  [&](uintE i) { finish_compress(degree_two[i]); });
    }
    round++;
    rt.next("Compress time");

    // Rake
    parallel_for(0, degree_one.size(), [&](uintE i) {
      auto outp = rake_f(degree_one[i]);
      degree_one[i] = outp ? rctree[degree_one[i]].parent : n;
    });
    rt.next("Rake time");

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
       [&](uintE i) { return (indices[i] < n) ? degree_one[indices[i]] : n; });
    auto new_degree_one = parlay::filter(
       unique_keys, [&](uintE id) { return (id < n) ? deg[id] == 1 : false; });
    auto new_degree_two = parlay::filter(unique_keys, [&](uintE id) {
      // Check if this node was symmetry broken to live.
      return (id < n) ? deg[id] == 2 : false;
    });
    degree_one = std::move(new_degree_one);

    auto rem_degree_two =
       parlay::filter(degree_two, [&](uintE id) { return deg[id] == 2; });
    auto combined_degree_two = parlay::append(rem_degree_two, new_degree_two);
    degree_two = std::move(combined_degree_two);
    round++;
    rt.next("Prepare next round time");
    // std::cout << rem << std::endl;
  }
  std::cout << "# Rounds = " << round << std::endl;
  t.next("RC Tree Time");

  return rctree;
}

}   // namespace gbbs
