#pragma once

#include "gbbs/gbbs.h"
#include "gbbs/julienne.h"
#include "concurrent_table.h"
#include "benchmarks/Connectivity/SimpleUnionAsync/Connectivity.h"

namespace gbbs {

// ~40 bytes / struct
template <typename weight_type>
struct RCtree_node {
  size_t parent;
  size_t round;
  size_t edge_index;
  weight_type wgh;  // TODO: get rid of edge_index in the weight and use the edge_index above
  // could get rid of this by checking if alt != UINT_E_MAX
  bool check = false;
  uintE alt;   // used in the compress case. (default = UINT_E_MAX)
};

// GA is the weighted input tree.
template <class Graph>
double DendrogramRCtreeTracing(Graph& GA) {
  using W = typename Graph::weight_type;
  using K = std::pair<uintE, uintE>;
  using V = size_t;

  timer t;
  t.start();
  // Preprocess
  auto n = GA.n;
  auto m = GA.m / 2;
  auto empty = std::make_tuple(std::make_pair(UINT_E_MAX, UINT_E_MAX), m);

  // Stores the degrees as nodes are peeled.
  auto deg = sequence<size_t>::from_function(
     n, [&](size_t u) { return GA.get_vertex(u).out_degree(); });

  // hash table needs to store ~n entries.
  auto S = make_concurrent_table<K, V>(m, empty, 1.05);

  auto edges = sequence<std::pair<uintE, uintE>>::uninitialized(2 * m);
  auto offsets = parlay::scan(deg).first;
  auto offset_f = [&](const uintE& src, const uintE& dst, const W& wgh,
                      const size_t& ind) {
    edges[offsets[src] + ind] = {std::min(src, dst), std::max(src, dst)};
  };
  parallel_for(0, n, [&](uintE u) {
    GA.get_vertex(u).out_neighbors().map_with_index(offset_f);
  });
  parlay::sort_inplace(edges);
  parallel_for(0, m, [&](size_t i) {
    auto key_value = std::make_tuple(edges[2 * i], i);
    S.insert(key_value);
  });
  t.next("Preprocess Time");

  // Step 1: Compute RC Tree
  auto rctree = sequence<RCtree_node<std::pair<W, size_t>>>(n);
  // The rc-tree node where this edge was contracted / merged. This
  // is used in the query-path when we traverse up the RC tree for an
  // edge. We may be able to get rid of this (e.g., by doing a
  // parallel map over all RC tree nodes, and using the edgeindex
  // stored at the node), but it's also fairly cheap and may not hurt
  // to keep around.
  auto edge2rcnode = sequence<size_t>::uninitialized(m);
  parallel_for(0, n, [&](size_t i) {
    rctree[i].parent = i;
    rctree[i].edge_index =
       m +
       i;   // to ensure unique buckets, even if the input tree is disconnected
  });
  size_t rem = m, round = 0;

  // Implements the rake operation.
  auto map_f = [&](const uintE& src, const uintE& dst, const W& wgh) {
    if (deg[dst] > 1) {
      deg[src]--;
      rctree[src].parent = dst;
      rctree[src].round = round;
      auto key = std::make_pair(std::min(src, dst), std::max(src, dst));
      auto edge_index = S.find(key);
      if (edge_index == m) {
        std::cout << src << " " << dst << " " << edge_index << std::endl;
      }
      rctree[src].edge_index = edge_index;
      edge2rcnode[edge_index] = src;
      rctree[src].wgh = {wgh, edge_index};
    } else if (deg[dst] == 1) {
      // do something; tiebreak based on id?
    }
  };

  auto degree_one =
       parlay::filter(parlay::iota(n), [&](size_t i) { return deg[i] == 1; });
  while (rem) {
    // Note that we need both the degree 1 and degree 2 buckets.
    // Only two buckets we care about: {1, 2, everything_else}
    // Dense iterations do help us here...

    // (1) filtering out the initial frontier (degree 1 vertices)
    // Peel:
    // (a) emit the id of the neighbor / peeled (degree1) vertex
    // (b) histogram the ids (do this using semisort / sort)
    // (c) update the degrees, and filter out those that become degree 1 / 2

    std::cout << "Degree_one.size = " << degree_one.size() << std::endl;
    if (degree_one.size() == 0) {
      std::cout << "Remaining: " << rem << " rounds = " << round << std::endl;
      exit(-1);
    }
    parallel_for(0, degree_one.size(), [&](size_t i) {
      GA.get_vertex(degree_one[i]).out_neighbors().map(map_f);
      degree_one[i] = rctree[degree_one[i]].parent;
    });

    parlay::sort_inplace(parlay::make_slice(degree_one));

    auto flags = parlay::delayed_seq<bool>(degree_one.size(), [&] (size_t i) {
      return (i == 0) || (degree_one[i] != degree_one[i-1]);
    });
    auto indices = parlay::pack_index(flags);

    parlay::parallel_for(0, indices.size(), [&] (size_t i) {
      auto start = indices[i];
      auto end = (i == indices.size() - 1) ? degree_one.size() : indices[i+1];
      size_t dst = degree_one[start];
      size_t degree_lost = end - start;
      deg[dst] -= degree_lost;
    });

    auto unique_keys = parlay::delayed_seq<size_t>(indices.size(), [&] (size_t i) {
      return degree_one[indices[i]];
    });
    auto new_degree_one = parlay::filter(unique_keys, [&] (size_t id) {
      return deg[id] == 1;
    });

    rem -= degree_one.size();
    round++;
    degree_one = std::move(new_degree_one);
  }
  t.next("RC Tree Time");

//  while (rem) {
//    auto cur =
//       parlay::filter(parlay::iota(n), [&](size_t i) { return deg[i] == 1; });
//    parallel_for(0, cur.size(), [&](size_t i) {
//      GA.get_vertex(cur[i]).out_neighbors().map(map_f);
//      cur[i] = rctree[cur[i]].parent;
//    });
//    // std::cout << "check22" << std::endl;
//    auto hist = parlay::histogram_by_key(cur);
//    // std::cout << "check221 " << hist.size() <<  std::endl;
//    parallel_for(0, hist.size(),
//                 [&](size_t i) { deg[hist[i].first] -= hist[i].second; });
//    // std::cout << "check23" << std::endl;
//    round++;
//    rem -= cur.size();
//  }
//  t.next("RC Tree Time");

  // Step 2: Compute bucket_id for each edge
  auto bkt_ids = sequence<std::pair<size_t, size_t>>::uninitialized(m);
  parallel_for(0, m, [&](size_t i) {
    auto cur = edge2rcnode[i];
    auto wgh = rctree[cur].wgh;
    cur = rctree[cur].parent;
    while (rctree[cur].parent != cur && rctree[cur].wgh < wgh) {
      cur = rctree[cur].parent;
    }
    bkt_ids[i] = {rctree[cur].edge_index, i};
  });
  t.next("Bucket Computation Time");

  // Step 3: Sort by bkt_ids
  auto parent = sequence<size_t>::from_function(m, [&](size_t i) { return i; });
  sort_inplace(bkt_ids);
  parallel_for(0, n - 1, [&](size_t i) {
    if (bkt_ids[i].first == bkt_ids[i + 1].first) {
      parent[bkt_ids[i].second] = parent[bkt_ids[i + 1].second];
    } else {
      if (bkt_ids[i].first < m) {
        parent[bkt_ids[i].second] = bkt_ids[i].first;
      }
    }
  });
  t.next("Bucket Sorting and Finish Time");
  double tt = t.total_time();
  // for (size_t i=0; i<m; i++){
  //     std::cout << parent[i] << " ";
  // }
  // std::cout << std::endl;
  return tt;
}

}   // namespace gbbs