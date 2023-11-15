#pragma once

#include "benchmarks/Connectivity/SimpleUnionAsync/Connectivity.h"
#include "gbbs/gbbs.h"
#include "gbbs/julienne.h"

#include "rctree_hashtable.h"
#include "rctree_crosslink.h"
#include "rctree_utils.h"

namespace gbbs {

// GA is the weighted input tree.
template <class Graph>
double DendrogramRCtreeTracing(Graph& GA) {
  using W = typename Graph::weight_type;
  using K = std::pair<uintE, uintE>;
  using V = size_t;
  using KA = uintE;
  using VA = std::pair<W, size_t>;

  timer t;
  t.start();
  // Preprocess
  auto n = GA.n;
  auto m = GA.m / 2;

  //auto rctree = build_rctree_ht(GA);
  auto rctree = build_rctree_crosslink(GA);
  t.next("Build RCTree Time");

      // Step 2: Compute bucket_id for each edge
  parallel_for(0, n, [&](size_t i) {
    if (rctree[i].alt != UINT_E_MAX) {
      auto alt = rctree[i].alt;
      if (rctree[alt].round < rctree[i].round) {
        rctree[i].parent = alt;
      }
    }
  });
  auto bkt_ids = sequence<std::pair<size_t, size_t>>::uninitialized(m);
  parallel_for(0, n, [&](size_t i) {
    auto edge_index = rctree[i].edge_index;
    if (edge_index < m) {
      auto wgh = std::make_pair(rctree[i].wgh, edge_index);
      auto cur = rctree[i].parent;
      auto new_wgh = std::make_pair(rctree[cur].wgh, rctree[cur].edge_index);
      while (rctree[cur].parent != cur && new_wgh < wgh) {
        cur = rctree[cur].parent;
        new_wgh = std::make_pair(rctree[cur].wgh, rctree[cur].edge_index);
      }
      bkt_ids[edge_index] = {rctree[cur].edge_index, edge_index};
    }
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
//  for (size_t i=0; i<m; i++){
//      std::cout << parent[i] << " ";
//  }
//  std::cout << std::endl;
  return tt;
}

}   // namespace gbbs