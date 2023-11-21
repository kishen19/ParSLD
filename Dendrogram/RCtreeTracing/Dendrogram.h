#pragma once

#include "gbbs/gbbs.h"
#include "gbbs/julienne.h"

#include "rctree_hashtable.h"
#include "rctree_crosslink.h"
#include "rctree_async.h"
#include "rctree_utils.h"

namespace gbbs {

// GA is the weighted input tree.
template <class Graph>
auto DendrogramRCtreeTracing(Graph& GA, bool debug = false) {

  timer t;
  t.start();
  // Preprocess
  auto n = GA.n;
  auto m = GA.m / 2;

  //auto rctree = build_rctree_ht(GA);
  // auto rctree = build_rctree_crosslink(GA);
  auto [rctree, offsets, neighbors] = build_rctree_async(GA);
  t.next("Build RCTree Time");

  // Only need to scan two neighbors since we prepare this
  // neighbor-list in finish_compress.
  auto get_both_neighbors = [&](uintE src) -> std::vector<uintE> {
    uintE offset = offsets[src];
    std::vector<uintE> ret;

    for (uintE i = 0; i < 2; ++i) {
      ret.push_back(std::get<0>(neighbors[offset + i]));
    }
    return ret;
  };

  for (uintE i=0; i<n; i++){
    std::cout << "i: " << i << ", round: " << rctree[i].round << ", parent: " << rctree[i].parent << ", edge: " << rctree[i].edge_index << std::endl;
  }

      // Step 2: Compute bucket_id for each edge
  parallel_for(0, n, [&](size_t i) {
    if (rctree[i].parent == UINT_E_MAX) {
      auto neighbors = get_both_neighbors(i);
      if (neighbors.size() != 2) {
        std::cerr << "wrong!" << std::endl;
        exit(-1);
      }

      auto p1 = neighbors[0];
      auto p2 = neighbors[1];
      std::cout << "edge: " << rctree[i].edge_index << ", i: " << i << ", p1: " << p1 << ", p2: " << p2 << std::endl;
      if (rctree[p2].round < rctree[p1].round) {
        std::swap(p1, p2);
      }
      // Parent cannot have round 0 (this indicates the root).
      if (rctree[p1].round == 0) {
        std::swap(p1, p2);
      }
      rctree[i].parent = p1;
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
        if (edge_index == 6) {
          std::cout << "cur = " << cur << " new_wgh = " << new_wgh.first << " " << new_wgh.second << std::endl;
        }
        cur = rctree[cur].parent;
        new_wgh = std::make_pair(rctree[cur].wgh, rctree[cur].edge_index);
      }
      if (edge_index == 6) {
        std::cout << "bucket = " << rctree[cur].edge_index << " " << edge_index << std::endl;
      }
      bkt_ids[edge_index] = {rctree[cur].edge_index, edge_index};
    }
  });
  t.next("Bucket Computation Time");

  // Step 3: Sort by bkt_ids
  auto parent = sequence<uintE>::from_function(m, [&](size_t i) { return i; });
  sort_inplace(bkt_ids);
  parallel_for(0, n - 1, [&](size_t i) {
    if (bkt_ids[i].first == bkt_ids[i + 1].first) {
      parent[bkt_ids[i].second] = bkt_ids[i + 1].second;
    } else {
      if (bkt_ids[i].first < m) {
        parent[bkt_ids[i].second] = bkt_ids[i].first;
      }
    }
  });
  t.next("Bucket Sorting and Finish Time");

  if (debug) {
    for (size_t i=0; i<m; i++){
        std::cout << parent[i] << " ";
    }
    std::cout << std::endl;
  }
  return parent;
}

}   // namespace gbbs