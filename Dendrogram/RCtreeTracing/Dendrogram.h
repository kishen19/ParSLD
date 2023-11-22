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
  using W = typename Graph::weight_type;

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

  auto bkt_ids = sequence<std::pair<uintE, std::pair<W, uintE>>>::uninitialized(m);
  parallel_for(0, n, [&](size_t i) {
    // Start at every RCtree node, and the edge associated with it.
    auto edge_index = rctree[i].edge_index;
    if (edge_index < m) { // Check for a non-root.
      // This is a valid non-root node.
      auto wgh = std::make_pair(rctree[i].wgh, edge_index);
      auto par = rctree[i].parent;
      auto par_wgh = std::make_pair(rctree[par].wgh, rctree[par].edge_index);
      // While the parent is not the root, and the parent weight is
      // smaller than our weight, go to the parent.
      while (rctree[par].parent != par && par_wgh < wgh) {
        par = rctree[par].parent;
        par_wgh = std::make_pair(rctree[par].wgh, rctree[par].edge_index);
      }
      // Sort by initial weight, breaking ties by edge-index.
      bkt_ids[edge_index] = {rctree[par].edge_index, wgh};
    }
  });
  t.next("Bucket Computation Time");

  // Step 3: Sort by bkt_ids
  auto parent = sequence<uintE>::from_function(m, [&](size_t i) { return i; });
  sort_inplace(bkt_ids);
  parallel_for(0, n - 1, [&](size_t i) {
    // There is a subsequent edge in this bucket---this is our parent.
    if (bkt_ids[i].first == bkt_ids[i + 1].first) {
      parent[bkt_ids[i].second.second] = bkt_ids[i + 1].second.second;
    } else {
      // If the bucket is a non-root node, assign it to the edge
      // associated with the bucket (this is the last edge in this
      // bucket).
      if (bkt_ids[i].first < m) {
        parent[bkt_ids[i].second.second] = bkt_ids[i].first;
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