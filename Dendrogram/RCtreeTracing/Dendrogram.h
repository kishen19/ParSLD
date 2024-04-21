#pragma once

#include "gbbs/gbbs.h"
#include "gbbs/julienne.h"

#include "rctree_hashtable.h"
#include "rctree_crosslink.h"
#include "rctree_async.h"
#include "rctree_utils.h"

namespace gbbs {

// GA is the weighted input tree.
template <class IdType, class Graph>
auto DendrogramRCtreeTracing_impl(Graph& GA, bool debug = false) {
  using W = typename Graph::weight_type;

  timer t;
  t.start();
  // Preprocess
  auto n = GA.n;
  auto m = GA.m / 2;

  //auto rctree = build_rctree_ht(GA);
  // auto rctree = build_rctree_crosslink(GA);
  auto [rctree, offsets, neighbors] = build_rctree_async<IdType>(GA);
  t.next("BuildRCTree");

  if (debug) {
    for (uintE i=0; i<n; i++){
      std::cout << "i: " << i << ", " << ", parent: " << rctree[i].parent << ", edge: " << rctree[i].edge_index << std::endl;
    }
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
  t.next("Tracing");

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
  t.next("Sort");

  if (debug) {
    for (size_t i=0; i<m; i++){
        std::cout << parent[i] << " ";
    }
    std::cout << std::endl;
  }
  return parent;
}

// GA is the weighted input tree.
template <class Graph>
auto DendrogramRCtreeTracing(Graph& GA, bool debug = false) {
  if (GA.n >= std::numeric_limits<int32_t>::max()) {
    return DendrogramRCtreeTracing_impl<size_t>(GA, debug);
  } else {
    return DendrogramRCtreeTracing_impl<uint32_t>(GA, debug);
  }
}

}   // namespace gbbs