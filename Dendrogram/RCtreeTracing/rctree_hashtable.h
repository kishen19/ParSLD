#pragma once

// A hash-table implementation of RC-tree construction

#include "concurrent_table.h"
#include "gbbs/gbbs.h"
#include "rctree_utils.h"

namespace gbbs {

template <class Graph>
auto build_rctree_ht(Graph& GA) {
  using W = typename Graph::weight_type;
  using K = std::pair<uintE, uintE>;
  using V = size_t;
  using KA = uintE;
  using VA = std::pair<W, size_t>;

  timer t;
  t.start();
  auto n = GA.n;
  auto m = GA.m / 2;

  // Stores the degrees as nodes are peeled.
  auto deg = sequence<size_t>::from_function(
     n, [&](size_t u) { return GA.get_vertex(u).out_degree(); });
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
  // hash table for edge -> indices
  //    needs to store ~n entries.
  auto emptyS = std::make_tuple(std::make_pair(UINT_E_MAX, UINT_E_MAX), m);
  auto S = make_concurrent_table<K, V>(m, emptyS, 1.05);
  parallel_for(0, m, [&](size_t i) {
    auto key_value = std::make_tuple(edges[2 * i], i);
    S.insert(key_value);
  });
  // Create Adjacency List: Dynamic, maintained using list of hash tables
  auto W_MAX = std::numeric_limits<W>::max();
  auto emptyAdj = std::make_tuple(UINT_E_MAX, std::make_pair(W_MAX, m));
  auto adj = sequence<concurrent_table<KA, VA>>::uninitialized(n);
  auto adj_f = [&](const uintE& src, const uintE& dst, const W& wgh) {
    auto index = S.find(std::make_pair(std::min(src, dst), std::max(src, dst)));
    auto key_value = std::make_tuple(dst, std::make_pair(wgh, index));
    adj[src].insert(key_value);
  };
  parallel_for(0, n, [&](uintE i) {
    adj[i] = make_concurrent_table<KA, VA>(deg[i], emptyAdj, 1.05);
    GA.get_vertex(i).out_neighbors().map(adj_f);
  });
  t.next("Preprocess Time");

  // Step 1: Compute RC Tree
  auto rctree = sequence<RCtree_node<W>>(n);
  // The rc-tree node where this edge was contracted / merged. This
  // is used in the query-path when we traverse up the RC tree for an
  // edge.
  parallel_for(0, n, [&](size_t i) {
    rctree[i].parent = i;
    // to ensure unique buckets, even if the input tree is disconnected
    rctree[i].edge_index = m + i;
  });
  parlay::random r(42);
  auto priorities = sequence<uint64_t>::from_function(n, [&](size_t i){
    return r.ith_rand(i);
  });
  // r = r.next();
  size_t rem = m, round = 0;

  // Implements the rake operation.
  // If deg 1 cluster merges into a deg>1 cluster:
  //    simply merge and unite the clusters
  //
  auto rake_f = [&](uintE src) -> bool {
    auto entries = adj[src].entries(UINT_E_MAX - 1);
    auto dst = std::get<0>(entries[0]);
    if (deg[dst] > 1 || (deg[dst] == 1 && src < dst)) {
      auto wgh = std::get<1>(entries[0]).first;
      auto edge_index = std::get<1>(entries[0]).second;
      deg[src]--;
      rctree[src].parent = dst;
      rctree[src].round = round;
      rctree[src].edge_index = edge_index;
      rctree[src].wgh = wgh;
      if (deg[dst] != 1) {
        adj[dst].remove(src);
      }
      return true;
    } else {
      return false;
    }
  };

  // // Implements the compress operation.
  //    Merge into the cluster with min-weight edge
  auto compress_f1 = [&](const uintE& src) {
    auto entries = adj[src].entries(UINT_E_MAX - 1);
    auto dst1 = std::get<0>(entries[0]);
    auto dst2 = std::get<0>(entries[1]);
    if (priorities[src] > priorities[dst1] && priorities[src] > priorities[dst2]) {
      deg[src] -= 2;
      rctree[src].round = round;
    }
  };

  // // Implements the compress operation.
  //    Merge into the cluster with min-weight edge
  auto compress_f2 = [&](const uintE& src) {
    // Restore Coin Toss values
    // tosses[src] = 0;
    if (rctree[src].round != UINT_E_MAX) {
      auto entries = adj[src].entries(UINT_E_MAX - 1);
      auto dst1 = std::get<0>(entries[0]);
      auto dst2 = std::get<0>(entries[1]);
      adj[dst1].remove(src);
      adj[dst2].remove(src);
    }
  };
  auto compress_f3 = [&](const uintE& src) {
    if (rctree[src].round != UINT_E_MAX) {
      auto entries = adj[src].entries(UINT_E_MAX - 1);
      auto dst1 = std::get<0>(entries[0]);
      auto dst2 = std::get<0>(entries[1]);

      auto edge_index1 = std::get<1>(entries[0]).second;
      auto wgh1 = std::get<1>(entries[0]).first;
      auto edge_index2 = std::get<1>(entries[1]).second;
      auto wgh2 = std::get<1>(entries[1]).first;
      if (std::make_pair(wgh1, edge_index1) <
          std::make_pair(wgh2, edge_index2)) {
        rctree[src].edge_index = edge_index1;
        rctree[src].wgh = wgh1;
        rctree[src].parent = dst1;
        rctree[src].alt = dst2;
        auto kv1 = std::make_tuple(dst2, std::make_pair(wgh2, edge_index2));
        auto kv2 = std::make_tuple(dst1, std::make_pair(wgh2, edge_index2));
        adj[dst1].insert(kv1, UINT_E_MAX - 1);
        adj[dst2].insert(kv2, UINT_E_MAX - 1);
      } else {
        rctree[src].edge_index = edge_index2;
        rctree[src].wgh = wgh2;
        rctree[src].parent = dst2;
        rctree[src].alt = dst1;
        auto kv1 = std::make_tuple(dst2, std::make_pair(wgh1, edge_index1));
        auto kv2 = std::make_tuple(dst1, std::make_pair(wgh1, edge_index1));
        adj[dst1].insert(kv1, UINT_E_MAX - 1);
        adj[dst2].insert(kv2, UINT_E_MAX - 1);
      }
    }
  };

  auto degree_one =
     parlay::filter(parlay::iota(n), [&](size_t i) { return deg[i] == 1; });
  auto degree_two =
     parlay::filter(parlay::iota(n), [&](size_t i) { return deg[i] == 2; });
  while (rem > 0) {
    // Note that we need both the degree 1 and degree 2 buckets.
    // Only two buckets we care about: {1, 2, everything_else}
    // Dense iterations do help us here...

    // (1) filtering out the initial frontier (degree 1 vertices)
    // Peel:
    // (a) emit the id of the neighbor / peeled (degree1) vertex
    // (b) histogram the ids (do this using semisort / sort)
    // (c) update the degrees, and filter out those that become degree 1 / 2

    // std::cout << "Degree_one.size = " << degree_one.size() << std::endl;
    // std::cout << "Degree_two.size = " << degree_two.size() << std::endl;

    // Compress
    if (degree_two.size() > 0) {
      // Apply compress
      parallel_for(0, degree_two.size(),
                   [&](size_t i) { compress_f1(degree_two[i]); });
      // Deletes edges
      parallel_for(0, degree_two.size(),
                   [&](size_t i) { compress_f2(degree_two[i]); });
      // Adds new edges
      parallel_for(0, degree_two.size(),
                   [&](size_t i) { compress_f3(degree_two[i]); });
    }
    rem -= parlay::count_if(
       degree_two, [&](size_t id) { return rctree[id].round != UINT_E_MAX; });
    round++;

    // Rake
    parallel_for(0, degree_one.size(), [&](size_t i) {
      auto outp = rake_f(degree_one[i]);
      degree_one[i] = outp ? rctree[degree_one[i]].parent : n;
    });
    parlay::sort_inplace(parlay::make_slice(degree_one));
    auto flags = parlay::delayed_seq<bool>(degree_one.size(), [&](size_t i) {
      return (i == 0) || (degree_one[i] != degree_one[i - 1]);
    });
    auto indices = parlay::pack_index(flags);

    rem -= degree_one.size();
    parlay::parallel_for(0, indices.size(), [&](size_t i) {
      auto start = indices[i];
      auto end = (i == indices.size() - 1) ? degree_one.size() : indices[i + 1];
      size_t dst = degree_one[start];
      size_t degree_lost = end - start;
      if (dst < n) {
        deg[dst] -= degree_lost;
      } else {
        rem += degree_lost;
      }
    });
    auto unique_keys = parlay::delayed_seq<size_t>(
       indices.size(),
       [&](size_t i) { return (indices[i] < n) ? degree_one[indices[i]] : n; });
    auto new_degree_one = parlay::filter(
       unique_keys, [&](size_t id) { return (id < n) ? deg[id] == 1 : false; });
    auto new_degree_two = parlay::filter(unique_keys, [&](size_t id) {
      return (id < n) ? deg[id] == 2 : false;
      ;
    });
    degree_one = std::move(new_degree_one);

    auto rem_degree_two =
       parlay::filter(degree_two, [&](size_t id) { return deg[id] == 2; });
    auto combined_degree_two = parlay::append(rem_degree_two, new_degree_two);
    degree_two = std::move(combined_degree_two);
    round++;
    // std::cout << rem << std::endl;
  }
  std::cout << "# Rounds = " << round << std::endl;
  t.next("RC Tree Time");

  return rctree;
}

}   // namespace gbbs
