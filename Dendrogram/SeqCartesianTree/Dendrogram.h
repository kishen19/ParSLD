#pragma once

#include "gbbs/gbbs.h"

namespace gbbs {

template <class W>
struct node { uintE parent; uintE index; W value; };

template <class W>
void cartesianTree(node<W>* Nodes, uintE n) {
  if (n < 1) return;
  for (size_t i=1; i<n; ++i) {
    uintE cur_idx = i-1;
    uintE cur_child = UINT_E_MAX;
    while (cur_idx != UINT_E_MAX) {
      if (Nodes[i].value > Nodes[cur_idx].value) {
        cur_child = cur_idx;
        cur_idx = Nodes[cur_idx].parent;
      } else {
        break;
      }
    }
    if (cur_child != UINT_E_MAX) {
      Nodes[cur_child].parent = i;
    }
    Nodes[i].parent = cur_idx;
  }
}

template <class Graph>
auto SeqCartesianTree_runner(Graph& GA, bool debug = false) {
	using W = typename Graph::weight_type;
	size_t n = GA.n; size_t m = GA.m/2;

  // Let us assume that the input is already list-ranked.
	auto parent = sequence<uintE>::from_function(m, [&](uintE i){ return i; });
	auto nodes = sequence<node<W>>::uninitialized(m+1);
  parlay::parallel_for(0, m+1, [&] (size_t i) {
    nodes[i].parent = UINT_E_MAX;
    nodes[i].value = i;
  });
  parlay::parallel_for(0, n, [&] (size_t i) {
    auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
      if (u < v) {
        nodes[i+1].value = wgh;
      }
    };
    GA.get_vertex(i).out_neighbors().map(map_f, false);
  });

  cartesianTree<W>(nodes.begin(), n);

  parlay::parallel_for(0, m, [&] (size_t i) {
    auto parent_i = nodes[i+1].parent;
    if (parent_i == UINT_E_MAX) { // this is the root
      parent[i] = i;
    } else {
      parent[i] = parent_i - 1;
    }
  });

//  std::cout << "Cartesian tree end: m = " << m << std::endl;
//  for (size_t i=0; i<m; ++i) {
//    std::cout << parent[i] << std::endl;
//  }
//  std::cerr << "Returning parent, size = " << parent.size() << std::endl;
  return parent;
}

}   // namespace gbbs
