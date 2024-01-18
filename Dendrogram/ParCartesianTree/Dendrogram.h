#pragma once

#include "gbbs/gbbs.h"

namespace gbbs {
namespace par_ct {

template <class W>
struct node { uintE parent; W value; };

template <class W>
void merge(node<W>* N, uintE left, uintE right) {
  uintE head;
  if (N[left].value < N[right].value) {
    head = left;
    left = N[left].parent;
  } else {
    head = right;
    right = N[right].parent;
  }

  while (1) {
    if (left == 0) {
      N[head].parent = right;
      break;
    }
    if (right == 0) {
      N[head].parent = left;
      break;
    }
    if (N[left].value < N[right].value) {
      N[head].parent = left;
      left = N[left].parent;
    } else {
      N[head].parent = right;
      right = N[right].parent;
    }
    head = N[head].parent;
  }
}

template <class W>
void cartesianTree(node<W>* Nodes, uintE s, uintE n) {
  if (n < 2) {
    return;
  }
  if (n == 2) {
    if (Nodes[s].value < Nodes[s + 1].value) {
      Nodes[s].parent = s + 1;
    } else {
      Nodes[s + 1].parent = s;
    }
    return;
  }
  if (n > 1000) {
    auto l = [&] () { cartesianTree(Nodes, s, n / 2); };
    auto r = [&] () { cartesianTree(Nodes, s + n / 2, n - n / 2); };
    parlay::par_do(l, r);
  } else {
    cartesianTree(Nodes, s, n / 2);
    cartesianTree(Nodes, s + n / 2, n - n / 2);
  }
  merge(Nodes, s + n / 2 - 1, s + n / 2);
}

}  // namespace par_ct

template <class Graph>
auto ParCartesianTree_runner(Graph& GA, bool debug = false) {
	using W = typename Graph::weight_type;
	size_t n = GA.n; size_t m = GA.m/2;

  // Let us assume that the input is already list-ranked.
	auto parent = sequence<uintE>::from_function(m, [&](uintE i){ return i; });
	auto nodes = sequence<par_ct::node<W>>::uninitialized(m+1);
  parlay::parallel_for(0, m+1, [&] (size_t i) {
    nodes[i].parent = 0;
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

  if (debug) {
    std::cout << "Cartesian tree start: m = " << m << std::endl;
    for (size_t i=0; i<m; ++i) {
      std::cout << nodes[i].parent << " " << nodes[i].value << std::endl;
    }
  }
  par_ct::cartesianTree<W>(nodes.begin(), 0, n);

  parlay::parallel_for(0, m, [&] (size_t i) {
    auto parent_i = nodes[i+1].parent;
    if (parent_i == 0) { // this is the root
      parent[i] = i;
    } else {
      parent[i] = parent_i - 1;
    }
  });

  if (debug) {
    std::cout << "Cartesian tree end: m = " << m << std::endl;
    for (size_t i=0; i<m; ++i) {
      std::cout << parent[i] << " " << nodes[i+1].value << std::endl;
    }
  }

  return parent;
}

}   // namespace gbbs
