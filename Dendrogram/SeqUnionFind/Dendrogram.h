#pragma once

#include "gbbs/gbbs.h"

#include "Dendrogram/common/union_find.h"

namespace gbbs {

template <class IdType, class Graph>
auto DendrogramSeqUF_impl(Graph& GA, bool debug = false) {
	using W = typename Graph::weight_type;

	timer t;
	t.start();

  parlay::parallel_for(0, GA.n, [&] (uintE u) {
    auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
      if (u == v) {
        std::cout << "Self loop: " << u << " " << v << std::endl;
        exit(-1);
      }
    };
    GA.get_vertex(u).out_neighbors().map(map_f, false);
  });

	// Step 1: Preprocess
	// Index the Edges
	size_t n = GA.n; size_t m = GA.m/2;
  auto deg = sequence<size_t>::from_function(n, [&](uintE u){
      return GA.get_vertex(u).out_degree();
  });
  auto edges = sequence<std::tuple<uintE, uintE, W>>::uninitialized(2*m);
  auto offsets = parlay::scan(deg).first;
  auto map_f = [&](const uintE& src, const uintE& dst, const W& wgh, const uintE& ind){
      edges[offsets[src] + ind] = {std::min(src,dst), std::max(src,dst), wgh};
  };
  parallel_for(0, n, [&](uintE u) {
      GA.get_vertex(u).out_neighbors().map_with_index(map_f);
  });
  parlay::sort_inplace(edges);
  for (size_t i=0; i<edges.size(); ++i) {
    auto [u, v, wgh] = edges[i];
    if (u == v) {
      std::cout << "Self loop in edges:" << " u = " << u << " v = " << v << std::endl;
      exit(-1);
    }
  }
	if (debug) {t.next("Preprocess Time");}

	// Step 2: Sorting the edges by (weight, index)
	auto proc_edges = sequence<std::pair<W,uintE>>::from_function(m, [&](size_t i){
		return std::make_pair(std::get<2>(edges[2*i]), i);
	});
//  for (size_t i=0; i<proc_edges.size(); ++i) {
//    std::cout << std::get<0>(edges[2*i]) << " " << std::get<1>(edges[2*i]) << " " << std::get<2>(edges[2*i]) << " i = " << i << std::endl;
//  }
	parlay::sort_inplace(proc_edges);
	if (debug) {t.next("Sorting Edges Time");}

	// Step 3: Applying Union Find to the sorted sequence of edges
	auto uf = union_find(n);
	auto parent = sequence<uintE>::from_function(m, [&](uintE i){ return i; });
	auto aux = sequence<uintE>(n, n); // extra info required for assigning parents
	// auto heights = sequence<uintE>(m,0); // Heights of every node in the dendrogram

  size_t total_work = 0;

	for(size_t ind = 0; ind < m; ind++) {
		size_t i = std::get<1>(proc_edges[ind]);
		auto [u, tw1] = uf.find_compress(std::get<0>(edges[2*i]));
		auto [v, tw2] = uf.find_compress(std::get<1>(edges[2*i]));
		// uintE height = 0;
		total_work += tw1 + tw2;
		if (aux[u] < n){
			parent[aux[u]] = i;
			// height = heights[aux[u]] + 1;
			aux[u] = n;
		}
		if (aux[v] < n){
			parent[aux[v]] = i;
			// height = std::max(height, heights[aux[v]] + 1);
			aux[v] = n;
		}
		auto [w, tw3] = uf.unite(u,v);
		total_work += tw3;
		aux[w] = i;
		// heights[i] = height;
	}
	if (debug) {
		t.next("Dendrogram Time");
		std::cout << "Total dendrogram work: " << total_work << std::endl;
		// std::cout << std::endl << "=> Dendrogram Height = " << parlay::reduce_max(heights) << std::endl;
	}
  return parent;
}

template <class Graph>
auto DendrogramSeqUF(Graph& GA, bool debug = false) {
  if (GA.n >= std::numeric_limits<int32_t>::max()) {
    return DendrogramSeqUF_impl<size_t>(GA, debug);
  } else {
    return DendrogramSeqUF_impl<uint32_t>(GA, debug);
  }
}



}  // namespace gbbs
