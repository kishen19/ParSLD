#pragma once

#include "gbbs/gbbs.h"

#include "Dendrogram/common/union_find.h"

namespace gbbs {

template <class IdType, class Graph>
auto DendrogramSeqUF_impl(Graph& GA, bool debug = false) {
	using W = typename Graph::weight_type;
	timer t;
	t.start();

	// Step 1: Preprocess
	// Index the Edges
	size_t n = GA.n; size_t m = GA.m/2;
  auto deg = sequence<size_t>::from_function(n, [&](uintE u){
      return GA.get_vertex(u).out_degree();
  });
  auto edges = sequence<std::tuple<W, uintE, uintE>>::uninitialized(2*m);
  auto offsets = parlay::scan(deg).first;
  auto map_f = [&](const uintE& src, const uintE& dst, const W& wgh, const uintE& ind){
      edges[offsets[src] + ind] = {wgh, std::min(src,dst), std::max(src,dst)};
  };
  parallel_for(0, n, [&](uintE u) {
      GA.get_vertex(u).out_neighbors().map_with_index(map_f);
  });
  parlay::sort_inplace(edges);
	if (debug) {t.next("Preprocess Time");}

	// Step 2: Applying Union Find to the sorted sequence of edges
	auto uf = union_find(n);
	auto dendrogram = sequence<std::pair<uintE,W>>::from_function(n+m, [&](uintE i){ 
		return std::make_pair(UINT_E_MAX, 0); });
	auto aux = sequence<uintE>(n, m); // extra info required for assigning parents
	// auto heights = sequence<uintE>(m,0); // Heights of every node in the dendrogram
	for(size_t i = 0; i < m; i++) {
		auto [u, tw1] = uf.find_compress(std::get<1>(edges[2*i]));
		auto [v, tw2] = uf.find_compress(std::get<2>(edges[2*i]));
		auto wgh = std::get<0>(edges[2*i]);
		// uintE height = 0;
		if (aux[u] < m){
			dendrogram[n+aux[u]] = std::make_pair(n+i, wgh);
			// height = heights[aux[u]] + 1;
			aux[u] = m;
		} else{
			dendrogram[u] = std::make_pair(n+i, wgh);
		}
		if (aux[v] < m){
			dendrogram[n+aux[v]] = std::make_pair(n+i, wgh);
			// height = std::max(height, heights[aux[v]] + 1);
			aux[v] = m;
		} else{
			dendrogram[v] = std::make_pair(n+i, wgh);
		}
		auto [w, tw3] = uf.unite(u,v);
		aux[w] = i;
		// heights[i] = height;
	}
	if (debug) {
		t.next("Dendrogram Time");
		// std::cout << std::endl << "=> Dendrogram Height = " << parlay::reduce_max(heights) << std::endl;
	}
  return dendrogram;
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
