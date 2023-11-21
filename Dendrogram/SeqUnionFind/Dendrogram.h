#pragma once

#include "gbbs/gbbs.h"

#include "Dendrogram/common/union_find.h"

namespace gbbs {

template <class Graph>
auto DendrogramSeqUF(Graph& GA){
	using W = typename Graph::weight_type;

	timer t;
	t.start();

	// Step 1: Preprocess
	// Index the Edges
	auto n = GA.n; auto m = GA.m/2;
    auto deg = sequence<uintE>::from_function(n, [&](uintE u){
        return GA.get_vertex(u).out_degree();
    });
    auto edges = sequence<std::tuple<uintE, uintE, W>>::uninitialized(2*m);
    auto offsets = parlay::scan(deg).first;
    auto map_f = [&](const uintE& src, const uintE& dst, const W& wgh, const uintE& ind){
        edges[offsets[src] + ind] = {std::min(src,dst), std::max(src,dst), wgh};
    };
    parallel_for(0, n, [&](uintE u){
        GA.get_vertex(u).out_neighbors().map_with_index(map_f);
    });
    parlay::sort_inplace(edges);
	t.next("Preprocess Time");

	// Step 2: Sorting the edges by (weight, index)
	auto proc_edges = sequence<std::pair<W,uintE>>::from_function(m, [&](uintE i){
		return std::make_pair(std::get<2>(edges[2*i]), i);
	});
	parlay::sort_inplace(proc_edges);
	t.next("Sorting Edges Time");

	// Step 3: Applying Union Find to the sorted sequence of edges
	auto uf = union_find(n);
	auto parent = sequence<uintE>::from_function(m, [&](uintE i){ return i; });
	auto aux = sequence<uintE>(n, n); // extra info required for assigning parents
	// auto heights = sequence<uintE>(m,0); // Heights of every node in the dendrogram

  size_t total_work = 0;

	for(uintE ind = 0; ind < m; ind++){
		auto i = std::get<1>(proc_edges[ind]);
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
	};
	t.next("Dendrogram Time");
  std::cout << "Total dendrogram work: " << total_work << std::endl;

	// for (size_t i=0; i<m; i++){
    //     std::cout << parent[i] << " ";
    // }
    // std::cout << std::endl;

	// std::cout << std::endl << "=> Dendrogram Height = " << parlay::reduce_max(heights) << std::endl;
  return parent;
}

}  // namespace gbbs
