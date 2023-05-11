#pragma once

#include "gbbs/gbbs.h"
#include "benchmarks/Connectivity/SimpleUnionAsync/Connectivity.h"
#include "benchmarks/MinimumSpanningForest/Boruvka/MinimumSpanningForest.h"
#include "Dendrogram/heaps/leftist_heap.h"
#include "Dendrogram/heaps/skew_heap.h"

namespace gbbs {

template <class Graph>
double DendrogramSeqUF(Graph& GA){
	size_t n = GA.n;

	timer t;
	t.start();
	// Step 1: Get MST (MSF) of the given graph
	auto mst_edges = MinimumSpanningForest_boruvka::MinimumSpanningForest(GA);
	size_t m = mst_edges.size();
	t.next("MST Time");

	// Step 2: Sorting the edges by (weight, index)
	auto indices = sequence<size_t>::from_function(m, [&](size_t i){return i;});
	auto comp = [&](const size_t& a, const size_t& b){
		auto w_a = std::get<2>(mst_edges[a]);
		auto w_b = std::get<2>(mst_edges[b]);
		if (w_a < w_b){
			return true;
		} else if (w_a == w_b){
			return (a < b);
		} else {
			return false;
		}
	};
	parlay::sort_inplace(indices, comp);
	t.next("Init Time");

	// Step 3: Applying Union Find to the sorted sequence of edges
	auto uf = simple_union_find::SimpleUnionAsyncStruct(n);
	auto aux = sequence<size_t>(n,n); // extra info required for assigning parents
	auto parents = sequence<uintE>::from_function(m, [&](size_t i){return i;}); //Output Dendrogram
	auto heights = sequence<uintE>(m,0); // Heights of every node in the dendrogram

	for(size_t ind = 0; ind < m; ind++){
		auto i = indices[ind];
		auto u = simple_union_find::find_compress(std::get<0>(mst_edges[i]), uf.parents);
		auto v = simple_union_find::find_compress(std::get<1>(mst_edges[i]), uf.parents);
		uintE height = 0;
		if (aux[u] < m){
			parents[aux[u]] = i;
			height = heights[aux[u]] + 1;
			aux[u] = n;
		}
		if (aux[v] < m){
			parents[aux[v]] = i;
			height = std::max(height, heights[aux[v]] + 1);
			aux[v] = n;
		}
		uf.unite(u,v);
		aux[simple_union_find::find_compress(u, uf.parents)] = i;
		heights[i] = height;
	};
	t.next("Dendrogram Time");
	double tt = t.total_time();
	std::cout << std::endl << "=> Dendrogram Height = " << parlay::reduce_max(heights) << std::endl;

	// return parents;
	return tt;
}

}  // namespace gbbs
