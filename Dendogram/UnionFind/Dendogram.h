#pragma once

#include "gbbs/gbbs.h"
#include "benchmarks/Connectivity/SimpleUnionAsync/Connectivity.h"
#include "benchmarks/MinimumSpanningForest/Boruvka/MinimumSpanningForest.h"
#include "Dendogram/heaps/leftist_heap.h"
#include "Dendogram/heaps/skew_heap.h"

namespace gbbs {

template <class Graph>
void Dendogram(Graph& GA){
	using W = typename Graph::weight_type;
	size_t n = GA.n;

	// Step 1: Get MST (MSF) of the given graph
	auto mst_edges = MinimumSpanningForest_boruvka::MinimumSpanningForest(GA);
	size_t m = mst_edges.size();

	// Step 2: Initialize (Leftist) Heaps and Union Find
	// Async initialization of heaps: O(nlogh) work and O(hlogh) depth
	auto heaps = sequence<leftist_heap::leftist_heap<W>*>(n);
	// auto heaps = sequence<skew_heap::skew_heap<W>*>(n);
	parallel_for(0, n, [&](size_t i){
		heaps[i] = new leftist_heap::leftist_heap<W>();
		// heaps[i] = new skew_heap::skew_heap<W>();
	});
	auto locks = sequence<bool>(n,0);
	auto func = [&](const uintE& u, const W& wgh, size_t ind){
		while (!gbbs::atomic_compare_and_swap(&locks[u],false,true)){}
		heaps[u]->insert(ind, wgh);
		locks[u] = false;
	};
	parallel_for(0, m, [&](size_t i){
		auto s = std::get<0>(mst_edges[i]);
		auto t = std::get<1>(mst_edges[i]);
		auto w = std::get<2>(mst_edges[i]);
		parlay::par_do(
			[&](){func(s,w,i);}, 
			[&](){func(t,w,i);}
		);
	});

	auto uf = simple_union_find::SimpleUnionAsyncStruct(n);
	auto parents = sequence<uintE>::from_function(m, [&](size_t i){return i;}); //Output Dendogram
	auto is_ready = sequence<bool>(m,false);
	parallel_for(0, m, [&](size_t i){
		auto u = std::get<0>(mst_edges[i]);
		auto v = std::get<1>(mst_edges[i]);
		auto min_elem_u = std::get<1>(heaps[u]->find_min());
		auto min_elem_v = std::get<1>(heaps[v]->find_min());
		if (min_elem_u == i && min_elem_v == i){
			is_ready[i] = true;
		}
	});

	parallel_for(0, m, [&](size_t i){
		auto edge = mst_edges[i];
		auto ind = i;
		while (1) {
			if (gbbs::atomic_compare_and_swap(&is_ready[ind], true, false)){
				auto u = simple_union_find::find_compress(std::get<0>(edge), uf.parents);
				auto v = simple_union_find::find_compress(std::get<1>(edge), uf.parents);
				parlay::par_do(
					[&](){heaps[u]->delete_min();},
					[&](){heaps[v]->delete_min();}
				);
				uf.unite(u,v);
				size_t temp;
				if (simple_union_find::find_compress(u, uf.parents) == v){
					heaps[v]->merge(heaps[u]);
					if (heaps[v]->is_empty()){ break;}
					temp = std::get<1>(heaps[v]->find_min());
				} else {
					heaps[u]->merge(heaps[v]);
					if (heaps[u]->is_empty()){ break;}
					temp = std::get<1>(heaps[u]->find_min());
				}
				parents[ind] = temp;
				ind = temp;
				edge = mst_edges[ind];

				u = simple_union_find::find_compress(std::get<0>(edge), uf.parents);
				v = simple_union_find::find_compress(std::get<1>(edge), uf.parents);
				auto min_elem_u = std::get<1>(heaps[u]->find_min());
				auto min_elem_v = std::get<1>(heaps[v]->find_min());
				if (min_elem_u == ind && min_elem_v == ind){
					gbbs::atomic_compare_and_swap(&is_ready[ind], false, true);
				} else{
					break;
				}
			} else{
				break;
			}
		}
	});

	for (size_t i=0; i<m; i++){
		std::cout << i << " -> " << parents[i] << std::endl;
	}
}

}  // namespace gbbs
