#pragma once

#include "gbbs/gbbs.h"
#include "benchmarks/Connectivity/SimpleUnionAsync/Connectivity.h"
#include "benchmarks/MinimumSpanningForest/Boruvka/MinimumSpanningForest.h"
#include "Dendogram/heaps/leftist_heap.h"
#include "Dendogram/heaps/skew_heap.h"

namespace gbbs {

template <class Graph>
double DendogramParUF(Graph& GA){
	using W = typename Graph::weight_type;
	using kv = std::pair<W, size_t>;
	size_t n = GA.n;

	// Step 1: Get MST (MSF) of the given graph
	auto mst_edges = MinimumSpanningForest_boruvka::MinimumSpanningForest(GA);
	size_t m = mst_edges.size();

	// Step 2: Initialize (Leftist/Skew/Pairing) Heaps and Union Find
	// Async initialization of heaps: O(nlogh) work and O(hlogh) depth
	// auto heaps = sequence<leftist_heap::leftist_heap<kv>*>(n);
	timer t;
	t.start();
	auto heaps = sequence<skew_heap::skew_heap<kv>*>(n);
	parallel_for(0, n, [&](size_t i){
		// heaps[i] = new leftist_heap::leftist_heap<kv>();
		heaps[i] = new skew_heap::skew_heap<kv>();
	});
	auto locks = sequence<bool>(n,0);
	auto func = [&](const uintE& u, const W& wgh, const size_t& ind){
		while (!gbbs::atomic_compare_and_swap(&locks[u],false,true)){}
		heaps[u]->insert({wgh, ind});
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

	// Step 3: Find initial "local minimum" edges
	// is_ready[ind] = 2 => edge is a local minimum
	auto is_ready = sequence<int>(m,0);
	parallel_for(0, n, [&](size_t u){
		if (heaps[u]->size > 0){
			auto min_elem = heaps[u]->find_min();
			gbbs::write_add(&is_ready[std::get<1>(min_elem)], 1);
		}
	});

	//Step 4: Apply Union-Find in (async) rounds, processing local minima edges in each round
	auto uf = simple_union_find::SimpleUnionAsyncStruct(n);
	auto parents = sequence<uintE>::from_function(m, [&](size_t i){return i;}); //Output Dendogram
	auto heights = sequence<uintE>(m,0);
	parallel_for(0, m, [&](size_t i){
		auto edge = mst_edges[i];
		auto ind = i;
		while (1) {
			if (gbbs::atomic_compare_and_swap(&is_ready[ind], 2, 0)){
				// All steps here are "contention-free" since other procs cannot process 
				// the edges incident on u and v until the edge (u,v) is processed
				auto u = simple_union_find::find_compress(std::get<0>(edge), uf.parents);
				auto v = simple_union_find::find_compress(std::get<1>(edge), uf.parents);
				parlay::par_do(
					[&](){heaps[u]->delete_min();},
					[&](){heaps[v]->delete_min();}
				);
				uf.unite(u,v);
				size_t temp;
				if (uf.parents[v] == v){
					heaps[v]->merge(heaps[u]);
					if (heaps[v]->is_empty()){ break;}
					temp = std::get<1>(heaps[v]->find_min());
				} else {
					heaps[u]->merge(heaps[v]);
					if (heaps[u]->is_empty()){ break;}
					temp = std::get<1>(heaps[u]->find_min());
				}
				gbbs::write_max(&heights[temp], heights[ind]+1);
				parents[ind] = temp;
				ind = temp;
				edge = mst_edges[ind];
				gbbs::write_add(&is_ready[ind], 1);
			} else{ // Edge is not ready, will be processed by neighbor in the future
				break;
			}
		}
	});

	double tt = t.stop();
	std::cout << std::endl << "=> Dendogram Height = " << parlay::reduce_max(heights) << std::endl;

	return tt;
}

}  // namespace gbbs
