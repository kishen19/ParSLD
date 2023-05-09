#pragma once

#include "gbbs/gbbs.h"
#include "benchmarks/Connectivity/SimpleUnionAsync/Connectivity.h"
#include "benchmarks/MinimumSpanningForest/Boruvka/MinimumSpanningForest.h"
#include "Dendrogram/heaps/leftist_heap.h"
#include "Dendrogram/heaps/skew_heap.h"
#include "Dendrogram/heaps/pairing_heap.h"

namespace gbbs {

template <class Graph>
double DendrogramParUF(Graph& GA){
	using W = typename Graph::weight_type;
	using kv = std::tuple<W, size_t>;
	using edge = std::tuple<uintE, uintE, W, size_t>;
	size_t n = GA.n;

	timer t;
	t.start();

	// Step 1: Get MST (MSF) of the given graph
	auto mst_edges = MinimumSpanningForest_boruvka::MinimumSpanningForest(GA);
	size_t m = mst_edges.size();
	t.next("MST Time");

	// Step 2: Initialize (Leftist/Skew/Pairing) Heaps and Union Find
	auto new_edges = sequence<edge>::uninitialized(2*m);
	parallel_for(0, m, [&](size_t i){
		uintE u,v;
		W w;
		std::tie(u,v,w) = mst_edges[i];
		new_edges[2*i] = std::make_tuple(u,v,w,i);
		new_edges[2*i+1] = std::make_tuple(v,u,w,i);
	});
	sort_inplace(new_edges);
	
	auto sizes = sequence<uintT>(n, 0);
	sizes[std::get<0>(new_edges[2*m-1])] = 2*m;
    parallel_for(0, 2*m-1, [&](size_t i) {
		if (std::get<0>(new_edges[i]) != std::get<0>(new_edges[i+1])){
			sizes[std::get<0>(new_edges[i])] = i+1;
		}
	});
	parallel_for(1, 2*m, [&](size_t i) {
		if (std::get<0>(new_edges[i]) != std::get<0>(new_edges[i-1])){
			sizes[std::get<0>(new_edges[i])] -= i;
		}
	});
	
    auto offsets = parlay::scan(make_slice(sizes)).first;
	auto heap_stuff = parlay::sequence<kv>::from_function(2*m, [&](size_t i){
		return std::make_tuple(std::get<2>(new_edges[i]), std::get<3>(new_edges[i]));
	});

    // auto heaps = sequence<leftist_heap::leftist_heap<kv>>(n);
	// auto heaps = sequence<skew_heap::skew_heap<kv>>(n);
	auto heaps = sequence<pairing_heap::pairing_heap<kv>>(n);

    parallel_for(0, n, [&](size_t i){
        heaps[i].create_heap(heap_stuff.cut(offsets[i], offsets[i]+sizes[i]));
    });
	t.next("Heap Init Time");

	// Step 3: Find initial "local minimum" edges
	// is_ready[ind] = 2 => edge is a local minimum
	auto is_ready = sequence<int>(m,0);
	parallel_for(0, n, [&](size_t u){
		if (!(heaps[u].is_empty())){
			auto min_elem = heaps[u].find_min();
			gbbs::write_add(&is_ready[std::get<1>(min_elem)], 1);
		}
	});
	t.next("Finding Local Min Edges Time");

	//Step 4: Apply Union-Find in (async) rounds, processing local minima edges in each round
	auto uf = simple_union_find::SimpleUnionAsyncStruct(n);
	auto parents = sequence<uintE>::from_function(m, [&](size_t i){return i;}); //Output Dendrogram
	auto heights = sequence<uintE>(m,0);
	parallel_for(0, m, [&](size_t i){
		auto edge = mst_edges[i];
		size_t ind = i, temp;
		uintE u,v,x,y;
		while (1) {
			if (gbbs::atomic_compare_and_swap(&is_ready[ind], 2, 0)){
				// All steps here are "contention-free" since other procs cannot process 
				// the edges incident on u and v until the edge (u,v) is processed
				u = simple_union_find::find_compress(std::get<0>(edge), uf.parents);
				v = simple_union_find::find_compress(std::get<1>(edge), uf.parents);
				parlay::par_do(
					[&](){heaps[u].delete_min();},
					[&](){heaps[v].delete_min();}
				);
				uf.unite(u,v);
				x = simple_union_find::find_compress(u, uf.parents);
				y = x^u^v;
				heaps[x].merge(&heaps[y]);
				if (heaps[x].is_empty()){
					break;
				}
				temp = std::get<1>(heaps[x].find_min());
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
	t.next("Dendrogram Time");

	std::cout << std::endl << "=> Dendrogram Height = " << parlay::reduce_max(heights) << std::endl;

	double tt = t.stop();

	// return parents;
	return tt;
}

}  // namespace gbbs
