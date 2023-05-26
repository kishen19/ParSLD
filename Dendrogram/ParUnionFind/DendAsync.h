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
	// using heap = leftist_heap::heap<kv>;
	// using heap = skew_heap::heap<kv>;
	using heap = pairing_heap::heap<kv>;
	size_t n = GA.n;

    // Step 1: Get MST (MSF) of the given graph
	auto mst_edges = MinimumSpanningForest_boruvka::MinimumSpanningForest(GA);
	size_t m = mst_edges.size();

	// Step 2: Initialize (Leftist/Skew/Pairing) Heaps and Union Find
	// Async initialization of heaps: O(nlogh) work and O(hlogh) depth
	timer t;
	t.start();

	auto heaps = gbbs::new_array<heap>(n);

	auto locks = sequence<bool>(n,0);
	auto func = [&](const uintE& u, const W& wgh, const size_t& ind){
		while (!gbbs::atomic_compare_and_swap(&locks[u],false,true)){}
		heaps[u].insert({wgh, ind});
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
		if (!(heaps[u].is_empty())){
			auto min_elem = heaps[u].get_min();
			gbbs::write_add(&is_ready[std::get<1>(min_elem)], 1);
		}
	});
    t.next("Init Time");

	//Step 4: Apply Union-Find in (async) rounds, processing local minima edges in each round
	auto uf = simple_union_find::SimpleUnionAsyncStruct(n);
	auto parents = sequence<uintE>::from_function(m, [&](size_t i){return i;}); //Output Dendrogram
	auto heights = sequence<uintE>(m,0);
	auto proc_edges = parlay::filter(parlay::iota(m), [&](size_t i){return is_ready[i]==2;});
	auto num = proc_edges.size();
	size_t num_iters = 0;
	std::cout << "num = " << num << std::endl;
	while (1){
		if (num == 1){
			is_ready[proc_edges[0]] = 2;
			break;
		}
		parallel_for(0, num, [&](size_t i){
			size_t ind = proc_edges[i];
			auto edge = mst_edges[ind];
			is_ready[ind]++;
			uintE u = simple_union_find::find_compress(std::get<0>(edge), uf.parents);
			uintE v = simple_union_find::find_compress(std::get<1>(edge), uf.parents);
			parlay::par_do(
				[&](){heaps[u].delete_min();},
				[&](){heaps[v].delete_min();}
			);
			uf.unite(u,v);
			uintE x = simple_union_find::find_compress(u, uf.parents);
			uintE y = x^u^v;
			heaps[x].merge(&heaps[y]);
			if (heaps[x].is_empty()){
				proc_edges[i] = m;
			} else{
				size_t temp = std::get<1>(heaps[x].get_min());
				heights[temp] = std::max(heights[temp], heights[ind]+1);
				parents[ind] = temp;
				gbbs::write_add(&is_ready[temp], 1);
				if (gbbs::atomic_compare_and_swap(&is_ready[temp], 2, 3)){
					proc_edges[i] = temp;
				} else{
					proc_edges[i] = m;
				}
			}
		});
		auto new_edges = parlay::filter(proc_edges, [&](auto i){return (i!=m);});
		proc_edges = new_edges;
		num = proc_edges.size();
		// std::cout << "num = " << num  << std::endl;
		num_iters++;
	}
	t.next("Dendrogram Stage 1 Time");
	auto rem_edges = parlay::filter(parlay::iota(m), [&](size_t i){return is_ready[i] < 3;});
	// auto rem_edges = parlay::pack(mst_edges, parlay::delayed_seq<bool>(m, [&](size_t i){return is_ready[i]<3;}));
	sort_inplace(rem_edges, [&](auto e1, auto e2){
		auto p1 = std::make_tuple(std::get<2>(mst_edges[e1]), e1);
		auto p2 = std::make_tuple(std::get<2>(mst_edges[e2]), e2);
		return (p1 < p2);
	});
	num = rem_edges.size();
	parallel_for(0, num-1, [&](size_t i){
		auto ind = rem_edges[i];
		auto temp = rem_edges[i+1];
		heights[temp] = num_iters + i + 1; //std::max(heights[temp], heights[i]+1);
		parents[ind] = temp;
	});
	t.next("Dendrogram Stage 2 Time");
	double tt = t.total_time();

	std::cout << "Remaining Edges = " << num << std::endl;
	std::cout << "Num Iters = " << num_iters << std::endl;
	std::cout << std::endl << "=> Dendrogram Height = " << parlay::reduce_max(heights) << std::endl;

	// return parents;
	return tt;
}

}  // namespace gbbs