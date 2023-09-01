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
	using edge = std::tuple<uintE, W, size_t>;
	// using heap = leftist_heap::heap<kv>;
	// using heap = skew_heap::heap<kv>;
	using heap = pairing_heap::heap<kv>;
	// using heap = pairing_heap::block_heap<kv>;
	size_t n = GA.n;

	// Step 1: Get MST (MSF) of the given graph
	auto mst_edges = MinimumSpanningForest_boruvka::MinimumSpanningForest(GA);
	size_t m = mst_edges.size();

	// Step 2: Initialize (Leftist/Skew/Pairing/Block) Heaps and Union Find
	timer t;
	t.start();
	
	auto new_edges = sequence<edge>::uninitialized(2*m);
	parallel_for(0, m, [&](size_t i){
		uintE u,v;
		W w;
		std::tie(u,v,w) = mst_edges[i];
		new_edges[2*i] = std::make_tuple(u,w,i);
		new_edges[2*i+1] = std::make_tuple(v,w,i);
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
	auto heap_stuff = gbbs::new_array_no_init<kv>(2*m);
	parallel_for(0, 2*m, [&](size_t i){
		heap_stuff[i] = std::make_tuple(std::get<1>(new_edges[i]), std::get<2>(new_edges[i]));
	});

	auto heaps = gbbs::new_array_no_init<heap>(n);
    parallel_for(0, n, [&](size_t i){
		heaps[i].init(heap_stuff + offsets[i], sizes[i]);
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
	t.next("Heap Init Time");

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
			uintE u = simple_union_find::find_compress(std::get<0>(edge), uf.parents.begin());
			uintE v = simple_union_find::find_compress(std::get<1>(edge), uf.parents.begin());
			parlay::par_do(
				[&](){heaps[u].delete_min();},
				[&](){heaps[v].delete_min();}
			);
			uf.unite(u,v);
			uintE x = simple_union_find::find_compress(u, uf.parents.begin());
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
