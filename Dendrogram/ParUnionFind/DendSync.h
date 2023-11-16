#pragma once

#include "gbbs/gbbs.h"
#include "Dendrogram/common/leftist_heap.h"
#include "Dendrogram/common/skew_heap.h"
#include "Dendrogram/common/pairing_heap.h"
#include "Dendrogram/common/union_find.h"

namespace gbbs {

template <class Graph>
double DendrogramParUF(Graph& GA){
	using W = typename Graph::weight_type;
	using kv = std::pair<W, uintE>;
	// using heap = leftist_heap::heap<kv>;
	using heap = skew_heap::heap<kv>;
	// using heap = pairing_heap::heap<kv>;
	// using heap = pairing_heap::block_heap<kv>;

	timer t;
	t.start();

	// Step 1: Preprocess
	// Index the Edges
	auto n = GA.n; auto m = GA.m/2;
	auto offsets = sequence<uintE>::from_function(
		n, [&](uintE u) { return GA.get_vertex(u).out_degree(); });
	parlay::scan_inplace(offsets);

	auto edges = sequence<std::tuple<uintE, uintE, uintE, W>>::uninitialized(2 * m);
	auto offset_f = [&](const uintE& src, const uintE& dst, const W& wgh,
						const uintE& ind) {
		edges[offsets[src] + ind] = {std::min(src, dst), std::max(src, dst), offsets[src] + ind, wgh};
		// ngh's offset in src is ind.
	};
	parallel_for(0, n, [&](uintE u) {
		GA.get_vertex(u).out_neighbors().map_with_index(offset_f);
	});
	parlay::sort_inplace(edges);
	t.next("Sort edges");

	// Neighbors of vertex i at offsets[i] -- offsets[i+1]:
	// We store:
	// - weight of the edge (redundant, but saves a cache miss)
	// - edge_index of this edge (the original identity of this edge)
	auto neighbors = parlay::sequence<kv>::uninitialized(2*m);

	parlay::parallel_for(0, m, [&](uintE i) {
		auto [u, v, ind1, wgh] = edges[2*i];
		auto [_u, _v, ind2, _wgh] = edges[2*i+1];
		// auto [_, wgh] = GA.get_vertex(u).out_neighbors().get_ith_neighbor(ind1-offsets[u]);
		// smaller of the two inds corresponds to the smaller vertex id.
		neighbors[ind1] = {wgh, i};
		neighbors[ind2] = {wgh, i};
	});
	t.next("Initialize Neighbors pairs");

	// Step 2: Initialize (Leftist/Skew/Pairing/Block) Heaps and Union Find
	auto heaps = gbbs::new_array_no_init<heap>(n);
    parallel_for(0, n, [&](size_t i){
		heaps[i].init(neighbors.cut(offsets[i], ((i == n-1)? 2*m : offsets[i+1])));
    });
	t.next("Heap Init Time");

	// Step 3: Find initial "local minimum" edges
	// is_ready[ind] = 2 => edge is a local minimum
	auto is_ready = sequence<int>(m,0);
	parallel_for(0, n, [&](size_t u){
		if (!(heaps[u].is_empty())){
			auto min_elem = heaps[u].get_min();
			gbbs::write_add(&is_ready[min_elem.second], 1);
		}
	});
	t.next("Is_Ready Time");

	//Step 4: Apply Union-Find in (sync) rounds, processing local minima edges in each round
	auto uf = union_find(n);
	auto parent = sequence<uintE>::from_function(m, 
		[&](uintE i){ return i; }); //Output Dendrogram
	// auto heights = sequence<uintE>(m,0);
	auto proc_edges = parlay::filter(parlay::iota<uintE>(m), 
		[&](size_t i){ return is_ready[i]==2; });
	auto num = proc_edges.size();
	uintE num_iters = 0;
	std::cout << "num = " << num << std::endl;
	while (1){
		if (num == 1){
			is_ready[proc_edges[0]] = 2;
			break;
		}
		parallel_for(0, num, [&](uintE ind){
			uintE i = proc_edges[ind];
			is_ready[i]++;
			uintE u = uf.find_compress(std::get<0>(edges[2*i]));
			uintE v = uf.find_compress(std::get<1>(edges[2*i]));
			parlay::par_do(
				[&](){heaps[u].delete_min();},
				[&](){heaps[v].delete_min();}
			);
			auto w = uf.unite(u,v);
			uintE t = w^u^v;
			heaps[w].merge(&heaps[t]);
			if (heaps[w].is_empty()){
				proc_edges[ind] = m;
			} else{
				uintE temp = heaps[w].get_min().second;
				// heights[temp] = std::max(heights[temp], heights[i]+1);
				parent[i] = temp;
				gbbs::write_add(&is_ready[temp], 1);
				if (gbbs::atomic_compare_and_swap(&is_ready[temp], 2, 3)){
					proc_edges[ind] = temp;
				} else{
					proc_edges[ind] = m;
				}
			}
		});
		proc_edges = parlay::filter(proc_edges, [&](auto i){return (i != m);});
		num = proc_edges.size();
		// std::cout << "num = " << num  << std::endl;
		num_iters++;
	}
	t.next("Dendrogram Stage 1 Time");
	auto rem_edges = parlay::filter(parlay::iota<uintE>(m), 
		[&](uintE i){ return is_ready[i] < 3; });
	sort_inplace(rem_edges, [&](auto e1, auto e2){
		auto p1 = std::make_pair(std::get<3>(edges[e1]), e1);
		auto p2 = std::make_pair(std::get<3>(edges[e2]), e2);
		return (p1 < p2);
	});
	num = rem_edges.size();
	parallel_for(0, num-1, [&](size_t i){
		auto ind = rem_edges[i];
		auto temp = rem_edges[i+1];
		// heights[temp] = num_iters + i + 1; //std::max(heights[temp], heights[i]+1);
		parent[ind] = temp;
	});
	t.next("Dendrogram Stage 2 Time");
	double tt = t.total_time();

	std::cout << "Remaining Edges = " << num << std::endl;
	std::cout << "Num Iters = " << num_iters << std::endl;
	// std::cout << std::endl << "=> Dendrogram Height = " << parlay::reduce_max(heights) << std::endl;

	// return parents;
	return tt;
}

}  // namespace gbbs
