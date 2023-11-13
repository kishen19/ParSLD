#pragma once

#include "gbbs/gbbs.h"
#include "benchmarks/Connectivity/SimpleUnionAsync/Connectivity.h"
#include "benchmarks/MinimumSpanningForest/Boruvka/MinimumSpanningForest.h"

namespace gbbs {

template <class Graph>
double DendrogramSeqUF(Graph& GA){
	using W = typename Graph::weight_type;

	// Init
	auto n = GA.n; auto m = GA.m/2;
    auto deg = sequence<size_t>::from_function(n, [&](size_t u){
        return GA.get_vertex(u).out_degree();
    });
    auto edges = sequence<std::tuple<uintE, uintE, W>>::uninitialized(2*m);
    auto offsets = parlay::scan(deg).first;
    auto map_f = [&](const uintE& src, const uintE& dst, const W& wgh, const size_t& ind){
        edges[offsets[src] + ind] = {std::min(src,dst), std::max(src,dst), wgh};
    };
    parallel_for(0, n, [&](uintE u){
        GA.get_vertex(u).out_neighbors().map_with_index(map_f);
    });
    parlay::sort_inplace(edges);
    // auto index = [&](uintE u, uintE v){
    //     auto it = parlay::find(edges,std::make_tuple(std::min(u,v), std::max(u,v)));
    //     return (it - edges.begin())/2;
    // };
    auto parent = sequence<size_t>::from_function(m, [&](size_t i){
        return i;
    });

	// Step 2: Sorting the edges by (weight, index)
	timer t;
	t.start();
	
	auto proc_edges = sequence<std::tuple<W,size_t>>::from_function(m, [&](size_t i){
		return std::make_tuple(std::get<2>(edges[2*i]), i);
	});
	parlay::sort_inplace(proc_edges);
	// t.next("Init Time");

	// Step 3: Applying Union Find to the sorted sequence of edges
	auto uf = simple_union_find::SimpleUnionAsyncStruct(n);
	auto aux = sequence<size_t>(n,n); // extra info required for assigning parents
	// auto heights = sequence<uintE>(m,0); // Heights of every node in the dendrogram

	for(size_t ind = 0; ind < m; ind++){
		auto wgh = std::get<0>(proc_edges[ind]);
		auto i = std::get<1>(proc_edges[ind]);
		auto u = simple_union_find::find_compress(std::get<0>(edges[i]), uf.parents.begin());
		auto v = simple_union_find::find_compress(std::get<1>(edges[i]), uf.parents.begin());
		// uintE height = 0;
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
		uf.unite(u,v);
		aux[simple_union_find::find_compress(u, uf.parents.begin())] = i;
		// heights[i] = height;
	};
	// t.next("Dendrogram Time");
	double tt = t.total_time();

	// for (size_t i=0; i<m; i++){
    //     std::cout << parent[i] << " ";
    // }
    // std::cout << std::endl;

	// std::cout << std::endl << "=> Dendrogram Height = " << parlay::reduce_max(heights) << std::endl;

	// return parents;
	return tt;
}

template <class Graph>
double DendrogramSeqUFMST(Graph& GA){
	size_t n = GA.n;

	// Step 1: Get MST (MSF) of the given graph
	auto mst_edges = MinimumSpanningForest_boruvka::MinimumSpanningForest(GA);
	size_t m = mst_edges.size();

	// Step 2: Sorting the edges by (weight, index)
	timer t;
	t.start();
	
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
		auto u = simple_union_find::find_compress(std::get<0>(mst_edges[i]), uf.parents.begin());
		auto v = simple_union_find::find_compress(std::get<1>(mst_edges[i]), uf.parents.begin());
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
		aux[simple_union_find::find_compress(u, uf.parents.begin())] = i;
		heights[i] = height;
	};
	t.next("Dendrogram Time");
	double tt = t.total_time();

	std::cout << std::endl << "=> Dendrogram Height = " << parlay::reduce_max(heights) << std::endl;

	// return parents;
	return tt;
}

}  // namespace gbbs
