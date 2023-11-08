#pragma once

#include "gbbs/gbbs.h"
// #include "gbbs/helpers/sequential_ht.h"
#include "benchmarks/Connectivity/SimpleUnionAsync/Connectivity.h"

namespace gbbs {

template <typename weight_type>
struct RCtree_node{
    size_t parent;
    size_t round;
    size_t edge_index;
    weight_type wgh;
    bool check = false;
    uintE alt;
};

template <class Graph>
double DendrogramRCtreeTracing(Graph& GA){
	using W = typename Graph::weight_type;
    using K = std::pair<uintE, uintE>;
    using V = size_t;

    timer t;
	t.start();
    // Preprocess
    auto n = GA.n; auto m = GA.m/2;
    auto ht_size = 1 << (parlay::log2_up((size_t)(m + 1)) + 7);
    auto empty = std::make_tuple(std::make_pair(UINT_E_MAX, UINT_E_MAX), m);

    auto deg = sequence<size_t>::from_function(n, [&](size_t u){ 
        return GA.get_vertex(u).out_degree(); });
    auto table = sequence<std::tuple<K,V>>(ht_size, empty);
    sequentialHT<K, V> S(table.begin(), ht_size, empty);
    
    auto edges = sequence<std::pair<uintE, uintE>>::uninitialized(2*m);
    auto offsets = parlay::scan(deg).first;
    auto offset_f = [&](const uintE& src, const uintE& dst, const W& wgh, const size_t& ind){
        edges[offsets[src] + ind] = {std::min(src,dst), std::max(src,dst)};
    };
    parallel_for(0, n, [&](uintE u){
        GA.get_vertex(u).out_neighbors().map_with_index(offset_f);
    });
    parlay::sort_inplace(edges);
    auto ht_f = [&](auto cur, auto v){ return std::get<1>(v); };
    parallel_for(0, m, [&](size_t i){
        auto key_value = std::make_tuple(edges[2*i], i);
        S.insertF(key_value, ht_f);
    });
    t.next("Preprocess Time");
    
    // Step 1: Compute RC Tree
    auto rctree = sequence<RCtree_node<std::pair<W, size_t>>>(n);
    auto edge2rcnode = sequence<size_t>::uninitialized(m);
    parallel_for(0, n, [&](size_t i){ 
        rctree[i].parent = i;
        rctree[i].edge_index = m + i;
    });
    size_t rem = m, round = 0;
    auto map_f = [&](const uintE& src, const uintE& dst, const W& wgh){
        if (deg[dst] > 1){
            deg[src]--;
            rctree[src].parent = dst;
            rctree[src].round = round;
            auto key = std::make_pair(std::min(src, dst), std::max(src, dst));
            auto edge_index = std::get<1>(S.find(key));
            if (edge_index == m){
                std::cout << src << " " << dst << " " << edge_index << std::endl;
            }
            rctree[src].edge_index = edge_index;
            edge2rcnode[edge_index] = src;
            rctree[src].wgh = {wgh,edge_index};
        } else if (deg[dst] == 1){
            // do something
        }
    };
    
    while(rem){
        auto cur = parlay::filter(parlay::iota(n), [&](size_t i){return deg[i]==1;});
        parallel_for(0, cur.size(),[&](size_t i){
            GA.get_vertex(cur[i]).out_neighbors().map(map_f);
            cur[i] = rctree[cur[i]].parent;
        });
        // std::cout << "check22" << std::endl;
        auto hist = parlay::histogram_by_key(cur);
        // std::cout << "check221 " << hist.size() <<  std::endl;
        parallel_for(0, hist.size(),[&](size_t i){
            deg[hist[i].first] -= hist[i].second;
        });
        // std::cout << "check23" << std::endl;
        round++;
        rem -= cur.size();
    }
    t.next("RC Tree Time");

    // Step 2: Compute bucket_id for each edge
    auto bkt_ids = sequence<std::pair<size_t, size_t>>::uninitialized(m);
    parallel_for(0, m, [&](size_t i){
        auto cur = edge2rcnode[i];
        auto wgh = rctree[cur].wgh;
        cur = rctree[cur].parent;
        while (rctree[cur].parent != cur && rctree[cur].wgh < wgh){
            cur = rctree[cur].parent;
        }
        bkt_ids[i] = {rctree[cur].edge_index, i};
    });
    t.next("Bucket Computation Time");

    // Step 3: Sort by bkt_ids 
    auto parent = sequence<size_t>::from_function(m, [&](size_t i){ return i;});
    sort_inplace(bkt_ids);
    parallel_for(0, n-1, [&](size_t i){
        if (bkt_ids[i].first == bkt_ids[i+1].first){
            parent[bkt_ids[i].second] = parent[bkt_ids[i+1].second];
        } else{
            if (bkt_ids[i].first < m){
                parent[bkt_ids[i].second] = bkt_ids[i].first;
            }
        }
    });
    t.next("Bucket Sorting and Finish Time");
    double tt = t.total_time();
    // for (size_t i=0; i<m; i++){
    //     std::cout << parent[i] << " ";
    // }
    // std::cout << std::endl;
    return tt;
}

}  // namespace gbbs