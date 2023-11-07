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
void DendrogramRCtreeTracing(Graph& GA){
	using W = typename Graph::weight_type;
    using K = std::tuple<uintE, uintE>;
    using V = size_t;

    // Preprocess
    auto n = GA.n; auto m = GA.m/2;
    auto deg = sequence<size_t>::from_function(n, [&](size_t u){
        return GA.get_vertex(u).out_degree();
    });
    auto edges = sequence<std::tuple<uintE, uintE>>::uninitialized(2*m);
    auto offsets = parlay::scan(deg).first;
    auto map_f = [&](const uintE& src, const uintE& dst, const W& wgh, const size_t& ind){
        edges[offsets[src] + ind] = {std::min(src,dst), std::max(src,dst)};
    };
    parallel_for(0, n, [&](uintE u){
        GA.get_vertex(u).out_neighbors().map_with_index(map_f);
    });
    parlay::sort_inplace(edges);
    auto index = [&](uintE u, uintE v){
        auto it = parlay::find(edges,std::make_tuple(std::min(u,v), std::max(u,v)));
        return (it - edges.begin())/2;
    };
    auto parent = sequence<size_t>::from_function(m, [&](size_t i){
        return i;
    });

    // Step 1: Compute RC Tree
    auto rctree = sequence<RCtree_node<std::pair<W, size_t>>>(n);
    parallel_for(0, n, [&](size_t i){ 
        rctree[i].parent = i;
        rctree[i].edge_index = m + i;
    });
    auto edge2rcnode = sequence<size_t>::uninitialized(m);

    size_t rem = m, round = 0;
    while (rem){
        auto cur = parlay::filter(parlay::iota(n), [&](size_t i){return deg[i]==1;});
        auto map_f = [&](const uintE& src, const uintE& dst, const W& wgh){
            if (deg[dst] > 1){
                deg[src]--;
                rctree[src].parent = dst;
                rctree[src].round = round;
                auto edge_index = index(src, dst);
                rctree[src].edge_index = edge_index;
                edge2rcnode[edge_index] = src;
                rctree[src].wgh = {wgh,edge_index};
            } else if (deg[dst] == 1){
                // do something
            }
        };
        parallel_for(0, cur.size(),[&](size_t i){
            GA.get_vertex(cur[i]).out_neighbors().map(map_f);
            cur[i] = rctree[cur[i]].parent;
        });
        auto hist = parlay::histogram_by_key(cur);
        parallel_for(0, hist.size(),[&](size_t i){
            deg[hist[i].first] -= hist[i].second;
        });
        round++;
        rem -= cur.size();
    }

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

    // Step 3: Sort by bkt_ids 
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
    // for (size_t i=0; i<m; i++){
    //     std::cout << parent[i] << " ";
    // }
    // std::cout << std::endl;
}

}  // namespace gbbs