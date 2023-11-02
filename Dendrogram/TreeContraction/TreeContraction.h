#pragma once
#include "gbbs/gbbs.h"
#include "benchmarks/Connectivity/SimpleUnionAsync/Connectivity.h"

namespace gbbs{
namespace tree_contraction{

template <class Graph, class heap>
struct tree_contraction{
    using W = typename Graph::weight_type;
    Graph G;
    size_t n,m;
    heap* heaps;
    // simple_union_find::SimpleUnionAsyncStruct* uf;
    sequence<size_t> parent;
    sequence<size_t> deg, deg_old;

    sequence<std::tuple<uintE, uintE>> edges;
    // parlay::random_generator gen;
    // sequence<bool> tosses;

    tree_contraction(Graph& G_): G(G_){
        n = G.n; m = G.m/2;
        heaps = gbbs::new_array<heap>(n);
        // tosses = sequence<bool>(n,0);
        // uf = new simple_union_find::SimpleUnionAsyncStruct(n);
        deg = sequence<size_t>::from_function(n, [&](size_t u){
            return G.get_vertex(u).out_degree();
        });
        deg_old = sequence<size_t>::uninitialized(n);
        parent = sequence<size_t>::from_function(m,[&](size_t i){return i; });

        // Compute indices for edges
        edges = sequence<std::tuple<uintE, uintE>>::uninitialized(2*m);
        auto offsets = parlay::scan(deg).first;
        auto map_f = [&](const uintE& src, const uintE& dst, const W& wgh, const size_t& ind){
            edges[offsets[src] + ind] = {std::min(src,dst), std::max(src,dst)};
        };
        parallel_for(0, n, [&](uintE u){
            G.get_vertex(u).out_neighbors().map_with_index(map_f);
        });
        parlay::sort_inplace(edges);
    }

    size_t get_index(uintE u, uintE v){
        auto it = parlay::find(edges,std::make_tuple(std::min(u,v), std::max(u,v)));
        return (it - edges.begin())/2;
    }

    template <class Seq>
    heap* reduce_heaps(const Seq& seq){
        auto seq_size = seq.size();
        if (seq_size <= 1000){
            auto u = std::get<0>(seq[0]);
            for (size_t i=1; i<seq_size; ++i){
                heaps[u].merge(&heaps[std::get<0>(seq[i])]);
            }
            return &heaps[u];
        }
        else {
            heap* h1;
            heap* h2;
            auto f1 = [&]() { h1 = reduce_heaps(seq.cut(0,seq_size / 2)); };
            auto f2 = [&]() { h2 = reduce_heaps(seq.cut(seq_size/2, seq_size)); };
            par_do(f1, f2);
            h1->merge(h2);
            return h1;
        }
    }

    void run(){
        std::uniform_int_distribution<> dis(0, 1);
        // for(uintE i=0; i<n; i++){
        //     std::cout << deg[i] << " ";
        // }
        // std::cout << std::endl;
        size_t rem = n;
        while (rem){
            parlay::copy(deg, deg_old);
            // Step 1: Rake
            parallel_for(0, n, [&](uintE i){
                if (deg_old[i]>1){
                    auto filter_f = [&](const uintE& src, const uintE& nbhr, const W& wgh){
                        return deg_old[nbhr] == 1;
                    };
                    auto num = G.get_vertex(i).out_neighbors().count(filter_f);
                    if (num>0){
                        auto tmp = sequence<std::tuple<uintE, W>>::uninitialized(num);
                        auto out_f = [&](size_t j, const std::tuple<uintE, W>& nw){
                            // Get the edge and index
                            auto v = std::get<0>(nw);
                            auto wgh = std::get<1>(nw);
                            auto ind = get_index(i,v);
                            // Filter out the protected edges and assign parent
                            auto filtered = sequence<std::tuple<W, size_t>>::uninitialized(heaps[v].size);
                            auto k = heaps[v].filter({wgh,ind}, filtered);
                            if (k < 1000){
                                for(size_t s = 0; s<k; s++){
                                    parent[std::get<1>(filtered[s])] = std::get<1>(filtered[s+1]);
                                }
                            } else{
                                parallel_for(0, k-1, [&](size_t s){
                                    parent[std::get<1>(filtered[s])] = std::get<1>(filtered[s+1]);
                                });
                            }
                            if (k>0){
                                parent[std::get<1>(filtered[k-1])] = ind;
                            }
                            // Update heaps[v], but don't merge. We will merge later in parallel
                            // Update deg[v]
                            heaps[v].insert({wgh,ind});
                            deg[v] -=1;
                            tmp[j] = {v, wgh};
                        };
                        G.get_vertex(i).out_neighbors().filter(filter_f, out_f, tmp.begin());
                        auto out = reduce_heaps(make_slice(tmp));
                        heaps[i].merge(out);
                        deg[i] -= num;
                    }
                }
            });
            rem = n-parlay::count(deg,0);
            std::cout << "n: " << rem << std::endl;
            // Step 2: Compress
            // auto ind = parlay::iota(n);
            // auto deg2 = parlay::filter(ind, [&](size_t u){
            // 	return deg[u] == 2;
            // });
            // parallel_for(0, deg2.size(),[&](size_t i){
            // 	auto r = gen[i];
            // 	tosses[deg2[i]] =  dis(r);
            // });
            // gen += n;
            // parallel_for(0, deg2.size(), [&](size_t i){
            // 	auto u = deg2[i];
            // });
        }

        // Update parent for edges in heaps
        parallel_for(0, n, [&](size_t i){
            if (heaps[i].size > 0){
                auto kv = heaps[i].get_min();
                heaps[i].delete_min();
                while(heaps[i].size){
                    auto kv1 = heaps[i].get_min();
                    heaps[i].delete_min();
                    parent[std::get<1>(kv)] = std::get<1>(kv1);
                    kv = kv1;
                }
            }
        });
        // for (auto i=0; i<m; i++){
        //     std::cout << parent[i] << " ";
        // }
        // std::cout << std::endl;
    }    
};

} // namespace tree_contraction
} // namespace gbbs