// Utilities for creating main functions
#pragma once

// #include "assert.h"
#include "caterpillar.h"
#include "full_binary.h"
#include "gbbs/gbbs.h"
#include "gbbs/graph_io.h"
#include "generate_MSF.h"
#include "path.h"
#include "star.h"
#include "uniform_hook.h"
#include "graph_io.h"

#include <iomanip>

#include "benchmarks/TriangleCounting/ShunTangwongsan15/Triangle.h"

#define run_app(G, APP, mutates, rounds)    \
  double total_time = 0.0;                  \
  for (size_t r = 0; r < rounds; r++) {     \
    if (mutates) {                          \
      auto G_copy = G;                      \
      total_time += APP(G_copy, P);         \
    } else {                                \
      total_time += APP(G, P);              \
    }                                       \
  }                                         \
  auto time_per_iter = total_time / rounds; \
  std::cout << "# time per iter: " << time_per_iter << "\n";

namespace gbbs {


template <class W, class P>
void apply_weights(gbbs::edge_array<W>& edges, P& cmd_line) {
  if (cmd_line.getOption("-perm_weights")) {
    std::cout << "Applying perm weights" << std::endl;
    auto A = parlay::random_permutation(edges.size());
    parlay::parallel_for(0, edges.size(), [&](size_t i) {
      auto [u, v, wgh] = edges.E[i];
      edges.E[i] = {u, v, A[i] + 1};
    });
  } else if (cmd_line.getOption("-low_parallelism")) {
    std::cout << "Applying low_parallelism" << std::endl;
    auto m = edges.size();
    parlay::parallel_for(0, m, [&](size_t i) {
      auto [u, v, wgh] = edges.E[i];
      if (i < m / 2) {
        edges.E[i] = {u, v, i};
      } else {
        edges.E[i] = {u, v, m - 1 - i + m / 2};
      }
    });
  }
}

template <class Graph>
void print_graph(Graph& GA) {
  using W = typename Graph::weight_type;
  size_t n = GA.n;
  for (size_t i=0; i<n; ++i) {
    auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
      std::cout << std::setprecision(15) << u << " " << v << " " << wgh << std::endl;
    };
    GA.get_vertex(i).out_neighbors().map(map_f, false);
  }
}

// Counts the number of triangles in the input graph.
//
// Implementation note: this converts the input graph to a directed graph in
// which we point edges from lower-degree vertices to higher-degree vertices,
// hence the function name.
//
// Arguments:
//   G
//     Graph on which we'll count triangles.
//   f: (uintE, uintE, uintE) -> void
//     Function that's run each triangle. On a triangle with vertices {u, v, w},
//     we run `f(u, v, w)`.
//
// Returns:
//   The number of triangles in `G`.
template <class Graph, class F>
auto Triangle_EmitWeightedEdges(Graph& G, const F& f) {
  using W = typename Graph::weight_type;
  timer gt;
  gt.start();
  uintT n = G.n;
  auto counts = sequence<size_t>::uninitialized(n);
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) { counts[i] = 0; });

  // 1. Rank vertices based on degree
  timer rt;
  rt.start();
  uintE* rank = rankNodes(G, G.n);
  rt.stop();
  rt.next("rank time");

  // 2. Direct edges to point from lower to higher rank vertices.
  // Note that we currently only store out-neighbors for this graph to save
  // memory.
  auto pack_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    return rank[u] < rank[v];
  };
  auto DG = filterGraph(G, pack_predicate);
  // auto DG = Graph::filterGraph(G, pack_predicate);
  gt.stop();
  gt.next("build graph time");

  // 3. Count triangles on the digraph
  timer ct;
  ct.start();

  size_t count = CountDirectedBalanced(DG, counts.begin(), f);
  std::cout << "### Num triangles = " << count << "\n";
  ct.stop();
  ct.next("count time");
  gbbs::free_array(rank, G.n);
  return count;
}

template <class P, class Graph>
void apply_weights_with_graph(gbbs::edge_array<float>& edges, Graph& GA,
                              P& cmd_line) {
  apply_weights(edges, cmd_line);
  if (cmd_line.getOption("-inv_deg_weights")) {
    parlay::parallel_for(0, edges.size(), [&](size_t i) {
      auto [u, v, wgh] = edges.E[i];
      uintE deg_u = GA.get_vertex(u).out_degree();
      uintE deg_v = GA.get_vertex(v).out_degree();
      edges.E[i] = {u, v, 1.0f / (deg_u + deg_v)};
    });
  } else if (cmd_line.getOption("-triangle_weights")) {

    using K = std::pair<uintE, uintE>;
    using V = size_t;
    size_t m = GA.m / 2;
    std::tuple<K, V> empty =
       std::make_tuple(std::make_pair(UINT_E_MAX, UINT_E_MAX), 0);
    auto hash_func = [](K k) -> size_t {
      auto [u, v] = k;
      return parlay::hash64((static_cast<size_t>(u) << 32UL) |
                            static_cast<size_t>(v));
    };
    std::cout << "Creating table." << std::endl;
    auto table = gbbs::make_sparse_table<K, V>(m, empty, hash_func, 1.1);
    std::cout << "Created table, inserting." << std::endl;

    size_t n = GA.n;
    using W = typename Graph::weight_type;
    timer t;
    t.start();
    parlay::parallel_for(
       0, n,
       [&](size_t i) {
         auto map_f = [&](uintE u, uintE v, W wgh) {
           if (u < v) {
             K k = {u, v};
             V v = 0;
             table.insert(std::make_tuple(k, v));
           }
         };
         GA.get_vertex(i).out_neighbors().map(map_f);
       },
       1);
    t.next("insert time");
    std::cout << "Finished building initial table." << std::endl;

    using KV = std::tuple<K, V>;
    auto ins_f = [&](V* v, KV& kv) {
      gbbs::write_add(v, 1);
    };

    auto f = [&](uintE u, uintE v, uintE w) {
      uintE s[3];
      s[0] = u;
      s[1] = v;
      s[2] = w;
      std::sort(s, s + 3);

      K k1 = {s[0], s[1]};
      K k2 = {s[0], s[2]};
      K k3 = {s[1], s[2]};
      V val = 1;
      table.insert_f(std::make_tuple(k1, val), ins_f);
      table.insert_f(std::make_tuple(k2, val), ins_f);
      table.insert_f(std::make_tuple(k3, val), ins_f);
    };
    size_t ct = Triangle_EmitWeightedEdges(GA, f);
    t.next("count time");
    std::cout << "ct = " << ct << std::endl;
    parlay::parallel_for(0, edges.size(), [&](size_t i) {
      auto [u, v, wgh] = edges.E[i];
      float cnt = table.find(std::make_pair(u, v), 1);
      edges.E[i] = {u, v, 1.0f / cnt};
    });
  }
}

template <class CmdLine>
auto run_triangle_weights(CmdLine& P) {
  char* iFile = P.getArgument(0);
  bool mmap = P.getOptionValue("-m");
  bool binary = P.getOptionValue("-b");
  auto GA = gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, mmap, binary);

  using K = std::pair<uintE, uintE>;
  using V = size_t;
  size_t m = GA.m / 2;
  std::tuple<K, V> empty =
     std::make_tuple(std::make_pair(UINT_E_MAX, UINT_E_MAX), 0);
  auto hash_func = [](K k) -> size_t {
    auto [u, v] = k;
    return parlay::hash64((static_cast<size_t>(u) << 32UL) |
                          static_cast<size_t>(v));
  };
  std::cout << "Creating table. m = " << m << std::endl;
  double mult = P.getOptionDoubleValue("-mult", 1.1);
  auto table = gbbs::make_sparse_table<K, V>(mult * m, empty, hash_func);
  std::cout << "Created table, inserting." << std::endl;

  using KV = std::tuple<K, V>;
  auto ins_f = [&](V* v, KV& kv) {
    gbbs::write_add(v, 1);
  };

  size_t ctr = 0;
  auto f = [&](uintE u, uintE v, uintE w) {
    uintE s[3];
    s[0] = u;
    s[1] = v;
    s[2] = w;
    std::sort(s, s + 3);

    K k1 = {s[0], s[1]};
    K k2 = {s[0], s[2]};
    K k3 = {s[1], s[2]};
    V val = 1;
    bool r1 = table.insert_f(std::make_tuple(k1, val), ins_f);
    bool r2 = table.insert_f(std::make_tuple(k2, val), ins_f);
    bool r3 = table.insert_f(std::make_tuple(k3, val), ins_f);

//    size_t failed = !r1;
//    failed += !r2;
//    failed += !r3;
//    if (failed) {
//      C.update_value(failed);
//    }
  };
  size_t ct = Triangle_EmitWeightedEdges(GA, f);
//  std::cout << "Num failed = " << C.get_value() << std::endl;

  auto ents = table.entries();
  std::cout << "Sizeof ents = " << ents.size() << std::endl;

  auto edges = gbbs::edge_array<float>();
  using edge = typename gbbs::edge_array<float>::edge;
  edges.E = parlay::sequence<edge>::from_function(ents.size(), [&](size_t i) {
    auto [k, wgh] = ents[i];
    auto [u, v] = k;
    return std::make_tuple(u, v, 1.0f / (float)wgh);
  });
  auto MSF = gbbs::gbbs_io::edge_list_to_symmetric_graph<float>(edges);
  return MSF;
}

template <class Graph, class CmdLine>
auto run_triangle_weights_2(Graph& GA, CmdLine& P) {
  char* iFile = P.getArgument(0);
  bool mmap = P.getOptionValue("-m");
  bool binary = P.getOptionValue("-b");

  size_t n = GA.n;
  auto counts = sequence<size_t>::uninitialized(GA.n);
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) { counts[i] = 0; });
  std::cout << "n = " << n << " m = " << GA.m << std::endl;

  // 1. Rank vertices based on degree
  timer rt;
  rt.start();
  uintE* rank = rankNodes(GA, GA.n);
  rt.stop();
  rt.next("rank time");

  // 2. Direct edges to point from lower to higher rank vertices.
  // Note that we currently only store out-neighbors for this graph to save
  // memory.
  using W = gbbs::empty;
  auto pack_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    return rank[u] < rank[v];
  };
  auto DG = filterGraph(GA, pack_predicate);
  // auto DG = Graph::filterGraph(G, pack_predicate);

  auto offsets = parlay::sequence<size_t>::uninitialized(GA.n);
  parlay::parallel_for(0, GA.n, [&](size_t i) {
    offsets[i] = std::min(DG.get_vertex(i).out_degree(), (uintE)5);
  });
  size_t m = parlay::scan_inplace(parlay::make_slice(offsets));

  // 3. Count triangles on the digraph
  timer ct;
  ct.start();
  auto edges = gbbs::edge_array<float>();
  using edge = typename gbbs::edge_array<float>::edge;
  std::cout << "Allocating m = " << m << " edges, size = " << sizeof(edge)
            << std::endl;
  edges.E = parlay::sequence<edge>::uninitialized(m);

  std::cout << "Counting!" << std::endl;
  using Wgh = gbbs::empty;
  auto f = [](const uintE& u, const uintE& v, const uintE& w) {
  };
  timer t;
  t.start();
  parlay::parallel_for(
     0, GA.n,
     [&](size_t i) {
       size_t offset_i = offsets[i];
       auto our_neighbors = DG.get_vertex(i).out_neighbors();
       std::vector<std::pair<uintE, uintE>> cts;
       auto map_f = [&](const uintE& u, const uintE& v, const Wgh& wgh,
                        size_t idx) {
         auto their_neighbors = DG.get_vertex(v).out_neighbors();
         size_t ct = our_neighbors.intersect_f_par(&their_neighbors, f);
         cts.push_back({ct, v});
       };
       auto gt = [&](std::pair<uintE, uintE> l, std::pair<uintE, uintE> r) {
         return l.first > r.first;
       };
       DG.get_vertex(i).out_neighbors().map_with_index(map_f, false);
       std::sort(cts.begin(), cts.end(), gt);
       for (size_t j = 0; j < std::min(cts.size(), (size_t)5); ++j) {
         edges.E[offset_i + j] = {i, cts[j].second,
                                  1.0f / (cts[j].first + 1.0f)};
       }
     },
     1);
  t.next("count time");

  std::cout << "Building ret graph" << std::endl;
  auto MSF = gbbs::gbbs_io::edge_list_to_symmetric_graph<float>(edges);
  std::cout << "Returning MSF edges!" << std::endl;
  std::cout << MSF.n << " " << MSF.m << std::endl;
  return MSF;
}

template <class W>
using Edge = gbbs_io::Edge<W>;

template <class P>
gbbs::edge_array<gbbs::intE> add_weights(
   std::vector<Edge<gbbs::empty>> edge_list, P& cmd_line) {
  auto edges = gbbs::edge_array<gbbs::intE>();
  using edge = typename gbbs::edge_array<gbbs::intE>::edge;
  edges.E = parlay::sequence<edge>::uninitialized(edge_list.size());
  parlay::parallel_for(0, edge_list.size(), [&](size_t i) {
    auto [u, v, _] = edge_list[i];
    edges.E[i] = {u, v, 0};
  });
  // Note: n is not set in edges.
  apply_weights(edges, cmd_line);
  return edges;
}

template <class W, class Graph>
gbbs::edge_array<W> to_edges(Graph& GA) {
  auto edges = gbbs::edge_array<W>();
  using edge = typename gbbs::edge_array<W>::edge;
  edges.E = parlay::sequence<edge>::uninitialized(GA.m);
  using Wgh = typename Graph::weight_type;

  size_t n = GA.n;
  auto deg = sequence<uintE>::from_function(
     n, [&](uintE u) { return GA.get_vertex(u).out_degree(); });

  auto offsets = parlay::scan(deg).first;

  auto offset_f = [&](const uintE& src, const uintE& dst, const Wgh& wgh,
                      const uintE& ind) {
    edges.E[offsets[src] + ind] = {src, dst, W()};
  };
  parallel_for(0, n, [&](uintE u) {
    GA.get_vertex(u).out_neighbors().map_with_index(offset_f);
  });
  return edges;
}

}   // namespace gbbs

/* Integer Weighted Graphs */
#define generate_int_weighted_main(APP, mutates)                               \
  int main(int argc, char* argv[]) {                                           \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                         \
    size_t rounds = P.getOptionLongValue("-rounds", 1);                        \
    bool is_path = P.getOption("-path");                                       \
    bool is_caterpillar = P.getOption("-caterpillar");                         \
    bool is_star = P.getOption("-star");                                       \
    bool is_fullb = P.getOption("-fullb");                                     \
    bool is_randomb = P.getOption("-randomb");                                 \
    bool is_randomk = P.getOption("-randomk");                                 \
    bool is_unifhook = P.getOption("-unifhook");                               \
    size_t build_graph = is_path + is_caterpillar + is_star + is_fullb +       \
                         is_randomb + is_randomk + is_unifhook;                \
    bool is_edge_array = P.getOption("-EA");                                   \
    if (build_graph) {                                                         \
      gbbs::edge_array<gbbs::intE> edge_list;                                  \
      if (is_path) {                                                           \
        size_t n = P.getOptionLongValue("-n", 10);                             \
        edge_list = gbbs::generate_path_graph<gbbs::intE>(n);                  \
      } else if (is_star) {                                                    \
        size_t n = P.getOptionLongValue("-n", 10);                             \
        edge_list = gbbs::generate_star_graph<gbbs::intE>(n);                  \
      } else if (is_caterpillar) {                                             \
        size_t k = P.getOptionLongValue("-k", 3);                              \
        edge_list = gbbs::generate_caterpillar_graph<gbbs::intE>(k);           \
      } else if (is_fullb) {                                                   \
        size_t k = P.getOptionLongValue("-k", 3);                              \
        edge_list = gbbs::generate_full_binary_tree<gbbs::intE>(k);            \
      } else if (is_unifhook) {                                                \
        std::cout << "Generating unifhook" << std::endl;                       \
        size_t n = P.getOptionLongValue("-n", 10);                             \
        edge_list = gbbs::generate_uniform_hook<gbbs::intE>(n);                \
      } else if (is_randomb) {                                                 \
      } else if (is_randomk) {                                                 \
      }                                                                        \
      gbbs::apply_weights(edge_list, P);                                       \
      auto G =                                                                 \
         gbbs::gbbs_io::edge_list_to_symmetric_graph<gbbs::intE>(edge_list);   \
      run_app(G, APP, mutates, rounds)                                         \
    } else if (is_edge_array) {                                                \
      char* iFile = P.getArgument(0);                                          \
      auto edge_list{gbbs::gbbs_io::read_unweighted_edge_list(iFile)};         \
      auto weighted_edge_list = gbbs::add_weights(edge_list, P);               \
      auto G = gbbs::gbbs_io::edge_list_to_symmetric_graph<gbbs::intE>(        \
         weighted_edge_list);                                                  \
      run_app(G, APP, mutates, rounds)                                         \
    } else if (P.getOption("-fw")) {                                           \
      char* iFile = P.getArgument(0);                                          \
      bool compressed = P.getOptionValue("-c");                                \
      bool mmap = P.getOptionValue("-m");                                      \
      bool binary = P.getOptionValue("-b");                                    \
      bool is_forest = P.getOption("-f");                                      \
      if (compressed) {                                                        \
        auto G =                                                               \
           gbbs::gbbs_io::read_compressed_symmetric_graph<float>(iFile, mmap); \
        if (is_forest) {                                                       \
          run_app(G, APP, mutates, rounds)                                     \
        } else {                                                               \
          auto MSF = generate_MSF(G, P);                                       \
          run_app(MSF, APP, mutates, rounds)                                   \
        }                                                                      \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_weighted_symmetric_graph<float>(          \
           iFile, mmap, binary);                                               \
        if (is_forest) {                                                       \
          run_app(G, APP, mutates, rounds)                                     \
        } else {                                                               \
          auto MSF = generate_MSF(G, P);                                       \
          run_app(MSF, APP, mutates, rounds)                                   \
        }                                                                      \
      }                                                                        \
    } else if (P.getOption("-iw")) {                                           \
      char* iFile = P.getArgument(0);                                          \
      bool compressed = P.getOptionValue("-c");                                \
      bool mmap = P.getOptionValue("-m");                                      \
      bool binary = P.getOptionValue("-b");                                    \
      bool is_forest = P.getOption("-f");                                      \
      if (compressed) {                                                        \
        auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::intE>(   \
           iFile, mmap);                                                       \
        if (is_forest) {                                                       \
          run_app(G, APP, mutates, rounds)                                     \
        } else {                                                               \
          auto MSF = generate_MSF(G, P);                                       \
          run_app(MSF, APP, mutates, rounds)                                   \
        }                                                                      \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_weighted_symmetric_graph<gbbs::intE>(     \
           iFile, mmap, binary);                                               \
        if (is_forest) {                                                       \
          run_app(G, APP, mutates, rounds)                                     \
        } else {                                                               \
          auto MSF = generate_MSF(G, P);                                       \
          run_app(MSF, APP, mutates, rounds)                                   \
        }                                                                      \
      }                                                                        \
    } else {                                                                   \
      char* iFile = P.getArgument(0);                                          \
      bool compressed = P.getOptionValue("-c");                                \
      bool mmap = P.getOptionValue("-m");                                      \
      if (P.getOption("-triangle_weights") && !compressed) {                   \
        bool binary = P.getOptionValue("-b");                                  \
        auto GA = gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, mmap,  \
                                                                 binary);      \
        auto G = gbbs::run_triangle_weights_2(GA, P);                          \
        if (P.getOption("-write_graph")) {                                     \
          std::cout << "Writing graph" << std::endl;                           \
          gbbs::write_sym_graph_binary(G, P.getOptionValue("-write_graph"));   \
          std::cout << "Finished writing graph" << std::endl;                  \
        }                                                                      \
        auto MSF = generate_MSF(G, P);                                         \
        if (P.getOption("-write_tree")) {                                      \
          std::cout << "Writing graph" << std::endl;                           \
          gbbs::gbbs_io::write_graph_to_file(P.getOptionValue("-write_tree"),  \
                                             MSF);                             \
          std::cout << "Finished writing graph" << std::endl;                  \
        }                                                                      \
        run_app(MSF, APP, mutates, rounds) return 1;                           \
      }                                                                        \
      if (P.getOption("-triangle_weights") && compressed) {                    \
        auto GA = gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::empty>( \
           iFile, mmap);                                                       \
        std::cout << "read graph" << std::endl;                                \
        auto G = gbbs::run_triangle_weights_2(GA, P);                          \
        if (P.getOption("-write_graph")) {                                     \
          std::cout << "Writing graph" << std::endl;                           \
          gbbs::write_sym_graph_binary(G, P.getOptionValue("-write_graph"));   \
          std::cout << "Finished writing graph" << std::endl;                  \
        }                                                                      \
        auto MSF = generate_MSF(G, P);                                         \
        if (P.getOption("-write_tree")) {                                      \
          std::cout << "Writing graph" << std::endl;                           \
          gbbs::gbbs_io::write_graph_to_file(P.getOptionValue("-write_tree"),  \
                                             MSF);                             \
          std::cout << "Finished writing graph" << std::endl;                  \
        }                                                                      \
        run_app(MSF, APP, mutates, rounds) return 1;                           \
      }                                                                        \
      bool binary = P.getOptionValue("-b");                                    \
      bool is_forest = P.getOption("-f");                                      \
      if (compressed) {                                                        \
        auto G_o =                                                             \
           gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::empty>(iFile,  \
                                                                       mmap);  \
        auto edges = gbbs::to_edges<float>(G_o);                               \
        apply_weights_with_graph(edges, G_o, P);                               \
        auto G = gbbs::gbbs_io::edge_list_to_symmetric_graph<float>(edges);    \
        if (is_forest) {                                                       \
          run_app(G, APP, mutates, rounds)                                     \
        } else {                                                               \
          auto MSF = generate_MSF(G, P);                                       \
          run_app(MSF, APP, mutates, rounds)                                   \
        }                                                                      \
      } else {                                                                 \
        auto G_o = gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, mmap, \
                                                                  binary);     \
        auto edges = gbbs::to_edges<float>(G_o);                               \
        apply_weights_with_graph(edges, G_o, P);                               \
        auto G = gbbs::gbbs_io::edge_list_to_symmetric_graph<float>(edges);    \
        if (is_forest) {                                                       \
          run_app(G, APP, mutates, rounds)                                     \
        } else {                                                               \
          auto MSF = generate_MSF(G, P);                                       \
          run_app(MSF, APP, mutates, rounds)                                   \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  }
