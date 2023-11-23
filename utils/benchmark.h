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

template <class P, class Graph>
void apply_weights_with_graph(gbbs::edge_array<float>& edges, Graph& GA,
                              P& cmd_line) {
  apply_weights(edges, cmd_line);
  if (cmd_line.getOption("-inv_deg_weights")) {
     parlay::parallel_for(0, edges.size(), [&](size_t i) {
      auto [u, v, wgh] = edges.E[i];
      uintE deg_u = GA.get_vertex(u).out_degree();
      uintE deg_v = GA.get_vertex(v).out_degree();
      edges.E[i] = {u, v, 1.0f/(deg_u + deg_v)};
    });
  } else if (cmd_line.getOption("-triangle_weights")) {
     // TODO: compute jaccard instead of doing the intersection.
     parlay::parallel_for(0, edges.size(), [&](size_t i) {
      auto [u, v, wgh] = edges.E[i];
      auto v_neighbors = GA.get_vertex(v).out_neighbors();
      size_t ct;
      using vertex = typename Graph::vertex;
      if constexpr (std::is_same<vertex, gbbs::symmetric_vertex<gbbs::empty>>()) {
        auto f = [&] () {};
        ct = std::max(GA.get_vertex(u).out_neighbors().intersect(&v_neighbors), (size_t)1);
      } else {
        ct = 1;
      }
      edges.E[i] = {u, v, 1.0f/ct};
    });
  }
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
    } else if (P.getOption("-w")) {                                            \
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
          auto MSF = generate_MSF(G);                                          \
          run_app(MSF, APP, mutates, rounds)                                   \
        }                                                                      \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_weighted_symmetric_graph<gbbs::intE>(     \
           iFile, mmap, binary);                                               \
        if (is_forest) {                                                       \
          run_app(G, APP, mutates, rounds)                                     \
        } else {                                                               \
          auto MSF = generate_MSF(G);                                          \
          run_app(MSF, APP, mutates, rounds)                                   \
        }                                                                      \
      }                                                                        \
    } else {                                                                   \
      char* iFile = P.getArgument(0);                                          \
      bool compressed = P.getOptionValue("-c");                                \
      bool mmap = P.getOptionValue("-m");                                      \
      bool binary = P.getOptionValue("-b");                                    \
      bool is_forest = P.getOption("-f");                                      \
      if (compressed) {                                                        \
        auto G_o =                                                             \
           gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::empty>(iFile,  \
                                                                       mmap);  \
        auto edges = gbbs::to_edges<float>(G_o);                                            \
        apply_weights_with_graph(edges, G_o, P);                               \
        auto G =                                                               \
           gbbs::gbbs_io::edge_list_to_symmetric_graph<float>(edges);     \
        if (is_forest) {                                                       \
          run_app(G, APP, mutates, rounds)                                     \
        } else {                                                               \
          auto MSF = generate_MSF(G);                                          \
          run_app(MSF, APP, mutates, rounds)                                   \
        }                                                                      \
      } else {                                                                 \
        auto G_o = gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, mmap, \
                                                                  binary);     \
        auto edges = gbbs::to_edges<float>(G_o);                                            \
        apply_weights_with_graph(edges, G_o, P);                               \
        auto G =                                                               \
           gbbs::gbbs_io::edge_list_to_symmetric_graph<float>(edges);     \
        if (is_forest) {                                                       \
          run_app(G, APP, mutates, rounds)                                     \
        } else {                                                               \
          auto MSF = generate_MSF(G);                                          \
          run_app(MSF, APP, mutates, rounds)                                   \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  }
