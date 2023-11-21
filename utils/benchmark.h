// Utilities for creating main functions
#pragma once

// #include "assert.h"
#include "caterpillar.h"
#include "full_binary.h"
#include "gbbs/graph_io.h"
#include "gbbs/gbbs.h"
#include "generate_MSF.h"
#include "path.h"
#include "star.h"
#include "uniform_hook.h"

namespace gbbs {
template <class P>
void apply_weights(
    gbbs::edge_array<gbbs::intE>& edges, P& cmd_line) {
  if (cmd_line.getOption("-perm_weights")) {
    auto A = parlay::random_permutation(edges.size());
    parlay::parallel_for(0, edges.size(), [&] (size_t i) {
      auto [u, v, wgh] = edges.E[i];
      edges.E[i] = {u, v, A[i] + 1};
    });
  }
}
}  // namespace gbbs

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
#define generate_int_weighted_main(APP, mutates)                             \
  int main(int argc, char* argv[]) {                                         \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                       \
    size_t rounds = P.getOptionLongValue("-rounds", 1);                      \
    bool is_path = P.getOption("-path");                                     \
    bool is_caterpillar = P.getOption("-caterpillar");                       \
    bool is_star = P.getOption("-star");                                     \
    bool is_fullb = P.getOption("-fullb");                                   \
    bool is_randomb = P.getOption("-randomb");                               \
    bool is_randomk = P.getOption("-randomk");                               \
    bool is_unifhook = P.getOption("-unifhook");                             \
    size_t build_graph = is_path + is_caterpillar + is_star + is_fullb +     \
                         is_randomb + is_randomk + is_unifhook;              \
    if (build_graph) {                                                       \
      gbbs::edge_array<gbbs::intE> edge_list;                  \
      if (is_path) {                                                         \
        size_t n = P.getOptionLongValue("-n", 10);                           \
        edge_list = gbbs::generate_path_graph<gbbs::intE>(n);             \
      } else if (is_star) {                                                  \
        size_t n = P.getOptionLongValue("-n", 10);                           \
        edge_list = gbbs::generate_star_graph<gbbs::intE>(n);                \
      } else if (is_caterpillar) {                                           \
        size_t k = P.getOptionLongValue("-k", 3);                            \
        edge_list = gbbs::generate_caterpillar_graph<gbbs::intE>(k);      \
      } else if (is_fullb) {                                                 \
        size_t k = P.getOptionLongValue("-k", 3);                            \
        edge_list = gbbs::generate_full_binary_tree<gbbs::intE>(k);          \
      } else if (is_unifhook) {                                              \
        std::cout << "Generating unifhook" << std::endl;                     \
        size_t n = P.getOptionLongValue("-n", 10);                           \
        edge_list = gbbs::generate_uniform_hook<gbbs::intE>(n);              \
      } else if (is_randomb) {                                               \
      } else if (is_randomk) {                                               \
      }                                                                      \
      gbbs::apply_weights(edge_list, P);                                               \
      auto G = gbbs::gbbs_io::edge_list_to_symmetric_graph<gbbs::intE>(edge_list);          \
      run_app(G, APP, mutates, rounds)                                       \
    } else {                                                                 \
      char* iFile = P.getArgument(0);                                        \
      bool compressed = P.getOptionValue("-c");                              \
      bool mmap = P.getOptionValue("-m");                                    \
      bool binary = P.getOptionValue("-b");                                  \
      bool is_forest = P.getOption("-f");                                    \
      if (compressed) {                                                      \
        auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::intE>( \
            iFile, mmap);                                                    \
        if (is_forest) {                                                     \
          run_app(G, APP, mutates, rounds)                                   \
        } else {                                                             \
          auto MSF = generate_MSF(G);                                        \
          run_app(MSF, APP, mutates, rounds)                                 \
        }                                                                    \
      } else {                                                               \
        auto G = gbbs::gbbs_io::read_weighted_symmetric_graph<gbbs::intE>(   \
            iFile, mmap, binary);                                            \
        if (is_forest) {                                                     \
          run_app(G, APP, mutates, rounds)                                   \
        } else {                                                             \
          auto MSF = generate_MSF(G);                                        \
          run_app(MSF, APP, mutates, rounds)                                 \
        }                                                                    \
      }                                                                      \
    }                                                                        \
  }
