// Utilities for creating main functions
#pragma once

// #include "assert.h"
#include "generate_MSF.h"
#include "full_binary.h"
#include "gbbs/graph_io.h"
#include "path.h"
#include "star.h"
#include "uniform_hook.h"

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
    bool is_star = P.getOption("-star");                                     \
    bool is_fullb = P.getOption("-fullb");                                   \
    bool is_randomb = P.getOption("-randomb");                               \
    bool is_randomk = P.getOption("-randomk");                               \
    bool is_unifhook = P.getOption("-unifhook");                             \
    if (is_path) {                                                           \
      size_t n = P.getOptionLongValue("-n", 10);                             \
      auto G = gbbs::generate_path_graph<gbbs::intE>(n);                     \
      run_app(G, APP, mutates, rounds)                                       \
    } else if (is_star) {                                                    \
      size_t n = P.getOptionLongValue("-n", 10);                             \
      auto G = gbbs::generate_star_graph<gbbs::intE>(n);                     \
      run_app(G, APP, mutates, rounds)                                       \
    } else if (is_fullb) {                                                   \
      size_t k = P.getOptionLongValue("-k", 3);                              \
      auto G = gbbs::generate_full_binary_tree<gbbs::intE>(k);               \
      run_app(G, APP, mutates, rounds)                                       \
    } else if (is_unifhook) {                                                \
      std::cout << "Generating unifhook" << std::endl;                       \
      size_t n = P.getOptionLongValue("-n", 10);                             \
      auto G = gbbs::generate_uniform_hook<gbbs::intE>(n);                   \
      run_app(G, APP, mutates, rounds)                                       \
    } else if (is_randomb) {                                                 \
    } else if (is_randomk) {                                                 \
    } else {                                                                 \
      char* iFile = P.getArgument(0);                                        \
      bool compressed = P.getOptionValue("-c");                              \
      bool mmap = P.getOptionValue("-m");                                    \
      bool binary = P.getOptionValue("-b");                                  \
      bool is_forest = P.getOption("-f");                                    \
      if (compressed) {                                                      \
        auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::intE>( \
           iFile, mmap);                                                     \
        if (is_forest){                                                      \
          run_app(G, APP, mutates, rounds)                                   \
        } else {                                                             \
          auto MSF = generate_MSF(G);                                        \
          run_app(MSF, APP, mutates, rounds)                                 \
        }                                                                    \
      } else {                                                               \
        auto G = gbbs::gbbs_io::read_weighted_symmetric_graph<gbbs::intE>(   \
           iFile, mmap, binary);                                             \
        if (is_forest){                                                      \
          run_app(G, APP, mutates, rounds)                                   \
        } else {                                                             \
          auto MSF = generate_MSF(G);                                        \
          run_app(MSF, APP, mutates, rounds)                                 \
        }                                                                    \
      }                                                                      \
    }                                                                        \
  }
