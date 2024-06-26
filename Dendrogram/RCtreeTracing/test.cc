// Usage:
// ./Dendrogram -s -rounds 1 ../../inputs/rMatGraph_WJ_5_100
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm

#include "utils/benchmark.h"

namespace gbbs {
namespace {
}   // namespace
}   // namespace gbbs


void run_t1() {
  size_t n = 100000000;
  using uintE = gbbs::uintE;

  using K=std::pair<uintE, uintE>;
  using V = size_t;

  std::tuple<K, V> empty = std::make_tuple(std::make_pair(UINT_E_MAX, UINT_E_MAX), 0);
  std::cout << "Size: " << (sizeof(std::tuple<K, V>)) << std::endl;
  auto hash_func = [] (K k) -> size_t {
    size_t u, v;
    std::tie(u,v) = k;
    return parlay::hash64((static_cast<size_t>(u) << 32UL) | static_cast<size_t>(v));
  };
  gbbs::timer t; t.start();
  auto table = gbbs::make_sparse_table<K, V>(n, empty, hash_func);
  t.next("Init time");
  parlay::parallel_for(0, n/4, [&] (size_t i) {
    K k = {i, i+1};
    V v = 0;
    table.insert(std::make_tuple(k,v));
  });
  t.next("Insert time");
}

//void run_t2() {
//  size_t n = 100000000;
//  using uintE = gbbs::uintE;
//
//  using K=size_t;
//  using V = uintE;
//
//  std::tuple<K, V> empty = std::make_tuple(UINT_E_MAX, 0);
//  auto hash_func = [] (K k) -> size_t {
//    return parlay::hash64(k);
//  };
//  gbbs::timer t; t.start();
//  auto table = gbbs::make_sparse_table<K, V>(n, empty, hash_func);
//  t.next("Init time");
//  parlay::parallel_for(0, n/4, [&] (size_t i) {
//    K k = (i << 32UL) | (i+1); 
//    V v = 0;
//    table.insert(std::make_tuple(k,v));
//  });
//  t.next("Insert time");
//}


int main() {

  std::cout << "here!" << std::endl;

  run_t1();

//  run_t2();

  return 1;
}


