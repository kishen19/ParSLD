#pragma once

#include <cassert>
#include <tuple>

#include "gbbs/bridge.h"

namespace gbbs {

template <class K, class V>
class concurrent_table {
 public:
  using T = std::tuple<K, V>;

  size_t m;
  size_t mask;
  T empty;
  K empty_key;

  sequence<T> backing;
  slice<T> table;

  parlay::hash<K> hash;

  size_t size() const { return m; }

  inline size_t hashToRange(size_t h) const { return h & mask; }
  inline size_t firstIndex(K& k) const { return hashToRange(hash(k)); }
  inline size_t incrementIndex(size_t h) const { return hashToRange(h + 1); }

  concurrent_table()
      : m(0), mask(0), table(make_slice((T*)nullptr, (T*)nullptr)) {}

  // Size is the maximum number of values the hash table will hold.
  // Overfilling the table could put it into an infinite loop.
  concurrent_table(size_t _m, T _empty, long inp_space_mult = -1)
      : empty(_empty),
        empty_key(std::get<0>(empty)),
        table(make_slice((T*)nullptr, (T*)nullptr)) {
    double space_mult = 1.1;
    if (inp_space_mult != -1)
      space_mult = inp_space_mult;
    m = (size_t)1 << parlay::log2_up((size_t)(space_mult * _m) + 1);
    mask = m - 1;
    backing = sequence<T>::uninitialized(m);
    table = make_slice(backing);
    std::cout << "Table size = " << m << std::endl;
    clear_table();
  }

  // Size is the maximum number of values the hash table will hold.
  // Overfilling the table could put it into an infinite loop.
  concurrent_table(size_t _m, T _empty, T* _tab, bool clear = true)
      : m(_m),
        mask(m - 1),
        empty(_empty),
        empty_key(std::get<0>(empty)),
        backing(sequence<T>()),
        table(make_slice(_tab, _tab + m)) {
    if (clear) {
      clear_table();
    }
  }

  bool insert(std::tuple<K, V> kv) {
    K k = std::get<0>(kv);
    size_t h = firstIndex(k);
    while (true) {
      if (std::get<0>(table[h]) == empty_key) {
        if (gbbs::atomic_compare_and_swap(&std::get<0>(table[h]), empty_key,
                                          k)) {
          std::get<1>(table[h]) = std::get<1>(kv);
          return true;
        }
      }
      if (std::get<0>(table[h]) == k) {
        return false;
      }
      h = incrementIndex(h);
    }
    return false;
  }

  bool contains(K k) const {
    size_t h = firstIndex(k);
    while (true) {
      if (std::get<0>(table[h]) == k) {
        return true;
      } else if (std::get<0>(table[h]) == empty_key) {
        return false;
      }
      h = incrementIndex(h);
    }
    return false;
  }

  V find(K k) const {
    size_t h = firstIndex(k);
    while (true) {
      if (std::get<0>(table[h]) == k) {
        return std::get<1>(table[h]);
      } else if (std::get<0>(table[h]) == empty_key) {
        std::cerr << "Unexpected empty." << std::endl;
        exit(-1);
      }
      h = incrementIndex(h);
    }
  }

  sequence<T> entries() const {
    auto pred = [&](const T& t) {
      return std::get<0>(t) != empty_key;
    };
    return parlay::filter(table, pred);
  }

  void clear_table() {
    parallel_for(0, m, [&](size_t i) { table[i] = empty; });
  }
};

template <class K, class V>
inline concurrent_table<K, V> make_concurrent_table(size_t m,
                                                    std::tuple<K, V> empty,
                                                    long space_mult = -1) {
  return concurrent_table<K, V>(m, empty,space_mult);
}

template <class K, class V>
inline concurrent_table<K, V> make_concurrent_table(std::tuple<K, V>* tab,
                                                    size_t m,
                                                    std::tuple<K, V> empty,
                                                    bool clear = true) {
  return concurrent_table<K, V>(m, empty, tab, clear);
}

}   // namespace gbbs
