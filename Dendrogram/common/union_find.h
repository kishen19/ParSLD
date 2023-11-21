#pragma once

namespace gbbs {

struct union_find {
  size_t n;
  sequence<uintE> parents;
  union_find(size_t n) : n(n) {
    parents =
        sequence<uintE>::from_function(n, [&](uintE i) { return i; });
  }

  std::pair<uintE, uintE> find_compress(uintE i) {
    uintE j = i;
    uintE total_work = 0;
    if (parents[j] == j) return {j, 1};
    do {
      j = parents[j];
      total_work++;
    } while (parents[j] != j);
    uintE tmp;
    while ((tmp = parents[i]) > j) {
      parents[i] = j;
      i = tmp;
      total_work++;
    }
    return {j, total_work};
  }

  // Returns the parent
  std::pair<uintE, uintE> unite(uintE u_orig, uintE v_orig) {
    uintE u = u_orig;
    uintE v = v_orig;
    while (u != v) {
      uintE tw1, tw2;
      std::tie(u, tw1) = find_compress(u);
      std::tie(v, tw2) = find_compress(v);
      if (u > v && parents[u] == u &&
              gbbs::atomic_compare_and_swap(&parents[u], u, v)) {
        return {v, tw1 + tw2};
      } else if (v > u && parents[v] == v &&
                gbbs::atomic_compare_and_swap(&parents[v], v, u)) {
        return {u, tw1 + tw2};
      }
    }
    std::cerr << "Uniting two already merged clusters: " << u_orig 
              << ", " << v_orig << std::endl;
    assert(false);
    exit(-1);
  }
};

} // namespace gbbs