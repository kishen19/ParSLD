#pragma once

namespace gbbs {

struct union_find {
  size_t n;
  sequence<uintE> parents;
  union_find(size_t n) : n(n) {
    parents =
        sequence<uintE>::from_function(n, [&](uintE i) { return i; });
  }

  uintE find_compress(uintE i) {
    uintE j = i;
    if (parents[j] == j) return j;
    do {
      j = parents[j];
    } while (parents[j] != j);
    uintE tmp;
    while ((tmp = parents[i]) > j) {
      parents[i] = j;
      i = tmp;
    }
    return j;
  }

  // Returns the parent
  uintE unite(uintE u_orig, uintE v_orig) {
    uintE u = u_orig;
    uintE v = v_orig;
    while (u != v) {
      u = find_compress(u);
      v = find_compress(v);
      if (u > v && parents[u] == u &&
              gbbs::atomic_compare_and_swap(&parents[u], u, v)) {
        return v;
      } else if (v > u && parents[v] == v &&
                gbbs::atomic_compare_and_swap(&parents[v], v, u)) {
        return u;
      }
    }
    std::cerr << "Uniting two already merged clusters: " << u_orig 
              << ", " << v_orig << std::endl;
    assert(false);
    exit(-1);
  }
};

} // namespace gbbs