
#include "gbbs/gbbs.h"

namespace gbbs {

// Caterpillar on 2^k vertices
template <typename W>
auto generate_caterpillar_graph(size_t k) {
  using edge = std::tuple<uintE, uintE, W>;
  size_t n = 1 << k;
  auto m = n-1;
  auto edges = sequence<edge>::uninitialized(m);
  auto midp = m / 2;
  parallel_for(0, midp, [&](uintE i) {
    edges[i] = {i, i+1, 1};
  });
  parallel_for(midp, m, [&](uintE i) {
    edges[i] = {i-midp, i+1, 1};
  });
  return edge_array<W>(std::move(edges), n);
}

} // namespace gbbs
