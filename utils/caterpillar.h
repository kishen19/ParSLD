
#include "gbbs/gbbs.h"

namespace gbbs {

// Caterpillar on 2^k vertices
template <typename W, class P>
auto generate_caterpillar_graph(size_t k, P& cmd_line) {
  using edge = std::tuple<uintE, uintE, W>;
  bool permute = cmd_line.getOption("-permute");

  parlay::sequence<size_t> p;
  size_t n = 1 << k;
  if (permute) {
    p = parlay::random_permutation(n);
  } else {
    p = parlay::sequence<size_t>::uninitialized(n);
    parlay::copy(parlay::iota(n), p);
  }

  auto m = n-1;
  auto edges = sequence<edge>::uninitialized(m);
  auto midp = m / 2;
  parallel_for(0, midp, [&](uintE i) {
    edges[i] = {p[i], p[i+1], 1};
  });
  parallel_for(midp, m, [&](uintE i) {
    edges[i] = {p[i-midp], p[i+1], 1};
  });
  auto edge_list = edge_array<W>(std::move(edges), n);
  return gbbs_io::edge_list_to_symmetric_graph<W>(edge_list);
}

} // namespace gbbs
