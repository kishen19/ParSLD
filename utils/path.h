
#include "gbbs/gbbs.h"

namespace gbbs {

template <typename W>
auto generate_path_graph(size_t n) {
  using edge = std::tuple<uintE, uintE, W>;

  auto m = n-1;
  auto edges = sequence<edge>::uninitialized(m);
  parallel_for(0, m, [&](uintE i) {
    edges[i] = {i, i+1, 1};
  });
  return edge_array<W>(std::move(edges), n);
}

} // namespace gbbs