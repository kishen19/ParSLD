
#include "gbbs/gbbs.h"

namespace gbbs {

template <typename W>
auto generate_star_graph(size_t n) {
  using edge = std::tuple<uintE, uintE, W>;

  auto m = n-1;
  auto edges = sequence<edge>::uninitialized(m);
  parallel_for(1, n, [&](uintE i){
    edges[i-1] = {0, i, 1};
  });
  return edge_array<W>(std::move(edges), n);
}

} // namespace gbbs