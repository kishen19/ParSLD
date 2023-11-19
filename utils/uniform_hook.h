
#include "gbbs/gbbs.h"

namespace gbbs {

template <typename W>
auto generate_uniform_hook(size_t n) {
  // Each edge is picked uniformly at random among the previous i
  // nodes.
  using edge = std::tuple<uintE, uintE, W>;

  auto m = n-1;
  auto edges = sequence<edge>::uninitialized(m);
  parlay::random r;
  parallel_for(1, n, [&](uintE i) {
    // Generate the edge for node 1
    auto target = r.ith_rand(i) % i;
    edges[i-1] = {i, target, 1};
  });
  std::cout << "Generated edges" << std::endl;
  auto edge_list = edge_array<W>(std::move(edges), n);
  return gbbs_io::edge_list_to_symmetric_graph<W>(edge_list);
}

} // namespace gbbs