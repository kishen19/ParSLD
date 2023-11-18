#include "gbbs/gbbs.h"

namespace gbbs {

template <typename W>
auto generate_full_binary_tree(size_t k) {
  using edge = std::tuple<uintE, uintE, W>;

  auto n = (1UL << k) - 1;
  auto m = n-1;

  auto edges = sequence<edge>::uninitialized(m);
  parallel_for(0, n/2, [&](uintE i){
    edges[2*i] = {i, 2*i+1, 1};
    edges[2*i+1] = {i, 2*i+2, 1};
  });
  auto edge_list = edge_array<W>(std::move(edges), n);
  return gbbs_io::edge_list_to_symmetric_graph<W>(edge_list);
}

} // namespace gbbs