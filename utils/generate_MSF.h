#include "gbbs/gbbs.h"
#include "benchmarks/MinimumSpanningForest/Boruvka/MinimumSpanningForest.h"

namespace gbbs {

template <class Graph>
auto generate_MSF(Graph& G) {
  using W = typename Graph::weight_type;
  size_t n = G.n;
  auto mst_edges = MinimumSpanningForest_boruvka::MinimumSpanningForest(G);
  return gbbs_io::edge_list_to_symmetric_graph<W>(edge_array<W>(std::move(mst_edges), n));
}

} // namespace gbbs