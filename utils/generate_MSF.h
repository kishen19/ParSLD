#include "gbbs/gbbs.h"
#include "benchmarks/MinimumSpanningForest/Boruvka/MinimumSpanningForest.h"

namespace gbbs {

template <class Graph, class CmdLine>
auto generate_MSF(Graph& G, CmdLine& P) {
  using W = typename Graph::weight_type;
  size_t n = G.n;
  auto mst_edges = MinimumSpanningForest_boruvka::MinimumSpanningForest(G);
  auto GMST = gbbs_io::edge_list_to_symmetric_graph<W>(edge_array<W>(std::move(mst_edges), n));
  auto out_f = P.getOptionValue("-write-tree");
  if (out_f) {
    gbbs_io::write_graph_to_file(out_f, GMST);
  }
  return GMST;
}

} // namespace gbbs