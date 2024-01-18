#include "gbbs/gbbs.h"
#include "benchmarks/MinimumSpanningForest/Boruvka/MinimumSpanningForest.h"
#include "benchmarks/MinimumSpanningForest/Kruskal/MinimumSpanningForest.h"
#include "graph_io.h"

namespace gbbs {

template <class Graph, class CmdLine>
auto generate_MSF(Graph& G, CmdLine& P) {
  using W = typename Graph::weight_type;
  size_t n = G.n;
  //auto mst_edges = MinimumSpanningForest_boruvka::MinimumSpanningForest(G);
  std::cout << "Generating MSF." << std::endl;
  auto mst_edges = MinimumSpanningForest_kruskal::MinimumSpanningForest(G);
  auto GMST = gbbs_io::edge_list_to_symmetric_graph<W>(edge_array<W>(std::move(mst_edges), n));
  std::cout << "GMST size = " << GMST.m << std::endl;
  auto out_f = P.getOptionValue("-write-tree");
  if (out_f) {
    // gbbs_io::write_graph_to_file(out_f, GMST);
    gbbs::write_sym_graph_binary(GMST, out_f);
  }
  return GMST;
}

} // namespace gbbs