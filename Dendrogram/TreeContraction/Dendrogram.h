#pragma once

#include "gbbs/gbbs.h"
#include "benchmarks/Connectivity/SimpleUnionAsync/Connectivity.h"
#include "benchmarks/MinimumSpanningForest/Boruvka/MinimumSpanningForest.h"
#include "Dendrogram/heaps/leftist_heap.h"
#include "Dendrogram/heaps/skew_heap.h"
#include "Dendrogram/heaps/pairing_heap.h"
#include "Dendrogram/TreeContraction/TreeContraction.h"

namespace gbbs {

template <class Graph>
void DendrogramContract(Graph& GA){
	using W = typename Graph::weight_type;
	using kv = std::tuple<W, size_t>;
	using heap = pairing_heap::heap<kv>;
	auto tc = tree_contraction::tree_contraction<Graph, heap>(GA);
	tc.run();
}

}  // namespace gbbs