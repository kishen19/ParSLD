// Usage:
// ./Dendrogram -s -rounds 1 ../../inputs/rMatGraph_WJ_5_100
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm

// #include "DendAsync.h"
#include "DendSync.h"

namespace gbbs {
namespace {

template <class Graph>
double Dendrogram_runner(Graph& G, commandLine P) {
	std::cout << "### Application: Dendrogram Construction" << std::endl;
	std::cout << "### Graph: " << P.getArgument(0) << std::endl;
	std::cout << "### Threads: " << num_workers() << std::endl;
	std::cout << "### n: " << G.n << std::endl;
	std::cout << "### m: " << G.m << std::endl;
	std::cout << "### ------------------------------------" << std::endl;
	assert(P.getOption("-s"));

	double tt = DendrogramParUF(G);

	std::cout << "### Running Time: " << tt << std::endl;
	return tt;
}

}  // namespace
}  // namespace gbbs

generate_symmetric_weighted_main(gbbs::Dendrogram_runner, true);
