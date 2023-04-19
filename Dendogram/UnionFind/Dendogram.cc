// Usage:
// ./Dendogram -s -rounds 1 ../../inputs/rMatGraph_WJ_5_100
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm

#include "Dendogram.h"

namespace gbbs {
namespace {

template <class Graph>
double Dendogram_runner(Graph& G, commandLine P) {
	std::cout << "### Application: Dendogram Construction" << std::endl;
	std::cout << "### Graph: " << P.getArgument(0) << std::endl;
	std::cout << "### Threads: " << num_workers() << std::endl;
	std::cout << "### n: " << G.n << std::endl;
	std::cout << "### m: " << G.m << std::endl;
	std::cout << "### ------------------------------------" << std::endl;
	assert(P.getOption("-s"));

	timer t;
	t.start();
	Dendogram(G);
	double tt = t.stop();

	std::cout << "### Running Time: " << tt << std::endl;
	return tt;
}

}  // namespace
}  // namespace gbbs

generate_symmetric_weighted_main(gbbs::Dendogram_runner, true);
