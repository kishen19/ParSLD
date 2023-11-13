#include "gbbs/gbbs.h"

#include <fstream>
#include <iostream>

namespace gbbs {

int BuildFullBinaryTree(int argc, char* argv[]) {
  commandLine P(argc, argv, "");

  size_t k = P.getOptionLongValue("-k", 20);
  auto out_f = P.getOptionValue("-outfile", "");

  if (out_f == "") {
    std::cout << "specify a valid outfile using -outfile" << std::endl;
    abort();
  }

  auto n = (1UL << k) - 1;
  auto m = 2*(n-1);
  std::cout << "WeightedAdjacencyGraph\n" << n << "\n" << m << std::endl;

  auto adj = sequence<sequence<uintE>>::uninitialized(n);
  auto wghs = sequence<sequence<uintE>>::uninitialized(n);
  parallel_for(0, n, [&](uintE i){
    if (i == 0){
      adj[i] = sequence<uintE>(2);
      wghs[i] = sequence<uintE>(2, 1);
      adj[i][0] = 2*i+1;
      adj[i][1] = 2*i+2;
    } else if (i >= n/2){
      adj[i] = sequence<uintE>(1);
      wghs[i] = sequence<uintE>(1, 1);
      adj[i][0] = (i-1)/2;
    } else{
      adj[i] = sequence<uintE>(3);
      wghs[i] = sequence<uintE>(3, 1);
      adj[i][0] = (i-1)/2;
      adj[i][1] = 2*i+1;
      adj[i][2] = 2*i+2;
    }
  });

  auto deg = sequence<uintE>::from_function(n, [&](uintE i){
    if (i==0) return 2;
    else if (i<n/2) return 3;
    else return 1;
  });
  parlay::scan_inplace(deg);
  auto temp = parlay::append(parlay::flatten(adj), parlay::flatten(wghs));
  auto output = parlay::append(deg, temp);
  auto C = parlay::sequence_to_string(output);

  size_t nn = C.size();
  std::ofstream file(out_f.c_str(), std::ios::out | std::ios::binary);
  if (!file.is_open()) {
    std::cout << "Unable to open file for writing: " << out_f << std::endl;
    return -1;
  }
  //  file << "# COO Format" << std::endl;
  //  file << "# n = " << n << std::endl;
  //  file << "# m = " << m << std::endl;

  file << "WeightedAdjacencyGraph\n" << n << "\n" << m << "\n";
  file.write(C.begin(), nn);
  file.close();

  std::cout << "done" << std::endl;
  return 0;
}

}  // namespace gbbs

int main(int argc, char* argv[]) { return gbbs::BuildFullBinaryTree(argc, argv); }