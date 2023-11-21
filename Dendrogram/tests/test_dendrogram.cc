#include "Dendrogram/SeqUnionFind/Dendrogram.h"
#include "Dendrogram/ParUnionFind/Dendrogram.h"
#include "Dendrogram/RCtreeTracing/Dendrogram.h"

#include <unordered_set>

#include "gbbs/graph.h"
#include "gbbs/macros.h"
#include "gbbs/gbbs.h"
#include "gbbs/unit_tests/graph_test_utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "utils/caterpillar.h"
#include "utils/full_binary.h"
#include "utils/path.h"
#include "utils/star.h"
#include "utils/uniform_hook.h"

using ::testing::AnyOf;
using ::testing::ElementsAre;

namespace gbbs {

auto GetGraph(std::string s) {
  size_t n = 10000;
  size_t k = 14;
  if (s == "path") {
    return gbbs::generate_path_graph<gbbs::intE>(n);
  } else if (s == "star") {
    return gbbs::generate_star_graph<gbbs::intE>(n);
  } else if (s == "caterpillar") {
    return gbbs::generate_caterpillar_graph<gbbs::intE>(k);
  } else if (s == "fullb") {
    return gbbs::generate_full_binary_tree<gbbs::intE>(k);
  } else if (s == "unifhook") {
    return gbbs::generate_uniform_hook<gbbs::intE>(n);
  }
}

class MyTestSuite : public testing::TestWithParam<std::tuple<std::string, std::string>> {
};

TEST_P(MyTestSuite, RunAlgorithm) {
  auto E = GetGraph(std::get<0>(GetParam()));
  auto G = gbbs_io::edge_list_to_symmetric_graph<gbbs::intE>(E);
  auto parents_seq = DendrogramSeqUF(G);
  std::cout << "Generated parents!" << std::endl;
  parlay::sequence<uintE> parents_par;
  if (std::get<1>(GetParam()) == "ParUF") {
    parents_par = DendrogramParUF(G, 10);
  } else if (std::get<1>(GetParam()) == "RCtree") {
    parents_par = DendrogramRCtreeTracing(G);
  }
  std::cout << "Generated parents_par!" << std::endl;
  EXPECT_EQ(parents_seq.size(), parents_par.size());
  for (size_t i=0; i<parents_seq.size(); ++i) {
    EXPECT_EQ(parents_seq[i], parents_par[i]);
  }
}

INSTANTIATE_TEST_SUITE_P(
    MyGroup, MyTestSuite,
    testing::Combine(
        testing::Values("path", "star", "caterpillar", "fullb", "unifhook"),
        testing::Values("ParUF", "RCtree")),
    [](const testing::TestParamInfo<MyTestSuite::ParamType>& info) {
      std::string name = std::get<0>(info.param)+std::get<1>(info.param);
      return name;
    });

}  // namespace gbbs