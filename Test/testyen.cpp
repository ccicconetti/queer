/*
              __ __ __
             |__|__|  | __
             |  |  |  ||__|
  ___ ___ __ |  |  |  |
 |   |   |  ||  |  |  |    Ubiquitous Internet @ IIT-CNR
 |   |   |  ||  |  |  |    C++ quantum routing libraries and tools
 |_______|__||__|__|__|    https://github.com/ccicconetti/quantum-routing

Licensed under the MIT License <http://opensource.org/licenses/MIT>
Copyright (c) 2022 C. Cicconetti <https://ccicconetti.github.io/>

Permission is hereby  granted, free of charge, to any  person obtaining a copy
of this software and associated  documentation files (the "Software"), to deal
in the Software  without restriction, including without  limitation the rights
to  use, copy,  modify, merge,  publish, distribute,  sublicense, and/or  sell
copies  of  the Software,  and  to  permit persons  to  whom  the Software  is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE  IS PROVIDED "AS  IS", WITHOUT WARRANTY  OF ANY KIND,  EXPRESS OR
IMPLIED,  INCLUDING BUT  NOT  LIMITED TO  THE  WARRANTIES OF  MERCHANTABILITY,
FITNESS FOR  A PARTICULAR PURPOSE AND  NONINFRINGEMENT. IN NO EVENT  SHALL THE
AUTHORS  OR COPYRIGHT  HOLDERS  BE  LIABLE FOR  ANY  CLAIM,  DAMAGES OR  OTHER
LIABILITY, WHETHER IN AN ACTION OF  CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE  OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "yen/yen_ksp.hpp"

#include "gtest/gtest.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/detail/adjacency_list.hpp>
#include <boost/graph/graph_selectors.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/property_map.hpp>

#include <glog/logging.h>

#include <cassert>
#include <optional>

namespace uiiit {
namespace qr {

struct TestYen : public ::testing::Test {
  template <class Graph>
  struct EdgeExtractor {
    EdgeExtractor(const Graph& aGraph)
        : theGraph(aGraph) {
      // noop
    }
    typename Graph::edge_descriptor operator()(const unsigned long aSrc,
                                               const unsigned long aDst) const {
      assert(boost::edge(aSrc, aDst, theGraph).second);
      return boost::edge(aSrc, aDst, theGraph).first;
    }
    const Graph& theGraph;
  };
};

TEST_F(TestYen, test_ksp) {
  using Graph =
      boost::adjacency_list<boost::vecS,
                            boost::vecS,
                            boost::bidirectionalS,
                            boost::no_property,
                            boost::property<boost::edge_weight_t, int>>;

  const std::list<std::tuple<unsigned long, unsigned long, int>> myEdges({
      {0, 1, 1},
      {0, 2, 1},
      {0, 3, 1},
      {1, 5, 1},
      {2, 5, 2},
      {3, 5, 3},
      {4, 5, 4},
      {3, 4, 1},
  });

  Graph         myGraph;
  EdgeExtractor edge(myGraph);
  auto          myWeights = boost::get(boost::edge_weight, myGraph);
  for (const auto& myEdge : myEdges) {
    const auto mySrc    = std::get<0>(myEdge);
    const auto myDst    = std::get<1>(myEdge);
    const auto myWeight = std::get<2>(myEdge);
    boost::add_edge(mySrc, myDst, myGraph);
    myWeights[edge(mySrc, myDst)] = myWeight;
  }

  using Path = std::list<Graph::edge_descriptor>;
  std::optional<unsigned> myK;
  auto                    myIndexMap = boost::get(boost::vertex_index, myGraph);

  // weighted graph, no K
  auto myResults = boost::yen_ksp(myGraph, 0, 5, myWeights, myIndexMap, myK);
  ASSERT_EQ(4, myResults.size());
  auto it = myResults.begin();
  EXPECT_EQ(2, it->first);
  EXPECT_EQ(Path({edge(0, 1), edge(1, 5)}), it->second);
  ++it;
  EXPECT_EQ(3, it->first);
  EXPECT_EQ(Path({edge(0, 2), edge(2, 5)}), it->second);
  ++it;
  EXPECT_EQ(4, it->first);
  EXPECT_EQ(Path({edge(0, 3), edge(3, 5)}), it->second);
  ++it;
  EXPECT_EQ(6, it->first);
  EXPECT_EQ(Path({edge(0, 3), edge(3, 4), edge(4, 5)}), it->second);
  ++it;
  ASSERT_EQ(myResults.end(), it);

  // weighted graph, K = 2
  myK       = 2;
  myResults = boost::yen_ksp(myGraph, 0, 5, myWeights, myIndexMap, myK);
  ASSERT_EQ(2, myResults.size());
  it = myResults.begin();
  EXPECT_EQ(2, it->first);
  EXPECT_EQ(Path({edge(0, 1), edge(1, 5)}), it->second);
  ++it;
  EXPECT_EQ(3, it->first);
  EXPECT_EQ(Path({edge(0, 2), edge(2, 5)}), it->second);
  ++it;
  ASSERT_EQ(myResults.end(), it);

  // unweighted graph, K = 2
  myResults =
      boost::yen_ksp(myGraph,
                     0,
                     5,
                     boost::make_static_property_map<Graph::edge_descriptor>(1),
                     myIndexMap,
                     myK);
  ASSERT_EQ(2, myResults.size());
  it = myResults.begin();
  EXPECT_EQ(2, it->first);
  EXPECT_EQ(Path({edge(0, 1), edge(1, 5)}), it->second);
  ++it;
  EXPECT_EQ(2, it->first);
  EXPECT_EQ(Path({edge(0, 2), edge(2, 5)}), it->second);
  ++it;
  ASSERT_EQ(myResults.end(), it);

  // unfeasible route
  myResults = boost::yen_ksp(myGraph, 5, 0, myWeights, myIndexMap, myK);
  ASSERT_TRUE(myResults.empty());
}

} // namespace qr
} // namespace uiiit
