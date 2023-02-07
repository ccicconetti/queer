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

#include "QuantumRouting/capacitynetwork.h"
#include "Support/random.h"
#include "Support/tostring.h"

#include "gtest/gtest.h"

#include <glog/logging.h>

#include <ctime>
#include <set>
#include <stdexcept>

namespace uiiit {
namespace qr {

struct TestCapacityNetwork : public ::testing::Test {
  CapacityNetwork::EdgeVector exampleEdges() {
    return CapacityNetwork::EdgeVector({
        {0, 1},
        {1, 2},
        {2, 3},
        {0, 4},
        {4, 3},
    });
  }

  //   /--> 1 -- >2 -+
  //  /              v
  // 0               3   all weights are 4, except 0->4 which is 1
  //  \              ^
  //   \---> 4 ------+
  CapacityNetwork::WeightVector exampleEdgeWeights() {
    return CapacityNetwork::WeightVector({
        {0, 1, 4},
        {1, 2, 4},
        {2, 3, 4},
        {0, 4, 1},
        {4, 3, 4},
    });
  }

  //
  //  +----> 1 <----+ +---> 4 ----+
  //  |             | |           |
  //  |             v v           v
  //  0              3            6 all weights are 1
  //  |             ^ ^           ^
  //  |             | |           |
  //  +----> 2 <----+ +---> 5 ----+
  //
  CapacityNetwork::WeightVector anotherExampleEdgeWeights() {
    return CapacityNetwork::WeightVector({
        {0, 1, 1},
        {0, 2, 1},
        {1, 3, 1},
        {2, 3, 1},
        {3, 1, 1},
        {3, 2, 1},
        {3, 4, 1},
        {3, 5, 1},
        {4, 3, 1},
        {4, 6, 1},
        {5, 3, 1},
        {5, 6, 1},
    });
  }

  using Vec = std::vector<unsigned long>;
  using Set = std::set<unsigned long>;
  Set vecToSet(const Vec& aVec) {
    Set ret;
    for (const auto& elem : aVec) {
      ret.insert(elem);
    }
    return ret;
  }
  bool setInSet(const Set& aOuter, const Set& aInner) {
    for (const auto& elem : aInner) {
      if (aOuter.count(elem) > 0) {
        return true;
      }
    }
    return false;
  }
};

TEST_F(TestCapacityNetwork, test_random_weights) {
  support::UniformRv myRv(0, 100, std::time(nullptr), 0, 0);

  for (const auto myBidir : std::vector({true, false})) {
    LOG(INFO) << "bidir = " << myBidir;
    CapacityNetwork myNetwork(exampleEdges(), myRv, myBidir);

    const auto myWeights = myNetwork.weights();
    ASSERT_EQ(myBidir ? 10 : 5, myWeights.size());
    std::set<double> myWeightSet;
    for (const auto& elem : myWeights) {
      myWeightSet.insert(std::get<2>(elem));
      ASSERT_TRUE(std::get<2>(elem) >= 0 and std::get<2>(elem) <= 100)
          << "(" << std::get<0>(elem) << "," << std::get<1>(elem) << ") ["
          << std::get<2>(elem) << "]";
    }
    ASSERT_EQ(5, myWeightSet.size());

    if (VLOG_IS_ON(1)) {
      myNetwork.toDot("TestCapacityNetwork.test_random_weights-" +
                      std::to_string(myBidir) + ".dot");
    }
  }
}

TEST_F(TestCapacityNetwork, test_graph_properties) {
  CapacityNetwork myNetwork(exampleEdgeWeights());
  ASSERT_EQ(5, myNetwork.numNodes());
  ASSERT_EQ(5, myNetwork.numEdges());
  ASSERT_FLOAT_EQ(17, myNetwork.totalCapacity());
  ASSERT_EQ(0, myNetwork.inDegree().first);
  ASSERT_EQ(2, myNetwork.inDegree().second);
  ASSERT_EQ(0, myNetwork.outDegree().first);
  ASSERT_EQ(2, myNetwork.outDegree().second);

  ASSERT_EQ(std::vector<double>({5, 4, 4, 0, 4}), myNetwork.nodeCapacities());
}

TEST_F(TestCapacityNetwork, test_reachable_nodes) {
  CapacityNetwork myNetwork(anotherExampleEdgeWeights());

  std::size_t myDiameter;

  const auto myAll = myNetwork.reachableNodes(0, 99, myDiameter);
  ASSERT_EQ(4, myDiameter);
  ASSERT_EQ(7, myAll.size());
  ASSERT_EQ(std::set<unsigned long>({1, 2, 3, 4, 5, 6}), myAll.find(0)->second);
  ASSERT_EQ(std::set<unsigned long>({2, 3, 4, 5, 6}), myAll.find(1)->second);
  ASSERT_EQ(std::set<unsigned long>({1, 3, 4, 5, 6}), myAll.find(2)->second);
  ASSERT_EQ(std::set<unsigned long>({1, 2, 4, 5, 6}), myAll.find(3)->second);
  ASSERT_EQ(std::set<unsigned long>({1, 2, 3, 5, 6}), myAll.find(4)->second);
  ASSERT_EQ(std::set<unsigned long>({1, 2, 3, 4, 6}), myAll.find(5)->second);
  ASSERT_EQ(std::set<unsigned long>({}), myAll.find(6)->second);

  const auto mySome = myNetwork.reachableNodes(0, 2, myDiameter);
  ASSERT_EQ(7, myAll.size());
  ASSERT_EQ(std::set<unsigned long>({1, 2, 3}), mySome.find(0)->second);
  ASSERT_EQ(std::set<unsigned long>({2, 3, 4, 5}), mySome.find(1)->second);
  ASSERT_EQ(std::set<unsigned long>({1, 3, 4, 5}), mySome.find(2)->second);
  ASSERT_EQ(std::set<unsigned long>({1, 2, 4, 5, 6}), mySome.find(3)->second);
  ASSERT_EQ(std::set<unsigned long>({1, 2, 3, 5, 6}), mySome.find(4)->second);
  ASSERT_EQ(std::set<unsigned long>({1, 2, 3, 4, 6}), mySome.find(5)->second);
  ASSERT_EQ(std::set<unsigned long>({}), mySome.find(6)->second);

  const auto myTwo = myNetwork.reachableNodes(2, 2, myDiameter);
  ASSERT_EQ(7, myTwo.size());
  ASSERT_EQ(std::set<unsigned long>({3}), myTwo.find(0)->second);
  ASSERT_EQ(std::set<unsigned long>({2, 4, 5}), myTwo.find(1)->second);
  ASSERT_EQ(std::set<unsigned long>({1, 4, 5}), myTwo.find(2)->second);
  ASSERT_EQ(std::set<unsigned long>({6}), myTwo.find(3)->second);
  ASSERT_EQ(std::set<unsigned long>({1, 2, 5}), myTwo.find(4)->second);
  ASSERT_EQ(std::set<unsigned long>({1, 2, 4}), myTwo.find(5)->second);
  ASSERT_EQ(std::set<unsigned long>({}), myTwo.find(6)->second);

  const auto myNone = myNetwork.reachableNodes(99, 99, myDiameter);
  ASSERT_EQ(7, myNone.size());
  for (const auto& elem : myNone) {
    ASSERT_TRUE(elem.second.empty());
  }
}

TEST_F(TestCapacityNetwork, test_closest_nodes) {
  support::UniformRv myRv(0, 1, 42, 0, 0);
  CapacityNetwork    myNetwork(exampleEdgeWeights());

  ASSERT_EQ(Set(), vecToSet(myNetwork.closestNodes(0, 0, myRv)));

  ASSERT_TRUE(setInSet(vecToSet(myNetwork.closestNodes(0, 1, myRv)), {1, 4}));
  ASSERT_EQ(Set({1, 4}), vecToSet(myNetwork.closestNodes(0, 2, myRv)));
  ASSERT_EQ(Set({1, 2, 4}), vecToSet(myNetwork.closestNodes(0, 3, myRv)));
  ASSERT_EQ(Set({1, 2, 3, 4}), vecToSet(myNetwork.closestNodes(0, 4, myRv)));
  ASSERT_EQ(Set({1, 2, 3, 4}), vecToSet(myNetwork.closestNodes(0, 99, myRv)));

  ASSERT_EQ(Set({2}), vecToSet(myNetwork.closestNodes(1, 1, myRv)));
  ASSERT_EQ(Set({2, 3}), vecToSet(myNetwork.closestNodes(1, 2, myRv)));

  ASSERT_EQ(Set(), vecToSet(myNetwork.closestNodes(3, 1, myRv)));
}

TEST_F(TestCapacityNetwork, test_closest_nodes_whitelist) {
  support::UniformRv myRv(0, 1, 42, 0, 0);
  CapacityNetwork    myNetwork(exampleEdgeWeights());

  ASSERT_EQ(Set(), vecToSet(myNetwork.closestNodes(0, 0, myRv, {1, 2, 3, 4})));
  ASSERT_EQ(Set({1}), vecToSet(myNetwork.closestNodes(0, 1, myRv, {1, 2, 3})));
  ASSERT_EQ(Set({4}), vecToSet(myNetwork.closestNodes(0, 1, myRv, {2, 3, 4})));
  ASSERT_EQ(Set({2}), vecToSet(myNetwork.closestNodes(0, 1, myRv, {2, 3})));
  ASSERT_EQ(Set({3}), vecToSet(myNetwork.closestNodes(0, 1, myRv, {3})));
  ASSERT_EQ(Set({1, 4}), vecToSet(myNetwork.closestNodes(0, 2, myRv)));
  ASSERT_EQ(Set({2, 3}), vecToSet(myNetwork.closestNodes(0, 2, myRv, {2, 3})));
  ASSERT_EQ(Set({3}), vecToSet(myNetwork.closestNodes(0, 99, myRv, {3})));
  ASSERT_EQ(Set(), vecToSet(myNetwork.closestNodes(0, 99, myRv, {0})));
  ASSERT_EQ(Set(), vecToSet(myNetwork.closestNodes(0, 99, myRv, {42})));
}

TEST_F(TestCapacityNetwork, test_closest_nodes_another) {
  support::UniformRv myRv(0, 1, 42, 0, 0);
  CapacityNetwork    myNetwork(anotherExampleEdgeWeights());

  ASSERT_EQ(Set(), vecToSet(myNetwork.closestNodes(0, 0, myRv)));
  ASSERT_TRUE(setInSet(vecToSet(myNetwork.closestNodes(0, 1, myRv)), {1, 2}));
  ASSERT_EQ(Set({1, 2}), vecToSet(myNetwork.closestNodes(0, 2, myRv)));
  ASSERT_EQ(Set({1, 2, 3}), vecToSet(myNetwork.closestNodes(0, 3, myRv)));
  auto myNodes = vecToSet(myNetwork.closestNodes(0, 4, myRv));
  ASSERT_TRUE(myNodes == Set({1, 2, 3, 4}) or myNodes == Set({1, 2, 3, 5}));
  ASSERT_EQ(Set({1, 2, 3, 4, 5}), vecToSet(myNetwork.closestNodes(0, 5, myRv)));
  ASSERT_EQ(Set({1, 2, 3, 4, 5, 6}),
            vecToSet(myNetwork.closestNodes(0, 6, myRv)));
  ASSERT_EQ(Set({1, 2, 3, 4, 5, 6}),
            vecToSet(myNetwork.closestNodes(0, 99, myRv)));

  ASSERT_EQ(Set({1, 2, 4, 5, 6}),
            vecToSet(myNetwork.closestNodes(3, 99, myRv)));
}

TEST_F(TestCapacityNetwork, test_min_capacity_nodes) {
  CapacityNetwork myNetwork(exampleEdgeWeights());

  ASSERT_FLOAT_EQ(4, myNetwork.minCapacity(0, {1, 2, 3}));
  ASSERT_FLOAT_EQ(4, myNetwork.minCapacity(1, {2, 3}));
  ASSERT_FLOAT_EQ(1, myNetwork.minCapacity(0, {4}));
  ASSERT_FLOAT_EQ(1, myNetwork.minCapacity(0, {4, 3}));

  ASSERT_FLOAT_EQ(std::numeric_limits<double>::max(),
                  myNetwork.minCapacity(0, {}));

  ASSERT_NO_THROW(myNetwork.minCapacity(0, {1}));
  ASSERT_THROW(myNetwork.minCapacity(1, {0}), std::runtime_error);
  ASSERT_THROW(myNetwork.minCapacity(0, {3}), std::runtime_error);
  ASSERT_THROW(myNetwork.minCapacity(0, {99}), std::runtime_error);
}

TEST_F(TestCapacityNetwork, test_min_capacity_edges) {
  CapacityNetwork myNetwork(exampleEdgeWeights());
  const auto&     G = myNetwork.theGraph;
  using P           = CapacityNetwork::Path;

  ASSERT_FLOAT_EQ(4,
                  CapacityNetwork::minCapacity(P({
                                                   boost::edge(0, 1, G).first,
                                                   boost::edge(1, 2, G).first,
                                                   boost::edge(2, 3, G).first,
                                               }),
                                               G));
  ASSERT_FLOAT_EQ(1,
                  CapacityNetwork::minCapacity(P({
                                                   boost::edge(0, 4, G).first,
                                                   boost::edge(4, 3, G).first,
                                               }),
                                               G));
  ASSERT_FLOAT_EQ(4,
                  CapacityNetwork::minCapacity(P({
                                                   boost::edge(4, 3, G).first,
                                               }),
                                               G));
  ASSERT_FLOAT_EQ(std::numeric_limits<double>::max(),
                  CapacityNetwork::minCapacity(P(), G));
}

} // namespace qr
} // namespace uiiit
