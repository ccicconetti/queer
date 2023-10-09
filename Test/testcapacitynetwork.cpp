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
#include "QuantumRouting/networkfactory.h"
#include "Support/random.h"
#include "Support/tostring.h"

#include "gtest/gtest.h"

#include <glog/logging.h>

#include <ctime>
#include <optional>
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

TEST_F(TestCapacityNetwork, test_make_capacity_network_waxman) {
  const std::size_t       myNodes       = 50;
  const double            myMaxDistance = 100;
  std::vector<Coordinate> myCoordinates;
  const auto myNetwork = makeCapacityNetworkWaxman<CapacityNetwork>(
      [myMaxDistance](const double d) {
        return 100e3 * std::exp(-d / myMaxDistance);
      },
      42,
      myNodes,
      myMaxDistance,
      0.5,
      0.5,
      myCoordinates);

  ASSERT_EQ(myNodes, myNetwork->numNodes());
  ASSERT_GT(myNetwork->numEdges(), 0);
  ASSERT_GT(myNetwork->totalCapacity(), 0);
  ASSERT_LT(myNetwork->totalCapacity() / myNetwork->numEdges(), 100e3);
  if (VLOG_IS_ON(1)) {
    myNetwork->toDot("network.dot");
  }
}

TEST_F(TestCapacityNetwork, test_cspf) {
  CapacityNetwork         myNetwork(exampleEdgeWeights());
  std::set<unsigned long> myDestinations({3, 4});
  using CspfRes = std::map<unsigned long, std::vector<unsigned long>>;

  ASSERT_EQ(CspfRes({{3, {4, 3}}, {4, {4}}}),
            myNetwork.cspf(0, 0, myDestinations));
  ASSERT_EQ(CspfRes({{3, {4, 3}}, {4, {4}}}),
            myNetwork.cspf(0, 1, myDestinations));
  ASSERT_EQ(CspfRes({{3, {1, 2, 3}}, {4, {}}}),
            myNetwork.cspf(0, 2, myDestinations));
  ASSERT_EQ(CspfRes({{3, {}}, {4, {}}}), myNetwork.cspf(0, 99, myDestinations));
  myDestinations = std::set<unsigned long>({0, 1, 2, 4});
  ASSERT_EQ(CspfRes({{0, {}}, {1, {}}, {2, {}}, {4, {}}}),
            myNetwork.cspf(3, 0, myDestinations));
}

TEST_F(TestCapacityNetwork, test_remove_capacity_from_path) {
  std::size_t     myDiameter;
  CapacityNetwork myNetwork(exampleEdgeWeights());
  ASSERT_FLOAT_EQ(17.0, myNetwork.totalCapacity());
  ASSERT_EQ(Set({1, 2, 3, 4}), myNetwork.reachableNodes(0, 99, myDiameter)[0]);

  // cannot remove that much capacity from 0->1
  ASSERT_THROW(CapacityNetwork::removeCapacityFromPath(
                   0, Vec({1}), 99.0, std::nullopt, myNetwork.theGraph),
               std::runtime_error);

  // edge does not exist
  ASSERT_THROW(CapacityNetwork::removeCapacityFromPath(
                   0, Vec({2}), 1.0, std::nullopt, myNetwork.theGraph),
               std::runtime_error);

  // node does not exist
  ASSERT_THROW(CapacityNetwork::removeCapacityFromPath(
                   99, Vec({1}), 1.0, std::nullopt, myNetwork.theGraph),
               std::runtime_error);
  ASSERT_THROW(CapacityNetwork::removeCapacityFromPath(
                   0, Vec({99}), 1.0, std::nullopt, myNetwork.theGraph),
               std::runtime_error);
  ASSERT_THROW(CapacityNetwork::removeCapacityFromPath(
                   0, Vec({1, 99}), 1.0, std::nullopt, myNetwork.theGraph),
               std::runtime_error);

  // remove completely the capacity from 0->1, but leave the edge
  CapacityNetwork::removeCapacityFromPath(
      0, Vec({1}), 4.0, std::nullopt, myNetwork.theGraph);
  ASSERT_FLOAT_EQ(13.0, myNetwork.totalCapacity());
  ASSERT_EQ(Set({1, 2, 3, 4}), myNetwork.reachableNodes(0, 99, myDiameter)[0]);

  // add some capacity to 0->1
  CapacityNetwork::removeCapacityFromPath(
      0, Vec({1}), -1.0, std::nullopt, myNetwork.theGraph);
  ASSERT_FLOAT_EQ(14.0, myNetwork.totalCapacity());

  // remove completely the capacity from 0->1 and prune the edge, too
  CapacityNetwork::removeCapacityFromPath(
      0, Vec({1}), 1.0, 0.1, myNetwork.theGraph);
  ASSERT_FLOAT_EQ(13.0, myNetwork.totalCapacity());
  ASSERT_EQ(Set({3, 4}), myNetwork.reachableNodes(0, 99, myDiameter)[0]);

  // remove capacity from 1->4->3 and prune the edges, as needed
  CapacityNetwork::removeCapacityFromPath(
      0, Vec({4, 3}), 1.0, 0.1, myNetwork.theGraph);
  ASSERT_FLOAT_EQ(11.0, myNetwork.totalCapacity());
  ASSERT_EQ(Set({}), myNetwork.reachableNodes(0, 99, myDiameter)[0]);

  // remove all edges, the graph is now empty
  CapacityNetwork::removeCapacityFromPath(
      4, Vec({3}), 0.0, 99, myNetwork.theGraph);
  CapacityNetwork::removeCapacityFromPath(
      1, Vec({2, 3}), 0.0, 99, myNetwork.theGraph);
  ASSERT_FLOAT_EQ(0.0, myNetwork.totalCapacity());
  ASSERT_EQ(0, myNetwork.numEdges());
}

TEST_F(TestCapacityNetwork, test_remove_capacity_from_path_alt) {
  std::size_t     myDiameter;
  CapacityNetwork myNetwork(exampleEdgeWeights());
  ASSERT_FLOAT_EQ(17.0, myNetwork.totalCapacity());
  ASSERT_EQ(Set({1, 2, 3, 4}), myNetwork.reachableNodes(0, 99, myDiameter)[0]);
  const auto& myGraph = myNetwork.theGraph;

  // cannot remove that much capacity from 0->1
  ASSERT_THROW(CapacityNetwork::removeCapacityFromPath(
                   {boost::edge(0, 1, myGraph).first},
                   99.0,
                   std::nullopt,
                   myNetwork.theGraph),
               std::runtime_error);

  // remove completely the capacity from 0->1, but leave the edge
  CapacityNetwork::removeCapacityFromPath({boost::edge(0, 1, myGraph).first},
                                          4.0,
                                          std::nullopt,
                                          myNetwork.theGraph);
  ASSERT_FLOAT_EQ(13.0, myNetwork.totalCapacity());
  ASSERT_EQ(Set({1, 2, 3, 4}), myNetwork.reachableNodes(0, 99, myDiameter)[0]);

  // add some capacity to 0->1
  CapacityNetwork::removeCapacityFromPath({boost::edge(0, 1, myGraph).first},
                                          -1.0,
                                          std::nullopt,
                                          myNetwork.theGraph);
  ASSERT_FLOAT_EQ(14.0, myNetwork.totalCapacity());

  // remove completely the capacity from 0->1 and prune the edge, too
  CapacityNetwork::removeCapacityFromPath(
      0, Vec({1}), 1.0, 0.1, myNetwork.theGraph);
  ASSERT_FLOAT_EQ(13.0, myNetwork.totalCapacity());
  ASSERT_EQ(Set({3, 4}), myNetwork.reachableNodes(0, 99, myDiameter)[0]);

  // remove capacity from 1->4->3 and prune the edges, as needed
  CapacityNetwork::removeCapacityFromPath(
      {boost::edge(0, 4, myGraph).first, boost::edge(4, 3, myGraph).first},
      1.0,
      0.1,
      myNetwork.theGraph);
  ASSERT_FLOAT_EQ(11.0, myNetwork.totalCapacity());
  ASSERT_EQ(Set({}), myNetwork.reachableNodes(0, 99, myDiameter)[0]);

  // remove all edges, the graph is now empty
  CapacityNetwork::removeCapacityFromPath(
      {boost::edge(4, 3, myGraph).first}, 0.0, 99, myNetwork.theGraph);
  CapacityNetwork::removeCapacityFromPath(
      {boost::edge(1, 2, myGraph).first, boost::edge(2, 3, myGraph).first},
      0.0,
      99,
      myNetwork.theGraph);
  ASSERT_FLOAT_EQ(0.0, myNetwork.totalCapacity());
  ASSERT_EQ(0, myNetwork.numEdges());
}

} // namespace qr
} // namespace uiiit
