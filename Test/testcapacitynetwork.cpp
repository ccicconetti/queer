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

#define ROUTE_DRR(q, k)                                                        \
  myNetwork.route(myApps, AppRouteAlgo::Drr, q, myRouteRv, k)
#define ROUTE_RND(k)                                                           \
  myNetwork.route(myApps, AppRouteAlgo::Random, -1, myRouteRv, k)
#define ROUTE_BFT(k)                                                           \
  myNetwork.route(myApps, AppRouteAlgo::BestFit, -1, myRouteRv, k)

namespace uiiit {
namespace qr {

struct TestCapacityNetwork : public ::testing::Test {
  using Apps = std::vector<CapacityNetwork::AppDescriptor>;

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

TEST_F(TestCapacityNetwork, test_app_route_algo) {
  for (const auto& myAlgo : allAppRouteAlgos()) {
    ASSERT_EQ(myAlgo, appRouteAlgofromString(toString(myAlgo)));
  }
  ASSERT_EQ("unknown", toString(static_cast<AppRouteAlgo>(999)));
  ASSERT_THROW(appRouteAlgofromString("not-existing-algo"), std::runtime_error);
}

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

TEST_F(TestCapacityNetwork, test_measurement_probability) {
  CapacityNetwork myNetwork(exampleEdgeWeights());
  ASSERT_FLOAT_EQ(1, myNetwork.measurementProbability());
  myNetwork.measurementProbability(0.314);
  ASSERT_FLOAT_EQ(0.314, myNetwork.measurementProbability());
  ASSERT_THROW(myNetwork.measurementProbability(-0.5), std::runtime_error);
  ASSERT_THROW(myNetwork.measurementProbability(2), std::runtime_error);
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

TEST_F(TestCapacityNetwork, test_route_flows) {
  CapacityNetwork myNetwork(exampleEdgeWeights());
  myNetwork.measurementProbability(0.5);

  // no route existing
  std::vector<CapacityNetwork::FlowDescriptor> myFlows({
      {3, 0, 1.0},
  });
  myNetwork.route(myFlows);
  ASSERT_EQ(1, myFlows.size());
  ASSERT_TRUE(myFlows[0].thePath.empty());
  ASSERT_EQ(1, myFlows[0].theDijsktra);

  // add an unfeasible and a feasible route
  myFlows.clear();
  myFlows.emplace_back(3, 0, 1.0);
  myFlows.emplace_back(0, 3, 1.0);
  myNetwork.route(myFlows);
  ASSERT_EQ(2, myFlows.size());
  ASSERT_TRUE(myFlows[0].thePath.empty());
  ASSERT_FLOAT_EQ(0, myFlows[0].theGrossRate);
  ASSERT_EQ(1, myFlows[0].theDijsktra);
  ASSERT_EQ(std::vector<unsigned long>({1, 2, 3}), myFlows[1].thePath);
  ASSERT_FLOAT_EQ(4, myFlows[1].theGrossRate);
  ASSERT_EQ(2, myFlows[1].theDijsktra);
  ASSERT_EQ(CapacityNetwork::WeightVector({
                {0, 1, 0},
                {1, 2, 0},
                {2, 3, 0},
                {0, 4, 1},
                {4, 3, 4},
            }),
            myNetwork.weights());

  // the same route is not feasible anymore
  myFlows.clear();
  myFlows.emplace_back(0, 3, 1.0);
  myNetwork.route(myFlows);
  ASSERT_EQ(1, myFlows.size());
  ASSERT_TRUE(myFlows[0].thePath.empty());

  // request with smaller capacity, but cannot be admitted due to constraint
  myFlows.clear();
  myFlows.emplace_back(0, 3, 0.5);
  myNetwork.route(myFlows,
                  [](const auto& aFlow) { return aFlow.thePath.size() == 1; });
  ASSERT_EQ(1, myFlows.size());
  ASSERT_TRUE(myFlows[0].thePath.empty());

  // same request without constrating can be admitted
  myFlows.clear();
  myFlows.emplace_back(0, 3, 0.5);
  myNetwork.route(myFlows);
  ASSERT_EQ(1, myFlows.size());
  ASSERT_EQ(std::vector<unsigned long>({4, 3}), myFlows[0].thePath);
  ASSERT_FLOAT_EQ(1, myFlows[0].theGrossRate);
  ASSERT_EQ(CapacityNetwork::WeightVector({
                {0, 1, 0},
                {1, 2, 0},
                {2, 3, 0},
                {0, 4, 0},
                {4, 3, 3},
            }),
            myNetwork.weights());

  // add a request for an adjacent node
  myFlows.clear();
  myFlows.emplace_back(4, 3, 3);
  myNetwork.route(myFlows);
  ASSERT_EQ(1, myFlows.size());
  ASSERT_EQ(std::vector<unsigned long>({3}), myFlows[0].thePath);
  ASSERT_FLOAT_EQ(3, myFlows[0].theGrossRate);
  ASSERT_EQ(CapacityNetwork::WeightVector({
                {0, 1, 0},
                {1, 2, 0},
                {2, 3, 0},
                {0, 4, 0},
                {4, 3, 0},
            }),
            myNetwork.weights());
  ASSERT_FLOAT_EQ(0, myNetwork.totalCapacity());

  // no request can be served now
  myFlows.clear();
  for (size_t i = 0; i < 5; i++) {
    for (size_t j = 0; j < 5; j++) {
      if (i != j) {
        myFlows.emplace_back(i, j, 0.001);
      }
    }
  }
  myNetwork.route(myFlows);
  for (const auto& myFlow : myFlows) {
    ASSERT_TRUE(myFlow.thePath.empty());
    ASSERT_FLOAT_EQ(0, myFlow.theGrossRate);
  }

  // add ill-formed requests
  myFlows.clear();
  myFlows.emplace_back(0, 0, 1);
  ASSERT_THROW(myNetwork.route(myFlows), std::runtime_error);
  myFlows.clear();
  myFlows.emplace_back(0, 1, 0);
  ASSERT_THROW(myNetwork.route(myFlows), std::runtime_error);
  myFlows.clear();
  myFlows.emplace_back(0, 1, -1);
  ASSERT_THROW(myNetwork.route(myFlows), std::runtime_error);
  myFlows.clear();
  myFlows.emplace_back(0, 99, 1);
  ASSERT_THROW(myNetwork.route(myFlows), std::runtime_error);
  myFlows.clear();
  myFlows.emplace_back(99, 0, 1);
  ASSERT_THROW(myNetwork.route(myFlows), std::runtime_error);
}

TEST_F(TestCapacityNetwork, test_route_flows_another) {
  // swap weights
  auto myWeights = exampleEdgeWeights();
  for (auto& elem : myWeights) {
    auto& myWeight = std::get<2>(elem);
    if (myWeight == 1) {
      myWeight = 4;
    } else {
      myWeight = 1;
    }
  }

  CapacityNetwork myNetwork(myWeights);
  myNetwork.measurementProbability(0.5);

  std::vector<CapacityNetwork::FlowDescriptor> myFlows({
      {0, 3, 0.1},
  });
  myNetwork.route(myFlows);
  ASSERT_EQ(1, myFlows.size());
  ASSERT_EQ(1, myFlows[0].theDijsktra);
  ASSERT_EQ(std::vector<unsigned long>({4, 3}), myFlows[0].thePath);
}

TEST_F(TestCapacityNetwork, test_route_apps_drr) {
  CapacityNetwork myNetwork(exampleEdgeWeights());
  myNetwork.measurementProbability(0.5);
  support::UniformRv myRouteRv(0, 1, 42, 0, 0);
  Apps               myApps;

  // ill-formed requests
  myApps = Apps({{0, {0}, 1, 0}});
  ASSERT_THROW(ROUTE_DRR(1, 1), std::runtime_error);
  myApps = Apps({{0, {42}, 1, 0}});
  ASSERT_THROW(ROUTE_DRR(1, 1), std::runtime_error);
  myApps = Apps({{0, {1}, 0, 0}});
  ASSERT_THROW(ROUTE_DRR(1, 1), std::runtime_error);
  myApps = Apps({{0, {1}, -1, 0}});
  ASSERT_THROW(ROUTE_DRR(1, 1), std::runtime_error);
  myApps = Apps({{0, {1}, 1, 0}});
  ASSERT_THROW(ROUTE_DRR(0, 1), std::runtime_error);
  myApps = Apps({{0, {1}, 1, 0}});
  ASSERT_THROW(ROUTE_DRR(-1, 1), std::runtime_error);
  myApps = Apps({{0, {1}, 1, 0}});
  ASSERT_THROW(ROUTE_DRR(1, 0), std::runtime_error);

  // no route existing
  myApps = Apps({
      {3, {2, 0}, 1, 0},
      {2, {1}, 1, 0},
  });
  ROUTE_DRR(1.4, 99);
  ASSERT_EQ(2, myApps.size());
  ASSERT_EQ(0, myApps[0].theAllocated.size());
  ASSERT_FLOAT_EQ(0, myApps[0].grossRate());
  ASSERT_EQ(0, myApps[1].theAllocated.size());
  ASSERT_FLOAT_EQ(0, myApps[1].grossRate());

  // existing routes
  myApps = Apps({
      {0, {2, 3}, 1, 0},
      {1, {3}, 1, 0},
  });
  ROUTE_DRR(1.4, 99);
  ASSERT_EQ(2, myApps.size());
  ASSERT_TRUE(myApps[0].theRemainingPaths.empty());
  ASSERT_EQ(8, myApps[0].theVisits);
  ASSERT_EQ(2, myApps[0].theAllocated.size());
  ASSERT_EQ(1, myApps[0].theAllocated[2].size());
  ASSERT_EQ(CapacityNetwork::AppDescriptor::Hops({1, 2}),
            myApps[0].theAllocated[2][0].theHops);
  ASSERT_EQ(1, myApps[0].theAllocated[3].size());
  ASSERT_EQ(CapacityNetwork::AppDescriptor::Hops({4, 3}),
            myApps[0].theAllocated[3][0].theHops);
  ASSERT_TRUE(myApps[1].theRemainingPaths.empty());
  ASSERT_EQ(4, myApps[1].theVisits);
  ASSERT_EQ(1, myApps[1].theAllocated.size());
  ASSERT_EQ(1, myApps[1].theAllocated[3].size());
  ASSERT_EQ(CapacityNetwork::AppDescriptor::Hops({2, 3}),
            myApps[1].theAllocated[3][0].theHops);

  double myGrossRate = 0;
  double myNetRate   = 0;
  for (const auto& myApp : myApps) {
    myGrossRate += myApp.grossRate();
    myNetRate += myApp.netRate();
  }
  ASSERT_FLOAT_EQ(5, myGrossRate);
  ASSERT_FLOAT_EQ(2.5, myNetRate);
  ASSERT_FLOAT_EQ(7, myNetwork.totalCapacity());
  const auto myWeights = myNetwork.weights();
  ASSERT_EQ(3, myWeights.size());
  ASSERT_EQ(0, std::get<0>(myWeights[0]));
  ASSERT_EQ(1, std::get<1>(myWeights[0]));
  ASSERT_FLOAT_EQ(1.9, std::get<2>(myWeights[0]));
  ASSERT_EQ(2, std::get<0>(myWeights[1]));
  ASSERT_EQ(3, std::get<1>(myWeights[1]));
  ASSERT_FLOAT_EQ(2.1, std::get<2>(myWeights[1]));
  ASSERT_EQ(4, std::get<0>(myWeights[2]));
  ASSERT_EQ(3, std::get<1>(myWeights[2]));
  ASSERT_FLOAT_EQ(3, std::get<2>(myWeights[2]));

  // consume the remaining capacity
  myApps = Apps({
      {0, {1, 2, 3, 4}, 1, 0}, // only 0->1 is still available
      {2, {0, 1, 3, 4}, 1, 0}, // same for 2->3
      {4, {0, 1, 2, 3}, 1, 0}, // same for 4->3
  });
  ROUTE_DRR(0.1, 99);
  ASSERT_EQ(3, myApps.size());
  ASSERT_EQ(1, myApps[0].theAllocated.size());
  ASSERT_EQ(1, myApps[1].theAllocated.size());
  ASSERT_EQ(1, myApps[2].theAllocated.size());
  ASSERT_EQ(1, myApps[0].theAllocated[1].size());
  ASSERT_EQ(1, myApps[1].theAllocated[3].size());
  ASSERT_EQ(1, myApps[2].theAllocated[3].size());
  ASSERT_EQ(58, myApps[0].theVisits);
  ASSERT_EQ(64, myApps[1].theVisits);
  ASSERT_EQ(92, myApps[2].theVisits);
  ASSERT_FLOAT_EQ(0, myNetwork.totalCapacity());
}

TEST_F(TestCapacityNetwork, test_route_apps_rnd) {
  CapacityNetwork myNetwork(exampleEdgeWeights());
  myNetwork.measurementProbability(0.5);
  support::UniformRv myRouteRv(0, 1, 42, 0, 0);

  Apps myApps = Apps({
      {0, {2, 3}, 1, 0},
      {1, {3}, 1, 0},
  });
  ROUTE_RND(99);

  ASSERT_EQ(2, myApps.size());
  for (const auto& myApp : myApps) {
    ASSERT_TRUE(myApp.theRemainingPaths.empty());
  }
  LOG(INFO) << myApps[0].toString();
  LOG(INFO) << myApps[1].toString();

  ASSERT_EQ(4, myApps[0].theVisits);
  ASSERT_EQ(1, myApps[0].theAllocated.size());
  ASSERT_EQ(1, myApps[0].theAllocated[3].size());
  ASSERT_EQ(CapacityNetwork::AppDescriptor::Hops({4, 3}),
            myApps[0].theAllocated[3][0].theHops);
  ASSERT_FLOAT_EQ(0.5, myApps[0].netRate());
  ASSERT_FLOAT_EQ(1, myApps[0].grossRate());

  ASSERT_EQ(2, myApps[1].theVisits);
  ASSERT_EQ(1, myApps[1].theAllocated.size());
  ASSERT_EQ(1, myApps[1].theAllocated[3].size());
  ASSERT_EQ(CapacityNetwork::AppDescriptor::Hops({2, 3}),
            myApps[1].theAllocated[3][0].theHops);
  ASSERT_FLOAT_EQ(2, myApps[1].netRate());
  ASSERT_FLOAT_EQ(4, myApps[1].grossRate());

  // look for a different allocation
  std::map<unsigned long, std::set<std::string>> myAllocations;
  for (std::size_t mySeed = 0; mySeed < 100; mySeed++) {
    CapacityNetwork myAnotherNetwork(exampleEdgeWeights());
    myAnotherNetwork.measurementProbability(0.5);
    support::UniformRv myAnotherRouteRv(0, 1, mySeed, 0, 0);
    myApps = Apps({
        {0, {2, 3}, 1, 0},
        {1, {3}, 1, 0},
    });
    myAnotherNetwork.route(
        myApps, AppRouteAlgo::Random, -1, myAnotherRouteRv, 99);
    for (unsigned long i = 0; i < myApps.size(); i++) {
      myAllocations[i].insert(myApps[i].toString());
    }
  }

  for (const auto& elem : myAllocations) {
    VLOG(1) << "app " << elem.first;
    for (const auto& alloc : elem.second) {
      VLOG(2) << "\t" << alloc;
    }
  }
  ASSERT_EQ(2, myAllocations[0].size());
  ASSERT_EQ(2, myAllocations[1].size());
}

TEST_F(TestCapacityNetwork, test_route_apps_bft) {
  CapacityNetwork myNetwork(exampleEdgeWeights());
  myNetwork.measurementProbability(0.5);
  support::UniformRv myRouteRv(0, 1, 42, 0, 0);

  Apps myApps = Apps({
      {0, {3}, 1, 0},
      {1, {2, 3}, 1, 0},
      {2, {3}, 1, 0},
  });
  ROUTE_BFT(99);

  ASSERT_EQ(3, myApps.size());
  for (const auto& myApp : myApps) {
    ASSERT_TRUE(myApp.theRemainingPaths.empty());
  }
  LOG(INFO) << myApps[0].toString();
  LOG(INFO) << myApps[1].toString();
  LOG(INFO) << myApps[2].toString();

  ASSERT_EQ(3, myApps[0].theVisits);
  ASSERT_EQ(1, myApps[0].theAllocated.size());
  ASSERT_EQ(1, myApps[0].theAllocated[3].size());
  ASSERT_EQ(CapacityNetwork::AppDescriptor::Hops({4, 3}),
            myApps[0].theAllocated[3][0].theHops);
  ASSERT_FLOAT_EQ(0.5, myApps[0].netRate());
  ASSERT_FLOAT_EQ(1, myApps[0].grossRate());

  ASSERT_EQ(3, myApps[1].theVisits);
  ASSERT_EQ(1, myApps[1].theAllocated.size());
  ASSERT_EQ(1, myApps[1].theAllocated[2].size());
  ASSERT_EQ(CapacityNetwork::AppDescriptor::Hops({2}),
            myApps[1].theAllocated[2][0].theHops);
  ASSERT_FLOAT_EQ(4, myApps[1].netRate());
  ASSERT_FLOAT_EQ(4, myApps[1].grossRate());

  ASSERT_EQ(2, myApps[2].theVisits);
  ASSERT_EQ(1, myApps[2].theAllocated.size());
  ASSERT_EQ(1, myApps[2].theAllocated[3].size());
  ASSERT_EQ(CapacityNetwork::AppDescriptor::Hops({3}),
            myApps[2].theAllocated[3][0].theHops);
  ASSERT_FLOAT_EQ(4, myApps[2].netRate());
  ASSERT_FLOAT_EQ(4, myApps[2].grossRate());
}

TEST_F(TestCapacityNetwork, test_add_capacity_to_edge) {
  CapacityNetwork myNetwork(exampleEdgeWeights());
  myNetwork.measurementProbability(0.5);

  // add one (admissible) flow
  const auto myCapacityTot = myNetwork.totalCapacity();
  std::vector<CapacityNetwork::FlowDescriptor> myFlows({
      {0, 3, 1.0},
  });
  myNetwork.route(myFlows);
  ASSERT_EQ(1, myFlows.size());
  ASSERT_EQ(std::vector<unsigned long>({1, 2, 3}), myFlows[0].thePath);
  ASSERT_FLOAT_EQ(4, myFlows[0].theGrossRate);
  ASSERT_EQ(myCapacityTot - myFlows[0].thePath.size() * myFlows[0].theGrossRate,
            myNetwork.totalCapacity());

  // re-add the capacity along the path
  myNetwork.addCapacityToPath(0, {1, 2, 3}, myFlows[0].theGrossRate);
  ASSERT_EQ(myCapacityTot, myNetwork.totalCapacity());

  // re-add an identical flow
  std::vector<CapacityNetwork::FlowDescriptor> myOtherFlows({
      {0, 3, 1.0},
  });
  myNetwork.route(myOtherFlows);

  // add capacity partially
  ASSERT_EQ(myFlows[0].thePath, myOtherFlows[0].thePath);
  myNetwork.addCapacityToPath(2, {3}, myOtherFlows[0].theGrossRate);
  ASSERT_EQ(myCapacityTot - 2 * myOtherFlows[0].theGrossRate,
            myNetwork.totalCapacity());

  // remove too much capacity
  ASSERT_THROW(myNetwork.addCapacityToPath(2, {3}, -10), std::runtime_error);

  // non-existing edge
  ASSERT_THROW(myNetwork.addCapacityToPath(1, {0}, 1), std::runtime_error);

  ASSERT_NO_THROW(myNetwork.addCapacityToPath(0, {1}, 1));
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
  using P           = CapacityNetwork::AppDescriptor::Path;

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

TEST_F(TestCapacityNetwork, test_max_net_rate) {
  CapacityNetwork myNetwork(exampleEdgeWeights());
  myNetwork.measurementProbability(0.5);
  const std::vector<unsigned long> myNoPeers;

  // path 0-1-2-3, min capacity = 4 and 2 swaps
  ASSERT_FLOAT_EQ(1,
                  myNetwork.maxNetRate(
                      CapacityNetwork::AppDescriptor(0, myNoPeers, 1, 0.5), 3));

  // path 0-4, min capacity = 1 and 0 swaps
  ASSERT_FLOAT_EQ(1,
                  myNetwork.maxNetRate(
                      CapacityNetwork::AppDescriptor(0, myNoPeers, 1, 0.5), 3));

  // path 0-4-3, min capacity = 1 and 1 swap
  ASSERT_FLOAT_EQ(
      0.5,
      myNetwork.maxNetRate(CapacityNetwork::AppDescriptor(0, myNoPeers, 1, 0.5),
                           3,
                           [](const auto&, const auto& aPath) {
                             for (const auto& elem : aPath) {
                               if (elem.m_target == 2) {
                                 return false;
                               }
                             }
                             return true;
                           }));

  // non-existing path 3-0
  ASSERT_FLOAT_EQ(0,
                  myNetwork.maxNetRate(
                      CapacityNetwork::AppDescriptor(3, myNoPeers, 1, 0.5), 0));
}

} // namespace qr
} // namespace uiiit
