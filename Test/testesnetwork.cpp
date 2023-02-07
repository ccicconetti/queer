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

#include "QuantumRouting/esnetwork.h"
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

struct TestEsNetwork : public ::testing::Test {
  using Apps = std::vector<EsNetwork::AppDescriptor>;

  //   /--> 1 -- >2 -+
  //  /              v
  // 0               3   all weights are 4, except 0->4 which is 1
  //  \              ^
  //   \---> 4 ------+
  EsNetwork::WeightVector exampleEdgeWeights() {
    return EsNetwork::WeightVector({
        {0, 1, 4},
        {1, 2, 4},
        {2, 3, 4},
        {0, 4, 1},
        {4, 3, 4},
    });
  }
};

TEST_F(TestEsNetwork, test_app_route_algo) {
  for (const auto& myAlgo : allAppRouteAlgos()) {
    ASSERT_EQ(myAlgo, appRouteAlgofromString(toString(myAlgo)));
  }
  ASSERT_EQ("unknown", toString(static_cast<AppRouteAlgo>(999)));
  ASSERT_THROW(appRouteAlgofromString("not-existing-algo"), std::runtime_error);
}

TEST_F(TestEsNetwork, test_measurement_probability) {
  EsNetwork myNetwork(exampleEdgeWeights());
  ASSERT_FLOAT_EQ(1, myNetwork.measurementProbability());
  myNetwork.measurementProbability(0.314);
  ASSERT_FLOAT_EQ(0.314, myNetwork.measurementProbability());
  ASSERT_THROW(myNetwork.measurementProbability(-0.5), std::runtime_error);
  ASSERT_THROW(myNetwork.measurementProbability(2), std::runtime_error);
}

TEST_F(TestEsNetwork, test_route_flows) {
  EsNetwork myNetwork(exampleEdgeWeights());
  myNetwork.measurementProbability(0.5);

  // no route existing
  std::vector<EsNetwork::FlowDescriptor> myFlows({
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
  ASSERT_EQ(EsNetwork::WeightVector({
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
  ASSERT_EQ(EsNetwork::WeightVector({
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
  ASSERT_EQ(EsNetwork::WeightVector({
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

TEST_F(TestEsNetwork, test_add_capacity_to_edge) {
  EsNetwork myNetwork(exampleEdgeWeights());
  myNetwork.measurementProbability(0.5);

  // add one (admissible) flow
  const auto myCapacityTot = myNetwork.totalCapacity();
  std::vector<EsNetwork::FlowDescriptor> myFlows({
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
  std::vector<EsNetwork::FlowDescriptor> myOtherFlows({
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

TEST_F(TestEsNetwork, test_route_flows_another) {
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

  EsNetwork myNetwork(myWeights);
  myNetwork.measurementProbability(0.5);

  std::vector<EsNetwork::FlowDescriptor> myFlows({
      {0, 3, 0.1},
  });
  myNetwork.route(myFlows);
  ASSERT_EQ(1, myFlows.size());
  ASSERT_EQ(1, myFlows[0].theDijsktra);
  ASSERT_EQ(std::vector<unsigned long>({4, 3}), myFlows[0].thePath);
}

TEST_F(TestEsNetwork, test_route_apps_drr) {
  EsNetwork myNetwork(exampleEdgeWeights());
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
  ASSERT_EQ(EsNetwork::AppDescriptor::Hops({1, 2}),
            myApps[0].theAllocated[2][0].theHops);
  ASSERT_EQ(1, myApps[0].theAllocated[3].size());
  ASSERT_EQ(EsNetwork::AppDescriptor::Hops({4, 3}),
            myApps[0].theAllocated[3][0].theHops);
  ASSERT_TRUE(myApps[1].theRemainingPaths.empty());
  ASSERT_EQ(4, myApps[1].theVisits);
  ASSERT_EQ(1, myApps[1].theAllocated.size());
  ASSERT_EQ(1, myApps[1].theAllocated[3].size());
  ASSERT_EQ(EsNetwork::AppDescriptor::Hops({2, 3}),
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

TEST_F(TestEsNetwork, test_route_apps_rnd) {
  EsNetwork myNetwork(exampleEdgeWeights());
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
  ASSERT_EQ(EsNetwork::AppDescriptor::Hops({4, 3}),
            myApps[0].theAllocated[3][0].theHops);
  ASSERT_FLOAT_EQ(0.5, myApps[0].netRate());
  ASSERT_FLOAT_EQ(1, myApps[0].grossRate());

  ASSERT_EQ(2, myApps[1].theVisits);
  ASSERT_EQ(1, myApps[1].theAllocated.size());
  ASSERT_EQ(1, myApps[1].theAllocated[3].size());
  ASSERT_EQ(EsNetwork::AppDescriptor::Hops({2, 3}),
            myApps[1].theAllocated[3][0].theHops);
  ASSERT_FLOAT_EQ(2, myApps[1].netRate());
  ASSERT_FLOAT_EQ(4, myApps[1].grossRate());

  // look for a different allocation
  std::map<unsigned long, std::set<std::string>> myAllocations;
  for (std::size_t mySeed = 0; mySeed < 100; mySeed++) {
    EsNetwork myAnotherNetwork(exampleEdgeWeights());
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

TEST_F(TestEsNetwork, test_route_apps_bft) {
  EsNetwork myNetwork(exampleEdgeWeights());
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
  ASSERT_EQ(EsNetwork::AppDescriptor::Hops({4, 3}),
            myApps[0].theAllocated[3][0].theHops);
  ASSERT_FLOAT_EQ(0.5, myApps[0].netRate());
  ASSERT_FLOAT_EQ(1, myApps[0].grossRate());

  ASSERT_EQ(3, myApps[1].theVisits);
  ASSERT_EQ(1, myApps[1].theAllocated.size());
  ASSERT_EQ(1, myApps[1].theAllocated[2].size());
  ASSERT_EQ(EsNetwork::AppDescriptor::Hops({2}),
            myApps[1].theAllocated[2][0].theHops);
  ASSERT_FLOAT_EQ(4, myApps[1].netRate());
  ASSERT_FLOAT_EQ(4, myApps[1].grossRate());

  ASSERT_EQ(2, myApps[2].theVisits);
  ASSERT_EQ(1, myApps[2].theAllocated.size());
  ASSERT_EQ(1, myApps[2].theAllocated[3].size());
  ASSERT_EQ(EsNetwork::AppDescriptor::Hops({3}),
            myApps[2].theAllocated[3][0].theHops);
  ASSERT_FLOAT_EQ(4, myApps[2].netRate());
  ASSERT_FLOAT_EQ(4, myApps[2].grossRate());
}

TEST_F(TestEsNetwork, test_max_net_rate) {
  EsNetwork myNetwork(exampleEdgeWeights());
  myNetwork.measurementProbability(0.5);
  const std::vector<unsigned long> myNoPeers;

  // path 0-1-2-3, min capacity = 4 and 2 swaps
  ASSERT_FLOAT_EQ(
      1,
      myNetwork.maxNetRate(EsNetwork::AppDescriptor(0, myNoPeers, 1, 0.5), 3));

  // path 0-4, min capacity = 1 and 0 swaps
  ASSERT_FLOAT_EQ(
      1,
      myNetwork.maxNetRate(EsNetwork::AppDescriptor(0, myNoPeers, 1, 0.5), 3));

  // path 0-4-3, min capacity = 1 and 1 swap
  ASSERT_FLOAT_EQ(
      0.5,
      myNetwork.maxNetRate(EsNetwork::AppDescriptor(0, myNoPeers, 1, 0.5),
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
  ASSERT_FLOAT_EQ(
      0,
      myNetwork.maxNetRate(EsNetwork::AppDescriptor(3, myNoPeers, 1, 0.5), 0));
}

} // namespace qr
} // namespace uiiit
