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
}

TEST_F(TestCapacityNetwork, test_route) {
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

} // namespace qr
} // namespace uiiit
