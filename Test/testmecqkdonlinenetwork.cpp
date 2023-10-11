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

#include "QuantumRouting/mecqkdonlinenetwork.h"
#include "Support/random.h"

#include "gtest/gtest.h"

#include <glog/logging.h>

#include <vector>

namespace uiiit {
namespace qr {

struct TestMecQkdOnlineNetwork : public ::testing::Test {
  TestMecQkdOnlineNetwork() = default;

  // aFanIn == false
  //
  //   +--> 1 --> 2 --> 3 --> 4 [8]
  //  /
  // 0 --> 5 --> 6 --> 7 [4]
  // | \                      all weights are 2, except 3->4 which is 1
  // |  +---> 7 --> 8 [2]     the processing capacities are in square brackets
  // |
  //  \+--> 9 [1]
  //
  // aFanIn == true
  //
  //   +--> 1 --> 2 --> 3 --> 4 ---+
  //  /                            v
  // 0 --> 5 --> 6 --> 7 --------> 11 <-+
  // | \                           ^    |
  // |  +---> 7 --> 8 -------------|    |
  // |                                  |
  //  \+--> 9 --------------------------+
  //
  // all weights on the first branch are 8, on the second 6, on the third 4,
  // and on the last 2; the only edge node is 11, with a capacity of 8
  //
  static std::unique_ptr<MecQkdOnlineNetwork>
  makeNetwork(const MecQkdOnlineAlgo aAlgo, const bool aFanIn) {
    CapacityNetwork::WeightVector   myWeights;
    std::map<unsigned long, double> myEdgeProcessing;
    if (aFanIn) {
      myWeights        = CapacityNetwork::WeightVector({
          {0, 1, 2},
          {1, 2, 8},
          {2, 3, 8},
          {3, 4, 8},
          {4, 11, 8},
          {0, 5, 6},
          {5, 6, 6},
          {6, 7, 6},
          {7, 11, 6},
          {0, 8, 4},
          {8, 9, 4},
          {9, 11, 4},
          {0, 10, 2},
          {10, 11, 2},
      });
      myEdgeProcessing = std::map<unsigned long, double>({
          {11, 8},
      });
    } else {
      myWeights        = CapacityNetwork::WeightVector({
          {0, 1, 2},
          {1, 2, 2},
          {2, 3, 2},
          {3, 4, 1},
          {0, 5, 2},
          {5, 6, 2},
          {6, 7, 2},
          {0, 8, 2},
          {8, 9, 2},
          {0, 10, 2},
      });
      myEdgeProcessing = std::map<unsigned long, double>({
          {4, 8},
          {7, 4},
          {9, 2},
          {10, 1},
      });
    }
    auto ret  = std::make_unique<MecQkdOnlineNetwork>(myWeights);
    auto myRv = std::make_unique<support::DeterministicRv<std::vector<double>>>(
        std::vector<double>({0.1, 0.4, 0.25, 0.9, 0.6, 0.7, 0.7, 0.3}));
    ret->configure(
        aAlgo, std::forward<decltype(myRv)>(myRv), {0}, myEdgeProcessing);
    return ret;
  }
};

TEST_F(TestMecQkdOnlineNetwork, test_algorithms) {
  for (const auto& myAlgo : allMecQkdOnlineAlgos()) {
    const auto myAlgoString = toString(myAlgo);
    ASSERT_EQ(myAlgo, mecQkdOnlineAlgofromString(myAlgoString));
  }
}

TEST_F(TestMecQkdOnlineNetwork, test_policy014_k_1) {
  auto myNetwork = makeNetwork(MecQkdOnlineAlgo::Policy014_k1, false);
  ASSERT_FLOAT_EQ(19, myNetwork->totalCapacity());
  ASSERT_FLOAT_EQ(15, myNetwork->totProcessing());
  ASSERT_FLOAT_EQ(0, myNetwork->totNetRate());

  std::vector<uint64_t> myAllocated;

  // rate is too high
  MecQkdOnlineNetwork::Allocation myApp(0, 2.1, 1);
  myNetwork->add(myApp);
  ASSERT_TRUE(not myApp.allocated());

  // load is too high
  myApp.theRate = 1;
  myApp.theLoad = 8.1;
  myNetwork->add(myApp);
  ASSERT_TRUE(not myApp.allocated());

  // allocate one app on the closest node
  myApp.theRate = 1;
  myApp.theLoad = 0.5;
  myNetwork->add(myApp);
  ASSERT_TRUE(myApp.allocated());
  myAllocated.emplace_back(myApp.theId);
  ASSERT_EQ(10, myApp.edgeNode());
  ASSERT_EQ(1, myApp.thePath.size());
  ASSERT_FLOAT_EQ(18, myNetwork->totalCapacity());
  ASSERT_FLOAT_EQ(14.5, myNetwork->totProcessing());
  ASSERT_FLOAT_EQ(1, myNetwork->totNetRate());

  // allocate another app on the same node
  myApp.theRate = 1;
  myApp.theLoad = 0.1;
  myNetwork->add(myApp);
  ASSERT_TRUE(myApp.allocated());
  myAllocated.emplace_back(myApp.theId);
  ASSERT_EQ(10, myApp.edgeNode());
  ASSERT_EQ(1, myApp.thePath.size());
  ASSERT_FLOAT_EQ(17, myNetwork->totalCapacity());
  ASSERT_FLOAT_EQ(14.4, myNetwork->totProcessing());
  ASSERT_FLOAT_EQ(2, myNetwork->totNetRate());

  // now allocate an app on the three-hop path
  myApp.theRate = 2;
  myApp.theLoad = 4;
  myNetwork->add(myApp);
  ASSERT_TRUE(myApp.allocated());
  myAllocated.emplace_back(myApp.theId);
  ASSERT_EQ(7, myApp.edgeNode());
  ASSERT_EQ(3, myApp.thePath.size());
  ASSERT_FLOAT_EQ(11, myNetwork->totalCapacity());
  ASSERT_FLOAT_EQ(10.4, myNetwork->totProcessing());
  ASSERT_FLOAT_EQ(4, myNetwork->totNetRate());

  // now allocate an app on the four-hop path
  myApp.theRate = 1;
  myApp.theLoad = 3;
  myNetwork->add(myApp);
  ASSERT_TRUE(myApp.allocated());
  myAllocated.emplace_back(myApp.theId);
  ASSERT_EQ(4, myApp.edgeNode());
  ASSERT_EQ(4, myApp.thePath.size());
  ASSERT_FLOAT_EQ(7, myNetwork->totalCapacity());
  ASSERT_FLOAT_EQ(7.4, myNetwork->totProcessing());
  ASSERT_FLOAT_EQ(5, myNetwork->totNetRate());

  // remove all apps
  for (const auto& myId : myAllocated) {
    myNetwork->del(myId);
  }
  ASSERT_FLOAT_EQ(19, myNetwork->totalCapacity());
  ASSERT_FLOAT_EQ(15, myNetwork->totProcessing());
  ASSERT_FLOAT_EQ(0, myNetwork->totNetRate());

  ASSERT_EQ(0, myNetwork->signalling());
}

TEST_F(TestMecQkdOnlineNetwork, test_policy014_k_1_another) {
  auto myNetwork = makeNetwork(MecQkdOnlineAlgo::Policy014_k1, true);

  // allocate one app on the closest node
  MecQkdOnlineNetwork::Allocation myApp(0, 2, 1);
  myNetwork->add(myApp);
  ASSERT_TRUE(myApp.allocated());
  ASSERT_EQ(11, myApp.edgeNode());
  ASSERT_EQ(2, myApp.thePath.size());

  // no more room for other apps
  myNetwork->add(myApp);
  ASSERT_TRUE(not myApp.allocated());

  ASSERT_EQ(0, myNetwork->signalling());
}

TEST_F(TestMecQkdOnlineNetwork, test_policy014_k_3) {
  auto myNetwork = makeNetwork(MecQkdOnlineAlgo::Policy014_k3, true);
  ASSERT_FLOAT_EQ(74, myNetwork->totalCapacity());
  ASSERT_FLOAT_EQ(8, myNetwork->totProcessing());
  ASSERT_FLOAT_EQ(0, myNetwork->totNetRate());

  std::vector<uint64_t> myAllocated;

  // rate is too high
  MecQkdOnlineNetwork::Allocation myApp(0, 8.1, 1);
  myNetwork->add(myApp);
  ASSERT_TRUE(not myApp.allocated());

  // load is too high
  myApp.theRate = 1;
  myApp.theLoad = 8.1;
  myNetwork->add(myApp);
  ASSERT_TRUE(not myApp.allocated());

  // allocate one app on the closest node
  myApp.theRate = 2;
  myApp.theLoad = 1;
  myNetwork->add(myApp);
  ASSERT_TRUE(myApp.allocated());
  myAllocated.emplace_back(myApp.theId);
  ASSERT_EQ(11, myApp.edgeNode());
  ASSERT_EQ(2, myApp.thePath.size());
  ASSERT_FLOAT_EQ(70, myNetwork->totalCapacity());
  ASSERT_FLOAT_EQ(7, myNetwork->totProcessing());
  ASSERT_FLOAT_EQ(2, myNetwork->totNetRate());

  // allocate two more apps on secondary paths
  myNetwork->add(myApp);
  ASSERT_TRUE(myApp.allocated());
  myAllocated.emplace_back(myApp.theId);
  ASSERT_EQ(11, myApp.edgeNode());
  ASSERT_EQ(3, myApp.thePath.size());
  ASSERT_EQ(5, myNetwork->signalling());

  myNetwork->add(myApp);
  ASSERT_TRUE(myApp.allocated());
  myAllocated.emplace_back(myApp.theId);
  ASSERT_EQ(11, myApp.edgeNode());
  ASSERT_EQ(4, myApp.thePath.size());
  ASSERT_EQ(5 + 7, myNetwork->signalling());

  ASSERT_FLOAT_EQ(50, myNetwork->totalCapacity());
  ASSERT_FLOAT_EQ(5, myNetwork->totProcessing());
  ASSERT_FLOAT_EQ(6, myNetwork->totNetRate());
}

TEST_F(TestMecQkdOnlineNetwork, test_policy015) {
  auto myNetwork = makeNetwork(MecQkdOnlineAlgo::Policy015, false);
  ASSERT_FLOAT_EQ(19, myNetwork->totalCapacity());
  ASSERT_FLOAT_EQ(15, myNetwork->totProcessing());
  ASSERT_FLOAT_EQ(0, myNetwork->totNetRate());

  // rate is too high
  MecQkdOnlineNetwork::Allocation myApp(0, 8.1, 1);
  myNetwork->add(myApp);
  ASSERT_TRUE(not myApp.allocated());

  // load is too high
  myApp.theRate = 1;
  myApp.theLoad = 8.1;
  myNetwork->add(myApp);
  ASSERT_TRUE(not myApp.allocated());

  // allocate one app on the closest node
  myApp.theRate = 2;
  myApp.theLoad = 1;
  myNetwork->add(myApp);
  ASSERT_TRUE(myApp.allocated());
  ASSERT_EQ(10, myApp.edgeNode());
  ASSERT_EQ(1, myApp.thePath.size());
  ASSERT_FLOAT_EQ(17, myNetwork->totalCapacity());
  ASSERT_FLOAT_EQ(14, myNetwork->totProcessing());
  ASSERT_FLOAT_EQ(2, myNetwork->totNetRate());

  // allocate another one on the farthest node
  myApp.theRate = 0.5;
  myApp.theLoad = 5;
  myNetwork->add(myApp);
  ASSERT_TRUE(myApp.allocated());
  ASSERT_EQ(4, myApp.edgeNode());
  ASSERT_EQ(4, myApp.thePath.size());
  ASSERT_FLOAT_EQ(15, myNetwork->totalCapacity());
  ASSERT_FLOAT_EQ(9, myNetwork->totProcessing());
  ASSERT_FLOAT_EQ(2.5, myNetwork->totNetRate());

  // allocate another one on the same node
  myApp.theRate = 0.5;
  myApp.theLoad = 3;
  myNetwork->add(myApp);
  ASSERT_TRUE(myApp.allocated());
  ASSERT_EQ(4, myApp.edgeNode());
  ASSERT_EQ(4, myApp.thePath.size());
  ASSERT_FLOAT_EQ(13, myNetwork->totalCapacity());
  ASSERT_FLOAT_EQ(6, myNetwork->totProcessing());
  ASSERT_FLOAT_EQ(3, myNetwork->totNetRate());

  // allocate one on the two-hop path
  myApp.theRate = 1;
  myApp.theLoad = 1;
  myNetwork->add(myApp);
  ASSERT_TRUE(myApp.allocated());
  ASSERT_EQ(9, myApp.edgeNode());
  ASSERT_EQ(2, myApp.thePath.size());
  ASSERT_FLOAT_EQ(11, myNetwork->totalCapacity());
  ASSERT_FLOAT_EQ(5, myNetwork->totProcessing());
  ASSERT_FLOAT_EQ(4, myNetwork->totNetRate());
}

TEST_F(TestMecQkdOnlineNetwork, test_policy015_alt) {
  auto myNetwork = makeNetwork(MecQkdOnlineAlgo::Policy015, false);
  ASSERT_FLOAT_EQ(19, myNetwork->totalCapacity());
  ASSERT_FLOAT_EQ(15, myNetwork->totProcessing());
  ASSERT_FLOAT_EQ(0, myNetwork->totNetRate());

  // allocate one app on the three-hop node
  MecQkdOnlineNetwork::Allocation myApp(0, 1.1, 3);
  myNetwork->add(myApp);
  ASSERT_TRUE(myApp.allocated());
  ASSERT_EQ(7, myApp.edgeNode());

  // subsequent allocations can be done on one of two nodes
  std::set<unsigned long> myEdgeNodes;
  myApp.theRate = 0.001;
  myApp.theLoad = 0.1;
  for (size_t i = 0; i < 10; i++) {
    myNetwork->add(myApp);
    ASSERT_TRUE(myApp.allocated());
    myEdgeNodes.emplace(myApp.edgeNode());
    myNetwork->del(myApp.theId);
  }
  ASSERT_EQ(std::set<unsigned long>({7, 10}), myEdgeNodes);
}

TEST_F(TestMecQkdOnlineNetwork, test_policy015_reuse) {
  auto myNetwork = makeNetwork(MecQkdOnlineAlgo::Policy015_reuse, false);
  ASSERT_FLOAT_EQ(19, myNetwork->totalCapacity());
  ASSERT_FLOAT_EQ(15, myNetwork->totProcessing());
  ASSERT_FLOAT_EQ(0, myNetwork->totNetRate());

  // allocate one app on the three-hop node
  MecQkdOnlineNetwork::Allocation myApp(0, 1.1, 3);
  myNetwork->add(myApp);
  ASSERT_TRUE(myApp.allocated());
  ASSERT_EQ(7, myApp.edgeNode());

  // all subsequent allocations are done on the same node
  myApp.theRate = 0.001;
  myApp.theLoad = 0.1;
  for (size_t i = 0; i < 10; i++) {
    myNetwork->add(myApp);
    ASSERT_TRUE(myApp.allocated());
    ASSERT_EQ(7, myApp.edgeNode()) << i;
    myNetwork->del(myApp.theId);
  }
}

} // namespace qr
} // namespace uiiit
