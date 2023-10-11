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

#include "QuantumRouting/mecqkdnetwork.h"
#include "Support/random.h"
#include "Support/tostring.h"

#include "gtest/gtest.h"

#include <boost/filesystem/operations.hpp>
#include <glog/logging.h>

#include <boost/filesystem.hpp>

#include <ctime>
#include <set>
#include <stdexcept>

namespace uiiit {
namespace qr {

struct TestMecQkdNetwork : public ::testing::Test {
  TestMecQkdNetwork()
      : theRv({0.1, 0.4, 0.25, 0.9, 0.6, 0.7, 0.7, 0.3}) {
    // noop
  }

  //   +--> 1 --> 2 --> 3 --> 4
  //  /
  // 0 --> 5 --> 6
  //  \                   all weights are 2, except 3->4 which is 1
  //   +---> 7
  static CapacityNetwork::WeightVector exampleEdgeWeights() {
    return CapacityNetwork::WeightVector({
        {0, 1, 2},
        {1, 2, 2},
        {2, 3, 2},
        {3, 4, 1},
        {0, 5, 2},
        {5, 6, 2},
        {0, 7, 2},
    });
  }

  static std::unique_ptr<MecQkdNetwork> makeNetwork() {
    auto ret = std::make_unique<MecQkdNetwork>(exampleEdgeWeights());
    ret->userNodes({0});
    ret->edgeNodes({{3, 5}, {4, 2}, {6, 10}, {7, 1}});
    return ret;
  }

  support::DeterministicRv<std::vector<double>> theRv;
};

TEST_F(TestMecQkdNetwork, test_user_edge_nodes) {
  MecQkdNetwork myNetwork(exampleEdgeWeights());

  ASSERT_NO_THROW(myNetwork.userNodes({}));
  ASSERT_NO_THROW(myNetwork.userNodes({1}));
  ASSERT_NO_THROW(myNetwork.userNodes({0, 1, 2, 3, 4}));
  ASSERT_THROW(myNetwork.userNodes({99}), std::runtime_error);

  ASSERT_NO_THROW(myNetwork.edgeNodes({}));
  ASSERT_NO_THROW(myNetwork.edgeNodes({{1, 3.14}, {4, 2.0}}));
  ASSERT_THROW(myNetwork.edgeNodes({{99, 2.0}}), std::runtime_error);
  ASSERT_THROW(myNetwork.edgeNodes({{4, -2.0}}), std::runtime_error);
}

TEST_F(TestMecQkdNetwork, test_allocation_single) {
  const std::vector<std::tuple<MecQkdAlgo, bool, unsigned long, std::size_t>>
      myExpectedAllocs({
          {MecQkdAlgo::Random, true, 6, 2},
          {MecQkdAlgo::Spf, true, 6, 2},
          {MecQkdAlgo::BestFit, true, 3, 3},
          {MecQkdAlgo::RandomBlind, true, 6, 2},
          {MecQkdAlgo::SpfBlind, false, 7, 1},
          {MecQkdAlgo::BestFitBlind, false, 4, 4},
          {MecQkdAlgo::SpfStatic, false, 4, 4},
      });

  for (const auto& myExpected : myExpectedAllocs) {
    MecQkdAlgo    myAlgo;
    bool          myAllocated;
    unsigned long myEdgeNode;
    std::size_t   myPathLength;
    std::tie(myAlgo, myAllocated, myEdgeNode, myPathLength) = myExpected;

    auto myNetwork = makeNetwork();
    ASSERT_FLOAT_EQ(18.0, myNetwork->totProcessing());
    ASSERT_FLOAT_EQ(13.0, myNetwork->totalCapacity());

    std::vector<MecQkdNetwork::Allocation> myOutput;
    myOutput.emplace_back(0, 1.5, 2.0);

    myNetwork->allocate(myOutput, myAlgo, theRv);
    ASSERT_EQ(1, myOutput.size());
    const auto& myAlloc = myOutput[0];
    ASSERT_EQ(myAllocated, myAlloc.theAllocated);
    if (myAllocated) {
      ASSERT_EQ(myEdgeNode, myAlloc.theEdgeNode);
      ASSERT_EQ(myPathLength, myAlloc.thePathLength);
      ASSERT_FLOAT_EQ(13.0 - myAlloc.totRate(), myNetwork->totalCapacity());
      ASSERT_FLOAT_EQ(18.0 - 2.0, myNetwork->totProcessing());
    }

    myOutput.front() = MecQkdNetwork::Allocation{0, 0.1, 99.0};
    myNetwork->allocate(myOutput, myAlgo, theRv);
    ASSERT_FALSE(myOutput[0].theAllocated);

    myOutput.front() = MecQkdNetwork::Allocation{0, 99.0, 0.1};
    myNetwork->allocate(myOutput, myAlgo, theRv);
    ASSERT_FALSE(myOutput[0].theAllocated);

    myOutput.front() = MecQkdNetwork::Allocation{0, 0.1, 0.1};
    myNetwork->allocate(myOutput, myAlgo, theRv);
    ASSERT_TRUE(myOutput[0].theAllocated);
  }
}

TEST_F(TestMecQkdNetwork, test_allocation_multi_spf) {
  auto myNetwork = makeNetwork();

  std::vector<MecQkdNetwork::Allocation> myOutput;
  myOutput.emplace_back(0, 1.0, 0.1);
  myOutput.emplace_back(0, 1.0, 0.1);
  myOutput.emplace_back(0, 2.0, 0.1);
  myOutput.emplace_back(0, 1.0, 0.1);
  myOutput.emplace_back(0, 1.0, 0.1);
  myOutput.emplace_back(0, 99.0, 0.1); // unfeas

  myNetwork->allocate(myOutput, MecQkdAlgo::SpfBlind, theRv);
  ASSERT_EQ(6, myOutput.size());
  ASSERT_EQ(7, myOutput[0].theEdgeNode);
  ASSERT_EQ(7, myOutput[1].theEdgeNode);
  ASSERT_EQ(6, myOutput[2].theEdgeNode);
  ASSERT_EQ(3, myOutput[3].theEdgeNode);
  ASSERT_EQ(3, myOutput[4].theEdgeNode);
  ASSERT_FALSE(myOutput[5].theAllocated);

  ASSERT_FLOAT_EQ(18.0 - 0.5, myNetwork->totProcessing());
  ASSERT_FLOAT_EQ(1.0, myNetwork->totalCapacity());
}

TEST_F(TestMecQkdNetwork, test_allocation_multi_bf) {
  auto myNetwork = makeNetwork();

  std::vector<MecQkdNetwork::Allocation> myOutput;
  myOutput.emplace_back(0, 0.1, 1.0);
  myOutput.emplace_back(0, 0.1, 1.0);
  myOutput.emplace_back(0, 0.1, 2.0);
  myOutput.emplace_back(0, 0.1, 1.0);
  myOutput.emplace_back(0, 0.1, 1.0);
  myOutput.emplace_back(0, 0.1, 99.0);

  myNetwork->allocate(myOutput, MecQkdAlgo::BestFitBlind, theRv);
  ASSERT_EQ(6, myOutput.size());

  ASSERT_FLOAT_EQ(18.0 - 6.0, myNetwork->totProcessing());
  ASSERT_FLOAT_EQ(11.5, myNetwork->totalCapacity());
  ASSERT_EQ(7, myOutput[0].theEdgeNode);
  ASSERT_EQ(4, myOutput[1].theEdgeNode);
  ASSERT_EQ(3, myOutput[2].theEdgeNode);
  ASSERT_EQ(4, myOutput[3].theEdgeNode);
  ASSERT_EQ(3, myOutput[4].theEdgeNode);
  ASSERT_FALSE(myOutput[5].theAllocated);
}

TEST_F(TestMecQkdNetwork, test_allocation_spf_static) {
  auto myNetwork = makeNetwork();

  std::vector<MecQkdNetwork::Allocation> myOutput;
  myOutput.emplace_back(0, 3.0, 0.1);
  myOutput.emplace_back(0, 0.1, 8.0);
  myOutput.emplace_back(0, 1.0, 0.1);
  myOutput.emplace_back(0, 1.0, 0.1);
  myOutput.emplace_back(0, 0.1, 0.1);

  myNetwork->allocate(myOutput, MecQkdAlgo::SpfStatic, theRv);
  ASSERT_EQ(5, myOutput.size());

  ASSERT_FLOAT_EQ(18.0 - 0.2, myNetwork->totProcessing());
  ASSERT_FLOAT_EQ(13.0 - 2.0, myNetwork->totalCapacity());
  ASSERT_FALSE(myOutput[0].theAllocated);
  ASSERT_FALSE(myOutput[1].theAllocated);
  ASSERT_EQ(7, myOutput[2].theEdgeNode);
  ASSERT_EQ(7, myOutput[3].theEdgeNode);
  ASSERT_FALSE(myOutput[4].theAllocated);
}

} // namespace qr
} // namespace uiiit
