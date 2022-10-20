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
#include "QuantumRouting/peerassignment.h"

#include "gtest/gtest.h"
#include <glog/logging.h>

#include <ctime>

namespace uiiit {
namespace qr {

struct TestPeerAssignment : public ::testing::Test {
  /*
  0 -------+   +-----> 4      0-3: 10     3-4: 10
           |  / +----> 5      1-3: 20     3-5: 20
           v / /              2-3: 40     3-6: 20
  1 -----> 3 --------> 6                  3-7: 40   7-8: 40
           ^ \
           |  \
  2 -------+   +-----> 7 ------> 8
  */
  CapacityNetwork::WeightVector exampleEdgeWeights() {
    return CapacityNetwork::WeightVector({
        {0, 3, 10},
        {1, 3, 20},
        {2, 3, 40},
        {3, 4, 10},
        {3, 5, 20},
        {3, 6, 20},
        {3, 7, 40},
        {7, 8, 40},
    });
  }
};

TEST_F(TestPeerAssignment, test_algorithms) {
  for (const auto& myAlgo : allPeerAssignmentAlgos()) {
    ASSERT_EQ(myAlgo, peerAssignmentAlgofromString(toString(myAlgo)));
  }
  ASSERT_EQ("unknown", toString(static_cast<PeerAssignmentAlgo>(999)));
  ASSERT_THROW(peerAssignmentAlgofromString("not-existing-algo"),
               std::runtime_error);
}

TEST_F(TestPeerAssignment, test_random) {
  using Nodes                 = std::vector<unsigned long>;
  const auto         NUM_RUNS = 100;
  support::UniformRv myRv(0, 1, std::time(nullptr), 0, 0);

  CapacityNetwork myNetwork(exampleEdgeWeights());

  auto myAssignment =
      makePeerAssignment(myNetwork, PeerAssignmentAlgo::Random, myRv);

  const std::vector<PeerAssignment::AppDescriptor> myApps({
      {0, 1, 0.5},
      {1, 1, 0.5},
      {2, 1, 0.5},
  });

  const Nodes myDataCenters({4, 5, 6, 8});

  const std::vector<unsigned long> W({1, 2}); // num data centers per app

  for (const auto w : W) {
    std::map<unsigned long, std::set<unsigned long>> myAllAssigned;
    for (auto i = 0u; i < NUM_RUNS; i++) {
      const auto myAssigned = myAssignment->assign(myApps, w, myDataCenters);
      ASSERT_EQ(myApps.size(), myAssigned.size());
      for (const auto& elem : myAssigned) {
        ASSERT_EQ(w, elem.thePeers.size());
        for (const auto myPeer : elem.thePeers) {
          myAllAssigned[elem.theHost].insert(myPeer);
        }
      }
    }
    for (const auto& elem : myAllAssigned) {
      for (const auto& node : myDataCenters) {
        EXPECT_EQ(1u, elem.second.count(node))
            << "w " << w << ", app " << elem.first << ", node " << node;
      }
    }
  }
}

} // namespace qr
} // namespace uiiit
