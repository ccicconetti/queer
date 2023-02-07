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
#include "QuantumRouting/peerassignment.h"
#include "Support/tostring.h"

#include "gtest/gtest.h"
#include <glog/logging.h>

#include <ctime>

namespace uiiit {
namespace qr {

struct TestPeerAssignment : public ::testing::Test {
  /*
  end users: 0, 1, 2
  data centers: 4, 5, 6, 8
  network-only nodes: 3, 7

                                 ┌─────┐
                                 │     │
                      ┌───10────►│  4  │
                      │          │     │
  ┌─────┐             │          └─────┘
  │     │             │
  │  0  ├────10─────┐ │          ┌─────┐
  │     │           │ │          │     │
  └─────┘           │ │     ┌───►│  5  │
                    ▼ │     │    │     │
  ┌─────┐         ┌───┴─┐   20   └─────┘
  │     │         │     ├────
  │  1  ├───20───►│  3  │             ┌─────┐
  │     │         │     ├───┐         │     │
  └─────┘         └───┬─┘   └──20────►│  6  │
                    ▲ │               │     │
  ┌─────┐           │ │               └─────┘
  │     │           │ │
  │  2  ├────40─────┘ │          ┌─────┐       ┌─────┐
  │     │             │          │     │       │     │
  └─────┘             └───40────►│  7  ├──40──►│  8  │
                                 │     │       │     │
                                 └─────┘       └─────┘
  */
  EsNetwork::WeightVector exampleEdgeWeights() {
    return EsNetwork::WeightVector({
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

  TestPeerAssignment()
      : theRv(0, 1, std::time(nullptr), 0, 0)
      , theApps({
            {0, 1, 0.5},
            {1, 1, 0.5},
            {2, 1, 0.5},
        })
      , theDataCenters({4, 5, 6, 8}) {
    // noop
  }

  using Nodes = std::vector<unsigned long>;

  const unsigned long NUM_RUNS = 100;

  support::UniformRv                               theRv;
  const std::vector<PeerAssignment::AppDescriptor> theApps;
  const Nodes                                      theDataCenters;
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
  EsNetwork myNetwork(exampleEdgeWeights());

  auto myAssignment =
      makePeerAssignment(myNetwork,
                         PeerAssignmentAlgo::Random,
                         theRv,
                         [](const auto&, const auto&) { return true; });

  const std::vector<unsigned long> W({1, 2}); // num data centers per app

  for (const auto w : W) {
    std::map<unsigned long, std::set<unsigned long>> myAllAssigned;
    for (auto i = 0u; i < NUM_RUNS; i++) {
      const auto myAssigned = myAssignment->assign(theApps, w, theDataCenters);
      ASSERT_EQ(theApps.size(), myAssigned.size());
      for (const auto& elem : myAssigned) {
        ASSERT_EQ(w, elem.thePeers.size());
        for (const auto myPeer : elem.thePeers) {
          myAllAssigned[elem.theHost].insert(myPeer);
        }
      }
    }
    for (const auto& elem : myAllAssigned) {
      for (const auto& node : theDataCenters) {
        ASSERT_EQ(1u, elem.second.count(node))
            << "w " << w << ", app " << elem.first << ", node " << node;
      }
    }
  }
}

TEST_F(TestPeerAssignment, test_shortest_path) {
  using Set = std::set<unsigned long>;
  EsNetwork myNetwork(exampleEdgeWeights());

  auto myAssignment =
      makePeerAssignment(myNetwork,
                         PeerAssignmentAlgo::ShortestPath,
                         theRv,
                         [](const auto&, const auto&) { return true; });

  // W = 1: all apps pick the same node
  auto myAssigned = myAssignment->assign(theApps, 1, theDataCenters);
  ASSERT_EQ(theApps.size(), myAssigned.size());
  for (const auto& elem : myAssigned) {
    ASSERT_EQ(1u, elem.thePeers.size());
    ASSERT_TRUE(Set({4, 5, 6}).count(elem.thePeers[0]) > 0);
  }

  // W = 2: all apps pick the same two nodes
  myAssigned = myAssignment->assign(theApps, 2, theDataCenters);
  ASSERT_EQ(theApps.size(), myAssigned.size());
  for (const auto& elem : myAssigned) {
    ASSERT_EQ(2u, elem.thePeers.size());
    ASSERT_TRUE(Set({4, 5, 6}).count(elem.thePeers[0]) > 0);
    ASSERT_TRUE(Set({4, 5, 6}).count(elem.thePeers[1]) > 0);
  }

  // W >= 4: all apps pick all nodes
  myAssigned = myAssignment->assign(theApps, 4, theDataCenters);
  ASSERT_EQ(theApps.size(), myAssigned.size());
  for (const auto& elem : myAssigned) {
    ASSERT_EQ(Set({4, 5, 6, 8}),
              Set(elem.thePeers.begin(), elem.thePeers.end()));
  }
  myAssigned = myAssignment->assign(theApps, 99, theDataCenters);
  ASSERT_EQ(theApps.size(), myAssigned.size());
  for (const auto& elem : myAssigned) {
    ASSERT_EQ(Set({4, 5, 6, 8}),
              Set(elem.thePeers.begin(), elem.thePeers.end()));
  }
}

TEST_F(TestPeerAssignment, test_load_balancing) {
  EsNetwork myNetwork(exampleEdgeWeights());
  myNetwork.measurementProbability(0.5);

  auto myAssignment =
      makePeerAssignment(myNetwork,
                         PeerAssignmentAlgo::LoadBalancing,
                         theRv,
                         [](const auto&, const auto&) { return true; });

  struct Checker {
    bool
    operator()(const std::vector<EsNetwork::AppDescriptor>&   aAssigned,
               const std::vector<std::vector<unsigned long>>& aExpected) const {
      auto ret = true;
      if (aAssigned.size() != aExpected.size()) {
        ret = false;
      } else {
        for (unsigned long a = 0; a < aAssigned.size(); a++) {
          if (aAssigned[a].thePeers != aExpected[a]) {
            ret = false;
            break;
          }
        }
      }
      if (ret == false) {
        LOG(INFO) << "assigned:";
        for (unsigned long a = 0; a < aAssigned.size(); a++) {
          LOG(INFO) << "\t#" << a << '\t' << aAssigned[a].toString();
        }
        LOG(INFO) << "expected:";
        for (unsigned long e = 0; e < aExpected.size(); e++) {
          LOG(INFO) << "\t#" << e << '\t' << ::toStringStd(aExpected[e], ",");
        }
      }
      return ret;
    }
  };
  Checker myChecker;

  // W = 1
  auto myAssigned = myAssignment->assign(theApps, 1, theDataCenters);
  ASSERT_TRUE(myChecker(
      myAssigned, std::vector<std::vector<unsigned long>>({{4}, {5}, {6}})));

  // W = 2
  myAssigned = myAssignment->assign(theApps, 2, theDataCenters);
  ASSERT_TRUE(myChecker(
      myAssigned,
      std::vector<std::vector<unsigned long>>({{4, 6}, {5, 6}, {5, 8}})));

  // W = 3
  myAssigned = myAssignment->assign(theApps, 3, theDataCenters);
  ASSERT_TRUE(myChecker(myAssigned,
                        std::vector<std::vector<unsigned long>>(
                            {{4, 5, 6}, {5, 6, 4}, {5, 6, 8}})));

  // W = 4: each app is assigned to all data centers
  myAssigned = myAssignment->assign(theApps, 4, theDataCenters);
  ASSERT_TRUE(myChecker(myAssigned,
                        std::vector<std::vector<unsigned long>>(
                            {{4, 5, 6, 8}, {5, 6, 4, 8}, {5, 6, 8, 4}})));

  // W = 10 > number of data centers
  myAssigned = myAssignment->assign(theApps, 10, theDataCenters);
  ASSERT_TRUE(myChecker(myAssigned,
                        std::vector<std::vector<unsigned long>>(
                            {{4, 5, 6, 8}, {5, 6, 4, 8}, {5, 6, 8, 4}})));
}

} // namespace qr
} // namespace uiiit
