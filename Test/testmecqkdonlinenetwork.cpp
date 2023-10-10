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
  TestMecQkdOnlineNetwork() {
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

  static std::unique_ptr<MecQkdOnlineNetwork>
  makeNetwork(const MecQkdOnlineAlgo aAlgo) {
    auto ret = std::make_unique<MecQkdOnlineNetwork>(exampleEdgeWeights());
    ret->configure(
        aAlgo,
        std::make_unique<support::DeterministicRv<std::vector<double>>>(
            std::vector<double>({0.1, 0.4, 0.25, 0.9, 0.6, 0.7, 0.7, 0.3})),
        {0},
        {{3, 5}, {4, 2}, {6, 10}, {7, 1}});
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
  auto myNetwork = makeNetwork(MecQkdOnlineAlgo::Policy014_k1);
}

} // namespace qr
} // namespace uiiit
