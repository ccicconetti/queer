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

struct TestMecQkdWorkload : public ::testing::Test {
  using AppInfo = MecQkdWorkload::AppInfo;

  TestMecQkdWorkload()
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

TEST_F(TestMecQkdWorkload, test_mec_qkd_workload) {
  std::vector<AppInfo> myAppInfo;
  ASSERT_THROW(MecQkdWorkload(myAppInfo, theRv, true), std::runtime_error);

  myAppInfo.emplace_back(AppInfo{0, 1, 0.5, 10});

  {
    MecQkdWorkload myGenerator(myAppInfo, theRv, true);
    for (auto i = 0; i < 10; i++) {
      const auto myInfo = myGenerator();
      ASSERT_EQ(0, myInfo.theRegion);
      ASSERT_EQ(1, myInfo.theWeight);
      ASSERT_EQ(0.5, myInfo.theLoad);
      ASSERT_EQ(10, myInfo.theRate);
    }
    ASSERT_EQ(std::set<unsigned long>({0}), myGenerator.regions());
  }

  myAppInfo.emplace_back(AppInfo{1, 2, 0.5, 10});

  {
    std::vector<unsigned long> myRegions;
    MecQkdWorkload             myGenerator(myAppInfo, theRv, true);
    for (auto i = 0; i < 10; i++) {
      myRegions.emplace_back(myGenerator().theRegion);
    }
    ASSERT_EQ(std::vector<unsigned long>({0, 1, 0, 1, 1, 1, 1, 0, 0, 1}),
              myRegions);
    ASSERT_EQ(std::set<unsigned long>({0, 1}), myGenerator.regions());
  }
}

TEST_F(TestMecQkdWorkload, test_mec_qkd_workload_from_file) {
  {
    std::ofstream myCsvFile("removeme.csv");
    myCsvFile << "0,1,0.02499083669321249,3\n"
                 "1,1,0.028765301871538224,5\n"
                 "2,1,0.03630966699496657,2\n"
                 "# comment line\n"
                 "3,1,0.027700677558573084,7\n"
                 "4,1,0.9177905729976219,5\n";
  }

  auto myGenerator = MecQkdWorkload::fromCsvFile("removeme.csv", theRv, true);

  boost::filesystem::remove(boost::filesystem::current_path() / "removeme.csv");

  std::set<std::string> myAllValues;
  for (auto i = 0; i < 10; i++) {
    myAllValues.emplace(myGenerator().toString());
  }
  ASSERT_EQ(5, myAllValues.size());
  ASSERT_EQ(std::set<unsigned long>({0, 1, 2, 3, 4}), myGenerator.regions());
}

} // namespace qr
} // namespace uiiit
