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

#include "QuantumRouting/qrutils.h"

#include "Details/examplenetwork.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <glog/logging.h>

#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>

namespace uiiit {
namespace qr {

struct TestQrUtils : public ::testing::Test {};

TEST_F(TestQrUtils, test_distance) {
  ASSERT_FLOAT_EQ(5, distance({0, 0, 0}, {3, 4, 0}));
  ASSERT_FLOAT_EQ(5, distance({0 - 5, 0 - 5, 0}, {3 - 5, 4 - 5, 0}));
  ASSERT_FLOAT_EQ(::sqrt(3.0), distance({1, 2, 3}, {2, 3, 4}));
}

TEST_F(TestQrUtils, test_find_links) {
  // put items on the four cornes of a square
  std::vector<Coordinate> myItems({
      {0, 0, 0},
      {0, 1, 0},
      {1, 0, 0},
      {1, 1, 0},
  });

  // items are too far
  ASSERT_EQ(0, findLinks(myItems, 0.5).size());

  // only edges are found
  const auto myEdges = findLinks(myItems, 1.1);
  ASSERT_EQ(4, myEdges.size());
  ASSERT_EQ((std::vector<std::pair<unsigned long, unsigned long>>({
                {1, 0},
                {2, 0},
                {3, 1},
                {3, 2},
            })),
            myEdges);

  // both edges and diagonals are found
  ASSERT_EQ(6, findLinks(myItems, 1.5).size());

  // random link generation
  ASSERT_EQ(0, findLinks(myItems, 1.5, 0, 0).size());
  ASSERT_EQ(4, findLinks(myItems, 1.5, 0.8, 0).size());
}

TEST_F(TestQrUtils, test_find_links_graphml) {
  std::stringstream myStream;
  myStream << exampleNetwork();

  std::vector<Coordinate> myCoordinates;
  const auto              myLinks = findLinks(myStream, myCoordinates);
  ASSERT_EQ(75, myLinks.size());
  ASSERT_EQ(61, myCoordinates.size());

  unsigned long myMin = std::numeric_limits<unsigned long>::max();
  unsigned long myMax = 0;
  for (const auto& myLink : myLinks) {
    VLOG(2) << "(" << myLink.first << "," << myLink.second << ")";
    myMin = std::min(myMin, myLink.first);
    myMin = std::min(myMin, myLink.second);
    myMax = std::max(myMax, myLink.first);
    myMax = std::max(myMax, myLink.second);
  }
  EXPECT_EQ(0, myMin);
  EXPECT_EQ(60, myMax);

  std::set<std::string> myFound;
  for (const auto& myLink : myLinks) {
    ASSERT_TRUE(myFound
                    .emplace(std::to_string(myLink.first) + "-" +
                             std::to_string(myLink.second))
                    .second)
        << "(" << myLink.first << "," << myLink.second << ")";
  }
}

TEST_F(TestQrUtils, test_bigraph_connected) {
  using Graph = std::vector<std::pair<unsigned long, unsigned long>>;

  // empty graph
  ASSERT_TRUE(bigraphConnected(Graph()));

  // connected graph
  Graph myEdges({
      {0, 1},
      {1, 2},
      {2, 3},
      {0, 3},
  });

  ASSERT_TRUE(bigraphConnected(myEdges));

  // add a disconnected clique of nodes
  myEdges.push_back({4, 5});
  myEdges.push_back({6, 5});
  myEdges.push_back({4, 6});

  ASSERT_FALSE(bigraphConnected(myEdges));

  // connect the clique to the other sub-graph
  myEdges.push_back({3, 4});
  ASSERT_TRUE(bigraphConnected(myEdges));
}

TEST_F(TestQrUtils, test_fidelity_swapping) {
  ASSERT_FLOAT_EQ(0.9925, fidelitySwapping(1, 1, 1, 0, 0.9925));
  ASSERT_FLOAT_EQ(0.985075, fidelitySwapping(1, 1, 1, 2, 0.9925));
  ASSERT_FLOAT_EQ(0.9704470075, fidelitySwapping(1, 1, 1, 4, 0.9925));
  ASSERT_FLOAT_EQ(0.92314349037037, fidelitySwapping(1, 1, 1, 4, 0.98));

  //
  // computed with gnuplot using:
  //   f(x,y)=1.0/4+3.0/4*((4*y-1)/3)**x *
  //     (0.9**2 * 0.5 * (4 * 0.95 * 0.95 - 1 ) / 3.0)**(x-1)
  // and then
  //   print f(4,0.98)
  //
  ASSERT_FLOAT_EQ(0.279446282739145, fidelitySwapping(0.9, 0.5, 0.95, 4, 0.98));
}

TEST_F(TestQrUtils, DISABLED_print_fidelity_per_hops) {
  constexpr double p1  = 1.0;
  constexpr double p2  = 1.0;
  constexpr double eta = 1.0;

  const std::string myFilename("fidelity.dat");
  std::ofstream     myOut(myFilename);
  LOG(INFO) << "saving to " << myFilename;
  for (std::size_t L = 1; L < 10; L++) {
    myOut << L;
    for (double myFidelityInit = 0.95; myFidelityInit < 0.991;
         myFidelityInit += 0.01) {
      myOut << ' ' << fidelitySwapping(p1, p2, eta, L, myFidelityInit);
    }
    myOut << '\n';
  }
  LOG(INFO) << "using Gnuplot plot with:\n"
               "plot 'fidelity.dat' u 1:2 w lp title \"F=0.95\", '' u 1:3 w lp "
               "title \"F=0.96\", '' u 1:4 w lp title \"F=0.97\", '' u 1:5 w "
               "lp title \"F=0.98\", '' u 1:6 w lp title \"F=0.99\"";
}

} // namespace qr
} // namespace uiiit
