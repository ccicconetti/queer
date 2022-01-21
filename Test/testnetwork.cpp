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

#include "QuantumRouting/network.h"

#include "Details/examplenetwork.h"

#include <boost/filesystem.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/detail/adjacency_list.hpp>
#include <boost/graph/graph_selectors.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include <boost/property_map/property_map.hpp>

#include <glog/logging.h>

#include <fstream>

#include "gtest/gtest.h"

namespace uiiit {
namespace qr {

struct TestNetwork : public ::testing::Test {
  TestNetwork()
      : theFilename(
            (boost::filesystem::current_path() / "example.gml").string()) {
    // noop
  }

  void createExampleNetworkFile() const {
    std::ofstream myOutfile(theFilename);
    myOutfile << exampleNetwork();
  }

  void SetUp() {
    boost::filesystem::remove(theFilename);
  }

  void TearDown() {
    if (not VLOG_IS_ON(1)) {
      boost::filesystem::remove(theFilename);
    }
  }

  const std::string theFilename;
};

TEST_F(TestNetwork, test_ctor) {
  ASSERT_NO_THROW(Network());
}

TEST_F(TestNetwork, test_read_graphml) {
  createExampleNetworkFile();

  struct GraphData {
    std::string theLabel;
    double      theLongitude;
    double      theLatitude;
  };

  using Graph =
      boost::adjacency_list<boost::listS,
                            boost::vecS,
                            boost::undirectedS,
                            GraphData,
                            boost::property<boost::edge_weight_t, float>,
                            boost::no_property,
                            boost::listS>;

  Graph myGraph;

  boost::dynamic_properties myDp(boost::ignore_other_properties);
  myDp.property("label", boost::get(&GraphData::theLabel, myGraph));
  myDp.property("Latitude", boost::get(&GraphData::theLatitude, myGraph));
  myDp.property("Longitude", boost::get(&GraphData::theLongitude, myGraph));

  std::ifstream myInfile(theFilename);
  boost::read_graphml(myInfile, myGraph, myDp, 0);

  EXPECT_EQ(61, boost::num_vertices(myGraph));
  EXPECT_EQ(89, boost::num_edges(myGraph));

  const auto myVertices = boost::vertices(myGraph);
  for (auto it = myVertices.first; it != myVertices.second; ++it) {
    LOG(INFO) << "vertex " << *it << "["
              << boost::get(&GraphData::theLabel, myGraph, *it) << "]"
              << ", out degree " << boost::out_degree(*it, myGraph);
  }

  boost::print_graph(myGraph, boost::get(&GraphData::theLabel, myGraph));
}

} // namespace qr
} // namespace uiiit
