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
#include "Support/tostring.h"

#include "Details/examplenetwork.h"

#include <boost/filesystem.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include <boost/graph/detail/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_selectors.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include <boost/property_map/property_map.hpp>

#include <cmath>
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

  struct VertexData {
    std::string theLabel;
    double      theLongitude;
    double      theLatitude;
  };

  struct EdgeData {
    double theLinkSpeed;
    int    theInvLinkSpeed;
  };

  using Graph = boost::adjacency_list<boost::listS,
                                      boost::vecS,
                                      boost::undirectedS,
                                      VertexData,
                                      EdgeData,
                                      boost::no_property,
                                      boost::listS>;

  Graph myGraph;

  boost::dynamic_properties myDp(boost::ignore_other_properties);
  myDp.property("label", boost::get(&VertexData::theLabel, myGraph));
  myDp.property("Latitude", boost::get(&VertexData::theLatitude, myGraph));
  myDp.property("Longitude", boost::get(&VertexData::theLongitude, myGraph));
  myDp.property("LinkSpeedRaw", boost::get(&EdgeData::theLinkSpeed, myGraph));

  std::ifstream myInfile(theFilename);
  boost::read_graphml(myInfile, myGraph, myDp, 0);

  if (VLOG_IS_ON(1)) {
    boost::print_graph(myGraph, boost::get(&VertexData::theLabel, myGraph));
  }

  EXPECT_EQ(61, boost::num_vertices(myGraph));
  EXPECT_EQ(89, boost::num_edges(myGraph));

  const auto  myVertices   = boost::vertices(myGraph);
  const auto& myLabels     = boost::get(&VertexData::theLabel, myGraph);
  const auto& myLatitudes  = boost::get(&VertexData::theLatitude, myGraph);
  const auto& myLongitudes = boost::get(&VertexData::theLongitude, myGraph);
  for (auto it = myVertices.first; it != myVertices.second; ++it) {
    LOG(INFO) << "vertex " << *it << " [" << myLabels[*it] << "] ("
              << myLatitudes[*it] << "," << myLongitudes[*it]
              << "), out degree " << boost::out_degree(*it, myGraph);
  }

  const auto myEdges = boost::edges(myGraph);
  for (auto it = myEdges.first; it != myEdges.second; ++it) {
    const auto myLinkSpeed = boost::get(&EdgeData::theLinkSpeed, myGraph, *it);
    const auto myInvLinkSpeed =
        ::isnormal(myLinkSpeed) ? static_cast<int>(1e12 / myLinkSpeed) : 0;
    ASSERT_GE(myLinkSpeed, 0.0);
    boost::put(&EdgeData::theInvLinkSpeed, myGraph, *it, myInvLinkSpeed);
    LOG(INFO) << "edge " << *it << ", link speed "
              << static_cast<size_t>(myLinkSpeed / 1e6) << " Mb/s"
              << ", inv link speed "
              << boost::get(&EdgeData::theInvLinkSpeed, myGraph, *it)
              << " ps/b";
  }

  using VertexDescriptor = boost::graph_traits<Graph>::vertex_descriptor;
  const VertexDescriptor        mySource = 0;
  std::vector<int>              myDistances(boost::num_vertices(myGraph));
  std::vector<VertexDescriptor> myPredecessors(boost::num_vertices(myGraph));
  boost::dijkstra_shortest_paths(
      myGraph,
      mySource,
      boost::predecessor_map(myPredecessors.data())
          .distance_map(boost::make_iterator_property_map(
              myDistances.data(), get(boost::vertex_index, myGraph)))
          .weight_map(boost::get(&EdgeData::theInvLinkSpeed, myGraph)));

  struct HopsFinder {
    HopsFinder(const std::vector<VertexDescriptor>& aPredecessors,
               const VertexDescriptor               aSource)
        : thePredecessors(aPredecessors)
        , theSource(aSource) {
      // noop
    }
    void operator()(std::list<VertexDescriptor>& aHops,
                    const VertexDescriptor       aNext) {
      if (aNext == theSource) {
        return;
      }
      aHops.push_front(aNext);
      (*this)(aHops, thePredecessors[aNext]);
    }
    const std::vector<VertexDescriptor>& thePredecessors;
    const VertexDescriptor               theSource;
  };

  HopsFinder myHopsFinder(myPredecessors, mySource);
  for (auto it = myVertices.first; it != myVertices.second; ++it) {
    std::list<VertexDescriptor> myHops;
    myHopsFinder(myHops, *it);
    LOG(INFO) << "from " << myLabels[mySource] << " to " << myLabels[*it]
              << ", distance " << myDistances[*it] << ", via "
              << toString(myHops, ",", [&myLabels](const auto& aValue) {
                   return myLabels[aValue];
                 });
  }
}

} // namespace qr
} // namespace uiiit
