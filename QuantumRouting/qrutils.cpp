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
#include "Support/random.h"

#include <glog/logging.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/graphml.hpp>

#include <cassert>
#include <cmath>
#include <memory>
#include <set>

namespace uiiit {
namespace qr {

double distance(const Coordinate& aLhs, const Coordinate& aRhs) {
  return std::sqrt(std::pow(std::get<0>(aLhs) - std::get<0>(aRhs), 2.0) +
                   std::pow(std::get<1>(aLhs) - std::get<1>(aRhs), 2.0) +
                   std::pow(std::get<2>(aLhs) - std::get<2>(aRhs), 2.0));
}

std::vector<std::pair<unsigned long, unsigned long>>
findLinks(const std::vector<Coordinate>& aItems,
          const double                   aThreshold,
          const double                   aProbability,
          const unsigned long            aSeed) {
  assert(aProbability >= 0 and aProbability <= 1);
  assert(aThreshold >= 0);

  const auto myRv =
      aProbability < 1 ?
          std::make_unique<support::UniformRv>(0, 1, aSeed, 0, 0) :
          nullptr;

  std::vector<std::pair<unsigned long, unsigned long>> ret;
  for (unsigned long i = 0; i < aItems.size(); i++) {
    for (unsigned long j = 0; j < i; j++) {
      if (distance(aItems[i], aItems[j]) < aThreshold and
          (myRv.get() == nullptr or (*myRv)() < aProbability)) {
        ret.push_back({i, j});
      }
    }
  }
  return ret;
}

std::vector<std::pair<unsigned long, unsigned long>>
findLinks(std::istream& aGraphMl) {
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

  boost::read_graphml(aGraphMl, myGraph, myDp, 0);

  if (VLOG_IS_ON(2)) {
    boost::print_graph(myGraph, boost::get(&VertexData::theLabel, myGraph));
  }

  std::vector<std::pair<unsigned long, unsigned long>> ret;
  const auto            myEdges = boost::edges(myGraph);
  std::set<std::string> myFound;
  for (auto it = myEdges.first; it != myEdges.second; ++it) {
    if (myFound
            .emplace(std::to_string(it->m_source) +
                     std::to_string(it->m_target))
            .second) {
      ret.push_back({it->m_source, it->m_target});
    }
  }
  return ret;
}

bool bigraphConnected(
    const std::vector<std::pair<unsigned long, unsigned long>>& aEdges) {

  // count the number of connected vertices
  using Graph =
      boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS>;
  using Vertex = boost::graph_traits<Graph>::vertex_descriptor;

  struct CountVisitor final : public boost::default_dfs_visitor {
    explicit CountVisitor(unsigned long& aCount)
        : theCount(aCount) {
      // noop
    }
    void discover_vertex([[maybe_unused]] Vertex       aVertex,
                         [[maybe_unused]] const Graph& aGraph) const {
      theCount++;
    }
    unsigned long& theCount;
  };

  Graph myGraph;
  for (const auto& myEdge : aEdges) {
    boost::add_edge(myEdge.first, myEdge.second, myGraph);
  }
  if (boost::num_vertices(myGraph) == 0) {
    return true;
  }

  unsigned long                          myCount = 0;
  CountVisitor                           myVisitor(myCount);
  std::vector<boost::default_color_type> myColorMap(
      boost::num_vertices(myGraph));
  boost::depth_first_visit(myGraph,
                           *boost::vertices(myGraph).first,
                           myVisitor,
                           boost::make_iterator_property_map(
                               myColorMap.begin(),
                               boost::get(boost::vertex_index, myGraph),
                               myColorMap[0]));

  // check if the graph is connected
  return myCount == boost::num_vertices(myGraph);
}

double fidelitySwapping(const double        p1,
                        const double        p2,
                        const double        eta,
                        const unsigned long L,
                        const double        F) {
  assert(p1 > 0 and p1 <= 1);
  assert(p2 > 0 and p2 <= 1);
  assert(eta >= 0.5 and eta <= 1);
  assert(F >= 0 and F <= 1);

  if (L == 0) {
    return F;
  }

  return 1.0 / 4.0 +
         3.0 / 4.0 * std::pow(p1 * p1 * p2 * (4 * eta * eta - 1) / 3.0, L - 1) *
             std::pow((4.0 * F - 1.0) / 3.0, L);
}

} // namespace qr
} // namespace uiiit