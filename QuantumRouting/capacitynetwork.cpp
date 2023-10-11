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

#include "Support/tostring.h"

#include <algorithm>
#include <boost/graph/detail/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/properties.hpp>

#include <boost/property_map/property_map.hpp>
#include <glog/logging.h>
#include <glog/vlog_is_on.h>

#include <limits>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>

namespace uiiit {
namespace qr {

CapacityNetwork::CapacityNetwork(
    const std::vector<std::pair<unsigned long, unsigned long>>& aEdges,
    support::RealRvInterface&                                   aWeightRv,
    const bool aMakeBidirectional)
    : Network()
    , theGraph() {
  std::set<std::string> myFound;
  for (const auto& myEdge : aEdges) {
    if (not myFound
                .emplace(std::to_string(myEdge.first) + "-" +
                         std::to_string(myEdge.second))
                .second) {
      VLOG(2) << "duplicate edge found: (" << myEdge.first << ','
              << myEdge.second;
      continue;
    }
    const auto myWeight = aWeightRv();
    Utils<Graph>::addEdge(theGraph, myEdge.first, myEdge.second, myWeight);
    if (aMakeBidirectional) {
      Utils<Graph>::addEdge(theGraph, myEdge.second, myEdge.first, myWeight);
    }
  }
}

CapacityNetwork::CapacityNetwork(const WeightVector& aEdgeWeights)
    : Network()
    , theGraph() {
  for (const auto& elem : aEdgeWeights) {
    Utils<Graph>::addEdge(
        theGraph, std::get<0>(elem), std::get<1>(elem), std::get<2>(elem));
  }
}

void CapacityNetwork::toDot(const std::string& aFilename) const {
  Utils<Graph>::toDot(theGraph, aFilename);
}

CapacityNetwork::WeightVector CapacityNetwork::weights() const {
  WeightVector ret;
  const auto   myEdges   = boost::edges(theGraph);
  const auto   myWeights = boost::get(boost::edge_weight, theGraph);
  for (auto it = myEdges.first; it != myEdges.second; ++it) {
    ret.push_back({it->m_source, it->m_target, myWeights[*it]});
  }
  return ret;
}

std::map<unsigned long, std::vector<unsigned long>>
CapacityNetwork::cspf(const unsigned long            aSource,
                      const double                   aCapacity,
                      const std::set<unsigned long>& aDestinations) const {
  const auto V = boost::num_vertices(theGraph);
#ifndef NDEBUG
  assert(aSource < V);
  for (const auto& v : aDestinations) {
    assert(v < V);
  }
#endif

  // make a working copy of the graph
  auto myGraph = theGraph;

  // remove all edges without enough capacity to satisfy the requirement
  const auto myWeights = boost::get(boost::edge_weight, myGraph);
  {
    boost::graph_traits<Graph>::edge_iterator it, end, next;
    std::tie(it, end) = boost::edges(myGraph);
    for (next = it; it != end; it = next) {
      next++;
      if (myWeights[*it] < aCapacity) {
        boost::remove_edge(*it, myGraph);
      }
    }
  }

  // run Dijkstra with distance measured in hops
  std::vector<VertexDescriptor> myPredecessors(V);
  boost::dijkstra_shortest_paths(
      myGraph,
      aSource,
      boost::predecessor_map(myPredecessors.data())
          .weight_map(
              boost::make_static_property_map<Graph::edge_descriptor>(1)));

  // save the shortest paths found
  std::map<unsigned long, std::vector<unsigned long>> ret;
  for (const auto myDestination : aDestinations) {
    auto it = ret.emplace(myDestination, std::vector<unsigned long>()).first;
    if (myPredecessors[myDestination] != myDestination) {
      HopsFinder myHopsFinder(myPredecessors, aSource);
      myHopsFinder(it->second, myDestination);
    }
  }

  return ret;
}

std::size_t CapacityNetwork::numNodes() const {
  return boost::num_vertices(theGraph);
}

std::size_t CapacityNetwork::numEdges() const {
  return boost::num_edges(theGraph);
}

std::pair<std::size_t, std::size_t> CapacityNetwork::inDegree() const {
  return minMaxVertexProp(
      [](Graph::vertex_descriptor aVertex, const Graph& aGraph) {
        return boost::in_degree(aVertex, aGraph);
      });
}

std::pair<std::size_t, std::size_t> CapacityNetwork::outDegree() const {
  return minMaxVertexProp(
      [](Graph::vertex_descriptor aVertex, const Graph& aGraph) {
        return boost::out_degree(aVertex, aGraph);
      });
}

double CapacityNetwork::totalCapacity() const {
  const auto myEdges   = boost::edges(theGraph);
  const auto myWeights = boost::get(boost::edge_weight, theGraph);
  double     ret       = 0;
  for (auto it = myEdges.first; it != myEdges.second; ++it) {
    ret += myWeights[*it];
  }
  return ret;
}

std::map<unsigned long, std::set<unsigned long>>
CapacityNetwork::reachableNodes(const std::size_t aMinHops,
                                const std::size_t aMaxHops,
                                std::size_t&      aDiameter) const {
  if (aMinHops > aMaxHops) {
    throw std::runtime_error(
        "Invalid min distance (" + std::to_string(aMinHops) +
        ") larger than max distance (" + std::to_string(aMaxHops) + ")");
  }

  const auto                    V = boost::num_vertices(theGraph);
  std::vector<VertexDescriptor> myDistances(V);

  aDiameter = 0;
  std::map<unsigned long, std::set<unsigned long>> ret;
  boost::graph_traits<Graph>::vertex_iterator      it, end;
  for (std::tie(it, end) = vertices(theGraph); it != end; ++it) {
    boost::dijkstra_shortest_paths(
        theGraph,
        *it,
        boost::weight_map(
            boost::make_static_property_map<Graph::edge_descriptor>(1))
            .distance_map(boost::make_iterator_property_map(
                myDistances.data(), get(boost::vertex_index, theGraph))));

    auto myEmplaceRet = ret.emplace(*it, std::set<unsigned long>());
    assert(myEmplaceRet.second);
    for (unsigned long i = 0; i < V; i++) {
      if (*it != i and
          myDistances[i] != std::numeric_limits<unsigned long>::max()) {
        aDiameter = std::max(aDiameter, myDistances[i]);
        VLOG(2) << *it << "->" << i << ": " << myDistances[i];
        if (myDistances[i] >= aMinHops and myDistances[i] <= aMaxHops) {
          myEmplaceRet.first->second.emplace(i);
        }
      }
    }
  }
  return ret;
}

std::vector<unsigned long>
CapacityNetwork::closestNodes(const unsigned long            aSrc,
                              const unsigned long            aNum,
                              support::RealRvInterface&      aRv,
                              const std::set<unsigned long>& aWhiteList) const {
  std::vector<unsigned long> ret;
  if (aNum == 0) {
    return ret;
  }

  // find the shortest paths to reach any node from aSrc
  const auto                    V = boost::num_vertices(theGraph);
  std::vector<VertexDescriptor> myDistances(V);
  boost::dijkstra_shortest_paths(
      theGraph,
      boost::vertex(aSrc, theGraph),
      boost::weight_map(
          boost::make_static_property_map<Graph::edge_descriptor>(1))
          .distance_map(boost::make_iterator_property_map(
              myDistances.data(), get(boost::vertex_index, theGraph))));

  // sort the distances
  // - do not include the nodes that are unreachable
  // - add a tiny floating point variation to break ties
  std::list<std::pair<double, unsigned long>> myDestinations;
  for (unsigned long i = 0; i < V; i++) {
    if (aSrc != i and
        myDistances[i] != std::numeric_limits<unsigned long>::max()) {
      myDestinations.emplace_back(myDistances[i] + aRv() * .1, i);
    }
  }
  myDestinations.sort([](const auto& aLhs, const auto& aRhs) {
    return aLhs.first < aRhs.first;
  });

  // take the first aNum instances
  auto it = myDestinations.begin();
  for (unsigned long i = 0; i < aNum and it != myDestinations.end(); ++it) {
    if (aWhiteList.empty() or aWhiteList.count(it->second) > 0) {
      ret.emplace_back(it->second);
      ++i;
    }
  }
  return ret;
}

void CapacityNetwork::addCapacityToPath(
    const VertexDescriptor               aSrc,
    const std::vector<VertexDescriptor>& aPath,
    const double                         aCapacity) {
  removeCapacityFromPath(aSrc, aPath, -aCapacity, std::nullopt, theGraph);
}

std::vector<double> CapacityNetwork::nodeCapacities() const {
  std::vector<double> ret(boost::num_vertices(theGraph), 0);
  for (const auto& myNode :
       boost::make_iterator_range(boost::vertices(theGraph))) {
    double myCapacity = 0;
    for (const auto& myEdge :
         boost::make_iterator_range(boost::out_edges(myNode, theGraph))) {
      myCapacity += boost::get(boost::edge_weight, theGraph, myEdge);
    }
    assert(myNode < ret.size());
    ret[myNode] = myCapacity;
  }
  return ret;
}

void CapacityNetwork::toGnuplot(
    const std::string&             aFilename,
    const std::vector<Coordinate>& aCoordinates) const {
  if (boost::num_vertices(theGraph) == 0 or aFilename.empty()) {
    return;
  }

  if (boost::num_vertices(theGraph) != aCoordinates.size()) {
    throw std::runtime_error("Invalid number of coordinates: expected " +
                             std::to_string(boost::num_vertices(theGraph)) +
                             ", got " + std::to_string(aCoordinates.size()));
  }

  const auto myCapacities = nodeCapacities();
  assert(myCapacities.size() == boost::num_vertices(theGraph));

  std::ofstream myOutVertices(aFilename + "-vertices.dat");
  for (std::size_t i = 0; i < aCoordinates.size(); i++) {
    myOutVertices << i << ',' << std::get<0>(aCoordinates[i]) << ','
                  << std::get<1>(aCoordinates[i]) << ','
                  << std::get<2>(aCoordinates[i]) << ',' << myCapacities[i]
                  << '\n';
  }

  std::ofstream myOutEdges(aFilename + "-edges.dat");
  for (const auto& myEdge :
       boost::make_iterator_range(boost::edges(theGraph))) {
    myOutEdges << std::get<0>(aCoordinates[myEdge.m_source]) << ','
               << std::get<1>(aCoordinates[myEdge.m_source]) << '\n'
               << std::get<0>(aCoordinates[myEdge.m_target]) << ','
               << std::get<1>(aCoordinates[myEdge.m_target]) << '\n'
               << '\n';
  }
}

bool CapacityNetwork::checkCapacity(const VertexDescriptor               aSrc,
                                    const std::vector<VertexDescriptor>& aPath,
                                    const double aCapacity,
                                    const Graph& aGraph) {
  auto mySrc = aSrc;
  for (std::size_t i = 0; i < aPath.size(); i++) {
    auto myDst = aPath[i];

    EdgeDescriptor        myEdge;
    [[maybe_unused]] auto myFound = false;
    std::tie(myEdge, myFound)     = boost::edge(mySrc, myDst, aGraph);
    assert(myFound);
    if (boost::get(boost::edge_weight, aGraph, myEdge) < aCapacity) {
      return false;
    }

    // move to the next edge
    mySrc = myDst;
  }
  return true;
}

double
CapacityNetwork::minCapacity(const VertexDescriptor               aSrc,
                             const std::vector<VertexDescriptor>& aPath) const {
  double ret   = std::numeric_limits<double>::max();
  auto   mySrc = aSrc;
  for (std::size_t i = 0; i < aPath.size(); i++) {
    auto myDst = aPath[i];

    EdgeDescriptor myEdge;
    auto           myFound    = false;
    std::tie(myEdge, myFound) = boost::edge(mySrc, myDst, theGraph);
    if (not myFound) {
      throw std::runtime_error("non-existing path for src node " +
                               std::to_string(aSrc) + ": " +
                               ::toStringStd(aPath, ","));
    }
    ret = std::min(ret, boost::get(boost::edge_weight, theGraph, myEdge));

    // move to the next edge
    mySrc = myDst;
  }
  return ret;
}

double CapacityNetwork::minCapacity(const Path& aPath, const Graph& aGraph) {
  double ret = std::numeric_limits<double>::max();
  for (const auto& edge : aPath) {
    ret = std::min(ret, boost::get(boost::edge_weight, aGraph, edge));
  }
  return ret;
}

void CapacityNetwork::removeSmallestCapacityEdge(
    const VertexDescriptor               aSrc,
    const std::vector<VertexDescriptor>& aPath,
    Graph&                               aGraph) {
  auto           mySrc = aSrc;
  EdgeDescriptor mySmallestCapacityEdge;
  double         mySmallestCapacity = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i < aPath.size(); i++) {
    auto myDst = aPath[i];

    EdgeDescriptor        myEdge;
    [[maybe_unused]] auto myFound = false;
    std::tie(myEdge, myFound)     = boost::edge(mySrc, myDst, aGraph);
    assert(myFound);
    const auto myCapacity = boost::get(boost::edge_weight, aGraph, myEdge);
    if (myCapacity < mySmallestCapacity) {
      mySmallestCapacity     = myCapacity;
      mySmallestCapacityEdge = myEdge;
    }

    // move to the next edge
    mySrc = myDst;
  }
  boost::remove_edge(mySmallestCapacityEdge, aGraph);
}

void CapacityNetwork::removeCapacityFromPath(
    const VertexDescriptor               aSrc,
    const std::vector<VertexDescriptor>& aPath,
    const double                         aCapacity,
    const std::optional<double>          aMinCapacity,
    Graph&                               aGraph) {
  const auto V = boost::num_vertices(aGraph);
  if (aSrc >= V) {
    throw std::runtime_error("source node does not exist: " +
                             std::to_string(aSrc));
  }
  auto myWeights = boost::get(boost::edge_weight, aGraph);

  // first pass: only checks, throw if needed
  auto mySrc = aSrc;
  for (std::size_t i = 0; i < aPath.size(); i++) {
    auto myDst = aPath[i];
    if (myDst >= V) {
      throw std::runtime_error("intermediate node does not exist: " +
                               std::to_string(myDst));
    }

    EdgeDescriptor myEdge;
    auto           myFound    = false;
    std::tie(myEdge, myFound) = boost::edge(mySrc, myDst, aGraph);
    if (not myFound) {
      throw std::runtime_error("edge not in the graph: " + ::toString(myEdge));
    }
    auto& myWeight = myWeights[myEdge];
    if (myWeight < aCapacity) {
      throw std::runtime_error(
          "cannot remove capacity " + std::to_string(aCapacity) + " > " +
          std::to_string(myWeight) + " for edge " + ::toString(myEdge));
    }

    // move to the next edge
    mySrc = myDst;
  }

  // second pass: no checks, just to the job
  mySrc = aSrc;
  for (std::size_t i = 0; i < aPath.size(); i++) {
    auto myDst = aPath[i];

    auto  myEdge   = boost::edge(mySrc, myDst, aGraph).first;
    auto& myWeight = myWeights[myEdge];
    myWeight -= aCapacity;

    if (aMinCapacity.has_value() and myWeight < aMinCapacity.value()) {
      boost::remove_edge(myEdge, aGraph);
    }

    // move to the next edge
    mySrc = myDst;
  }
}

void CapacityNetwork::removeCapacityFromPath(
    const Path&                 aPath,
    const double                aCapacity,
    const std::optional<double> aMinCapacity,
    Graph&                      aGraph) {
  const auto V            = boost::num_vertices(aGraph);
  auto       myCapacities = boost::get(boost::edge_weight, aGraph);

  // first pass: only checks, throw if needed
  for (const auto& edge : aPath) {
    if (edge.m_source >= V) {
      throw std::runtime_error("vertex does not exist: " +
                               std::to_string(edge.m_source));
    }
    if (edge.m_target >= V) {
      throw std::runtime_error("vertex does not exist: " +
                               std::to_string(edge.m_target));
    }
    if (not boost::edge(edge.m_source, edge.m_target, aGraph).second) {
      throw std::runtime_error("edge (" + std::to_string(edge.m_source) + "," +
                               std::to_string(edge.m_target) +
                               ") does not exist");
    }
    const auto& myCapacity = myCapacities[edge];
    if (myCapacity < aCapacity) {
      throw std::runtime_error(
          "cannot remove capacity " + std::to_string(aCapacity) + " > " +
          std::to_string(myCapacity) + " for edge " + ::toString(edge));
    }
  }

  // second pass: no checks, just to the job
  for (const auto& edge : aPath) {
    auto& myCapacity = myCapacities[edge];
    myCapacity -= aCapacity;
    if (aMinCapacity.has_value() and myCapacity < aMinCapacity.value()) {
      boost::remove_edge(edge, aGraph);
    }
  }
}

CapacityNetwork::Path CapacityNetwork::intersect(const Path& aLhsPath,
                                                 const Path& aRhsPath) {
  std::set<EdgeDescriptor> myCommonEdges;
  for (const auto& edge : aLhsPath) {
    if (std::find(aRhsPath.begin(), aRhsPath.end(), edge) != aRhsPath.end()) {
      myCommonEdges.emplace(edge);
    }
  }
  for (const auto& edge : aRhsPath) {
    if (std::find(aLhsPath.begin(), aLhsPath.end(), edge) != aLhsPath.end()) {
      myCommonEdges.emplace(edge);
    }
  }
  return {myCommonEdges.begin(), myCommonEdges.end()};
}

CapacityNetwork::Path CapacityNetwork::difference(const Path& aLhsPath,
                                                  const Path& aRhsPath) {
  Path ret;
  for (const auto& edge : aLhsPath) {
    if (std::find(aRhsPath.begin(), aRhsPath.end(), edge) == aRhsPath.end()) {
      ret.emplace_back(edge);
    }
  }
  return ret;
}

CapacityNetwork::Path
CapacityNetwork::toPath(const unsigned long               aSource,
                        const std::vector<unsigned long>& aHops,
                        const Graph&                      aGraph) {
  Path ret;
  for (size_t i = 0; i < aHops.size(); i++) {
    EdgeDescriptor myEdge;
    auto           myFound    = false;
    const auto     u          = i == 0 ? aSource : aHops[i - 1];
    const auto     v          = aHops[i];
    std::tie(myEdge, myFound) = boost::edge(u, v, aGraph);
    if (not myFound) {
      throw std::runtime_error("edge (" + std::to_string(u) + "," +
                               std::to_string(v) + ") does not exist");
    }
    ret.emplace_back(myEdge);
  }
  assert(ret.size() == aHops.size());
  return ret;
}

std::pair<std::size_t, std::size_t> CapacityNetwork::minMaxVertexProp(
    const std::function<std::size_t(Graph::vertex_descriptor, const Graph&)>&
        aPropFunctor) const {
  auto        myRange = boost::vertices(theGraph);
  std::size_t myMin   = std::numeric_limits<std::size_t>::max();
  std::size_t myMax   = 0;
  for (auto it = myRange.first; it != myRange.second; ++it) {
    const auto myCur = aPropFunctor(*it, theGraph);
    if (myCur < myMin) {
      myMin = myCur;
    }
    if (myCur > myMax) {
      myMax = myCur;
    }
  }
  return {myMin, myMax};
}

std::string toString(const CapacityNetwork::Path& aPath) {
  std::stringstream ret;
  ret << '[' << ::toString(aPath, ", ", [](const auto& elem) {
    return "(" + std::to_string(elem.m_source) + "," +
           std::to_string(elem.m_target) + ")";
  }) << ']';
  return ret.str();
}

} // namespace qr
} // namespace uiiit
