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

#include <boost/graph/detail/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/properties.hpp>

#include <boost/property_map/property_map.hpp>
#include <glog/logging.h>

#include <limits>
#include <sstream>
#include <stdexcept>

namespace uiiit {
namespace qr {

CapacityNetwork::FlowDescriptor::FlowDescriptor(const unsigned long aSrc,
                                                const unsigned long aDst,
                                                const double aNetRate) noexcept
    : theSrc(aSrc)
    , theDst(aDst)
    , theNetRate(aNetRate)
    , thePath()
    , theGrossRate(0)
    , theDijsktra(0) {
  // noop
}

void CapacityNetwork::FlowDescriptor::movePathRateFrom(
    FlowDescriptor& aAnother) {
  std::swap(thePath, aAnother.thePath);
  std::swap(theGrossRate, aAnother.theGrossRate);
}

std::string CapacityNetwork::FlowDescriptor::toString() const {
  std::stringstream myStream;
  myStream << "(" << theSrc << "," << theDst << ") [net rate = " << theNetRate
           << ", gross rate " << theGrossRate << "] path = [" << thePath.size()
           << "]{"
           << ::toString(
                  thePath,
                  ",",
                  [](const auto& aValue) { return std::to_string(aValue); })
           << "}, " << theDijsktra << " Dijkstra calls";
  return myStream.str();
}

CapacityNetwork::AppDescriptor::AppDescriptor(
    const unsigned long               aHost,
    const std::vector<unsigned long>& aPeers,
    const double                      aPriority) noexcept
    : theHost(aHost)
    , thePeers(aPeers)
    , thePriority(aPriority)
    , theAllPaths()
    , theRemainingPaths()
    , theAllocated()
    , theYen(0)
    , theVisits(0) {
  // noop
}

double CapacityNetwork::AppDescriptor::netRate() const {
  double ret = 0;
  for (const auto& elem : theAllocated) {
    ret += elem.second.theNetRate;
  }
  return ret;
}

double CapacityNetwork::AppDescriptor::grossRate() const {
  double ret = 0;
  for (const auto& elem : theAllocated) {
    ret += elem.second.theGrossRate;
  }
  return ret;
}

std::string CapacityNetwork::AppDescriptor::toString() const {
  std::stringstream myStream;
  myStream << "host " << theHost << ", peers {"
           << ::toString(
                  thePeers,
                  ",",
                  [](const auto& aPeer) { return std::to_string(aPeer); })
           << "}, prio " << thePriority << ", " << theRemainingPaths.size()
           << " remaining paths (" << theAllPaths.size() << " total), "
           << theAllocated.size() << " paths allocated with totale capacity "
           << netRate() << " (gross " << grossRate() << "), " << theYen
           << " Yen algorithm calls, " << theVisits << " visits made";
  return myStream.str();
}

CapacityNetwork::CapacityNetwork(
    const std::vector<std::pair<unsigned long, unsigned long>>& aEdges,
    support::RealRvInterface&                                   aWeightRv,
    const bool aMakeBidirectional)
    : Network()
    , theGraph()
    , theMeasurementProbability(1) {
  for (const auto& myEdge : aEdges) {
    const auto myWeight = aWeightRv();
    Utils<Graph>::addEdge(theGraph, myEdge.first, myEdge.second, myWeight);
    if (aMakeBidirectional) {
      Utils<Graph>::addEdge(theGraph, myEdge.second, myEdge.first, myWeight);
    }
  }
}

CapacityNetwork::CapacityNetwork(const WeightVector& aEdgeWeights)
    : Network()
    , theGraph()
    , theMeasurementProbability(1) {
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

void CapacityNetwork::measurementProbability(
    const double aMeasurementProbability) {
  if (aMeasurementProbability < 0 or aMeasurementProbability > 1) {
    throw std::runtime_error("Invalid measurement probability: " +
                             std::to_string(aMeasurementProbability));
  }
  VLOG(1) << "measurement probability set to " << aMeasurementProbability;
  theMeasurementProbability = aMeasurementProbability;
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

void CapacityNetwork::route(std::vector<FlowDescriptor>& aFlows,
                            const FlowCheckFunction&     aCheckFunction) {
  const auto V = boost::num_vertices(theGraph);

  // pre-condition checks
  for (const auto& myFlow : aFlows) {
    assert(myFlow.thePath.empty());
    assert(myFlow.theGrossRate == 0);
    assert(myFlow.theDijsktra == 0);

    if (boost::vertex(myFlow.theSrc, theGraph) >= V) {
      throw std::runtime_error("invalid source node in flow: " +
                               std::to_string(myFlow.theSrc));
    }
    if (boost::vertex(myFlow.theDst, theGraph) >= V) {
      throw std::runtime_error("invalid destination node in flow: " +
                               std::to_string(myFlow.theDst));
    }
    if (myFlow.theSrc == myFlow.theDst) {
      throw std::runtime_error("invalid flow: from " +
                               std::to_string(myFlow.theSrc) + " to itself");
    }
    if (myFlow.theNetRate <= 0) {
      throw std::runtime_error(
          "invalid nonpositive capacity request in flow: " +
          std::to_string(myFlow.theNetRate));
    }
  }

  std::vector<VertexDescriptor> myDistances(V);
  std::vector<VertexDescriptor> myPredecessors(V);

  for (auto& myFlow : aFlows) {
    assert(myFlow.thePath.empty());
    assert(myFlow.theGrossRate == 0);
    VLOG(1) << "flow " << myFlow.toString();

    auto myFoundOrDisconnected = false;
    auto myCopiedGraph         = theGraph;

    // loop until either there is no path from the source to the destination or
    // we find a candidate that can satisfy the flow requirements
    while (not myFoundOrDisconnected) {
      myFlow.theDijsktra++;
      boost::dijkstra_shortest_paths(
          myCopiedGraph,
          myFlow.theSrc,
          boost::predecessor_map(myPredecessors.data())
              .weight_map(
                  boost::make_static_property_map<Graph::edge_descriptor>(1))
              .distance_map(boost::make_iterator_property_map(
                  myDistances.data(),
                  get(boost::vertex_index, myCopiedGraph))));

      if (myPredecessors[myFlow.theDst] == myFlow.theDst) {
        myFoundOrDisconnected = true; // disconnected

      } else {
        // there is at least one path from source to destination
        HopsFinder     myHopsFinder(myPredecessors, myFlow.theSrc);
        FlowDescriptor myCandidate(myFlow);
        myHopsFinder(myCandidate.thePath, myFlow.theDst);
        assert(not myCandidate.thePath.empty());
        myCandidate.theGrossRate =
            toGrossRate(myCandidate.theNetRate, myCandidate.thePath.size());
        VLOG(1) << "candidate " << myCandidate.toString();

        // if the flow is not admissible because of the external function
        // we assume there is no need to continue the search, otherwise
        // we check that the gross EPR rate is feasible along the path selected
        if (not aCheckFunction(myCandidate)) {
          myFoundOrDisconnected = true;

        } else if (checkCapacity(myCandidate.theSrc,
                                 myCandidate.thePath,
                                 myCandidate.theGrossRate,
                                 myCopiedGraph)) {
          // flow is admissible on the shortest path, we can break from the
          // loop
          myFoundOrDisconnected = true;
          myFlow.movePathRateFrom(myCandidate);

        } else {
          // flow not admissible on the shortest path, remove the edge with
          // smallest capacity along the path and try again
          removeSmallestCapacityEdge(
              myCandidate.theSrc, myCandidate.thePath, myCopiedGraph);
        }
      }
    }

    if (myFlow.thePath.empty()) {
      VLOG(1) << "flow rejected " << myFlow.toString();

    } else {
      VLOG(1) << "flow admitted " << myFlow.toString();
      // flow admissible: remove the gross capacity from the edges along the
      // path and then move to the next flow in the list
      removeCapacityFromPath(
          myFlow.theSrc, myFlow.thePath, myFlow.theGrossRate, theGraph);
    }
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
    Graph&                               aGraph) {
  auto mySrc     = aSrc;
  auto myWeights = boost::get(boost::edge_weight, aGraph);
  for (std::size_t i = 0; i < aPath.size(); i++) {
    auto myDst = aPath[i];

    EdgeDescriptor        myEdge;
    [[maybe_unused]] auto myFound = false;
    std::tie(myEdge, myFound)     = boost::edge(mySrc, myDst, aGraph);
    assert(myFound);
    assert(myWeights[myEdge] >= aCapacity);
    myWeights[myEdge] -= aCapacity;

    // move to the next edge
    mySrc = myDst;
  }
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

double CapacityNetwork::toGrossRate(const double      aNetRate,
                                    const std::size_t aNumEdges) const {
  if (aNumEdges <= 1) {
    return aNetRate;
  }
  return aNetRate / std::pow(theMeasurementProbability, aNumEdges - 1);
}

} // namespace qr

} // end namespace uiiit
