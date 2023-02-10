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

#include "Support/tostring.h"

#include "yen/yen_ksp.hpp"

#include <glog/logging.h>
#include <glog/vlog_is_on.h>

namespace uiiit {
namespace qr {

std::vector<AppRouteAlgo> allAppRouteAlgos() {
  static const std::vector<AppRouteAlgo> myAlgos({
      AppRouteAlgo::Random,
      AppRouteAlgo::BestFit,
      AppRouteAlgo::Drr,
  });
  return myAlgos;
}

std::string toString(const AppRouteAlgo aAlgo) {
  switch (aAlgo) {
    case AppRouteAlgo::Random:
      return "random";
    case AppRouteAlgo::BestFit:
      return "bestfit";
    case AppRouteAlgo::Drr:
      return "drr";
    default:; /* fall-through */
  }
  return "unknown";
}

AppRouteAlgo appRouteAlgofromString(const std::string& aAlgo) {
  if (aAlgo == "random") {
    return AppRouteAlgo::Random;
  } else if (aAlgo == "bestfit") {
    return AppRouteAlgo::BestFit;
  } else if (aAlgo == "drr") {
    return AppRouteAlgo::Drr;
  }
  throw std::runtime_error(
      "invalid app route algorithm: " + aAlgo + " (valid options are: " +
      ::toString(allAppRouteAlgos(),
                 ",",
                 [](const auto& aAlgo) { return toString(aAlgo); }) +
      ")");
}

EsNetwork::FlowDescriptor::FlowDescriptor(const unsigned long aSrc,
                                          const unsigned long aDst,
                                          const double        aNetRate) noexcept
    : theSrc(aSrc)
    , theDst(aDst)
    , theNetRate(aNetRate)
    , thePath()
    , theGrossRate(0)
    , theDijsktra(0) {
  // noop
}

void EsNetwork::FlowDescriptor::movePathRateFrom(FlowDescriptor& aAnother) {
  std::swap(thePath, aAnother.thePath);
  std::swap(theGrossRate, aAnother.theGrossRate);
}

std::string EsNetwork::FlowDescriptor::toString() const {
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

EsNetwork::AppDescriptor::AppDescriptor(
    const unsigned long               aHost,
    const std::vector<unsigned long>& aPeers,
    const double                      aPriority,
    const double                      aFidelityThreshold) noexcept
    : theHost(aHost)
    , thePeers(aPeers)
    , thePriority(aPriority)
    , theFidelityThreshold(aFidelityThreshold)
    , theRemainingPaths()
    , theAllocated()
    , theVisits(0) {
  // noop
}

bool EsNetwork::AppDescriptor::operator<(const AppDescriptor& aOther) const {
  return not theRemainingPaths.empty() and
         (aOther.theRemainingPaths.empty() or
          theRemainingPaths.begin()->first <
              aOther.theRemainingPaths.begin()->first);
}

EsNetwork::AppDescriptor::AppDescriptor::Output::Output(const Path& aPath)
    : theHops() {
  for (const auto& myEdge : aPath) {
    theHops.emplace_back(myEdge.m_target);
  }
}

double EsNetwork::AppDescriptor::accumulate(
    const std::function<double(const AppDescriptor::Output&)>& aFn) const {
  return std::accumulate(
      theAllocated.begin(),
      theAllocated.end(),
      0.0,
      [&aFn](auto aLhs, const auto& aRhs) {
        return std::accumulate(
            aRhs.second.begin(),
            aRhs.second.end(),
            aLhs,
            [&aFn](const auto aInnerLhs, const auto& aInnerRhs) {
              return aInnerLhs + aFn(aInnerRhs);
            });
      });
}

double EsNetwork::AppDescriptor::netRate() const {
  return accumulate([](const auto& aDesc) { return aDesc.theNetRate; });
}

double EsNetwork::AppDescriptor::grossRate() const {
  return accumulate([](const auto& aDesc) { return aDesc.theGrossRate; });
}

std::string EsNetwork::AppDescriptor::toString() const {
  std::stringstream myStream;
  myStream << "host " << theHost << ", peers {"
           << ::toString(
                  thePeers,
                  ",",
                  [](const auto& aPeer) { return std::to_string(aPeer); })
           << "}, prio " << thePriority << ", " << theRemainingPaths.size()
           << " remaining paths, " << theAllocated.size()
           << " paths allocated with total capacity " << netRate() << " (gross "
           << grossRate() << "), " << theVisits << " visits made";
  return myStream.str();
}

EsNetwork::EsNetwork(
    const std::vector<std::pair<unsigned long, unsigned long>>& aEdges,
    support::RealRvInterface&                                   aWeightRv,
    const bool aMakeBidirectional)
    : CapacityNetwork(aEdges, aWeightRv, aMakeBidirectional)
    , theMeasurementProbability(1) {
  // noop
}

EsNetwork::EsNetwork(const WeightVector& aEdgeWeights)
    : CapacityNetwork(aEdgeWeights)
    , theMeasurementProbability(1) {
  // noop
}

void EsNetwork::measurementProbability(const double aMeasurementProbability) {
  if (aMeasurementProbability < 0 or aMeasurementProbability > 1) {
    throw std::runtime_error("Invalid measurement probability: " +
                             std::to_string(aMeasurementProbability));
  }
  VLOG(2) << "measurement probability set to " << aMeasurementProbability;
  theMeasurementProbability = aMeasurementProbability;
}

double EsNetwork::maxNetRate(const AppDescriptor&    aApp,
                             const unsigned long     aPeer,
                             const AppCheckFunction& aCheckFunction) const {
  const auto NUM_PATHS  = 10u;
  auto       myIndexMap = boost::get(boost::vertex_index, theGraph);
  auto       myResult =
      boost::yen_ksp(theGraph,
                     aApp.theHost,
                     aPeer,
                     boost::make_static_property_map<Graph::edge_descriptor>(1),
                     myIndexMap,
                     NUM_PATHS);
  assert(myResult.size() <= NUM_PATHS);

  std::list<double> myNetRates({0.0}); // 0 is always a possible output
  for (auto& elem : myResult) {
    assert(elem.first > 0);
    assert(static_cast<std::size_t>(elem.first) == elem.second.size());
    if (aCheckFunction(aApp, elem.second) == false) {
      // this path is invalid
      continue;
    }
    myNetRates.emplace_back(
        toNetRate(minCapacity(elem.second, theGraph), elem.first));
  }
  assert(myNetRates.size() <= (1 + NUM_PATHS));

  // return the maximum net rate found (can be 0)
  assert(not myNetRates.empty());
  myNetRates.sort();
  return myNetRates.back();
}

void EsNetwork::route(std::vector<FlowDescriptor>& aFlows,
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
    VLOG(2) << "flow " << myFlow.toString();

    auto myFoundOrDisconnected = false;
    auto myCopiedGraph         = theGraph;

    // loop until either there is no path from the source to the destination
    // or we find a candidate that can satisfy the flow requirements
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
        VLOG(2) << "candidate " << myCandidate.toString();

        // if the flow is not admissible because of the external function
        // we assume there is no need to continue the search, otherwise
        // we check that the gross EPR rate is feasible along the path
        // selected
        if (not aCheckFunction(myCandidate)) {
          myFoundOrDisconnected = true;

        } else if (checkCapacity(myCandidate.theSrc,
                                 myCandidate.thePath,
                                 myCandidate.theGrossRate,
                                 myCopiedGraph)) {
          // flow is admissible on the shortest path, break from loop
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
      VLOG(2) << "flow rejected " << myFlow.toString();

    } else {
      VLOG(2) << "flow admitted " << myFlow.toString();
      // flow admissible: remove the gross capacity from the edges along the
      // path and then move to the next flow in the list
      removeCapacityFromPath(myFlow.theSrc,
                             myFlow.thePath,
                             myFlow.theGrossRate,
                             std::nullopt,
                             theGraph);
    }
  }
}

void EsNetwork::route(std::vector<AppDescriptor>& aApps,
                      const AppRouteAlgo          aAlgo,
                      const double                aQuantum,
                      support::RealRvInterface&   aRv,
                      const std::size_t           aK,
                      const AppCheckFunction&     aCheckFunction) {
  // check arguments
  if (aK == 0) {
    throw std::runtime_error("invalid k: cannot be null");
  }

  // pre-condition checks
  const auto V = boost::num_vertices(theGraph);
  for (const auto& myApp : aApps) {
    assert(myApp.theRemainingPaths.empty());
    assert(myApp.theAllocated.empty());
    assert(myApp.theVisits == 0);

    if (boost::vertex(myApp.theHost, theGraph) >= V) {
      throw std::runtime_error("invalid host node in app: " +
                               std::to_string(myApp.theHost));
    }
    for (const auto& myPeer : myApp.thePeers) {
      if (boost::vertex(myPeer, theGraph) >= V) {
        throw std::runtime_error("invalid peer node in app: " +
                                 std::to_string(myPeer));
      }
      if (myApp.theHost == myPeer) {
        throw std::runtime_error("invalid app: from " +
                                 std::to_string(myApp.theHost) + " to itself");
      }
    }
    if (myApp.thePriority <= 0) {
      throw std::runtime_error("invalid nonpositive priority for app: " +
                               std::to_string(myApp.thePriority));
    }
  }

  // for each app, find k-shortest paths towards each peer with Yen's
  // algorithm
  auto myIndexMap = boost::get(boost::vertex_index, theGraph);
  for (auto& myApp : aApps) {
    for (const auto& myPeer : myApp.thePeers) {
      auto myResult = boost::yen_ksp(
          theGraph,
          myApp.theHost,
          myPeer,
          boost::make_static_property_map<Graph::edge_descriptor>(1),
          myIndexMap,
          aK);
      for (auto& elem : myResult) {
        const auto myValid = aCheckFunction(myApp, elem.second);
        VLOG(2) << myApp.theHost << " -> " << myPeer << ": "
                << (myValid ? "valid" : "invalid") << " path found ["
                << elem.first << "] {"
                << ::toString(elem.second,
                              ",",
                              [](const auto& aEdge) {
                                return std::to_string(aEdge.m_target);
                              })
                << "}";
        if (myValid) {
          auto myRes =
              myApp.theRemainingPaths.emplace(elem.first, std::list<Path>());
          myRes.first->second.emplace_back(std::move(elem.second));
        }
      }
    }
  }

  // allocate resources based on the specific algorithm used
  if (aAlgo == AppRouteAlgo::Random) {
    routeRandom(aApps, aRv);
  } else if (aAlgo == AppRouteAlgo::BestFit) {
    routeBestFit(aApps);
  } else if (aAlgo == AppRouteAlgo::Drr) {
    routeDrr(aApps, aQuantum);
  } else {
    throw std::runtime_error("allocation strategy not implemented: " +
                             toString(aAlgo));
  }
}

double EsNetwork::toGrossRate(const double      aNetRate,
                              const std::size_t aNumEdges) const {
  if (aNumEdges <= 1) {
    return aNetRate;
  }
  return aNetRate / std::pow(theMeasurementProbability, aNumEdges - 1);
}

double EsNetwork::toNetRate(const double      aGrossRate,
                            const std::size_t aNumEdges) const {
  if (aNumEdges <= 1) {
    return aGrossRate;
  }
  return aGrossRate * std::pow(theMeasurementProbability, aNumEdges - 1);
}

void EsNetwork::routeRandom(std::vector<AppDescriptor>& aApps,
                            support::RealRvInterface&   aRv) {
  // create a structure that contains all apps with remaining paths
  std::list<std::size_t> myAppIndices;
  for (std::size_t ndx = 0; ndx < aApps.size(); ndx++) {
    if (not aApps[ndx].theRemainingPaths.empty()) {
      myAppIndices.emplace_back(ndx);
    }
  }

  // iterate by choosing one application at random, which is then
  // removed from list of indices if there are no more remaining paths
  auto myInfinite = std::numeric_limits<double>::max();
  while (not myAppIndices.empty()) {
    const auto myRndNdx = support::choice(myAppIndices, aRv);
    auto&      myCurApp = aApps[myRndNdx];
    schedule(myCurApp, myInfinite);
    if (myCurApp.theRemainingPaths.empty()) {
      auto it = std::find(myAppIndices.begin(), myAppIndices.end(), myRndNdx);
      assert(it != myAppIndices.end());
      myAppIndices.erase(it);
    }
  }
}

void EsNetwork::routeBestFit(std::vector<AppDescriptor>& aApps) {
  // iterate until there are no more applications with remaining paths
  // at each iteration select an app with minimum path length
  auto myInfinite = std::numeric_limits<double>::max();
  while (true) {
    auto it = std::min_element(aApps.begin(), aApps.end());
    if (it->theRemainingPaths.empty()) {
      break;
    }
    schedule(*it, myInfinite);
  }
}

void EsNetwork::routeDrr(std::vector<AppDescriptor>& aApps,
                         const double                aQuantum) {
  if (aQuantum <= 0) {
    throw std::runtime_error("invalid non-positive quantum value: " +
                             std::to_string(aQuantum));
  }

  // determine the quantum per application
  std::vector<double> myQuanta(aApps.size());
  double              mySumPriorities = 0;
  for (size_t i = 0; i < aApps.size(); i++) {
    mySumPriorities += aApps[i].thePriority;
  }
  for (size_t i = 0; i < aApps.size(); i++) {
    myQuanta[i] = aQuantum * aApps[i].thePriority / mySumPriorities;
  }

  // create an active list initialized with all (feasible) applications
  std::list<std::size_t> myActiveApps;
  for (std::size_t i = 0; i < aApps.size(); i++) {
    // only add apps with at least one path
    if (not aApps[i].theRemainingPaths.empty()) {
      myActiveApps.emplace_back(i);
    }
  }

  // do the allocation using deficit round robin
  auto myCurAppIt = myActiveApps.begin();
  while (not myActiveApps.empty()) {
    auto  myResidualCapacity = myQuanta[*myCurAppIt];
    auto& myCurApp           = aApps[*myCurAppIt];

    // loop until there are valid paths and capacity to be allocated
    while (not myCurApp.theRemainingPaths.empty() and myResidualCapacity > 0) {
      schedule(myCurApp, myResidualCapacity);
    }

    // check if there are feasible paths remaining for this app:
    // - yes: move to the next app in the active list
    // - no: remove this app from active list
    // in both cases it is possible that myCurAppIt reaches the end of the
    // list
    if (not myCurApp.theRemainingPaths.empty()) {
      ++myCurAppIt;
    } else {
      VLOG(2) << "remove from active list app w/ host " << myCurApp.theHost;
      myCurAppIt = myActiveApps.erase(myCurAppIt);
    }

    //  wrap-around at the end of the active list
    if (myCurAppIt == myActiveApps.end()) {
      myCurAppIt = myActiveApps.begin();
    }
  }
}

bool EsNetwork::schedule(AppDescriptor& aApp, double& aResidualCapacity) {
  // one more visit to this application
  aApp.theVisits++;

  // select the first of the shortest paths of the current app
  const auto& myCandidate = aApp.theRemainingPaths.begin()->second.front();
  VLOG(2) << "host " << aApp.theHost << ", path {"
          << ::toString(myCandidate,
                        ",",
                        [](const auto& aEdge) {
                          return std::to_string(aEdge.m_target);
                        })
          << "}";

  // check that all the edges still exist and find that with less capacity
  auto   myValidPath   = true;
  double myMinCapacity = std::numeric_limits<double>::max();
  for (const auto& elem : myCandidate) {
    EdgeDescriptor myEdge;
    bool           myFound;
    std::tie(myEdge, myFound) =
        boost::edge(elem.m_source, elem.m_target, theGraph);
    if (not myFound) {
      myValidPath = false;
      break;
    }
    myMinCapacity = std::min(myMinCapacity,
                             boost::get(boost::edge_weight, theGraph, myEdge));
  }

  // remove the path if it is not available anymore and return early
  if (not myValidPath) {
    aApp.theRemainingPaths.begin()->second.pop_front();
    if (aApp.theRemainingPaths.begin()->second.empty()) {
      // no more same-length paths, remove entry from theRemainingPaths
      aApp.theRemainingPaths.erase(aApp.theRemainingPaths.begin());
    }
    return false;
  }

  // set rate allocated
  const auto myAllocatedGross = std::min(myMinCapacity, aResidualCapacity);
  aResidualCapacity -= myAllocatedGross;

  // remove the gross capacity from all edges along the path
  // if the capacity becomes zero, remove the edge, too
  for (const auto& elem : myCandidate) {
    assert(boost::get(boost::edge_weight, theGraph, elem) >= myAllocatedGross);
    auto& myWeight = boost::get(boost::edge_weight, theGraph, elem) -=
        myAllocatedGross;
    if (myWeight == 0) {
      VLOG(2) << "removing edge (" << elem.m_source << "," << elem.m_target
              << ")";
      boost::remove_edge(elem, theGraph);
    }
  }

  // add the allocation to the path
  AppDescriptor::Output myOutput(myCandidate);
  VLOG(2) << "allocated gross capacity " << myAllocatedGross
          << " EPR-pairs/s for host " << aApp.theHost << " towards "
          << myOutput.theHops.back() << " along path {"
          << ::toString(myOutput.theHops,
                        ",",
                        [](const auto& aHop) { return std::to_string(aHop); })
          << "}, residual capacity " << aResidualCapacity;
  assert(not myOutput.theHops.empty());
  auto it = aApp.theAllocated.emplace(
      myOutput.theHops.back(), std::vector<AppDescriptor::Output>({myOutput}));

  const auto jt = std::find_if(it.first->second.begin(),
                               it.first->second.end(),
                               [&myOutput](const auto& aElem) {
                                 return myOutput.theHops == aElem.theHops;
                               });
  if (jt != it.first->second.end()) {
    jt->theNetRate += toNetRate(myAllocatedGross, jt->theHops.size());
    jt->theGrossRate += myAllocatedGross;
  }

  return true;
}

} // namespace qr
} // namespace uiiit
