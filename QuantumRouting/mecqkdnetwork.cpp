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

#include "Support/tostring.h"

#include <algorithm>
#include <boost/graph/subgraph.hpp>
#include <glog/logging.h>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>

namespace uiiit {
namespace qr {

std::vector<MecQkdAlgo> allMecQkdAlgos() {
  static const std::vector<MecQkdAlgo> myAlgos({
      MecQkdAlgo::Random,
      MecQkdAlgo::Spf,
      MecQkdAlgo::BestFit,
      MecQkdAlgo::RandomBlind,
      MecQkdAlgo::SpfBlind,
      MecQkdAlgo::BestFitBlind,
      MecQkdAlgo::SpfStatic,
  });
  return myAlgos;
}

std::string toString(const MecQkdAlgo aAlgo) {
  switch (aAlgo) {
    case MecQkdAlgo::Random:
      return "random";
    case MecQkdAlgo::Spf:
      return "spf";
    case MecQkdAlgo::BestFit:
      return "bestfit";
    case MecQkdAlgo::RandomBlind:
      return "random-blind";
    case MecQkdAlgo::SpfBlind:
      return "spf-blind";
    case MecQkdAlgo::BestFitBlind:
      return "bestfit-blind";
    case MecQkdAlgo::SpfStatic:
      return "spf-static";
    default:; /* fall-through */
  }
  return "unknown";
}

MecQkdAlgo mecQkdAlgofromString(const std::string& aAlgo) {
  if (aAlgo == "random") {
    return MecQkdAlgo::Random;
  } else if (aAlgo == "bestfit") {
    return MecQkdAlgo::BestFit;
  } else if (aAlgo == "spf") {
    return MecQkdAlgo::Spf;
  } else if (aAlgo == "random-blind") {
    return MecQkdAlgo::RandomBlind;
  } else if (aAlgo == "bestfit-blind") {
    return MecQkdAlgo::BestFitBlind;
  } else if (aAlgo == "spf-blind") {
    return MecQkdAlgo::SpfBlind;
  } else if (aAlgo == "spf-static") {
    return MecQkdAlgo::SpfStatic;
  }
  throw std::runtime_error(
      "invalid edge QKD algorithm: " + aAlgo + " (valid options are: " +
      ::toString(allMecQkdAlgos(),
                 ",",
                 [](const auto& aAlgo) { return toString(aAlgo); }) +
      ")");
}

std::string MecQkdNetwork::EdgeNode::toString() const {
  std::stringstream ret;
  ret << "#" << theId << ", residual " << theResidual << " (tot "
      << theAvailable << "), path [" << thePathSize << "] {"
      << ::toStringStd(thePath, ",") << "} " << (feasible() ? "v" : "x");
  return ret.str();
}

MecQkdNetwork::Allocation::Allocation(const unsigned long aUserNode,
                                      const double        aRate,
                                      const double        aLoad)
    : theUserNode(aUserNode)
    , theRate(aRate)
    , theLoad(aLoad) {
  // noop
}

std::string MecQkdNetwork::Allocation::toString() const {
  std::stringstream ret;
  ret << "user " << theUserNode << ", rate " << theRate << " kb/s, load "
      << theLoad << ": ";
  if (theAllocated) {
    ret << "allocated to edge " << theEdgeNode << ", path length "
        << thePathLength << ", tot rate " << totRate() << " kb/s";
  } else {
    ret << "not allocated";
  }
  return ret.str();
}

MecQkdNetwork::MecQkdNetwork(
    const std::vector<std::pair<unsigned long, unsigned long>>& aEdges,
    support::RealRvInterface&                                   aWeightRv,
    const bool aMakeBidirectional)
    : CapacityNetwork(aEdges, aWeightRv, aMakeBidirectional) {
  // noop
}

MecQkdNetwork::MecQkdNetwork(const WeightVector& aEdgeWeights)
    : CapacityNetwork(aEdgeWeights) {
  // noop
}

void MecQkdNetwork::userNodes(const std::set<unsigned long>& aUserNodes) {
  const auto V = boost::num_vertices(theGraph);
  for (const auto& v : aUserNodes) {
    if (v >= V) {
      throw std::runtime_error("Invalid user node: " + std::to_string(v));
    }
  }

  theUserNodes = aUserNodes;
}

void MecQkdNetwork::edgeNodes(
    const std::map<unsigned long, double>& aEdgeProcessing) {
  const auto V = boost::num_vertices(theGraph);
  theEdgeNodes.clear();
  for (const auto& elem : aEdgeProcessing) {
    if (elem.second < 0) {
      throw std::runtime_error("Invalid negative processing capability: " +
                               std::to_string(elem.second));
    }
    if (elem.first >= V) {
      throw std::runtime_error("Invalid edge node: " +
                               std::to_string(elem.first));
    }
    theEdgeNodes.emplace(elem.first);
  }
  theEdgeProcessing = aEdgeProcessing;
}

double MecQkdNetwork::totProcessing() const {
  return std::accumulate(
      theEdgeProcessing.begin(),
      theEdgeProcessing.end(),
      0.0,
      [](auto aSum, const auto& aElem) { return aSum + aElem.second; });
}

void MecQkdNetwork::allocate(std::vector<Allocation>&  aApps,
                             const MecQkdAlgo          aAlgo,
                             support::RealRvInterface& aRv) {
  static const auto EPSILON = 1e-5;
  if (theUserNodes.empty()) {
    throw std::runtime_error("invalid empty set of user nodes");
  }
  if (theEdgeNodes.empty()) {
    throw std::runtime_error("invalid empty set of edge nodes");
  }

  if (VLOG_IS_ON(1)) {
    LOG(INFO) << "allocation is about to start with " << toString(aAlgo);
    LOG(INFO) << "network capacity   : " << totalCapacity() << " b/s";
    LOG(INFO) << "processing capacity: " << totProcessing();
    LOG(INFO) << "user nodes         : " << ::toStringStd(theUserNodes, ",");
    LOG(INFO) << "edge nodes         : "
              << ::toString(theEdgeProcessing, ",", [](const auto& elem) {
                   return std::to_string(elem.first) + " (" +
                          std::to_string(elem.second) + ")";
                 });
  }

  std::vector<EdgeNode> myCandidates;
  for (auto& elem : theEdgeProcessing) {
    myCandidates.emplace_back(
        EdgeNode{elem.first, elem.second, 0, 0, std::vector<unsigned long>()});
  }

  // special allocation with SpfStatic
  if (aAlgo == MecQkdAlgo::SpfStatic) {
    allocateSpfStatic(aApps);
    return;
  }

  for (auto& myApp : aApps) {
    // compute, for each candidate, the residual capacity, which can be
    // negative, and the constrained shortest path length, which can be empty
    const auto myPaths = cspf(myApp.theUserNode, myApp.theRate, theEdgeNodes);
    assert(myPaths.size() == theEdgeNodes.size());
    for (auto& myCandidate : myCandidates) {
      const auto it = myPaths.find(myCandidate.theId);
      assert(it != myPaths.end());
      myCandidate.thePath = it->second;
      myCandidate.thePathSize =
          it->second.empty() ? 0.0 : (it->second.size() + aRv() * 0.1);

      myCandidate.theResidual = myCandidate.theAvailable - myApp.theLoad;
    }

    if (VLOG_IS_ON(2)) {
      LOG(INFO) << "candidates for " << myApp.toString();
      for (const auto& myCandidate : myCandidates) {
        LOG(INFO) << myCandidate.toString();
      }
    }

    auto mySelected = selectCandidate(myCandidates, aAlgo, aRv);
    if (mySelected != myCandidates.end() and mySelected->feasible()) {
      // save the allocation data into the output
      myApp.theAllocated  = true;
      myApp.theEdgeNode   = mySelected->theId;
      myApp.thePathLength = static_cast<std::size_t>(mySelected->thePathSize);

      // update the available processing capacity on the node selected
      assert(mySelected->theAvailable >= myApp.theLoad);
      mySelected->theAvailable -= myApp.theLoad;

      // update the network capacities
      assert(mySelected->thePath.size() == myApp.thePathLength);
      removeCapacityFromPath(myApp.theUserNode,
                             mySelected->thePath,
                             myApp.theRate,
                             EPSILON,
                             theGraph);
    }

    VLOG(1) << myApp.toString();
  }

  // update the edge node available processing power after the allocation
  for (const auto& myCandidate : myCandidates) {
    auto it = theEdgeProcessing.find(myCandidate.theId);
    assert(it != theEdgeProcessing.end());
    it->second = myCandidate.theAvailable;
  }
}

MecQkdNetwork::Candidates::iterator
MecQkdNetwork::selectCandidate(std::vector<EdgeNode>&    aCandidates,
                               const MecQkdAlgo          aAlgo,
                               support::RealRvInterface& aRv) {
  // return immediately a pointer beyond the end if there are no candidates
  if (aCandidates.begin() == aCandidates.end()) {
    return aCandidates.end();
  }

  // random algorithms
  if (aAlgo == MecQkdAlgo::Random or aAlgo == MecQkdAlgo::RandomBlind) {
    // with non-blind: retrieve the intersection of edge nodes that feasible
    // according to both rate and processing constraints
    // with blind: just return a candidate truly at random
    std::set<Candidates::iterator> myFilteredCandidates;
    for (auto it = aCandidates.begin(); it != aCandidates.end(); ++it) {
      if (aAlgo == MecQkdAlgo::RandomBlind or it->feasible()) {
        myFilteredCandidates.emplace(it);
      }
    }
    if (myFilteredCandidates.empty()) {
      return aCandidates.end();
    }
    return support::choice(myFilteredCandidates, aRv);
  }

  // for the other algorithms, just loop through the vector and keep an iterator
  // to the current best solution found to be returned at the end
  assert(aAlgo == MecQkdAlgo::BestFit or aAlgo == MecQkdAlgo::BestFitBlind or
         aAlgo == MecQkdAlgo::Spf or aAlgo == MecQkdAlgo::SpfBlind);
  struct Compare {
    Compare(const MecQkdAlgo aAlgo)
        : theAlgo(aAlgo) {
      // noop
    }
    Candidates::iterator best(const Candidates::iterator aLhs,
                              const Candidates::iterator aRhs) const noexcept {
      assert(theAlgo != MecQkdAlgo::Random and
             theAlgo != MecQkdAlgo::RandomBlind);
      switch (theAlgo) {
        case MecQkdAlgo::Random:
          return aLhs;
        case MecQkdAlgo::Spf:
          return (aRhs->feasible() and
                  (not aLhs->feasible() or
                   aLhs->thePathSize > aRhs->thePathSize)) ?
                     aRhs :
                     aLhs;
        case MecQkdAlgo::BestFit:
          return (aRhs->feasible() and
                  (not aLhs->feasible() or
                   aLhs->theResidual > aRhs->theResidual)) ?
                     aRhs :
                     aLhs;
        case MecQkdAlgo::SpfBlind:
          return (aRhs->feasiblePath() and
                  (not aLhs->feasiblePath() or
                   aLhs->thePathSize > aRhs->thePathSize)) ?
                     aRhs :
                     aLhs;
        case MecQkdAlgo::BestFitBlind:
          return (aRhs->feasibleResidual() and
                  (not aLhs->feasibleResidual() or
                   aLhs->theResidual > aRhs->theResidual)) ?
                     aRhs :
                     aLhs;
        case MecQkdAlgo::RandomBlind:
          return aLhs;
        case MecQkdAlgo::SpfStatic:
          assert(false);
          return aLhs;
      }
      assert(false);
      return aLhs;
    }
    const MecQkdAlgo theAlgo;
  };
  Compare              myCmp(aAlgo);
  Candidates::iterator myBest = aCandidates.begin();
  for (auto myCurr = aCandidates.begin() + 1; myCurr != aCandidates.end();
       ++myCurr) {
    myBest = myCmp.best(myBest, myCurr);
  }
  return myBest;
}

void MecQkdNetwork::allocateSpfStatic(std::vector<Allocation>& aApps) {
  assert(not theEdgeNodes.empty());

  // find the static paths for the source nodes of the applications
  std::map<unsigned long, std::vector<unsigned long>> myStatic;
  for (const auto& myApp : aApps) {
    auto cur =
        myStatic.emplace(myApp.theUserNode, std::vector<unsigned long>());
    if (not cur.second) {
      // static path already found for this source node
      continue;
    }
    const auto myPaths = cspf(cur.first->first, 0, theEdgeNodes);
    const auto min     = std::min_element(
        myPaths.begin(), myPaths.end(), [](const auto& aLhs, const auto& aRhs) {
          return aLhs.second.size() < aRhs.second.size();
        });
    assert(min != myPaths.end());
    cur.first->second = min->second;
  }

#ifndef NDEBUG
  VLOG(2) << "static paths:";
  for (const auto& myPath : myStatic) {
    VLOG(2) << myPath.first << " -> " << ::toStringStd(myPath.second, ",");
  }
#endif

  // allocate the applications
  for (auto& myApp : aApps) {
    auto myCandidate = myStatic.find(myApp.theUserNode);
    assert(myCandidate != myStatic.end());

    // if there is no path available, then skip immediately this app
    if (myCandidate->second.empty()) {
      VLOG(1) << myApp.toString();
      continue;
    }
    const auto myEdgeNode = myCandidate->second.back();

    auto myProcIt = theEdgeProcessing.find(myEdgeNode);
    assert(myProcIt != theEdgeProcessing.end());

    const auto myMinCapacity =
        minCapacity(myApp.theUserNode, myCandidate->second);

    // if feasible for both capacity and load, allocate the app
    VLOG(2) << "user node " << myApp.theUserNode << ", edge node " << myEdgeNode
            << ", path " << ::toStringStd(myCandidate->second, ",")
            << ", capacity requested " << myApp.theRate << " vs. available "
            << myMinCapacity << ", load requested " << myApp.theLoad
            << " vs. available " << myProcIt->second;
    if (myApp.theRate <= myMinCapacity and myApp.theLoad < myProcIt->second) {
      myApp.theAllocated  = true;
      myApp.theEdgeNode   = myEdgeNode;
      myApp.thePathLength = myCandidate->second.size();

      // remove the load from the target edge node
      myProcIt->second -= myApp.theLoad;

      // remove the capacity along the path
      removeCapacityFromPath(myApp.theUserNode,
                             myCandidate->second,
                             myApp.theRate,
                             std::nullopt,
                             theGraph);
    }
    VLOG(1) << myApp.toString();
  }
}

} // namespace qr
} // namespace uiiit
