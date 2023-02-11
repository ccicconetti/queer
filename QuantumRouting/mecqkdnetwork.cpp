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

#include "Support/split.h"
#include "Support/tostring.h"

#include <boost/graph/subgraph.hpp>
#include <fstream>
#include <glog/logging.h>
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
  }
  throw std::runtime_error(
      "invalid edge QKD algorithm: " + aAlgo + " (valid options are: " +
      ::toString(allMecQkdAlgos(),
                 ",",
                 [](const auto& aAlgo) { return toString(aAlgo); }) +
      ")");
}

std::string MecQkdWorkload::AppInfo::toString() const {
  std::stringstream ret;
  ret << "region " << theRegion << ", weight " << theWeight << ", load "
      << theLoad << ", rate " << theRate;
  return ret.str();
}

MecQkdWorkload::MecQkdWorkload(const std::vector<AppInfo>& aAppInfo,
                               support::RealRvInterface&   aRv)
    : theAppInfo(aAppInfo)
    , theRv(aRv)
    , theRegions()
    , theWeights(aAppInfo.size()) {
  if (aAppInfo.empty()) {
    throw std::runtime_error("invalid empty MecQkd workload");
  }
  for (std::size_t i = 0; i < aAppInfo.size(); i++) {
    theRegions.emplace(aAppInfo[i].theRegion);
    theWeights[i] = aAppInfo[i].theWeight;
  }

  if (VLOG_IS_ON(1)) {
    for (std::size_t i = 0; i < aAppInfo.size(); i++) {
      LOG(INFO) << aAppInfo[i].toString();
    }
  }
}

MecQkdWorkload MecQkdWorkload::fromCsvFile(const std::string&        aFilename,
                                           support::RealRvInterface& aRv) {
  std::vector<AppInfo> myAppInfo;

  std::ifstream myInfile(aFilename);
  if (not myInfile) {
    throw std::runtime_error("could not open file for reading: " + aFilename);
  }

  std::string myLine;
  std::size_t myLineNo = 0;
  while (myInfile) {
    ++myLineNo;
    std::getline(myInfile, myLine);
    if (myLine.empty() or myLine[0] == '#') {
      continue;
    }
    const auto myTokens = support::split<std::vector<std::string>>(myLine, ",");
    if (myTokens.size() != 4) {
      throw std::runtime_error("invalid input at file '" + aFilename +
                               "' line: " + std::to_string(myLineNo));
    }
    myAppInfo.emplace_back(AppInfo{std::stoull(myTokens[0]),
                                   std::stod(myTokens[1]),
                                   std::stod(myTokens[2]),
                                   std::stod(myTokens[3])});
  }

  return MecQkdWorkload(myAppInfo, aRv);
}

MecQkdWorkload::AppInfo MecQkdWorkload::operator()() {
  const auto res = support::sampleWeighted(theAppInfo, theWeights, 1, theRv);
  assert(res.size() == 1);
  return res[0];
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

  for (auto& myApp : aApps) {
    // compute, for each candidate, the residual capacity, which can be
    // negative, and the constrained shortest path length, which can be empty
    const auto myPaths = cspf(myApp.theUserNode, myApp.theRate, theEdgeNodes);
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

} // namespace qr
} // namespace uiiit
