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

#include "QuantumRouting/mecqkdonlinenetwork.h"

#include "Support/tostring.h"

#include "yen/yen_ksp.hpp"

#include <algorithm>
#include <boost/graph/subgraph.hpp>
#include <glog/logging.h>
#include <optional>
#include <sstream>
#include <stdexcept>

namespace uiiit {
namespace qr {

std::vector<MecQkdOnlineAlgo> allMecQkdOnlineAlgos() {
  static const std::vector<MecQkdOnlineAlgo> myAlgos({
      MecQkdOnlineAlgo::Policy014_k1,
      MecQkdOnlineAlgo::Policy014_k3,
      MecQkdOnlineAlgo::Policy015,
      MecQkdOnlineAlgo::Policy015_reuse,
  });
  return myAlgos;
}

std::string toString(const MecQkdOnlineAlgo aAlgo) {
  switch (aAlgo) {
    case MecQkdOnlineAlgo::Policy014_k1:
      return "policy-014-k-1";
    case MecQkdOnlineAlgo::Policy014_k3:
      return "policy-014-k-3";
    case MecQkdOnlineAlgo::Policy015:
      return "policy-015";
    case MecQkdOnlineAlgo::Policy015_reuse:
      return "policy-015-reuse";
    default:; /* fall-through */
  }
  return "unknown";
}

MecQkdOnlineAlgo mecQkdOnlineAlgofromString(const std::string& aAlgo) {
  if (aAlgo == "policy-014-k-1") {
    return MecQkdOnlineAlgo::Policy014_k1;
  } else if (aAlgo == "policy-014-k-3") {
    return MecQkdOnlineAlgo::Policy014_k3;
  } else if (aAlgo == "policy-015") {
    return MecQkdOnlineAlgo::Policy015;
  } else if (aAlgo == "policy-015-reuse") {
    return MecQkdOnlineAlgo::Policy015_reuse;
  }
  throw std::runtime_error(
      "invalid edge QKD online algorithm: " + aAlgo + " (valid options are: " +
      ::toString(allMecQkdOnlineAlgos(),
                 ",",
                 [](const auto& aAlgo) { return toString(aAlgo); }) +
      ")");
}

std::string MecQkdOnlineNetwork::EdgeNode::toString() const {
  std::stringstream ret;
  ret << "#" << theId << ", residual " << theResidual << " (tot "
      << theAvailable << "), path [" << thePathSize << "] {"
      << ::toStringStd(thePath, ",") << "} " << (feasible() ? "v" : "x");
  return ret.str();
}

MecQkdOnlineNetwork::Allocation::Allocation(const unsigned long aUserNode,
                                            const double        aRate,
                                            const double        aLoad)
    : theUserNode(aUserNode)
    , theRate(aRate)
    , theLoad(aLoad) {
  // noop
}

std::string MecQkdOnlineNetwork::Allocation::toString() const {
  std::stringstream ret;
  ret << "user " << theUserNode << ", rate " << theRate << " kb/s, load "
      << theLoad;
  if (allocated()) {
    ret << ", #" << theId << ", allocated to edge " << edgeNode()
        << ", path length " << thePath.size() << ", tot rate " << totRate()
        << " kb/s";
  }
  return ret.str();
}

MecQkdOnlineNetwork::MecQkdOnlineNetwork(const WeightVector& aEdgeWeights)
    : CapacityNetwork(aEdgeWeights)
    , theAlgorithm(MecQkdOnlineAlgo::NotInitialized)
    , theRv()
    , theUserNodes()
    , theEdgeProcessing()
    , theEdgeNodes()
    , theActiveApps()
    , thePaths() {
  // noop
}

void MecQkdOnlineNetwork::configure(
    const MecQkdOnlineAlgo&                     aAlgorithm,
    std::unique_ptr<support::RealRvInterface>&& aRv,
    const std::set<unsigned long>&              aUserNodes,
    const std::map<unsigned long, double>&      aEdgeProcessing) {

  const auto V = boost::num_vertices(theGraph);

  // consistency checks
  if (aUserNodes.empty()) {
    throw std::runtime_error("empty set of user nodes");
  }
  if (aEdgeProcessing.empty()) {
    throw std::runtime_error("empty set of edge nodes");
  }
  if ((aUserNodes.size() + aEdgeProcessing.size()) > V) {
    throw std::runtime_error("too many user+edge nodes requested");
  }

  // throw if the caller is trying to reconfigure an operational network
  if (theNextId > 0) {
    throw std::runtime_error("cannot reconfigure a network while operating");
  }

  // set the algoritm and r.v.
  theAlgorithm = aAlgorithm;
  theRv        = std::move(aRv);

  // set the user and edge nodes
  for (const auto& v : aUserNodes) {
    if (v >= V) {
      throw std::runtime_error("Invalid user node: " + std::to_string(v));
    }
    if (aEdgeProcessing.find(v) != aEdgeProcessing.end()) {
      throw std::runtime_error("node " + std::to_string(v) +
                               " defined both as a user and as an edge");
    }
  }

  theUserNodes = aUserNodes;

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

  if (theAlgorithm == MecQkdOnlineAlgo::Policy014_k1 or
      theAlgorithm == MecQkdOnlineAlgo::Policy014_k3) {
    computeAllUserEdgePaths(theAlgorithm == MecQkdOnlineAlgo::Policy014_k1 ? 1 :
                                                                             3);
    assert(not thePaths.empty());
  }

  theTotProcessing = std::accumulate(
      theEdgeProcessing.begin(),
      theEdgeProcessing.end(),
      0.0,
      [](auto aSum, const auto& aElem) { return aSum + aElem.second; });
  theTotNetRate = 0;

  LOG(INFO) << "configured network with algorithm " << toString(theAlgorithm)
            << ", user nodes [" << ::toStringStd(theUserNodes, ",")
            << "], edge nodes ["
            << ::toString(theEdgeProcessing,
                          ",",
                          [](const auto& elem) {
                            std::stringstream ret;
                            ret << "(" << elem.first << "," << elem.second
                                << ")";
                            return ret.str();
                          })
            << "]";
}

void MecQkdOnlineNetwork::add(Allocation& aApp) {
  VLOG(2) << "request to add new app: " << aApp.toString();

  switch (theAlgorithm) {
    case MecQkdOnlineAlgo::NotInitialized:
      throw std::runtime_error("network not configured");
    case MecQkdOnlineAlgo::Policy014_k1:
    case MecQkdOnlineAlgo::Policy014_k3:
      addPolicy014(aApp);
      break;
    case MecQkdOnlineAlgo::Policy015:
    case MecQkdOnlineAlgo::Policy015_reuse:
      throw std::runtime_error("policy not yet implemented");
  }
}

void MecQkdOnlineNetwork::del(const uint64_t aAppId) {
  VLOG(2) << "request to del app #" << aAppId;

  auto it = theActiveApps.find(aAppId);
  if (it == theActiveApps.end()) {
    throw std::runtime_error("cannot remove non-active application #" +
                             std::to_string(aAppId));
  }

  // free the processing capacity on the edge node
  theEdgeProcessing[it->second.edgeNode()] += it->second.theLoad;
  theTotProcessing += it->second.theLoad;

  // free the used network rate in the network
  removeCapacityFromPath(
      it->second.thePath, -it->second.theRate, std::nullopt, theGraph);
  theTotNetRate -= it->second.theRate;

  // remove the app from the active ones
  theActiveApps.erase(it);
}

void MecQkdOnlineNetwork::computeAllUserEdgePaths(const std::size_t aK) {
  thePaths.clear(); // clear from previous assignment
  for (const auto& myUserNode : theUserNodes) {
    for (const auto& myEdgeNode : theEdgeNodes) {
      auto myIndexMap = boost::get(boost::vertex_index, theGraph);
      auto myResults  = boost::yen_ksp(
          theGraph,
          myUserNode,
          myEdgeNode,
          boost::make_static_property_map<Graph::edge_descriptor>(1),
          myIndexMap,
          aK);
      assert(myResults.size() <= aK);
      assert(thePaths.find(myUserNode) == thePaths.end() or
             thePaths[myUserNode].find(myEdgeNode) ==
                 thePaths[myUserNode].end());
      for (const auto& myResult : myResults) {
        assert(static_cast<std::size_t>(myResult.first) ==
               myResult.second.size()); // path length
        thePaths[myUserNode][myEdgeNode].emplace_back(myResult.second);
        assert(thePaths[myUserNode][myEdgeNode].size() <= aK);
      }
      assert(not myResults.empty() or
             thePaths.find(myUserNode) == thePaths.end() or
             thePaths[myUserNode].find(myEdgeNode) ==
                 thePaths[myUserNode].end());
    }
  }
  if (VLOG_IS_ON(2)) {
    LOG(INFO) << "initial static paths:\n" << pathsToString();
  }
}

void MecQkdOnlineNetwork::addPolicy014(Allocation& aApp) {
  auto it = thePaths.find(aApp.theUserNode);
  if (it == thePaths.end()) {
    throw std::runtime_error("invalid user node: " +
                             std::to_string(aApp.theUserNode));
  }

  // record all edge nodes that have sufficient network capacity in the primary
  // path and that have enough residual processing capacity as candidates:
  // key: residual capacity if the app is allocated to this edge node
  // value: the candidate edge node
  std::map<double, Path> myCandidates;
  for (const auto& elem : it->second) {
    assert(not elem.second.empty());
    auto kt = theEdgeProcessing.find(elem.first);
    assert(kt != theEdgeProcessing.end());
    if (minCapacity(elem.second.front(), theGraph) >= aApp.theRate and
        kt->second >= aApp.theLoad) {
      myCandidates.emplace(kt->second - aApp.theLoad + (*theRv)() * EPSILON,
                           elem.second.front());
    }
  }

  if (not myCandidates.empty()) {
    // select the candidate that minimizes the residual processing
    allocate(aApp, myCandidates.begin()->second);

  } else {
    // if not found, then search in the secondary paths
    for (const auto& elem : it->second) {
      assert(not elem.second.empty());

      // skip all the paths if the edge node cannot allocate the new app
      auto kt = theEdgeProcessing.find(elem.first);
      assert(kt != theEdgeProcessing.end());
      if (kt->second < aApp.theLoad) {
        continue;
      }

      // look into all the secondary paths, if any
      for (auto jt = std::next(elem.second.begin()); jt != elem.second.end();
           ++jt) {
        // compute the total rate from all other applications that would have
        // to be migrated to this path, if selected
        const auto myMigrationRate = std::accumulate(
            theActiveApps.begin(),
            theActiveApps.end(),
            0.0,
            [&aApp, &elem](auto aSum, const auto& aAnotherApp) {
              if (aAnotherApp.second.theUserNode == aApp.theUserNode and
                  aAnotherApp.second.edgeNode() == elem.first) {
                assert(aAnotherApp.second.thePath == elem.second.front());
                return aSum + aAnotherApp.second.theRate;
              }
              return aSum;
            });

        const auto myIntersection = intersect(*jt, elem.second.front());
        const auto myMinCapacityIntersection =
            minCapacity(myIntersection, theGraph);
        const auto myDifference = difference(*jt, elem.second.front());
        const auto myMinCapacityDifference =
            minCapacity(myDifference, theGraph);

        if (myMinCapacityIntersection >= aApp.theRate and
            myMinCapacityDifference >= (aApp.theRate + myMigrationRate)) {
          // add the current (secondary) path as a candidate iff the
          // following conditions are met:
          // 1. the capacity in the intersection between the primary and
          //    secondary path can withstand the new app
          // 2. the capacity in the edges of the secondary that are not already
          //    in the primary path (i.e., the difference path) can withstand
          //    the new app and all the active apps that would be migrated if
          //    the new app is allocated to this path
          // 3. the edge node has sufficient capacity to support the processing
          //    requirements of the new app (this condition is checked in the
          //    outer loop)
          myCandidates.emplace(kt->second - aApp.theLoad + (*theRv)() * EPSILON,
                               *jt);
          VLOG(2) << "candidate found in a secondary path, app "
                  << aApp.toString() << ", migration rate " << myMigrationRate
                  << ", primary path " << toString(elem.second.front())
                  << ", secondary path " << toString(*jt) << ", intersection "
                  << toString(myIntersection) << " ("
                  << myMinCapacityIntersection << "), difference "
                  << toString(myDifference) << " (" << myMinCapacityDifference
                  << "), residual processing " << kt->second;
        }
      }
    }

    // select the candidate that minimizes the residual processing, if any
    if (not myCandidates.empty()) {
      const auto myPath     = myCandidates.begin()->second;
      const auto myEdgeNode = myPath.crbegin()->m_target;
      auto       ht         = it->second.find(myEdgeNode);
      assert(ht != it->second.end());
      auto& myPaths = ht->second;

      // migrate all active apps with same user/edge node to the path selected
      for (auto& myApp : theActiveApps) {
        if (myApp.second.theUserNode == aApp.theUserNode and
            myApp.second.edgeNode() == myEdgeNode) {
          removeCapacityFromPath(myApp.second.thePath,
                                 -myApp.second.theRate,
                                 std::nullopt,
                                 theGraph);
          assert(minCapacity(myPath, theGraph) >= myApp.second.theRate);
          removeCapacityFromPath(
              myPath, myApp.second.theRate, std::nullopt, theGraph);
          myApp.second.thePath = myPath;
        }
      }

      // make the secondary path the new primary for this user-edge pair
      auto myPathIt = std::find(myPaths.begin(), myPaths.end(), myPath);
      assert(myPathIt != myPaths.begin());
      assert(myPathIt != myPaths.end());
      std::iter_swap(myPaths.begin(), myPathIt);

      // allocate the app on the selected path
      allocate(aApp, myPath);
    }
  }
}

std::string MecQkdOnlineNetwork::pathsToString() {
  std::stringstream ret;
  for (const auto& myOuter : thePaths) {
    for (const auto& myInner : myOuter.second) {
      for (const auto& myPath : myInner.second) {
        ret << myOuter.first << " -> " << myInner.first << ' '
            << toString(myPath) << '\n';
      }
    }
  }
  return ret.str();
}

void MecQkdOnlineNetwork::allocate(Allocation& aApp, const Path& aPath) {
  aApp.theId   = theNextId++;
  aApp.thePath = aPath;
  assert(theEdgeProcessing[aApp.edgeNode()] >= aApp.theLoad);
  theEdgeProcessing[aApp.edgeNode()] -= aApp.theLoad;
  theTotProcessing -= aApp.theLoad;
  assert(minCapacity(aApp.thePath, theGraph) >= aApp.theRate);
  removeCapacityFromPath(aApp.thePath, aApp.theRate, std::nullopt, theGraph);
  theTotNetRate += aApp.theRate;
  theActiveApps.emplace(aApp.theId,
                        aApp); // keep a copy in the active apps
}

} // namespace qr
} // namespace uiiit
