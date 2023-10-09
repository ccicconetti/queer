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

#include <boost/graph/subgraph.hpp>
#include <glog/logging.h>
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
    ret << ", #" << theId << ", allocated to edge " << theEdgeNode
        << ", path length " << thePathLength << ", tot rate " << totRate()
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
    , theEdgeNodes() {
  // noop
}

void MecQkdOnlineNetwork::configure(
    const MecQkdOnlineAlgo&                     aAlgorithm,
    std::unique_ptr<support::RealRvInterface>&& aRv,
    const std::set<unsigned long>&              aUserNodes,
    const std::map<unsigned long, double>&      aEdgeProcessing) {

  const auto V = boost::num_vertices(theGraph);

  // consistency checks
  if ((aUserNodes.size() + aEdgeProcessing.size()) > V) {
    throw std::runtime_error("too many user+edge nodes requested");
  }

  // throw if the caller is trying to reconfigure an operational network
  // XXX

  // set the algoritm and r.v.
  theAlgorithm = aAlgorithm;
  aRv          = std::move(theRv);

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

double MecQkdOnlineNetwork::totProcessing() const {
  return std::accumulate(
      theEdgeProcessing.begin(),
      theEdgeProcessing.end(),
      0.0,
      [](auto aSum, const auto& aElem) { return aSum + aElem.second; });
}

void MecQkdOnlineNetwork::add(Allocation& aApp) {
  VLOG(2) << "request to add new app: " << aApp.toString();
  aApp.theId = theNextId++;
}

void MecQkdOnlineNetwork::del(const uint64_t aAppId) {
  VLOG(2) << "request to del app #" << aAppId;
}

} // namespace qr
} // namespace uiiit
