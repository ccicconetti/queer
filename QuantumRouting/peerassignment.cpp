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

#include "QuantumRouting/peerassignment.h"

#include "Support/tostring.h"

#include <stdexcept>

namespace uiiit {
namespace qr {

std::vector<PeerAssignmentAlgo> allPeerAssignmentAlgos() {
  static const std::vector<PeerAssignmentAlgo> myAlgos({
      PeerAssignmentAlgo::Random,
      PeerAssignmentAlgo::ShortestPath,
      PeerAssignmentAlgo::LoadBalancing,
  });
  return myAlgos;
}

std::string toString(const PeerAssignmentAlgo aAlgo) {
  switch (aAlgo) {
    case PeerAssignmentAlgo::Random:
      return "random";
    case PeerAssignmentAlgo::ShortestPath:
      return "shortest-path";
    case PeerAssignmentAlgo::LoadBalancing:
      return "load-balancing";
    default:; /* fall-through */
  }
  return "unknown";
}

PeerAssignmentAlgo peerAssignmentAlgofromString(const std::string& aAlgo) {
  if (aAlgo == "random") {
    return PeerAssignmentAlgo::Random;
  } else if (aAlgo == "shortest-path") {
    return PeerAssignmentAlgo::ShortestPath;
  } else if (aAlgo == "load-balancing") {
    return PeerAssignmentAlgo::LoadBalancing;
  }
  throw std::runtime_error(
      "invalid peer assignment algorithm: " + aAlgo + " (valid options are: " +
      ::toString(allPeerAssignmentAlgos(),
                 ",",
                 [](const auto& aAlgo) { return toString(aAlgo); }) +
      ")");
}

std::unique_ptr<PeerAssignment>
makePeerAssignment(const CapacityNetwork&    aNetwork,
                   const PeerAssignmentAlgo  aAlgo,
                   support::RealRvInterface& aRv) {
  switch (aAlgo) {
    case PeerAssignmentAlgo::Random:
      return std::make_unique<PeerAssignmentRandom>(aNetwork, aRv);
    case PeerAssignmentAlgo::ShortestPath:
      return std::make_unique<PeerAssignmentShortestPath>(aNetwork, aRv);
    case PeerAssignmentAlgo::LoadBalancing:
      return std::make_unique<PeerAssignmentLoadBalancing>(aNetwork);
    default:; /* fall-through */
  }
  throw std::runtime_error("invalid peer assignment algorithm: " +
                           toString(aAlgo));
}

PeerAssignment::PeerAssignment(const CapacityNetwork&   aNetwork,
                               const PeerAssignmentAlgo aAlgo)
    : theNetwork(aNetwork)
    , theAlgo(aAlgo) {
  // noop
}

PeerAssignment::AppDescriptor::AppDescriptor(
    const unsigned long aHost,
    const double        aPriority,
    const double        aFidelityThreshold) noexcept
    : theHost(aHost)
    , thePriority(aPriority)
    , theFidelityThreshold(aFidelityThreshold) {
  // noop
}

PeerAssignmentRandom::PeerAssignmentRandom(const CapacityNetwork&    aNetwork,
                                           support::RealRvInterface& aRv)
    : PeerAssignment(aNetwork, PeerAssignmentAlgo::Random)
    , theRv(aRv) {
  // noop
}

std::vector<CapacityNetwork::AppDescriptor> PeerAssignmentRandom::assign(
    const std::vector<AppDescriptor>& aApps,
    const unsigned long               aNumPeers,
    const std::vector<unsigned long>& aCandidatePeers) {
  std::vector<CapacityNetwork::AppDescriptor> ret;
  for (const auto& myApp : aApps) {
    ret.emplace_back(myApp.theHost,
                     support::sample(aCandidatePeers, aNumPeers, theRv),
                     myApp.thePriority,
                     myApp.theFidelityThreshold);
  }
  return ret;
}

PeerAssignmentShortestPath::PeerAssignmentShortestPath(
    const CapacityNetwork& aNetwork, support::RealRvInterface& aRv)
    : PeerAssignment(aNetwork, PeerAssignmentAlgo::ShortestPath)
    , theRv(aRv) {
  // noop
}

std::vector<CapacityNetwork::AppDescriptor> PeerAssignmentShortestPath::assign(
    const std::vector<AppDescriptor>& aApps,
    const unsigned long               aNumPeers,
    const std::vector<unsigned long>& aCandidatePeers) {
  std::set<unsigned long> myDataCenters(aCandidatePeers.begin(),
                                        aCandidatePeers.end());
  std::vector<CapacityNetwork::AppDescriptor> ret;
  for (const auto& myApp : aApps) {
    ret.emplace_back(
        myApp.theHost,
        theNetwork.closestNodes(myApp.theHost, aNumPeers, theRv, myDataCenters),
        myApp.thePriority,
        myApp.theFidelityThreshold);
  }
  return ret;
}

PeerAssignmentLoadBalancing::PeerAssignmentLoadBalancing(
    const CapacityNetwork& aNetwork)
    : PeerAssignment(aNetwork, PeerAssignmentAlgo::LoadBalancing) {
  // noop
}

std::vector<CapacityNetwork::AppDescriptor> PeerAssignmentLoadBalancing::assign(
    [[maybe_unused]] const std::vector<AppDescriptor>& aApps,
    [[maybe_unused]] const unsigned long               aNumPeers,
    [[maybe_unused]] const std::vector<unsigned long>& aCandidatePeers) {
  std::vector<CapacityNetwork::AppDescriptor> ret;
  // XXX
  return ret;
}

} // namespace qr
} // namespace uiiit
