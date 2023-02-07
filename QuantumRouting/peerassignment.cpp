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
#include "hungarian-algorithm-cpp/Hungarian.h"

#include <algorithm>
#include <glog/logging.h>

#include <iterator>
#include <numeric>
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
makePeerAssignment(const EsNetwork&                   aNetwork,
                   const PeerAssignmentAlgo           aAlgo,
                   support::RealRvInterface&          aRv,
                   const EsNetwork::AppCheckFunction& aCheckFunction) {
  switch (aAlgo) {
    case PeerAssignmentAlgo::Random:
      return std::make_unique<PeerAssignmentRandom>(aNetwork, aRv);
    case PeerAssignmentAlgo::ShortestPath:
      return std::make_unique<PeerAssignmentShortestPath>(aNetwork, aRv);
    case PeerAssignmentAlgo::LoadBalancing:
      return std::make_unique<PeerAssignmentLoadBalancing>(aNetwork,
                                                           aCheckFunction);
    default:; /* fall-through */
  }
  throw std::runtime_error("invalid peer assignment algorithm: " +
                           toString(aAlgo));
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

PeerAssignment::PeerAssignment(const EsNetwork&         aNetwork,
                               const PeerAssignmentAlgo aAlgo)
    : theNetwork(aNetwork)
    , theAlgo(aAlgo) {
  // noop
}

void PeerAssignment::throwIfDuplicates(
    const std::vector<unsigned long>& aCandidatePeers) {
  std::set<unsigned long> myPeers(aCandidatePeers.begin(),
                                  aCandidatePeers.end());
  if (myPeers.size() != aCandidatePeers.size()) {
    assert(aCandidatePeers.size() > myPeers.size());
    throw std::runtime_error(
        "found " + std::to_string(aCandidatePeers.size() - myPeers.size()) +
        " duplicates among the candidate peers");
  }
}

PeerAssignmentRandom::PeerAssignmentRandom(const EsNetwork&          aNetwork,
                                           support::RealRvInterface& aRv)
    : PeerAssignment(aNetwork, PeerAssignmentAlgo::Random)
    , theRv(aRv) {
  // noop
}

std::vector<EsNetwork::AppDescriptor> PeerAssignmentRandom::assign(
    const std::vector<AppDescriptor>& aApps,
    const unsigned long               aNumPeers,
    const std::vector<unsigned long>& aCandidatePeers) {
  throwIfDuplicates(aCandidatePeers);
  std::vector<EsNetwork::AppDescriptor> ret;
  std::transform(aApps.cbegin(),
                 aApps.cend(),
                 std::back_inserter(ret),
                 [&aCandidatePeers, aNumPeers, this](const auto& aApp) {
                   return EsNetwork::AppDescriptor(
                       aApp.theHost,
                       support::sample(aCandidatePeers, aNumPeers, theRv),
                       aApp.thePriority,
                       aApp.theFidelityThreshold);
                 });
  return ret;
}

PeerAssignmentShortestPath::PeerAssignmentShortestPath(
    const EsNetwork& aNetwork, support::RealRvInterface& aRv)
    : PeerAssignment(aNetwork, PeerAssignmentAlgo::ShortestPath)
    , theRv(aRv) {
  // noop
}

std::vector<EsNetwork::AppDescriptor> PeerAssignmentShortestPath::assign(
    const std::vector<AppDescriptor>& aApps,
    const unsigned long               aNumPeers,
    const std::vector<unsigned long>& aCandidatePeers) {
  throwIfDuplicates(aCandidatePeers);
  std::set<unsigned long>               myDataCenters(aCandidatePeers.begin(),
                                        aCandidatePeers.end());
  std::vector<EsNetwork::AppDescriptor> ret;
  std::transform(aApps.cbegin(),
                 aApps.cend(),
                 std::back_inserter(ret),
                 [&myDataCenters, aNumPeers, this](const auto& aApp) {
                   return EsNetwork::AppDescriptor(
                       aApp.theHost,
                       theNetwork.closestNodes(
                           aApp.theHost, aNumPeers, theRv, myDataCenters),
                       aApp.thePriority,
                       aApp.theFidelityThreshold);
                 });
  return ret;
}

PeerAssignmentLoadBalancing::PeerAssignmentLoadBalancing(
    const EsNetwork&                   aNetwork,
    const EsNetwork::AppCheckFunction& aCheckFunction)
    : PeerAssignment(aNetwork, PeerAssignmentAlgo::LoadBalancing)
    , theCheckFunction(aCheckFunction) {
  // noop
}

std::vector<EsNetwork::AppDescriptor> PeerAssignmentLoadBalancing::assign(
    const std::vector<AppDescriptor>& aApps,
    const unsigned long               aNumPeers,
    const std::vector<unsigned long>& aCandidatePeers) {
  throwIfDuplicates(aCandidatePeers);

  // return immediately if there are no data centers, a.k.a. peer candidates
  if (aCandidatePeers.empty()) {
    return std::vector<EsNetwork::AppDescriptor>();
  }

  // keep the current peer assigned: the inner vectors will grow by
  // one unit at every iteration (#iterations = aNumPeers)
  std::vector<std::vector<unsigned long>> myCurAssign(aApps.size());

  // compute how many columns we need per data center
  const auto C = 1 + (aApps.size() * aNumPeers - 1) / aCandidatePeers.size();
  assert(C > 0);
  VLOG(2) << "num-apps " << aApps.size() << ", num-peers " << aNumPeers
          << ", num-data-centers " << aCandidatePeers.size() << ", C " << C;

  // assignment problem input matrix, with costs from end users to data centers
  hungarian::HungarianAlgorithm::DistMatrix myDistMatrix(
      aApps.size(), std::vector<double>(C * aCandidatePeers.size()));

  // the maximum profit in each row, which will be used to tranform the
  // problem into a cost-minimization
  std::vector<double> myPerRowMaxValues(aApps.size(), 0.0);
  for (unsigned long s = 0; s < aApps.size(); s++) {
    for (unsigned long d = 0; d < aCandidatePeers.size(); d++) {
      const auto myNetRate = theNetwork.maxNetRate(
          EsNetwork::AppDescriptor(aApps[s].theHost, {}, 1, 0.5),
          aCandidatePeers[d],
          theCheckFunction);
      myPerRowMaxValues[s] = std::max(myPerRowMaxValues[s], myNetRate);
      VLOG(2) << s << '/' << aApps.size() << " " << d << '/'
              << (aCandidatePeers.size()) << " net-rate " << myNetRate;
      for (unsigned long c = 0; c < C; c++) {
        myDistMatrix[s][d * C + c] = -myNetRate;
      }
    }
  }
  const auto mySumMax =
      std::accumulate(myPerRowMaxValues.begin(), myPerRowMaxValues.end(), 0.0);

  // transform profit-maximization profit into cost-minimization
  for (unsigned long s = 0; s < aApps.size(); s++) {
    for (unsigned long d = 0; d < (C * aCandidatePeers.size()); d++) {
      myDistMatrix[s][d] += mySumMax;
      assert(myDistMatrix[s][d] >= 0);
    }
  }

  // at every iteration add one peer to each app
  for (unsigned long myIteration = 0; myIteration < aNumPeers; myIteration++) {
    // solve assignment problem
    std::vector<int> myAssignment;
    const auto       myCost =
        hungarian::HungarianAlgorithm::Solve(myDistMatrix, myAssignment);
    const auto myProfit = aApps.size() * mySumMax - myCost;
    assert(myProfit >= 0);
    assert(myAssignment.size() == aApps.size());

    // print the allocation found
    if (VLOG_IS_ON(2)) {
      VLOG(2) << "[" << myIteration << "] found assignment solution with cost "
              << myCost << " (profit " << myProfit << ")";
      for (unsigned long a = 0; a < aApps.size(); a++) {
        VLOG(2) << "\tend user #" << a << "\tsrc " << aApps[a].theHost
                << "\tdst " << aCandidatePeers[myAssignment[a] / C] << "\tcost "
                << myDistMatrix[a][myAssignment[a]] << "\t(profit "
                << (mySumMax - myDistMatrix[a][myAssignment[a]]) << ")";
      }
    }

    for (unsigned long a = 0; a < aApps.size(); a++) {
      // copy the assigment of this iteration into myCurAssign if it is valid
      // and the app was not assigned a data center marked as unavailable
      if (myAssignment[a] >= 0 and
          myDistMatrix[a][myAssignment[a]] < mySumMax) {
        assert((myAssignment[a] / C) < aCandidatePeers.size());
        myCurAssign[a].emplace_back(aCandidatePeers[myAssignment[a] / C]);
      }

      // set maximum cost (mySumMax) in all cells that refer to the same data
      // center that has been assigned to current application a, so that a given
      // node is not assigned twice to the same app
      for (unsigned long c = 0; c < C; c++) {
        myDistMatrix[a][myAssignment[a] / C * C + c] = mySumMax;
      }

      // set maximum cost (mySumMax) in all columns already assigned
      for (unsigned long i = 0; i < aApps.size(); i++) {
        myDistMatrix[i][myAssignment[a]] = mySumMax;
      }
    }
  }

#ifndef NDEBUG
  for (unsigned long a = 0; a < aApps.size(); a++) {
    // no app is assigned more than data centers that requested
    assert(myCurAssign[a].size() <= aNumPeers);

    // no app is assigned twice the same data center
    std::set<unsigned long> myUniq;
    for (const auto d : myCurAssign[a]) {
      assert(myUniq.emplace(d).second);
    }
  }
#endif

  // copy the last assignment to the return value
  std::vector<EsNetwork::AppDescriptor> ret;
  assert(aApps.size() == myCurAssign.size());
  for (unsigned long a = 0; a < aApps.size(); a++) {
    ret.emplace_back(aApps[a].theHost,
                     myCurAssign[a],
                     aApps[a].thePriority,
                     aApps[a].theFidelityThreshold);
  }

  return ret;
}

} // namespace qr
} // namespace uiiit
