/*
              __ __ __
             |__|__|  | __
             |  |  |  ||__|
  ___ ___ __ |  |  |  |
 |   |   |  ||  |  |  |    Ubiquitous Internet @ IIT-CNR
 |   |   |  ||  |  |  |    C++ quantum routing libraries and tools
 |_______|__||__|__|__|    https://github.com/ccicconetti/serverlessonedge

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
#include "QuantumRouting/networkfactory.h"
#include "QuantumRouting/peerassignment.h"
#include "QuantumRouting/qrutils.h"
#include "Support/experimentdata.h"
#include "Support/fairness.h"
#include "Support/glograii.h"
#include "Support/parallelbatch.h"
#include "Support/queue.h"
#include "Support/random.h"
#include "Support/split.h"
#include "Support/stat.h"
#include "Support/tostring.h"
#include "Support/versionutils.h"

#include <boost/program_options.hpp>

#include <cmath>
#include <cstddef>
#include <glog/logging.h>

#include <cstdlib>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>

namespace po = boost::program_options;
namespace qr = uiiit::qr;
namespace us = uiiit::support;

struct Parameters {
  std::size_t theSeed;

  // scenario generation
  double theMu;
  double theGridLength;
  double theThreshold;
  double theLinkProbability;
  double theLinkMinEpr;
  double theLinkMaxEpr;

  // system
  qr::AppRouteAlgo theAlgo;
  double           theQ;
  double           theQuantum;
  std::size_t      theK;
  double           theFidelityInit;
  double           theFracEndUsers;
  double           theFracDataCenters;

  // applications
  std::size_t            theNumApps;
  std::size_t            theNumPeers;
  qr::PeerAssignmentAlgo thePeerAssignmentAlgo;
  std::vector<double>    thePriorities;
  std::vector<double>    theFidelityThresholds;
  double                 theTargetResidual;

  // not part of the experiment
  std::string theDotFile;

  static const std::vector<std::string>& names() {
    static std::vector<std::string> ret({
        "seed",
        "mu",
        "grid-length",
        "threshold",
        "link-prob",
        "link-min-epr",
        "link-max-epr",

        "algorithm",
        "q",
        "quantum",
        "k",
        "fidelity-init",
        "frac-end-users",
        "frac-data-centers",

        "num-apps",
        "num-peers",
        "peer-assign-algo",
        "priorities",
        "fidelity-thresholds",
        "target-residual",
    });
    return ret;
  }

  std::string toString() const {
    std::stringstream myStream;
    myStream << "num nodes drawn from PPP with mu " << theMu
             << " distributed on a flat square grid with edge size "
             << theGridLength << " m, a link is generated with probability "
             << theLinkProbability << " between any two nodes within "
             << theThreshold
             << " m apart, and the EPR generation rate of the list is drawn "
                "randomly from U["
             << theLinkMinEpr << ',' << theLinkMaxEpr
             << "]; resource allocation algorithm " << qr::toString(theAlgo)
             << ", probability of correct BSM " << theQ << ", quantum "
             << theQuantum << " EPR-pairs/s, max " << theK
             << " shortest paths per host/destination pair"
             << ", fidelity of freshly generated pairs " << theFidelityInit
             << "; " << static_cast<int>(std::round(theFracEndUsers * 100))
             << "\% nodes are end-users, "
             << static_cast<int>(std::round(theFracDataCenters * 100))
             << "\% nodes are data centers; there are " << theNumApps
             << " applications, with " << theNumPeers
             << " peers selected according using "
             << qr::toString(thePeerAssignmentAlgo) << ", with priority in {"
             << ::toStringStd(thePriorities, ",")
             << "} and minimum fidelity in {"
             << ::toStringStd(theFidelityThresholds, ",") << "}";
    if (theTargetResidual < 0) {
      myStream << ", no target residual capacity";
    } else {
      myStream << ", target residual capacity " << theTargetResidual
               << " EPR-pairs/s";
    }
    myStream << "; experiment seed " << theSeed;
    return myStream.str();
  }

  std::string toCsv() const {
    std::stringstream myStream;
    myStream << theSeed << ',' << theMu << ',' << theGridLength << ','
             << theThreshold << ',' << theLinkProbability << ','
             << theLinkMinEpr << ',' << theLinkMaxEpr << ','
             << qr::toString(theAlgo) << ',' << theQ << ',' << theQuantum << ','
             << theK << ',' << theFidelityInit << ',' << theFracEndUsers << ','
             << theFracDataCenters << ',' << theNumApps << ',' << theNumPeers
             << ',' << qr::toString(thePeerAssignmentAlgo) << ','
             << ::toStringStd(theFidelityThresholds, "@") << ','
             << ::toStringStd(thePriorities, "@") << ',' << theTargetResidual;
    return myStream.str();
  }
};

struct Output {
  // graph properties
  std::size_t theNumNodes      = 0;
  std::size_t theNumEdges      = 0;
  std::size_t theMinInDegree   = 0;
  std::size_t theMaxInDegree   = 0;
  std::size_t theMinOutDegree  = 0;
  std::size_t theMaxOutDegree  = 0;
  std::size_t theDiameter      = 0;
  double      theTotalCapacity = 0;

  // routing stats
  double theResidualCapacity = 0;
  double theFairnessWpf      = 0; // weighted proportional fairness index
  struct PerClass {
    // conf
    double thePriority          = 0;
    double theFidelityThreshold = 0;

    // stats
    std::size_t theNumApps        = 0;
    double      theAvgVisits      = 0;
    double      theSumGrossRate   = 0;
    double      theSumNetRate     = 0;
    double      theAvgPathSize    = 0;
    double      theAvgFidelity    = 0;
    double      theFairnessJain   = 0; // Jain's fairness index
    double      theFairnessJitter = 0; // max-min rate

    static const std::vector<std::string>& names() {
      static std::vector<std::string> ret({
          "priority",
          "fidelity-threshold",
          "num-apps",
          "avg-visits",
          "sum-gross-rate",
          "sum-net-rate",
          "avg-path-size",
          "avg-fidelity",
          "fairness-jain",
          "fairness-jitter",
      });
      return ret;
    }
  };
  std::vector<PerClass> thePerClassStats;

  static std::vector<std::string>
  names(const std::vector<double>& aPriorities,
        const std::vector<double>& aFidelityThresholds) {
    static std::vector<std::string> myStaticNames({
        // topology properties
        "num-nodes",
        "num-edges",
        "min-in-degree",
        "max-in-degree",
        "min-out-degree",
        "max-out-degree",
        "diameter",
        "capacity-tot",

        // routing stats
        "capacity-res",
        "fairness-wpf",
    });
    std::vector<std::string>        ret(myStaticNames);
    for (const auto& myPerClassName : PerClass::names()) {
      for (const auto& myPriority : aPriorities) {
        for (const auto& myFidelityThreshold : aFidelityThresholds) {
          ret.emplace_back(myPerClassName + "-" + std::to_string(myPriority) +
                           "-" + std::to_string(myFidelityThreshold));
        }
      }
    }
    return ret;
  }

  std::string toString() const {
    std::stringstream myStream;
    myStream << "G(" << theNumNodes << "," << theNumEdges << "), in-degree "
             << theMinInDegree << "-" << theMaxOutDegree << ", out-degree "
             << theMinOutDegree << "-" << theMaxOutDegree << ", diameter "
             << theDiameter << ", total capacity " << theTotalCapacity
             << " EPR-pairs/s; residual capacity " << theResidualCapacity
             << " EPR-pairs/s, sum-log-rates " << theFairnessWpf
             << " EPR-pairs/s; ";
    for (const auto& stat : thePerClassStats) {
      myStream << "app class priority " << stat.thePriority
               << ", fidelity threshold " << stat.theFidelityThreshold << ": "
               << stat.theNumApps << " applications served,"
               << " total EPR rate net " << stat.theSumNetRate << " (gross "
               << stat.theSumGrossRate << ") EPR-pairs/s, " << stat.theAvgVisits
               << " visits on average, average path size "
               << stat.theAvgPathSize
               << ", average fidelity of the end-to-end entangled pair "
               << stat.theAvgFidelity << ", Jain's fairness index "
               << stat.theFairnessJain << ", max rate - min rate "
               << stat.theFairnessJitter << " EPR-pairs/s";
    }
    return myStream.str();
  }

  std::string toCsv() const {
    std::stringstream myStream;
    myStream << theNumNodes << ',' << theNumEdges << ',' << theMinInDegree
             << ',' << theMaxInDegree << ',' << theMinOutDegree << ','
             << theMaxOutDegree << ',' << theDiameter << ',' << theTotalCapacity
             << ',' << theResidualCapacity << ',' << theFairnessWpf;
    struct Printer {
      Printer(const std::vector<PerClass>& aStats)
          : theStats(aStats) {
        // noop
      }
      void
      operator()(std::ostream&                                 aStream,
                 const std::function<double(const PerClass&)>& aGetter) const {
        for (const auto& stat : theStats) {
          aStream << ',' << aGetter(stat);
        }
      }
      const std::vector<PerClass>& theStats;
    };
    const Printer myPrinter(thePerClassStats);
    myPrinter(myStream, [](const auto& elem) { return elem.thePriority; });
    myPrinter(myStream,
              [](const auto& elem) { return elem.theFidelityThreshold; });
    myPrinter(myStream, [](const auto& elem) { return elem.theNumApps; });
    myPrinter(myStream, [](const auto& elem) { return elem.theAvgVisits; });
    myPrinter(myStream, [](const auto& elem) { return elem.theSumGrossRate; });
    myPrinter(myStream, [](const auto& elem) { return elem.theSumNetRate; });
    myPrinter(myStream, [](const auto& elem) { return elem.theAvgPathSize; });
    myPrinter(myStream, [](const auto& elem) { return elem.theAvgFidelity; });
    myPrinter(myStream, [](const auto& elem) { return elem.theFairnessJain; });
    myPrinter(myStream,
              [](const auto& elem) { return elem.theFairnessJitter; });

    return myStream.str();
  }
};

using Data = us::ExperimentData<Parameters, Output>;

void runExperiment(Data& aData, Parameters&& aParameters) {
  if (aParameters.theTargetResidual > 1) {
    throw std::runtime_error(
        "the target residual must be negative or in [0,1]");
  }
  if (aParameters.thePriorities.empty()) {
    throw std::runtime_error("no priority weights specified");
  }
  if (aParameters.theFidelityThresholds.empty()) {
    throw std::runtime_error("no fidelity thresholds specified");
  }

  // fidelity computation parameters
  constexpr double p1  = 1.0;
  constexpr double p2  = 1.0;
  constexpr double eta = 1.0;

  Data::Raii myRaii(aData, std::move(aParameters)); // experiment input
  Output     myOutput;                              // experiment output

  // create network
  us::UniformRv myLinkEprRv(myRaii.in().theLinkMinEpr,
                            myRaii.in().theLinkMaxEpr,
                            myRaii.in().theSeed,
                            1,
                            0);
  [[maybe_unused]] std::vector<qr::Coordinate> myCoordinates;
  auto myNetwork = qr::makeCapacityNetworkPpp(myLinkEprRv,
                                              myRaii.in().theSeed,
                                              myRaii.in().theMu,
                                              myRaii.in().theGridLength,
                                              myRaii.in().theThreshold,
                                              myRaii.in().theLinkProbability,
                                              myCoordinates);
  myNetwork->measurementProbability(myRaii.in().theQ);

  // define validity of a path based on the minimum fidelity
  const auto myCheckFunction = [&myRaii](const auto& aApp, const auto& aPath) {
    assert(not aPath.empty());
    return qr::fidelitySwapping(
               p1, p2, eta, aPath.size() - 1, myRaii.in().theFidelityInit) >=
           aApp.theFidelityThreshold;
  };

  // select the end users (candidate hosts) and data centers (candidate peers)
  struct Rounder {
    Rounder(const unsigned long aN)
        : theN(aN) {
      // noop
    }
    unsigned long operator()(const double aFrac) {
      return std::max((unsigned long)1,
                      static_cast<unsigned long>(std::round(theN * aFrac)));
    }
    const unsigned long theN;
  };
  Rounder    myRounder(myNetwork->numNodes());
  const auto myNumEndUsers    = myRounder(myRaii.in().theFracEndUsers);
  const auto myNumDataCenters = myRounder(myRaii.in().theFracDataCenters);
  if (myNumEndUsers + myNumDataCenters > myNetwork->numNodes()) {
    throw std::runtime_error("the network is too small: there are only " +
                             std::to_string(myNetwork->numNodes()) +
                             " nodes, so it is not possible to have " +
                             std::to_string(myNumEndUsers) + " end users and " +
                             std::to_string(myNumDataCenters) +
                             " data centers");
  }

  std::vector<unsigned long> myCurNodes(myNetwork->numNodes());
  std::iota(myCurNodes.begin(), myCurNodes.end(), 0);
  us::UniformRv myCandidateNodesRv(0, 1, myRaii.in().theSeed, 1, 1);
  auto          myCandidateHosts =
      us::sample(myCurNodes, myNumEndUsers, myCandidateNodesRv);
  std::vector<unsigned long> myCandidateDataCenters;
  for (unsigned long i = 0; i < myNetwork->numNodes(); i++) {
    if (std::find(myCandidateHosts.begin(), myCandidateHosts.end(), i) ==
        myCandidateHosts.end()) {
      myCandidateDataCenters.emplace_back(i);
    }
  }
  myCandidateDataCenters =
      us::sample(myCandidateDataCenters, myNumDataCenters, myCandidateNodesRv);

  // save the network properties
  assert(myNetwork.get() != nullptr);
  myOutput.theNumNodes      = myNetwork->numNodes();
  myOutput.theNumEdges      = myNetwork->numEdges();
  myOutput.theTotalCapacity = myNetwork->totalCapacity();
  std::tie(myOutput.theMinInDegree, myOutput.theMaxInDegree) =
      myNetwork->inDegree();
  std::tie(myOutput.theMinOutDegree, myOutput.theMaxOutDegree) =
      myNetwork->outDegree();

  // save to Graphviz, if needed
  if (not myRaii.in().theDotFile.empty()) {
    myNetwork->toDot(myRaii.in().theDotFile + "-" +
                     std::to_string(myRaii.in().theSeed) + ".dot");
  }

  // skip only if there are no applications at all (dry run)
  if (myRaii.in().theNumApps > 0) {
    // all the combinations of the priorities and the fidelity thresholds
    struct ClassParams {
      double thePriority          = 0;
      double theFidelityThreshold = 0;
    };
    const auto myNumClasses = myRaii.in().thePriorities.size() *
                              myRaii.in().theFidelityThresholds.size();
    std::vector<ClassParams> myClassParams;
    for (const auto& prio : myRaii.in().thePriorities) {
      for (const auto& thresh : myRaii.in().theFidelityThresholds) {
        myClassParams.emplace_back(ClassParams{prio, thresh});
      }
    }
    assert(myClassParams.size() == myNumClasses);

    // create all the random variables related to the applications
    us::UniformIntRv<unsigned long> myHostRv(
        0, myCandidateHosts.size() - 1, myRaii.in().theSeed, 2, 0);
    us::UniformRv myPeerAssignmentRv(0, 1, myRaii.in().theSeed, 2, 1);
    us::UniformIntRv<unsigned long> myClassRv(
        0, myClassParams.size() - 1, myRaii.in().theSeed, 2, 2);

    // create the object for peer assignment
    auto myPeerAssignment =
        qr::makePeerAssignment(*myNetwork,
                               myRaii.in().thePeerAssignmentAlgo,
                               myPeerAssignmentRv,
                               myCheckFunction);
    assert(myPeerAssignment.get() != nullptr);
    assert(myPeerAssignment->algo() == myRaii.in().thePeerAssignmentAlgo);

    // main loop:
    // at each iteration theNumApps applications are created and routed
    // if there is no target residual, the loop terminates after one
    // iteration, otherwise it continues until it is reached
    std::vector<qr::CapacityNetwork::AppDescriptor> myApps;
    do {
      // create applications, with assign QoS class parameters but no peers
      std::vector<qr::PeerAssignment::AppDescriptor> mySingleRunInApps;
      for (std::size_t i = 0; i < myRaii.in().theNumApps; i++) {
        const auto myHost       = myCandidateHosts.at(myHostRv());
        const auto myClassParam = myClassParams.at(myClassRv());
        mySingleRunInApps.emplace_back(myHost,
                                       myClassParam.thePriority,
                                       myClassParam.theFidelityThreshold);
      }

      // assign peers to the applications
      auto mySingleRunApps = myPeerAssignment->assign(
          mySingleRunInApps, myRaii.in().theNumPeers, myCandidateDataCenters);

      // route applications
      us::UniformRv myRouteRv(0, 1, myRaii.in().theSeed, 3, 0);
      myNetwork->route(mySingleRunApps,
                       myRaii.in().theAlgo,
                       myRaii.in().theQuantum * myRaii.in().theNumApps,
                       myRouteRv,
                       myRaii.in().theK,
                       myCheckFunction);

      std::move(mySingleRunApps.begin(),
                mySingleRunApps.end(),
                std::back_inserter(myApps));
    } while (myRaii.in().theTargetResidual >= 0 and
             myNetwork->totalCapacity() >
                 (myOutput.theTotalCapacity * myRaii.in().theTargetResidual));

    // traffic metrics
    struct PerClassSummaryStats {
      std::size_t         theNumApps = 0;
      us::SummaryStat     theVisits;
      us::SummaryStat     theGrossRate;
      us::SummaryStat     theNetRate;
      us::SummaryStat     theAdmissionRate;
      us::SummaryStat     thePathSize;
      us::SummaryStat     theFidelity;
      std::vector<double> theNetRates;
    };
    std::vector<PerClassSummaryStats> myPerClassSummaryStats(myNumClasses);

    // find class index by priority and fidelity threshold
    struct ClassFinder {
      ClassFinder(const std::vector<ClassParams>& aClassParams)
          : theClassParams(aClassParams) {
        // noop
      }
      std::size_t operator()(const double aPriority,
                             const double aFidelityThreshold) {
        for (std::size_t i = 0; i < theClassParams.size(); i++) {
          if (aPriority == theClassParams[i].thePriority and
              aFidelityThreshold == theClassParams[i].theFidelityThreshold) {
            return i;
          }
        }
        throw std::runtime_error(
            "could not find priority " + std::to_string(aPriority) +
            ", fidelity threshold " + std::to_string(aFidelityThreshold));
      }
      const std::vector<ClassParams>& theClassParams;
    };
    ClassFinder myClassFinder(myClassParams);

    // update statistics accumulators for all applications
    std::vector<double> myPerAppRates;
    std::vector<double> myPerAppWeights;
    myPerAppRates.reserve(myApps.size());
    myPerAppWeights.reserve(myApps.size());
    for (const auto& myApp : myApps) {
      auto& stat = myPerClassSummaryStats[myClassFinder(
          myApp.thePriority, myApp.theFidelityThreshold)];

      stat.theNumApps++;

      const auto myHostNetRate = myApp.netRate();
      stat.theNetRates.emplace_back(myHostNetRate);
      stat.theVisits(myApp.theVisits);
      stat.theGrossRate(myApp.grossRate());
      stat.theNetRate(myHostNetRate);
      myPerAppRates.emplace_back(myHostNetRate);
      myPerAppWeights.emplace_back(myApp.thePriority);

      us::SummaryStat myHostPathSize;
      us::SummaryStat myHostFidelity;
      if (myHostNetRate > 0) {
        assert(std::isnormal(myHostNetRate));
        for (const auto& myAllocation : myApp.theAllocated) {
          for (const auto& myPeer : myAllocation.second) {
            const auto myWeight = myPeer.theNetRate / myHostNetRate;
            myHostPathSize(myWeight * myPeer.theHops.size());
            myHostFidelity(myWeight *
                           qr::fidelitySwapping(p1,
                                                p2,
                                                eta,
                                                myPeer.theHops.size() - 1,
                                                myRaii.in().theFidelityInit));
          }
        }
        stat.thePathSize(myHostPathSize.mean() * myHostPathSize.count());
        stat.theFidelity(myHostFidelity.mean() * myHostFidelity.count());
      }
    }

    // consistency checks
    assert(myApps.size() == std::accumulate(myPerClassSummaryStats.begin(),
                                            myPerClassSummaryStats.end(),
                                            std::size_t(0),
                                            [](auto aSum, const auto& aElem) {
                                              return aSum + aElem.theNumApps;
                                            }));
    assert(myApps.size() == myPerAppRates.size());
    assert(myApps.size() == myPerAppWeights.size());
    assert(myOutput.thePerClassStats.empty());
    assert(myPerClassSummaryStats.size() == myClassParams.size());

    // write down output statistics about traffic
    myOutput.theResidualCapacity = myNetwork->totalCapacity();
    myOutput.theFairnessWpf =
        us::proportionalFairnessIndex(myPerAppRates, myPerAppWeights);
    for (std::size_t i = 0; i < myPerClassSummaryStats.size(); i++) {
      auto& stats = myPerClassSummaryStats[i];
      myOutput.thePerClassStats.emplace_back(Output::PerClass{
          myClassParams[i].thePriority,
          myClassParams[i].theFidelityThreshold,
          stats.theNumApps,
          stats.theVisits.mean(),
          stats.theGrossRate.count() * stats.theGrossRate.mean(),
          stats.theNetRate.count() * stats.theNetRate.mean(),
          stats.thePathSize.mean(),
          stats.theFidelity.mean(),
          stats.theNetRates.empty() ? -1 :
                                      us::jainFairnessIndex(stats.theNetRates),
          stats.theNetRates.empty() ?
              -1 :
              (stats.theNetRate.max() - stats.theNetRate.min())});
    }
  }

  // save data
  VLOG(1) << "experiment finished\n"
          << myRaii.in().toString() << '\n'
          << myOutput.toString();

  myRaii.finish(std::move(myOutput));
}

bool explainOrPrint(const po::variables_map&   aVarMap,
                    const std::vector<double>& aPriorities,
                    const std::vector<double>& aFidelityThresholds) {
  if (aVarMap.count("explain-output") == 1 and
      aVarMap.count("print-header") == 1) {
    throw std::runtime_error(
        "cannot specify both --explain-output and --print-header");
  }
  if (aVarMap.count("explain-output") == 1) {
    std::size_t myCol = 0;
    for (const auto& elem : Parameters::names()) {
      std::cout << '#' << ++myCol << '\t' << elem << '\n';
    }
    for (const auto& elem : Output::names(aPriorities, aFidelityThresholds)) {
      std::cout << '#' << ++myCol << '\t' << elem << '\n';
    }
    std::cout << '#' << ++myCol << "\tduration\n";
    return true;
  }

  if (aVarMap.count("print-header") == 1) {
    for (const auto& elem : Parameters::names()) {
      std::cout << elem << ',';
    }
    for (const auto& elem : Output::names(aPriorities, aFidelityThresholds)) {
      std::cout << elem << ',';
    }
    std::cout << "duration\n";
    return true;
  }
  return false;
}

int main(int argc, char* argv[]) {
  us::GlogRaii myGlogRaii(argv[0]);

  std::size_t myNumThreads;
  std::string myOutputFilename;
  std::size_t mySeedStart;
  std::size_t mySeedEnd;

  double      myMu;
  double      myLinkMinEpr;
  double      myLinkMaxEpr;
  std::size_t myNumApps;
  double      myGridSize;
  double      myThreshold;
  double      myLinkProbability;
  std::string myAlgorithm;
  double      myQ;
  double      myQuantum;
  std::size_t myK;
  double      myFidelityInit;
  double      myFracEndUsers;
  double      myFracDataCenters;
  std::string myPrioritiesStr;
  std::string myFidelityThresholdStr;
  std::size_t myNumPeers;
  std::string myPeerAssignmentAlgo;
  double      myTargetResidual;
  std::string myDotFile;

  po::options_description myDesc("Allowed options");
  // clang-format off
  myDesc.add_options()
    ("help,h", "produce help message")
    ("version,v", "print the version and quit")
    ("explain-output", "report the meaning of the columns in the output")
    ("print-header", "print the header of the CSV output file")
    ("dot-file",
     po::value<std::string>(&myDotFile)->default_value(""),
     "Save the network to this Graphviz file.")
    ("num-threads",
     po::value<std::size_t>(&myNumThreads)->default_value(1),
     "Number of threads used. If 0, then use the hardware concurrency value.")
    ("output",
     po::value<std::string>(&myOutputFilename)->default_value("output.csv"),
     "Output file name.")
    ("seed-start",
     po::value<std::size_t>(&mySeedStart)->default_value(0),
     "First seed used.")
    ("seed-end",
     po::value<std::size_t>(&mySeedEnd)->default_value(1),
    "Next seed after the last one to be used, i.e., the number of simulations is (seed-end - seed-start).")
    ("append", "Append to the output file.")
    ("mu",
     po::value<double>(&myMu)->default_value(100),
     "Average number of nodes.")
    ("link-min-epr",
     po::value<double>(&myLinkMinEpr)->default_value(1),
     "Min EPR rate of links.")
    ("link-max-epr",
     po::value<double>(&myLinkMaxEpr)->default_value(400),
     "Max EPR rate of links.")
    ("num-apps",
     po::value<std::size_t>(&myNumApps)->default_value(100),
     "Number of applications added at each iteration.")
    ("num-peers",
     po::value<std::size_t>(&myNumPeers)->default_value(1),
     "Number of peers per app.")
    ("peer-assignment-algo",
     po::value<std::string>(&myPeerAssignmentAlgo)->default_value(qr::toString(qr::PeerAssignmentAlgo::LoadBalancing)),
     (std::string("Algorithm to be used, one of: ") + toString(qr::allPeerAssignmentAlgos(), ", ", [](const auto& aAlgo) { return toString(aAlgo); })).c_str())
    ("grid-size",
     po::value<double>(&myGridSize)->default_value(60000),
     "Grid length, in km.")
    ("threshold",
     po::value<double>(&myThreshold)->default_value(15000),
     "Link creation threshold (Euclidean distance), in km.")
    ("link-probability",
     po::value<double>(&myLinkProbability)->default_value(1),
     "Link creation probability.")
    ("algorithm",
     po::value<std::string>(&myAlgorithm)->default_value(qr::toString(qr::AppRouteAlgo::Drr)),
     (std::string("Algorithm to be used, one of: ") + toString(qr::allAppRouteAlgos(), ", ", [](const auto& aAlgo) { return toString(aAlgo); })).c_str())
    ("q",
     po::value<double>(&myQ)->default_value(0.5),
     "Correct measurement probability.")
    ("quantum",
     po::value<double>(&myQuantum)->default_value(1),
     "The scheduler quantum, in EPR-pairs/s.")
    ("k",
     po::value<std::size_t>(&myK)->default_value(4),
     "Maximum number of shortest paths per host/peer.")
    ("fidelity-init",
     po::value<double>(&myFidelityInit)->default_value(0.99),
     "Fidelity of local entanglement between adjacent nodes.")
    ("frac-end-users",
     po::value<double>(&myFracEndUsers)->default_value(0.2),
     "Fraction of end users, i.e., nodes that are candidate hosts.")
    ("frac-data-centers",
     po::value<double>(&myFracDataCenters)->default_value(0.2),
     "Fraction of data centers, i.e., nodes that are candidate peers.")
    ("priorities",
     po::value<std::string>(&myPrioritiesStr)->default_value("4,2,1"),
     "Priorities, as comma-separated values.")
    ("fidelity-thresholds",
     po::value<std::string>(&myFidelityThresholdStr)->default_value("0.95,0.75"),
     "Fidelity thresholds, as comma-separated values.")
    ("target-residual",
     po::value<double>(&myTargetResidual)->default_value(-1),
     "Continue adding applications until this residual is reached, in fraction of the total capacity. A negative value means: single iteration.")
    ;
  // clang-format on

  try {
    po::variables_map myVarMap;
    po::store(po::parse_command_line(argc, argv, myDesc), myVarMap);
    po::notify(myVarMap);

    if (myVarMap.count("help")) {
      std::cout << myDesc << std::endl;
      return EXIT_FAILURE;
    }

    if (myVarMap.count("version")) {
      std::cout << us::version() << std::endl;
      return EXIT_SUCCESS;
    }

    const auto myPriorities =
        us::split<std::vector<double>>(myPrioritiesStr, ",");
    const auto myFidelityThresholds =
        us::split<std::vector<double>>(myFidelityThresholdStr, ",");

    if (explainOrPrint(myVarMap, myPriorities, myFidelityThresholds)) {
      return EXIT_SUCCESS;
    }

    if (myNumThreads == 0) {
      myNumThreads = std::thread::hardware_concurrency();
      VLOG(1) << "using " << myNumThreads << " threads";
    }

    std::ofstream myFile(myOutputFilename,
                         myVarMap.count("append") == 1 ? std::ios::app :
                                                         std::ios::trunc);
    if (not myFile) {
      throw std::runtime_error("could not open output file for writing: " +
                               myOutputFilename);
    }

    Data myData;

    us::Queue<Parameters> myParameters;
    for (auto mySeed = mySeedStart; mySeed < mySeedEnd; ++mySeed) {
      myParameters.push(
          Parameters{mySeed,
                     myMu,
                     myGridSize,
                     myThreshold,
                     myLinkProbability,
                     myLinkMinEpr,
                     myLinkMaxEpr,
                     qr::appRouteAlgofromString(myAlgorithm),
                     myQ,
                     myQuantum,
                     myK,
                     myFidelityInit,
                     myFracEndUsers,
                     myFracDataCenters,
                     myNumApps,
                     myNumPeers,
                     qr::peerAssignmentAlgofromString(myPeerAssignmentAlgo),
                     myPriorities,
                     myFidelityThresholds,
                     myTargetResidual,
                     myDotFile});
    }
    us::ParallelBatch<Parameters> myWorkers(
        myNumThreads, myParameters, [&myData](auto&& aParameters) {
          runExperiment(myData, std::move(aParameters));
        });
    const auto myExceptions = myWorkers.wait();
    LOG_IF(ERROR, not myExceptions.empty()) << "there were exceptions:";
    for (const auto& myException : myExceptions) {
      LOG(ERROR) << myException;
    }

    myData.toCsv(myFile);

    return EXIT_SUCCESS;
  } catch (const std::exception& aErr) {
    std::cerr << "Exception caught: " << aErr.what() << std::endl;
  } catch (...) {
    std::cerr << "Unknown exception caught" << std::endl;
  }

  return EXIT_FAILURE;
}
