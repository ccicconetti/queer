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

#include "QuantumRouting/esnetwork.h"
#include "QuantumRouting/networkfactory.h"
#include "QuantumRouting/qrutils.h"
#include "Support/experimentdata.h"
#include "Support/glograii.h"
#include "Support/parallelbatch.h"
#include "Support/queue.h"
#include "Support/random.h"
#include "Support/split.h"
#include "Support/stat.h"
#include "Support/tostring.h"
#include "Support/versionutils.h"

#include <algorithm>
#include <boost/program_options.hpp>

#include <fstream>
#include <glog/logging.h>

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace po = boost::program_options;
namespace qr = uiiit::qr;
namespace us = uiiit::support;

std::string v2s(const std::vector<double>& aValues) {
  return toString(
      aValues, "@", [](const double& x) { return std::to_string(x); });
}

struct Parameters {
  std::size_t theSeed;

  // scenario generation
  double      theMu;
  double      theGridLength;
  double      theThreshold;
  double      theLinkProbability;
  double      theLinkMinEpr;
  double      theLinkMaxEpr;
  std::string theSrcDstPolicy;
  std::string theGraphMlFilename;

  // system
  double theQ;
  double theFidelityInit;
  double theSimDuration;
  double theWarmup;

  // application flows
  double              theArrivalRate;
  double              theFlowDuration;
  std::vector<double> theNetRates;
  std::vector<double> theFidelityThresholds;

  // simulation
  std::string theTopoFilename;

  static const std::vector<std::string>& names() {
    static std::vector<std::string> ret({
        "seed",
        "mu",
        "grid-length",
        "threshold",
        "link-prob",
        "link-min-epr",
        "link-max-epr",
        "src-dst-policy",
        "graphml-filename",
        "q",
        "fidelity-init",
        "sim-duration",
        "warmup-duration",
        "arrival-rate",
        "flow-duration",
        "net-rate",
        "fidelity-thresh",
    });
    return ret;
  }

  std::string toString() const {
    std::stringstream myStream;
    if (theGraphMlFilename.empty()) {
      myStream << "num nodes drawn from PPP with mu " << theMu
               << " distributed on a flat square grid with edge size "
               << theGridLength << " m, a link is generated with probability "
               << theLinkProbability << " between any two nodes within "
               << theThreshold << " m apart";
    } else {
      myStream << "topology read from " << theGraphMlFilename;
    }

    myStream
        << ", src/dst nodes are selected with policy " << theSrcDstPolicy
        << ", and EPR generation rate of the list is drawn randomly from U["
        << theLinkMinEpr << ',' << theLinkMaxEpr << "]; simulation duration "
        << theSimDuration << " (warm-up " << theWarmup
        << "); probability of correct BSM " << theQ
        << " and fidelity of freshly generated pairs " << theFidelityInit
        << "; flows arrive with average rate " << theArrivalRate
        << " and they have an average duration " << theFlowDuration
        << ", minimum fidelity drawn randomly from {"
        << v2s(theFidelityThresholds) << "}"
        << ", and a net requested rate drawn randomly from {"
        << v2s(theNetRates) << "} EPR pairs/s, experiment seed " << theSeed;

    return myStream.str();
  }

  std::string toCsv() const {
    std::stringstream myStream;
    myStream << theSeed << ',' << theMu << ',' << theGridLength << ','
             << theThreshold << ',' << theLinkProbability << ','
             << theLinkMinEpr << ',' << theLinkMaxEpr << ',' << theSrcDstPolicy
             << ',' << theGraphMlFilename << ',' << theQ << ','
             << theFidelityInit << ',' << theSimDuration << ',' << theWarmup
             << ',' << theArrivalRate << ',' << theFlowDuration << ','
             << v2s(theNetRates) << ',' << v2s(theFidelityThresholds);
    return myStream.str();
  }
};

struct Output {
  explicit Output() = default;

  explicit Output(const std::vector<double>& aNetRates,
                  const std::vector<double>& aFidelityThresholds) {
    resize(aNetRates, aFidelityThresholds);
  }

  void resize(const std::vector<double>& aNetRates,
              const std::vector<double>& aFidelityThresholds) {
    if (aNetRates.empty() or aFidelityThresholds.empty()) {
      throw std::runtime_error("invalid empty set of rates or fidelities");
    }
    thePerClass.resize(
        aNetRates.size(),
        std::vector<PerClass>(aFidelityThresholds.size(), PerClass()));
    theNames = std::vector<std::string>({
        "num-nodes",
        "num-edges",
        "min-in-degree",
        "max-in-degree",
        "min-out-degree",
        "max-out-degree",
        "diameter",
        "capacity-tot",

        "capacity-res",
        "num-active-flows",
        "avg-dijkstra-calls",
        "avg-gross-rate",
        "avg-net-rate",
        "admission-rate",
        "avg-path-size",
        "avg-fidelity",
    });
    static const std::vector<std::string> myPerClassNames({
        "avg-gross-rate",
        "avg-net-rate",
        "admission-rate",
        "avg-path-size",
    });
    for (const auto& myPerClassName : myPerClassNames) {
      for (std::size_t i = 0; i < aNetRates.size(); i++) {
        for (std::size_t j = 0; j < aFidelityThresholds.size(); j++) {
          theNames.emplace_back(myPerClassName + "-" +
                                std::to_string(aNetRates[i]) + "-" +
                                std::to_string(aFidelityThresholds[j]));
        }
      }
    }
  }

  // graph properties
  std::size_t theNumNodes      = 0;
  std::size_t theNumEdges      = 0;
  std::size_t theMinInDegree   = 0;
  std::size_t theMaxInDegree   = 0;
  std::size_t theMinOutDegree  = 0;
  std::size_t theMaxOutDegree  = 0;
  std::size_t theDiameter      = 0;
  double      theTotalCapacity = 0;

  // routing properties
  double theResidualCapacity = 0;
  double theNumActiveFlows   = 0;
  double theAvgDijkstraCalls = 0;
  double theGrossRate        = 0;
  double theNetRate          = 0;
  double theAdmissionRate    = 0;
  double theAvgPathSize      = 0;
  double theAvgFidelity      = 0;

  struct PerClass {
    double theGrossRate     = 0;
    double theNetRate       = 0;
    double theAdmissionRate = 0;
    double theAvgPathSize   = 0;
  };
  std::vector<std::vector<PerClass>> thePerClass;

  std::vector<std::string> theNames;

  const std::vector<std::string> names() {
    return theNames;
  }

  std::string toString() const {
    std::stringstream myStream;
    myStream << "G(" << theNumNodes << "," << theNumEdges << "), in-degree "
             << theMinInDegree << "-" << theMaxOutDegree << ", out-degree "
             << theMinOutDegree << "-" << theMaxOutDegree << ", diameter "
             << theDiameter << ", total capacity " << theTotalCapacity
             << " EPR/s (residual " << theResidualCapacity
             << " EPR/s); number of active flows " << theNumActiveFlows
             << ", admission rate " << theAdmissionRate
             << ", average EPR rate net " << theNetRate << " (gross "
             << theGrossRate << "), with " << theAvgDijkstraCalls
             << " Dijkstra calls on average, average path size "
             << theAvgPathSize
             << ", average fidelity of the end-to-end entangled pair "
             << theAvgFidelity;
    return myStream.str();
  }

  std::string toCsv() const {
    std::stringstream myStream;
    myStream << theNumNodes << ',' << theNumEdges << ',' << theMinInDegree
             << ',' << theMaxInDegree << ',' << theMinOutDegree << ','
             << theMaxOutDegree << ',' << theDiameter << ',' << theTotalCapacity
             << ',' << theResidualCapacity << ',' << theNumActiveFlows << ','
             << theAvgDijkstraCalls << ',' << theGrossRate << ',' << theNetRate
             << ',' << theAdmissionRate << ',' << theAvgPathSize << ','
             << theAvgFidelity;
    static const std::vector<std::function<double(const PerClass&)>> myGetters =
        {
            [](const auto& aPerClass) { return aPerClass.theGrossRate; },
            [](const auto& aPerClass) { return aPerClass.theNetRate; },
            [](const auto& aPerClass) { return aPerClass.theAdmissionRate; },
            [](const auto& aPerClass) { return aPerClass.theAvgPathSize; },
        };
    for (const auto& myGetter : myGetters) {
      for (const auto& myPerClass : thePerClass) {
        for (const auto& myPerClassEntry : myPerClass) {
          myStream << ',' << myGetter(myPerClassEntry);
        }
      }
    }
    return myStream.str();
  }
};

using Data = us::ExperimentData<Parameters, Output>;

void runExperiment(Data& aData, Parameters&& aParameters) {
  // fidelity computation parameters
  constexpr double p1  = 1.0;
  constexpr double p2  = 1.0;
  constexpr double eta = 1.0;

  Data::Raii myRaii(aData, std::move(aParameters));

  Output myOutput(myRaii.in().theNetRates, myRaii.in().theFidelityThresholds);

  // consistency checks
  if (myRaii.in().theArrivalRate <= 0) {
    throw std::runtime_error("the arrival rate must be positive: " +
                             std::to_string(myRaii.in().theArrivalRate));
  }
  if (myRaii.in().theFlowDuration <= 0) {
    throw std::runtime_error("the flow duration must be positive: " +
                             std::to_string(myRaii.in().theFlowDuration));
  }

  // open the GraphML file with the network topology, if needed
  const auto myGraphMlStream =
      myRaii.in().theGraphMlFilename.empty() ?
          nullptr :
          std::make_unique<std::ifstream>(myRaii.in().theGraphMlFilename +
                                          ".graphml");
  if (myGraphMlStream.get() != nullptr and
      not static_cast<bool>(*myGraphMlStream)) {
    throw std::runtime_error("cannot read from file: " +
                             myRaii.in().theGraphMlFilename + ".graphml");
  }

  // create network
  us::UniformRv               myLinkEprRv(myRaii.in().theLinkMinEpr,
                            myRaii.in().theLinkMaxEpr,
                            myRaii.in().theSeed,
                            0,
                            0);
  std::vector<qr::Coordinate> myCoordinates;
  const auto myNetwork = myRaii.in().theGraphMlFilename.empty() ?
                             qr::makeCapacityNetworkPpp<qr::EsNetwork>(
                                 myLinkEprRv,
                                 myRaii.in().theSeed,
                                 myRaii.in().theMu,
                                 myRaii.in().theGridLength,
                                 myRaii.in().theThreshold,
                                 myRaii.in().theLinkProbability,
                                 myCoordinates) :
                             qr::makeCapacityNetworkGraphMl<qr::EsNetwork>(
                                 myLinkEprRv, *myGraphMlStream, myCoordinates);
  myNetwork->measurementProbability(myRaii.in().theQ);

  if (not myRaii.in().theTopoFilename.empty()) {
    myNetwork->toGnuplot(myRaii.in().theTopoFilename + "-" +
                             std::to_string(myRaii.in().theSeed),
                         myCoordinates);
  }

  if (myNetwork->numNodes() < 2) {
    throw std::runtime_error("the network must have at least "
                             "two nodes");
  }

  // simulated clock
  double myNow = 0;

  // network properties
  assert(myNetwork.get() != nullptr);
  myOutput.theNumNodes      = myNetwork->numNodes();
  myOutput.theNumEdges      = myNetwork->numEdges();
  myOutput.theTotalCapacity = myNetwork->totalCapacity();
  std::tie(myOutput.theMinInDegree, myOutput.theMaxInDegree) =
      myNetwork->inDegree();
  std::tie(myOutput.theMinOutDegree, myOutput.theMaxOutDegree) =
      myNetwork->outDegree();
  myNetwork->reachableNodes(0, 0, myOutput.theDiameter);

  // prepare statistics data structures
  struct PerClassStat {
    us::SummaryStat theGrossRate;
    us::SummaryStat theNetRate;
    us::SummaryStat theAdmissionRate;
    us::SummaryStat thePathSize;
  };
  us::SummaryWeightedStat myResidualCapacity(myNow, myRaii.in().theWarmup);
  us::SummaryWeightedStat myNumActiveFlows(myNow, myRaii.in().theWarmup);
  us::SummaryStat         myDijkstra;
  us::SummaryStat         myGrossRate;
  us::SummaryStat         myNetRate;
  us::SummaryStat         myAdmissionRate;
  us::SummaryStat         myPathSize;
  us::SummaryStat         myFidelity;
  std::vector<std::vector<std::shared_ptr<PerClassStat>>> myPerClassStats(
      myRaii.in().theNetRates.size(),
      std::vector<std::shared_ptr<PerClassStat>>(
          myRaii.in().theFidelityThresholds.size(), nullptr));
  for (auto& myElems : myPerClassStats) {
    for (auto& myStat : myElems) {
      myStat = std::make_shared<PerClassStat>();
    }
  }
  myResidualCapacity(myNetwork->totalCapacity());
  myNumActiveFlows(0);

  // flows admitted at any time
  struct AdmittedFlow {
    unsigned long              theSrc;
    std::vector<unsigned long> thePath;
    double                     theLeaveTime;
    double                     theGrossRate;
  };
  std::list<AdmittedFlow> myAdmittedFlows;

  // prepare random variables for flow
  // generation
  us::ExponentialRv myArrivalRv(
      myRaii.in().theArrivalRate, myRaii.in().theSeed, 0, 0);
  us::ExponentialRv myDurationRv(
      1.0 / myRaii.in().theFlowDuration, myRaii.in().theSeed, 1, 0);
  us::UniformIntRv<std::size_t> myNetRatesRv(
      0, myRaii.in().theNetRates.size() - 1, myRaii.in().theSeed, 2, 0);
  us::UniformIntRv<std::size_t> myFidelitiesRv(
      0,
      myRaii.in().theFidelityThresholds.size() - 1,
      myRaii.in().theSeed,
      3,
      0);
  us::UniformRv              mySrcDstRv(0, 1, myRaii.in().theSeed, 4, 0);
  std::vector<unsigned long> myNodes(myNetwork->numNodes());
  for (unsigned long i = 0; i < myNodes.size(); ++i) {
    myNodes[i] = i;
  }
  const auto myNodeCapacities = myNetwork->nodeCapacities();
  assert(myNodeCapacities.size() == myNodes.size());

  // run simulation
  double myNextArrival = 0;
  while (myNow <= myRaii.in().theSimDuration) {
    const auto myEarliestLeave =
        std::min_element(myAdmittedFlows.begin(),
                         myAdmittedFlows.end(),
                         [](const auto& aLhs, const auto& aRhs) {
                           return aLhs.theLeaveTime < aRhs.theLeaveTime;
                         });

    if (myEarliestLeave == myAdmittedFlows.end() or
        myNextArrival < myEarliestLeave->theLeaveTime) {
      myNow = myNextArrival;

      std::vector<unsigned long> mySrcDstNodes;
      if (myRaii.in().theSrcDstPolicy == "uniform") {
        mySrcDstNodes = us::sample(myNodes, 2, mySrcDstRv);
      } else if (myRaii.in().theSrcDstPolicy == "nodecapacities") {
        mySrcDstNodes =
            us::sampleWeighted(myNodes, myNodeCapacities, 2, mySrcDstRv);
      } else {
        throw std::runtime_error("unknown src/dst policy: " +
                                 myRaii.in().theSrcDstPolicy);
      }
      assert(mySrcDstNodes.size() == 2);
      assert(mySrcDstNodes[0] != mySrcDstNodes[1]);

      const auto myNetRateId = myNetRatesRv();
      assert(myNetRateId < myRaii.in().theNetRates.size());
      std::vector<qr::EsNetwork::FlowDescriptor> myFlows(
          {{mySrcDstNodes[0],
            mySrcDstNodes[1],
            myRaii.in().theNetRates[myNetRateId]}});

      // try to admit the new traffic flow
      const auto myFidelityThresholdId = myFidelitiesRv();
      assert(myFidelityThresholdId < myRaii.in().theFidelityThresholds.size());
      myNetwork->route(
          myFlows, [&myRaii, myFidelityThresholdId](const auto& aFlow) {
            assert(not aFlow.thePath.empty());
            return qr::fidelitySwapping(p1,
                                        p2,
                                        eta,
                                        aFlow.thePath.size() - 1,
                                        myRaii.in().theFidelityInit) >=
                   myRaii.in().theFidelityThresholds[myFidelityThresholdId];
          });
      assert(myFlows.size() == 1);

      // retrieve the per-class set of
      // statistics
      assert(myNetRateId < myPerClassStats.size());
      assert(myFidelityThresholdId < myPerClassStats[myNetRateId].size());
      auto& myPerClassStat =
          myPerClassStats[myNetRateId][myFidelityThresholdId];
      assert(myPerClassStat.get() != nullptr);

      if (myFlows[0].thePath.empty()) {
        VLOG(2) << "time " << myNow << " dropped  " << myFlows[0].toString()
                << ", fidelity threshold "
                << myRaii.in().theFidelityThresholds[myFidelityThresholdId];

        // record global and per-class
        // statistics
        if (myNow >= myRaii.in().theWarmup) {
          myDijkstra(myFlows[0].theDijsktra);
          myAdmissionRate(0.0);
          myPerClassStat->theAdmissionRate(0.0);
        }

      } else {
        const auto myLeaveTime = myNow + myDurationRv();
        VLOG(2) << "time " << myNow << " admitted " << myFlows[0].toString()
                << ", fidelity threshold "
                << myRaii.in().theFidelityThresholds[myFidelityThresholdId]
                << ", will leave at " << myLeaveTime;
        assert(myFlows[0].theGrossRate > 0);

        myAdmittedFlows.emplace_back(AdmittedFlow{myFlows[0].theSrc,
                                                  myFlows[0].thePath,
                                                  myLeaveTime,
                                                  myFlows[0].theGrossRate});

        // time-weighted statistics
        myResidualCapacity(myNetwork->totalCapacity());
        myNumActiveFlows(myAdmittedFlows.size());

        // time-independent statistics
        if (myNow >= myRaii.in().theWarmup) {
          // global statistics
          myDijkstra(myFlows[0].theDijsktra);
          myGrossRate(myFlows[0].theGrossRate);
          myNetRate(myFlows[0].theNetRate);
          myAdmissionRate(1.0);
          myPathSize(myFlows[0].thePath.size());
          myFidelity(qr::fidelitySwapping(p1,
                                          p2,
                                          eta,
                                          myFlows[0].thePath.size() - 1,
                                          myRaii.in().theFidelityInit));
          // per-class statistics
          myPerClassStat->theGrossRate(myFlows[0].theGrossRate);
          myPerClassStat->theNetRate(myFlows[0].theNetRate);
          myPerClassStat->theAdmissionRate(1.0);
          myPerClassStat->thePathSize(myFlows[0].thePath.size());
        }
      }

      myNextArrival = myNow + myArrivalRv();

    } else {
      assert(myEarliestLeave != myAdmittedFlows.end());
      myNow = myEarliestLeave->theLeaveTime;
      VLOG(3) << "time " << myNow << " an admitted flow leaves";

      // restore the capacity of the path
      myNetwork->addCapacityToPath(myEarliestLeave->theSrc,
                                   myEarliestLeave->thePath,
                                   myEarliestLeave->theGrossRate);

      // remove the flow from the list of
      // admitted/active ones
      myAdmittedFlows.erase(myEarliestLeave);

      // record statistics
      myResidualCapacity(myNetwork->totalCapacity());
      myNumActiveFlows(myAdmittedFlows.size());
    }
  }

  myOutput.theResidualCapacity = myResidualCapacity.mean();
  myOutput.theNumActiveFlows   = myNumActiveFlows.mean();
  myOutput.theAvgDijkstraCalls = myDijkstra.mean();
  myOutput.theGrossRate        = myGrossRate.mean();
  myOutput.theNetRate          = myNetRate.mean();
  myOutput.theAdmissionRate    = myAdmissionRate.mean();
  myOutput.theAvgPathSize      = myPathSize.mean();
  myOutput.theAvgFidelity      = myFidelity.mean();

  for (std::size_t i = 0; i < myPerClassStats.size(); i++) {
    for (std::size_t j = 0; j < myPerClassStats[i].size(); j++) {
      auto& myStat = myPerClassStats[i][j];
      assert(myStat.get() != nullptr);
      assert(i < myOutput.thePerClass.size());
      assert(j < myOutput.thePerClass[i].size());

      myOutput.thePerClass[i][j].theGrossRate = myStat->theGrossRate.mean();
      myOutput.thePerClass[i][j].theNetRate   = myStat->theNetRate.mean();
      myOutput.thePerClass[i][j].theAdmissionRate =
          myStat->theAdmissionRate.mean();
      myOutput.thePerClass[i][j].theAvgPathSize = myStat->thePathSize.mean();
    }
  }

  // save data
  VLOG(1) << "experiment finished\n"
          << myRaii.in().toString() << '\n'
          << myOutput.toString();

  myRaii.finish(std::move(myOutput));
}

bool explainOrPrint(const po::variables_map&   aVarMap,
                    const std::vector<double>& aNetRates,
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
    for (const auto& elem : Output(aNetRates, aFidelityThresholds).names()) {
      std::cout << '#' << ++myCol << '\t' << elem << '\n';
    }
    std::cout << '#' << ++myCol << "\tduration\n";
    return true;
  }

  if (aVarMap.count("print-header") == 1) {
    for (const auto& elem : Parameters::names()) {
      std::cout << elem << ',';
    }
    for (const auto& elem : Output(aNetRates, aFidelityThresholds).names()) {
      std::cout << elem << ',';
    }
    std::cout << "duration\n";
    return true;
  }
  return false;
}

int main(int argc, char* argv[]) {
  uiiit::support::GlogRaii myGlogRaii(argv[0]);

  std::size_t myNumThreads;
  std::string myOutputFilename;
  std::size_t mySeedStart;
  std::size_t mySeedEnd;

  double      myMu;
  double      myLinkMinEpr;
  double      myLinkMaxEpr;
  std::string mySrcDstPolicy;
  std::string myNetRatesStr;
  double      myGridSize;
  double      myThreshold;
  double      myLinkProbability;
  double      myQ;
  double      myFidelityInit;
  std::string myFidelityThresholdsStr;
  std::string myGraphMlFilename;
  double      mySimDuration;
  double      myWarmupDuration;
  double      myArrivalRate;
  double      myFlowDuration;
  std::string myTopoFilename;

  po::options_description myDesc("Allowed options");
  // clang-format off
  myDesc.add_options()
    ("help,h", "produce help message")
    ("version,v", "print the version and quit")
    ("explain-output", "report the meaning of the columns in the output")
    ("print-header", "print the header of the CSV output file")
    ("num-threads",
     po::value<std::size_t>(&myNumThreads)->default_value(1),
     "Number of threads used.")
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
    ("src-dst-policy",
     po::value<std::string>(&mySrcDstPolicy)->default_value("uniform"),
     "Policy to select the src/dst nodes. One of {uniform, nodecapacities}.")
    ("graphml-file",
     po::value<std::string>(&myGraphMlFilename)->default_value(""),
     "Read from a GraphML file.")
    ("sim-duration",
     po::value<double>(&mySimDuration)->default_value(100),
     "Simulation duration, in time units.")
    ("warmup-duration",
     po::value<double>(&myWarmupDuration)->default_value(10),
     "Simulation warm-up duration, in time units.")
    ("arrival-rate",
     po::value<double>(&myArrivalRate)->default_value(1),
     "Average arrival rates of new flows, in time units^-1.")
    ("flow-duration",
     po::value<double>(&myFlowDuration)->default_value(5),
     "Average duration of an admitted flow, in time units.")
    ("net-epr-rates",
     po::value<std::string>(&myNetRatesStr)->default_value("1"),
     "Set of possible net rates requested, in EPR-pairs/s: multiple values separated by @.")
    ("grid-size",
     po::value<double>(&myGridSize)->default_value(60000),
     "Grid length, in km (ignored with using a GraphML file).")
    ("threshold",
     po::value<double>(&myThreshold)->default_value(15000),
     "Link creation threshold (Euclidean distance), in km (ignored with using a GraphML file).")
    ("link-probability",
     po::value<double>(&myLinkProbability)->default_value(1),
     "Link creation probability (ignored with using a GraphML file).")
    ("q",
     po::value<double>(&myQ)->default_value(0.5),
     "Correct measurement probability.")
    ("fidelity-init",
     po::value<double>(&myFidelityInit)->default_value(0.99),
     "Fidelity of local entanglement between adjacent nodes.")
    ("fidelity-threshold",
     po::value<std::string>(&myFidelityThresholdsStr)->default_value("0.95"),
     "Set of possible fidelity thresholds: multiple values separated by @.")
    ("topo-filename",
     po::value<std::string>(&myTopoFilename)->default_value(""),
     "Save the topology to files with the given base name.")
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

    const auto myNetRates = us::split<std::vector<double>>(myNetRatesStr, "@");
    if (myNetRates.empty()) {
      throw std::runtime_error("invalid empty set of net EPR rates");
    }

    const auto myFidelityThresholds =
        us::split<std::vector<double>>(myFidelityThresholdsStr, "@");
    if (myFidelityThresholds.empty()) {
      throw std::runtime_error("invalid empty set of fidelity thresholds");
    }

    if (explainOrPrint(myVarMap, myNetRates, myFidelityThresholds)) {
      return EXIT_SUCCESS;
    }

    if (myGraphMlFilename.find(",") != std::string::npos) {
      throw std::runtime_error("the GraphML file name cannot contain ',': " +
                               myGraphMlFilename);
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
      myParameters.push(Parameters{mySeed,
                                   myMu,
                                   myGridSize,
                                   myThreshold,
                                   myLinkProbability,
                                   myLinkMinEpr,
                                   myLinkMaxEpr,
                                   mySrcDstPolicy,
                                   myGraphMlFilename,
                                   myQ,
                                   myFidelityInit,
                                   mySimDuration,
                                   myWarmupDuration,
                                   myArrivalRate,
                                   myFlowDuration,
                                   myNetRates,
                                   myFidelityThresholds,
                                   myTopoFilename});
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
