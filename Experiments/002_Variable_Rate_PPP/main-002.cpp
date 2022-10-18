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
#include "QuantumRouting/qrutils.h"
#include "Support/experimentdata.h"
#include "Support/fairness.h"
#include "Support/glograii.h"
#include "Support/parallelbatch.h"
#include "Support/queue.h"
#include "Support/random.h"
#include "Support/stat.h"
#include "Support/versionutils.h"

#include <boost/program_options.hpp>

#include <cmath>
#include <cstddef>
#include <glog/logging.h>

#include <cstdlib>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>

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
  double      theQ;
  double      theQuantum;
  std::size_t theK;
  double      theFidelityInit;

  // applications
  std::size_t theNumApps;
  std::size_t theNumPeersMin;
  std::size_t theNumPeersMax;
  std::size_t theDistanceMin;
  std::size_t theDistanceMax;
  double      theFidelityThreshold;
  double      theTargetResidual;

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

        "q",
        "quantum",
        "k",
        "fidelity-init",

        "num-apps",
        "num-peers-min",
        "num-peers-max",
        "distance-min",
        "distance-max",
        "fidelity-thresh",
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
             << "]; probability of correct BSM " << theQ << ", quantum "
             << theQuantum << " EPR-pairs/s, max " << theK
             << " shortest paths per host/destination pair"
             << ", fidelity of freshly generated pairs " << theFidelityInit
             << "; there are " << theNumApps
             << " applications, with number of peers drawn from U["
             << theNumPeersMin << "," << theNumPeersMax
             << "] from nodes with distance drawn from U[" << theDistanceMin
             << "," << theDistanceMax << "], with minimum fidelity "
             << theFidelityThreshold;
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
             << theLinkMinEpr << ',' << theLinkMaxEpr << ',' << theQ << ','
             << theQuantum << ',' << theK << ',' << theFidelityInit << ','
             << theNumApps << ',' << theNumPeersMin << ',' << theNumPeersMax
             << ',' << theDistanceMin << ',' << theDistanceMax << ','
             << theFidelityThreshold << ',' << theTargetResidual;
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

  // routing properties
  double      theResidualCapacity = 0;
  std::size_t theNumApps          = 0;
  double      theAvgVisits        = 0;
  double      theSumGrossRate     = 0;
  double      theSumNetRate       = 0;
  double      theAvgPathSize      = 0;
  double      theAvgFidelity      = 0;
  double      theFairnessJain     = 0;
  double      theFairnessJitter   = 0;

  static const std::vector<std::string>& names() {
    static std::vector<std::string> ret({
        "num-nodes",
        "num-edges",
        "min-in-degree",
        "max-in-degree",
        "min-out-degree",
        "max-out-degree",
        "diameter",
        "capacity-tot",

        "capacity-res",
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

  std::string toString() const {
    std::stringstream myStream;
    myStream << "G(" << theNumNodes << "," << theNumEdges << "), in-degree "
             << theMinInDegree << "-" << theMaxOutDegree << ", out-degree "
             << theMinOutDegree << "-" << theMaxOutDegree << ", diameter "
             << theDiameter << ", total capacity " << theTotalCapacity
             << " EPR-pairs/s (residual " << theResidualCapacity
             << " EPR-pairs/s); " << theNumApps
             << " applications served, total EPR rate net " << theSumNetRate
             << " (gross " << theSumGrossRate << ") EPR-pairs/s, "
             << theAvgVisits << " visits on average, average path size "
             << theAvgPathSize
             << ", average fidelity of the end-to-end entangled pair "
             << theAvgFidelity << ", Jain's fairness index " << theFairnessJain
             << ", max rate - min rate " << theFairnessJitter << " EPR-pairs/s";
    return myStream.str();
  }

  std::string toCsv() const {
    std::stringstream myStream;
    myStream << theNumNodes << ',' << theNumEdges << ',' << theMinInDegree
             << ',' << theMaxInDegree << ',' << theMinOutDegree << ','
             << theMaxOutDegree << ',' << theDiameter << ',' << theTotalCapacity
             << ',' << theResidualCapacity << ',' << theNumApps << ','
             << theAvgVisits << ',' << theSumGrossRate << ',' << theSumNetRate
             << ',' << theAvgPathSize << ',' << theAvgFidelity << ','
             << theFairnessJain << ',' << theFairnessJitter;
    return myStream.str();
  }
};

using Data = us::ExperimentData<Parameters, Output>;

void runExperiment(Data& aData, Parameters&& aParameters) {
  if (aParameters.theTargetResidual > 1) {
    throw std::runtime_error(
        "the target residual must be negative or in [0,1]");
  }

  // fidelity computation parameters
  constexpr double p1  = 1.0;
  constexpr double p2  = 1.0;
  constexpr double eta = 1.0;

  Data::Raii myRaii(aData, std::move(aParameters));

  Output myOutput;

  const auto                           MANY_TRIES = 1000000u;
  std::unique_ptr<qr::CapacityNetwork> myNetwork  = nullptr;
  qr::CapacityNetwork::ReachableNodes  myReachableNodes;
  us::UniformRv                        myLinkEprRv(myRaii.in().theLinkMinEpr,
                            myRaii.in().theLinkMaxEpr,
                            myRaii.in().theSeed,
                            0,
                            0);

  for (std::size_t mySeedOffset = 0;
       myNetwork.get() == nullptr and mySeedOffset < (MANY_TRIES * 10000);
       mySeedOffset += 10000) {
    const auto mySeed = myRaii.in().theSeed + mySeedOffset;
    // create network
    [[maybe_unused]] std::vector<qr::Coordinate> myCoordinates;
    myNetwork = qr::makeCapacityNetworkPpp(myLinkEprRv,
                                           mySeed,
                                           myRaii.in().theMu,
                                           myRaii.in().theGridLength,
                                           myRaii.in().theThreshold,
                                           myRaii.in().theLinkProbability,
                                           myCoordinates);

    // network properties
    assert(myNetwork.get() != nullptr);
    myOutput.theNumNodes      = myNetwork->numNodes();
    myOutput.theNumEdges      = myNetwork->numEdges();
    myOutput.theTotalCapacity = myNetwork->totalCapacity();
    std::tie(myOutput.theMinInDegree, myOutput.theMaxInDegree) =
        myNetwork->inDegree();
    std::tie(myOutput.theMinOutDegree, myOutput.theMaxOutDegree) =
        myNetwork->outDegree();

    // create applications
    myReachableNodes = myNetwork->reachableNodes(myRaii.in().theDistanceMin,
                                                 myRaii.in().theDistanceMax,
                                                 myOutput.theDiameter);
    const std::size_t myNumPossibleHosts =
        std::count_if(myReachableNodes.begin(),
                      myReachableNodes.end(),
                      [](const auto& elem) { return not elem.second.empty(); });
    if (myNumPossibleHosts == 0) {
      VLOG(1) << "graph does not have possible hosts (seed " << mySeed
              << "), trying again";
      myNetwork.reset();
    }
  }
  if (myNetwork.get() == nullptr) {
    throw std::runtime_error("Could not find a connected network after " +
                             std::to_string(MANY_TRIES) + " tries");
  }
  myNetwork->measurementProbability(myRaii.in().theQ);

  // save to Graphviz
  if (not myRaii.in().theDotFile.empty()) {
    myNetwork->toDot(myRaii.in().theDotFile + "-" +
                     std::to_string(myRaii.in().theSeed) + ".dot");
  }

  std::vector<qr::CapacityNetwork::AppDescriptor> myApps;

  if (myRaii.in().theNumApps > 0) {
    us::UniformIntRv<unsigned long> myHostRv(
        0, myNetwork->numNodes() - 1, myRaii.in().theSeed, 0, 0);
    us::UniformIntRv<unsigned long> myNumPeersRv(myRaii.in().theNumPeersMin,
                                                 myRaii.in().theNumPeersMax,
                                                 myRaii.in().theSeed,
                                                 0,
                                                 0);
    us::UniformRv myPeerSampleRv(0, 1, myRaii.in().theSeed, 0, 0);

    do {
      std::vector<qr::CapacityNetwork::AppDescriptor> mySingleRunApps;

      for (std::size_t i = 0; i < myRaii.in().theNumApps; i++) {
        const auto myHost = myHostRv();
        const auto it     = myReachableNodes.find(myHost);
        assert(it != myReachableNodes.end());
        const auto myPeersSet =
            us::sample(it->second, myNumPeersRv(), myPeerSampleRv);
        std::vector<unsigned long> myPeersVector;
        std::copy(myPeersSet.begin(),
                  myPeersSet.end(),
                  std::back_inserter(myPeersVector));
        const auto myPriority = 1.0;
        mySingleRunApps.emplace_back(myHost,
                                     myPeersVector,
                                     myPriority,
                                     myRaii.in().theFidelityThreshold);
      }

      // route applications
      us::UniformRv myRouteRv(0, 1, myRaii.in().theSeed, 3, 0); // unused
      myNetwork->route(mySingleRunApps,
                       qr::AppRouteAlgo::Drr,
                       myRaii.in().theQuantum * myRaii.in().theNumApps,
                       myRouteRv,
                       myRaii.in().theK,
                       [&myRaii](const auto&, const auto& aPath) {
                         assert(not aPath.empty());
                         return qr::fidelitySwapping(
                                    p1,
                                    p2,
                                    eta,
                                    aPath.size() - 1,
                                    myRaii.in().theFidelityInit) >=
                                myRaii.in().theFidelityThreshold;
                       });

      std::move(mySingleRunApps.begin(),
                mySingleRunApps.end(),
                std::back_inserter(myApps));
    } while (myRaii.in().theTargetResidual >= 0 and
             myNetwork->totalCapacity() >
                 (myOutput.theTotalCapacity * myRaii.in().theTargetResidual));

    // traffic metrics
    myOutput.theResidualCapacity = myNetwork->totalCapacity();
    myOutput.theNumApps          = myApps.size();
    us::SummaryStat     myVisits;
    us::SummaryStat     myGrossRate;
    us::SummaryStat     myNetRate;
    us::SummaryStat     myAdmissionRate;
    us::SummaryStat     myPathSize;
    us::SummaryStat     myFidelity;
    std::vector<double> myNetRates;
    for (const auto& myApp : myApps) {
      const auto myHostNetRate = myApp.netRate();
      myNetRates.emplace_back(myHostNetRate);
      myVisits(myApp.theVisits);
      myGrossRate(myApp.grossRate());
      myNetRate(myHostNetRate);

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
        myPathSize(myHostPathSize.mean() * myHostPathSize.count());
        myFidelity(myHostFidelity.mean() * myHostFidelity.count());
      }
    }
    myOutput.theAvgVisits      = myVisits.mean();
    myOutput.theSumGrossRate   = myGrossRate.count() * myGrossRate.mean();
    myOutput.theSumNetRate     = myNetRate.count() * myNetRate.mean();
    myOutput.theAvgPathSize    = myPathSize.mean();
    myOutput.theAvgFidelity    = myFidelity.mean();
    myOutput.theFairnessJain   = us::jainFairnessIndex(myNetRates);
    myOutput.theFairnessJitter = myNetRate.max() - myNetRate.min();
  }

  // save data
  VLOG(1) << "experiment finished\n"
          << myRaii.in().toString() << '\n'
          << myOutput.toString();

  myRaii.finish(std::move(myOutput));
}

bool explainOrPrint(const po::variables_map& aVarMap) {
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
    for (const auto& elem : Output::names()) {
      std::cout << '#' << ++myCol << '\t' << elem << '\n';
    }
    std::cout << '#' << ++myCol << "\tduration\n";
    return true;
  }

  if (aVarMap.count("print-header") == 1) {
    for (const auto& elem : Parameters::names()) {
      std::cout << elem << ',';
    }
    for (const auto& elem : Output::names()) {
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
  double      myQ;
  double      myQuantum;
  std::size_t myK;
  double      myFidelityInit;
  double      myFidelityThreshold;
  std::size_t myNumPeersMin;
  std::size_t myNumPeersMax;
  std::size_t myDistanceMin;
  std::size_t myDistanceMax;
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
    ("num-apps",
     po::value<std::size_t>(&myNumApps)->default_value(100),
     "Number of applications added at each iteration.")
    ("num-peers-min",
     po::value<std::size_t>(&myNumPeersMin)->default_value(1),
     "Minimum number of peers per app.")
    ("num-peers-max",
     po::value<std::size_t>(&myNumPeersMax)->default_value(1),
     "Maximum number of peers per app.")
    ("distance-min",
     po::value<std::size_t>(&myDistanceMin)->default_value(1),
     "Minimum distance of peer from host, in hops.")
    ("distance-max",
     po::value<std::size_t>(&myDistanceMax)->default_value(2),
     "Maximum distance of peer from host, in hops.")
    ("grid-size",
     po::value<double>(&myGridSize)->default_value(60000),
     "Grid length, in km.")
    ("threshold",
     po::value<double>(&myThreshold)->default_value(15000),
     "Link creation threshold (Euclidean distance), in km.")
    ("link-probability",
     po::value<double>(&myLinkProbability)->default_value(1),
     "Link creation probability.")
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
    ("fidelity-threshold",
     po::value<double>(&myFidelityThreshold)->default_value(0.95),
     "Fidelity threshold.")
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

    if (explainOrPrint(myVarMap)) {
      return EXIT_SUCCESS;
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
                                   myQ,
                                   myQuantum,
                                   myK,
                                   myFidelityInit,
                                   myNumApps,
                                   myNumPeersMin,
                                   myNumPeersMax,
                                   myDistanceMin,
                                   myDistanceMax,
                                   myFidelityThreshold,
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
