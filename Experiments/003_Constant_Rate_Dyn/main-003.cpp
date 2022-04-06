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
#include "Support/glograii.h"
#include "Support/parallelbatch.h"
#include "Support/queue.h"
#include "Support/random.h"
#include "Support/split.h"
#include "Support/stat.h"
#include "Support/tostring.h"
#include "Support/versionutils.h"

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
      aValues, "|", [](const double& x) { return std::to_string(x); });
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
  std::string theGraphMlFilename;

  // system
  double theQ;
  double theFidelityInit;
  double theSimDuration;

  // application flows
  double              theArrivalRate;
  double              theFlowDuration;
  std::vector<double> theNetRates;
  std::vector<double> theFidelityThresholds;

  static const std::vector<std::string>& names() {
    static std::vector<std::string> ret({
        "seed",
        "mu",
        "grid-length",
        "threshold",
        "link-prob",
        "link-min-epr",
        "link-max-epr",
        "graphml-filename",
        "q",
        "fidelity-init",
        "sim-duration",
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
        << ", and EPR generation rate of the list is drawn randomly from U["
        << theLinkMinEpr << ',' << theLinkMaxEpr << "]; simulation duration "
        << theSimDuration << "; probability of correct BSM " << theQ
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
             << theLinkMinEpr << ',' << theLinkMaxEpr << ','
             << theGraphMlFilename << ',' << theQ << ',' << theFidelityInit
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
        "capacity-tot",
        "capacity-res",
        "avg-dijkstra-calls",
        "sum-gross-rate",
        "sum-net-rate",
        "admission-rate",
        "avg-path-size",
        "avg-fidelity",
    });
    static const std::vector<std::string> myPerClassNames({
        "sum-gross-rate",
        "sum-net-rate",
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
  double      theTotalCapacity = 0;

  // routing properties
  double theResidualCapacity = 0;
  double theAvgDijkstraCalls = 0;
  double theSumGrossRate     = 0;
  double theSumNetRate       = 0;
  double theAdmissionRate    = 0;
  double theAvgPathSize      = 0;
  double theAvgFidelity      = 0;

  struct PerClass {
    double theSumGrossRate  = 0;
    double theSumNetRate    = 0;
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
             << theMinOutDegree << "-" << theMaxOutDegree
             << " with total capacity " << theTotalCapacity
             << " EPR/s (residual " << theResidualCapacity
             << " EPR/s); admission rate " << theAdmissionRate
             << ", total EPR rate net " << theSumNetRate << " (gross "
             << theSumGrossRate << "), with " << theAvgDijkstraCalls
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
             << theMaxOutDegree << ',' << theTotalCapacity << ','
             << theResidualCapacity << ',' << theAvgDijkstraCalls << ','
             << theSumGrossRate << ',' << theSumNetRate << ','
             << theAdmissionRate << ',' << theAvgPathSize << ','
             << theAvgFidelity;
    return myStream.str();
  }
};

using Data = us::ExperimentData<Parameters, Output>;

void runExperiment(Data& aData, Parameters&& aParameters) {
  // fidelity computation parameters
  // constexpr double p1  = 1.0;
  // constexpr double p2  = 1.0;
  // constexpr double eta = 1.0;

  Data::Raii myRaii(aData, std::move(aParameters));

  Output myOutput(myRaii.in().theNetRates, myRaii.in().theFidelityThresholds);

  // open the GraphML file with the network topology, if needed
  const auto myGraphMlStream =
      myRaii.in().theGraphMlFilename.empty() ?
          nullptr :
          std::make_unique<std::ifstream>(myRaii.in().theGraphMlFilename);
  if (myGraphMlStream.get() != nullptr and
      not static_cast<bool>(*myGraphMlStream)) {
    throw std::runtime_error("cannot read from file: " +
                             myRaii.in().theGraphMlFilename);
  }

  // create network
  us::UniformRv myLinkEprRv(myRaii.in().theLinkMinEpr,
                            myRaii.in().theLinkMaxEpr,
                            myRaii.in().theSeed,
                            0,
                            0);
  const auto    myNetwork =
      myRaii.in().theGraphMlFilename.empty() ?
          qr::makeCapacityNetworkPpp(myLinkEprRv,
                                     myRaii.in().theSeed,
                                     myRaii.in().theMu,
                                     myRaii.in().theGridLength,
                                     myRaii.in().theThreshold,
                                     myRaii.in().theLinkProbability) :
          qr::makeCapacityNetworkGraphMl(myLinkEprRv, *myGraphMlStream);
  myNetwork->measurementProbability(myRaii.in().theQ);

  // network properties
  assert(myNetwork.get() != nullptr);
  myOutput.theNumNodes      = myNetwork->numNodes();
  myOutput.theNumEdges      = myNetwork->numEdges();
  myOutput.theTotalCapacity = myNetwork->totalCapacity();
  std::tie(myOutput.theMinInDegree, myOutput.theMaxInDegree) =
      myNetwork->inDegree();
  std::tie(myOutput.theMinOutDegree, myOutput.theMaxOutDegree) =
      myNetwork->outDegree();

  // run simulation
  us::SummaryStat myDijkstra;
  us::SummaryStat myGrossRate;
  us::SummaryStat myNetRate;
  us::SummaryStat myAdmissionRate;
  us::SummaryStat myPathSize;
  us::SummaryStat myFidelity;
  double          myNow = 0;

  struct AdmittedFlow {
    qr::CapacityNetwork::FlowDescriptor theDesc;
    double                              theLeaveTime;
    double                              theFidelityThresh;
  };

  std::list<AdmittedFlow> myFlows;
  while (myNow <= myRaii.in().theSimDuration) {
    // us::UniformRv                   myNetRateRv(myRaii.in().theMinNetRate,
    //                           myRaii.in().theMaxNetRate,
    //                           myRaii.in().theSeed,
    //                           0,
    //                           0);
    // us::UniformIntRv<unsigned long> mySrcDstRv(
    //     0, myNetwork->numNodes() - 1, myRaii.in().theSeed, 0, 0);
    // for (std::size_t i = 0; i < myRaii.in().theNumFlows; i++) {
    //   unsigned long mySrc = 0;
    //   unsigned long myDst = 0;
    //   while (mySrc == myDst) {
    //     mySrc = mySrcDstRv();
    //     myDst = mySrcDstRv();
    //   }
    //   assert(mySrc != myDst);

    //   myFlows.emplace_back(mySrc, myDst, myNetRateRv());
    // }

    // route traffic flows
    // myNetwork->route(myFlows, [&myRaii](const auto& aFlow) {
    //   assert(not aFlow.thePath.empty());
    //   return qr::fidelitySwapping(p1,
    //                               p2,
    //                               eta,
    //                               aFlow.thePath.size() - 1,
    //                               myRaii.in().theFidelityInit) >= 0; // XXX
    // });

    // traffic metrics
    // myOutput.theResidualCapacity = myNetwork->totalCapacity();
    // for (const auto& myFlow : myFlows) {
    //   myDijkstra(myFlow.theDijsktra);
    //   myGrossRate(myFlow.theGrossRate);
    //   if (not myFlow.thePath.empty()) {
    //     myNetRate(myFlow.theNetRate);
    //     myAdmissionRate(1);
    //     myPathSize(myFlow.thePath.size());
    //     myFidelity(qr::fidelitySwapping(p1,
    //                                     p2,
    //                                     eta,
    //                                     myFlow.thePath.size() - 1,
    //                                     myRaii.in().theFidelityInit));
    //   } else {
    //     myAdmissionRate(0);
    //   }
    // }
  }

  myOutput.theAvgDijkstraCalls = myDijkstra.mean();
  myOutput.theSumGrossRate     = myGrossRate.count() * myGrossRate.mean();
  myOutput.theSumNetRate       = myNetRate.count() * myNetRate.mean();
  myOutput.theAdmissionRate    = myAdmissionRate.mean();
  myOutput.theAvgPathSize      = myPathSize.mean();
  myOutput.theAvgFidelity      = myFidelity.mean();

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
  std::string myNetRatesStr;
  double      myGridSize;
  double      myThreshold;
  double      myLinkProbability;
  double      myQ;
  double      myFidelityInit;
  std::string myFidelityThresholdsStr;
  std::string myGraphMlFilename;
  double      mySimDuration;
  double      myArrivalRate;
  double      myFlowDuration;

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
    ("graphml-file",
     po::value<std::string>(&myGraphMlFilename)->default_value(""),
     "Read from a GraphML file.")
    ("sim-duration",
     po::value<double>(&mySimDuration)->default_value(100),
     "Simulation duration, in time units.")
    ("arrival-rate",
     po::value<double>(&myArrivalRate)->default_value(1),
     "Average arrival rates of new flows, in time units^-1.")
    ("flow-duration",
     po::value<double>(&myFlowDuration)->default_value(5),
     "Average duration of an admitted flow, in time units.")
    ("net-epr-rates",
     po::value<std::string>(&myNetRatesStr)->default_value("1"),
     "Set of possible net rates requested, in EPR-pairs/s: multiple values separated by |.")
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
     "Set of possible fidelity thresholds: multiple values separated by |.")
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

    const auto myNetRates = us::split<std::vector<double>>(myNetRatesStr, "|");
    if (myNetRates.empty()) {
      throw std::runtime_error("invalid empty set of net EPR rates");
    }

    const auto myFidelityThresholds =
        us::split<std::vector<double>>(myFidelityThresholdsStr, "|");
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
                                   myGraphMlFilename,
                                   myQ,
                                   myFidelityInit,
                                   mySimDuration,
                                   myArrivalRate,
                                   myFlowDuration,
                                   myNetRates,
                                   myFidelityThresholds});
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
