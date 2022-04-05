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
#include "Support/stat.h"
#include "Support/versionutils.h"

#include <boost/program_options.hpp>

#include <glog/logging.h>

#include <cstdlib>
#include <iostream>
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
  double theQ;
  double theFidelityInit;

  // application flows
  std::size_t theNumFlows;
  double      theMinNetRate;
  double      theMaxNetRate;
  double      theFidelityThreshold;

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
        "fidelity-init",
        "num-flows",
        "min-net-rate",
        "max-net-rate",
        "fidelity-thresh",
    });
    return ret;
  }

  std::string toString() const {
    std::stringstream myStream;
    myStream
        << "num nodes drawn from PPP with mu " << theMu
        << " distributed on a flat square grid with edge size " << theGridLength
        << " m, a link is generated with probability " << theLinkProbability
        << " between any two nodes within " << theThreshold
        << " m apart, and the EPR generation rate of the list is drawn "
           "randomly from U["
        << theLinkMinEpr << ',' << theLinkMaxEpr
        << "]; probability of correct BSM " << theQ
        << " and fidelity of freshly generated pairs " << theFidelityInit
        << "; there are " << theNumFlows
        << " application flows requesting admission, with minimum fidelity "
        << theFidelityThreshold
        << " and a net EPR requested rate drawn randomly from U["
        << theMinNetRate << ',' << theMaxNetRate << "]"
        << ", experiment seed " << theSeed;
    ;
    return myStream.str();
  }

  std::string toCsv() const {
    std::stringstream myStream;
    myStream << theSeed << ',' << theMu << ',' << theGridLength << ','
             << theThreshold << ',' << theLinkProbability << ','
             << theLinkMinEpr << ',' << theLinkMaxEpr << ',' << theQ << ','
             << theFidelityInit << ',' << theNumFlows << ',' << theMinNetRate
             << ',' << theMaxNetRate << ',' << theFidelityThreshold;
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
  double      theTotalCapacity = 0;

  // routing properties
  double      theResidualCapacity = 0;
  double      theAvgDijkstraCalls = 0;
  double      theSumGrossRate     = 0;
  double      theSumNetRate       = 0;
  double      theAdmissionRate    = 0;
  std::size_t theAdmittedFlows    = 0;
  double      theAvgPathSize      = 0;
  double      theAvgFidelity      = 0;

  static const std::vector<std::string>& names() {
    static std::vector<std::string> ret({
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
        "admitted-flows",
        "avg-path-size",
        "avg-fidelity",
    });
    return ret;
  }

  std::string toString() const {
    std::stringstream myStream;
    myStream << "G(" << theNumNodes << "," << theNumEdges << "), in-degree "
             << theMinInDegree << "-" << theMaxOutDegree << ", out-degree "
             << theMinOutDegree << "-" << theMaxOutDegree
             << " with total capacity " << theTotalCapacity
             << " EPR/s (residual " << theResidualCapacity
             << " EPR/s); admitted " << theAdmittedFlows << " (rate "
             << theAdmissionRate << "), total EPR rate net " << theSumNetRate
             << " (gross " << theSumGrossRate << "), with "
             << theAvgDijkstraCalls
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
             << theAdmissionRate << ',' << theAdmittedFlows << ','
             << theAvgPathSize << ',' << theAvgFidelity;
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

  Output myOutput;

  // create network
  const auto myNetwork =
      qr::makeCapacityNetworkPpp(myRaii.in().theLinkMinEpr,
                                 myRaii.in().theLinkMaxEpr,
                                 myRaii.in().theSeed,
                                 myRaii.in().theMu,
                                 myRaii.in().theGridLength,
                                 myRaii.in().theThreshold,
                                 myRaii.in().theLinkProbability);
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

  // create traffic flows
  std::vector<qr::CapacityNetwork::FlowDescriptor> myFlows;
  us::UniformRv                   myNetRateRv(myRaii.in().theMinNetRate,
                            myRaii.in().theMaxNetRate,
                            myRaii.in().theSeed,
                            0,
                            0);
  us::UniformIntRv<unsigned long> mySrcDstRv(
      0, myNetwork->numNodes() - 1, myRaii.in().theSeed, 0, 0);
  for (std::size_t i = 0; i < myRaii.in().theNumFlows; i++) {
    unsigned long mySrc = 0;
    unsigned long myDst = 0;
    while (mySrc == myDst) {
      mySrc = mySrcDstRv();
      myDst = mySrcDstRv();
    }
    assert(mySrc != myDst);

    myFlows.emplace_back(mySrc, myDst, myNetRateRv());
  }

  // route traffic flows
  myNetwork->route(myFlows, [&myRaii](const auto& aFlow) {
    assert(not aFlow.thePath.empty());
    return qr::fidelitySwapping(p1,
                                p2,
                                eta,
                                aFlow.thePath.size() - 1,
                                myRaii.in().theFidelityInit) >=
           myRaii.in().theFidelityThreshold;
  });

  // traffic metrics
  myOutput.theResidualCapacity = myNetwork->totalCapacity();
  us::SummaryStat myDijkstra;
  us::SummaryStat myGrossRate;
  us::SummaryStat myNetRate;
  us::SummaryStat myAdmissionRate;
  us::SummaryStat myPathSize;
  us::SummaryStat myFidelity;
  for (const auto& myFlow : myFlows) {
    myDijkstra(myFlow.theDijsktra);
    myGrossRate(myFlow.theGrossRate);
    if (not myFlow.thePath.empty()) {
      myNetRate(myFlow.theNetRate);
      myAdmissionRate(1);
      myPathSize(myFlow.thePath.size());
      myFidelity(qr::fidelitySwapping(
          p1, p2, eta, myFlow.thePath.size() - 1, myRaii.in().theFidelityInit));
    } else {
      myAdmissionRate(0);
    }
  }
  myOutput.theAvgDijkstraCalls = myDijkstra.mean();
  myOutput.theSumGrossRate     = myGrossRate.count() * myGrossRate.mean();
  myOutput.theSumNetRate       = myNetRate.count() * myNetRate.mean();
  myOutput.theAdmissionRate    = myAdmissionRate.mean();
  myOutput.theAdmittedFlows = myAdmissionRate.count() * myAdmissionRate.mean();
  myOutput.theAvgPathSize   = myPathSize.mean();
  myOutput.theAvgFidelity   = myFidelity.mean();

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
  uiiit::support::GlogRaii myGlogRaii(argv[0]);

  std::size_t myNumThreads;
  std::string myOutputFilename;
  std::size_t mySeedStart;
  std::size_t mySeedEnd;

  double      myMu;
  double      myLinkMinEpr;
  double      myLinkMaxEpr;
  std::size_t myNumFlows;
  double      myMinNetRate;
  double      myMaxNetRate;
  double      myGridSize;
  double      myThreshold;
  double      myLinkProbability;
  double      myQ;
  double      myFidelityInit;
  double      myFidelityThreshold;

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
    ("num-flows",
     po::value<std::size_t>(&myNumFlows)->default_value(100),
     "Number of flows.")
    ("rate-min",
     po::value<double>(&myMinNetRate)->default_value(1),
     "Min net rate requested, in EPR-pairs/s.")
    ("rate-max",
     po::value<double>(&myMaxNetRate)->default_value(10),
     "Max net rate requested, in EPR-pairs/s.")
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
    ("fidelity-init",
     po::value<double>(&myFidelityInit)->default_value(0.99),
     "Fidelity of local entanglement between adjacent nodes.")
    ("fidelity-threshold",
     po::value<double>(&myFidelityThreshold)->default_value(0.95),
     "Fidelity threshold.")
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
                                   myFidelityInit,
                                   myNumFlows,
                                   myMinNetRate,
                                   myMaxNetRate,
                                   myFidelityThreshold});
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
