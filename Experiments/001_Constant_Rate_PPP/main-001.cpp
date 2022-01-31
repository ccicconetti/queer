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
#include "QuantumRouting/poissonpointprocess.h"
#include "QuantumRouting/qrutils.h"
#include "Support/experimentdata.h"
#include "Support/glograii.h"
#include "Support/parallelbatch.h"
#include "Support/queue.h"
#include "Support/random.h"
#include "Support/stat.h"

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
  double theEdgeProbability;
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

  std::string toString() const {
    std::stringstream myStream;
    myStream
        << "num nodes drawn from PPP with mu " << theMu
        << " distributed on a flat square grid with edge size " << theGridLength
        << " m, a link is generated with probability " << theEdgeProbability
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
    myStream << theMu << ',' << theGridLength << ',' << theThreshold << ','
             << theEdgeProbability << ',' << theLinkMinEpr << ','
             << theLinkMaxEpr << ',' << theQ << ',' << theFidelityInit << ','
             << theNumFlows << ',' << theMinNetRate << ',' << theMaxNetRate
             << ',' << theFidelityThreshold << ',' << theSeed;
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
             << theAvgPathSize;
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
             << theAvgPathSize;
    return myStream.str();
  }
};

using Data = us::ExperimentData<Parameters, Output>;

void runExperiment(Data& aData, Parameters&& aParameters) {
  Data::Raii myRaii(aData, std::move(aParameters));

  Output myOutput;

  // create network
  us::UniformRv                        myLinkEprRv(myRaii.in().theLinkMinEpr,
                            myRaii.in().theLinkMaxEpr,
                            myRaii.in().theSeed,
                            0,
                            0);
  auto                                 myPppSeed = myRaii.in().theSeed;
  std::unique_ptr<qr::CapacityNetwork> myNetwork = nullptr;
  while (not myNetwork) {
    const auto myCoordinates =
        qr::PoissonPointProcessGrid(myRaii.in().theMu,
                                    myPppSeed,
                                    myRaii.in().theGridLength,
                                    myRaii.in().theGridLength)();
    const auto myEdges = qr::findLinks(myCoordinates,
                                       myRaii.in().theThreshold,
                                       myRaii.in().theEdgeProbability,
                                       myRaii.in().theSeed);
    if (qr::bigraphConnected(myEdges)) {
      myNetwork =
          std::make_unique<qr::CapacityNetwork>(myEdges, myLinkEprRv, true);
      myNetwork->measurementProbability(myRaii.in().theQ);
    } else {
      VLOG(1) << "graph with seed " << myPppSeed
              << " is not connected, trying again";
      myPppSeed += 1000000;
    }
  }

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
  myNetwork->route(myFlows);

  // traffic metrics
  myOutput.theResidualCapacity = myNetwork->totalCapacity();
  us::SummaryStat myDijkstra;
  us::SummaryStat myGrossRate;
  us::SummaryStat myNetRate;
  us::SummaryStat myAdmissionRate;
  us::SummaryStat myPathSize;
  for (const auto& myFlow : myFlows) {
    myDijkstra(myFlow.theDijsktra);
    myGrossRate(myFlow.theGrossRate);
    if (not myFlow.thePath.empty()) {
      myNetRate(myFlow.theNetRate);
      myAdmissionRate(1);
      myPathSize(myFlow.thePath.size());
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

  // save data
  VLOG(1) << "experiment finished\n"
          << myRaii.in().toString() << '\n'
          << myOutput.toString();

  myRaii.finish(std::move(myOutput));
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

  po::options_description myDesc("Allowed options");
  // clang-format off
  myDesc.add_options()
    ("help,h", "produce help message")
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
     po::value<std::size_t>(&mySeedEnd)->default_value(0),
     "Last seed used.")
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
     "Min net rate requested, in EPR/s.")
    ("rate-max",
     po::value<double>(&myMaxNetRate)->default_value(10),
     "Max net rate requested, in EPR/s.")
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

    std::ofstream myFile(myOutputFilename,
                         myVarMap.count("append") == 1 ? std::ios::app :
                                                         std::ios::trunc);
    if (not myFile) {
      throw std::runtime_error("could not open output file for writing: " +
                               myOutputFilename);
    }

    Data myData;

    const double myGridSize          = 60000;
    const double myThreshold         = 10000;
    const double myLinkProbability   = 1;
    const double myQ                 = 0.5;
    const double myFidelityInit      = 0.9925;
    const double myFidelityThreshold = 0.7;

    us::Queue<Parameters> myParameters;
    for (auto mySeed = mySeedStart; mySeed <= mySeedEnd; ++mySeed) {
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
