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

#include "QuantumRouting/mecqkdnetwork.h"
#include "QuantumRouting/networkfactory.h"
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

#include <boost/algorithm/string/replace.hpp>
#include <boost/program_options.hpp>

#include <cmath>
#include <cstddef>
#include <fstream>
#include <glog/logging.h>

#include <cstdlib>
#include <iostream>
#include <iterator>
#include <limits>
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

  // topology generation
  std::size_t theNodes;
  double      theAlpha;
  double      theBeta;
  double      theMaxDistance; // in km
  double      theMaxCapacity; // in secret b/s

  // workload generation
  std::string theAppSpec;        // file containing the apps' specifications
  std::size_t theApplications;   // num of applications
  std::size_t theEdgeNodes;      // num of edge nodes
  std::string theEdgeProcessing; // r.v. to draw edge node processing power

  // system
  qr::MecQkdAlgo theAlgo;

  // not part of the experiment
  std::string theDotFile;

  static const std::vector<std::string>& names() {
    static std::vector<std::string> ret({"seed",

                                         "nodes",
                                         "alpha",
                                         "beta",
                                         "max-distance",
                                         "max-capacity",

                                         "app-spec",
                                         "applications",
                                         "edge-nodes",
                                         "edge-processing",

                                         "algorithm"});
    return ret;
  }

  std::string toString() const {
    std::stringstream myStream;
    myStream << "topology randomly generated according to the Waxman model "
                "with parameters: alpha "
             << theAlpha << ", beta " << theBeta << ", grid length "
             << theMaxDistance << " km; number of nodes: " << theNodes
             << "; maximum capacity of the QKD links: " << theMaxCapacity
             << " b/s; " << theApplications
             << " applications with specifications in " << theAppSpec << ", "
             << theEdgeNodes << " edge nodes with processing power "
             << theEdgeProcessing
             << "; allocation algorithm: " << qr::toString(theAlgo);
    myStream << "; experiment seed " << theSeed;
    return myStream.str();
  }

  std::string toCsv() const {
    std::stringstream myStream;
    myStream << theSeed << ',' << theNodes << ',' << theAlpha << ',' << theBeta
             << ',' << theMaxDistance << ',' << theMaxCapacity << ','
             << theAppSpec << ',' << theApplications << ',' << theEdgeNodes
             << ',' << boost::replace_all_copy(theEdgeProcessing, ",", ";")
             << ',' << qr::toString(theAlgo);
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

  static std::vector<std::string> names() {
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
    });
    return myStaticNames;
  }

  std::string toString() const {
    std::stringstream myStream;
    myStream << "G(" << theNumNodes << "," << theNumEdges << "), in-degree "
             << theMinInDegree << "-" << theMaxOutDegree << ", out-degree "
             << theMinOutDegree << "-" << theMaxOutDegree << ", diameter "
             << theDiameter << ", total capacity " << theTotalCapacity
             << " EPR-pairs/s";
    return myStream.str();
  }

  std::string toCsv() const {
    std::stringstream myStream;
    myStream << theNumNodes << ',' << theNumEdges << ',' << theMinInDegree
             << ',' << theMaxInDegree << ',' << theMinOutDegree << ','
             << theMaxOutDegree << ',' << theDiameter << ','
             << theTotalCapacity;
    return myStream.str();
  }
};

using Data = us::ExperimentData<Parameters, Output>;

void runExperiment(Data& aData, Parameters&& aParameters) {
  Data::Raii myRaii(aData, std::move(aParameters)); // experiment input
  Output     myOutput;                              // experiment output

  // create network
  std::vector<qr::Coordinate> myCoordinates;
  const auto                  myMaxDistance = myRaii.in().theMaxDistance;
  const auto                  myMaxCapacity = myRaii.in().theMaxCapacity;
  auto myNetwork = qr::makeCapacityNetworkWaxman<qr::MecQkdNetwork>(
      [myMaxDistance, myMaxCapacity](const double d) {
        return myMaxCapacity * std::exp(-d / myMaxDistance);
      },
      myRaii.in().theSeed,
      myRaii.in().theNodes,
      myRaii.in().theMaxDistance,
      myRaii.in().theAlpha,
      myRaii.in().theBeta,
      myCoordinates);
  assert(myNetwork->numNodes() == myRaii.in().theNodes);

  // load the workload generator parameters from file
  us::UniformRv myAppRv(0, 1, myRaii.in().theSeed, 3, 0);
  auto          myWorkload =
      qr::MecQkdWorkload::fromCsvFile(myRaii.in().theAppSpec, myAppRv);

  // select the user nodes
  if ((myRaii.in().theEdgeNodes + myWorkload.regions().size()) >
      myRaii.in().theNodes) {
    throw std::runtime_error(
        "cannot have " + std::to_string(myRaii.in().theEdgeNodes) +
        " edge nodes + " + std::to_string(myWorkload.regions().size()) +
        " user nodes in a network with " +
        std::to_string(myRaii.in().theNodes) + " nodes");
  }
  std::set<unsigned long> myAllNodes;
  for (unsigned long i = 0; i < myRaii.in().theNodes; i++) {
    myAllNodes.emplace(i);
  }
  assert(myAllNodes.size() == myRaii.in().theNodes);
  us::UniformRv myNodesRv(0, 1, myRaii.in().theSeed, 4, 0);
  const auto    myUserNodes =
      us::sample(myAllNodes, myWorkload.regions().size(), myNodesRv);
  for (const auto& myNode : myUserNodes) {
    myAllNodes.erase(myNode);
  }
  myNetwork->userNodes(myUserNodes);

  // select the edge nodes from the set of nodes _not_ selected as users
  // and assign them their processing power
  assert(myRaii.in().theEdgeNodes <= myAllNodes.size());
  auto myProcessingRv = us::RealRvInterface::fromString(
      myRaii.in().theEdgeProcessing, myRaii.in().theSeed, 5, 0);
  const auto myEdgeNodes =
      us::sample(myAllNodes, myRaii.in().theEdgeNodes, myNodesRv);
  std::map<unsigned long, double> myEdgeProcessing;
  for (const auto& myNode : myEdgeNodes) {
    myEdgeProcessing.emplace(myNode, (*myProcessingRv)());
  }
  myNetwork->edgeNodes(myEdgeProcessing);

  // draw the applications
  std::vector<qr::MecQkdNetwork::Allocation> myApps;
  for (std::size_t i = 0; i < myRaii.in().theApplications; i++) {
    const auto myAppInfo = myWorkload();
    myApps.emplace_back(qr::MecQkdNetwork::Allocation(
        myAppInfo.theRegion, myAppInfo.theRate, myAppInfo.theLoad));
  }

  // allocate the apps to edge nodes
  myNetwork->allocate(myApps, myRaii.in().theAlgo);

  // save the network properties
  assert(myNetwork.get() != nullptr);
  myOutput.theNumNodes      = myNetwork->numNodes();
  myOutput.theNumEdges      = myNetwork->numEdges();
  myOutput.theTotalCapacity = myNetwork->totalCapacity();
  std::tie(myOutput.theMinInDegree, myOutput.theMaxInDegree) =
      myNetwork->inDegree();
  std::tie(myOutput.theMinOutDegree, myOutput.theMaxOutDegree) =
      myNetwork->outDegree();
  myNetwork->reachableNodes(
      0, std::numeric_limits<std::size_t>::max(), myOutput.theDiameter);

  // save to Graphviz, if needed
  if (not myRaii.in().theDotFile.empty()) {
    myNetwork->toDot(myRaii.in().theDotFile + "-" +
                     std::to_string(myRaii.in().theSeed) + ".dot");

    std::ofstream myStream(myRaii.in().theDotFile + "-" +
                           std::to_string(myRaii.in().theSeed) + ".nodes");
    for (std::size_t i = 0; i < myCoordinates.size(); i++) {
      myStream << i << ' ' << std::get<0>(myCoordinates[i]) << ' '
               << std::get<1>(myCoordinates[i]) << ' '
               << std::get<2>(myCoordinates[i]) << '\n';
    }
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

  std::size_t myNodes;
  double      myAlpha;
  double      myBeta;
  double      myMaxDistance;
  double      myMaxCapacity;
  std::string myAppSpec;
  std::size_t myApplications;
  std::size_t myEdgeNodes;
  std::string myEdgeProcessing;
  std::string myAlgo;

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
    ("nodes",
     po::value<std::size_t>(&myNodes)->default_value(50),
     "Number of nodes.")
    ("alpha",
     po::value<double>(&myAlpha)->default_value(0.4),
     "Waxman' model alpha value.")
    ("beta",
     po::value<double>(&myBeta)->default_value(0.4),
     "Waxman' model beta value.")
    ("max-distance",
     po::value<double>(&myMaxDistance)->default_value(100),
     "Waxman' model max distance, in km.")
    ("max-capacity",
     po::value<double>(&myMaxCapacity)->default_value(100e3),
     "Max QKD capacity, in b/s.")
    ("app-spec",
     po::value<std::string>(&myAppSpec)->default_value("applications.dat"),
     "Name of the CSV file containing the specifications of the applications.")
    ("applications",
     po::value<std::size_t>(&myApplications)->default_value(100),
     "Number of applications.")
    ("edge-nodes",
     po::value<std::size_t>(&myEdgeNodes)->default_value(5),
     "Number of edge nodes.")
    ("edge-processing",
     po::value<std::string>(&myEdgeProcessing)->default_value("U(1.0,3.0)"),
     "Description of the r.v. to be used to determine the edge node processing power.")
    ("algo",
     po::value<std::string>(&myAlgo)->default_value(qr::toString(qr::MecQkdAlgo::Random)),
     (std::string("Algorithm to be used, one of: ") + toString(qr::allMecQkdAlgos(), ", ", [](const auto& aAlgo) { return toString(aAlgo); })).c_str())
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
      myParameters.push(Parameters{mySeed,
                                   myNodes,
                                   myAlpha,
                                   myBeta,
                                   myMaxDistance,
                                   myMaxCapacity,
                                   myAppSpec,
                                   myApplications,
                                   myEdgeNodes,
                                   myEdgeProcessing,
                                   qr::mecQkdAlgofromString(myAlgo),
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
