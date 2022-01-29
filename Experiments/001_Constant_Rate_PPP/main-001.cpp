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
#include "Support/experimentdata.h"
#include "Support/glograii.h"

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
  double      theAvgDijkstraCalls;
  double      theSumGrossRate;
  double      theSumNetRate;
  double      theAdmissionRate;
  std::size_t theAdmittedFlows;

  Output()
      : theAvgDijkstraCalls(0)
      , theSumGrossRate(0)
      , theSumNetRate(0)
      , theAdmissionRate(0)
      , theAdmittedFlows(0) {
    // noop
  }

  std::string toString() const {
    std::stringstream myStream;
    myStream << "admitted " << theAdmittedFlows << " (rate " << theAdmissionRate
             << "), total EPR rate net " << theSumNetRate << " (gross "
             << theSumGrossRate << "), with " << theAvgDijkstraCalls
             << " Dijkstra calls on average";
    return myStream.str();
  }

  std::string toCsv() const {
    std::stringstream myStream;
    myStream << theAvgDijkstraCalls << ',' << theSumGrossRate << ','
             << theSumNetRate << ',' << theAdmissionRate << ','
             << theAdmittedFlows;
    return myStream.str();
  }
};

int main(int argc, char* argv[]) {
  uiiit::support::GlogRaii myGlogRaii(argv[0]);

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

    using Data = us::ExperimentData<Parameters, Output>;
    Data myData;

    for (auto mySeed = mySeedStart; mySeed <= mySeedEnd; ++mySeed) {
      {
        Data::Raii myRaii(myData,
                          Parameters{mySeed,
                                     myMu,
                                     60000,
                                     10000,
                                     0.5,
                                     myLinkMinEpr,
                                     myLinkMaxEpr,
                                     0.5,
                                     0.9925,
                                     myNumFlows,
                                     myMinNetRate,
                                     myMaxNetRate,
                                     0.7});

        Output myOutput;
        myRaii.finish(std::move(myOutput));
      }

      VLOG(1) << "experiment finished\n"
              << myData.lastIn() << '\n'
              << myData.lastOut();
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
