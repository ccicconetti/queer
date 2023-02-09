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

#include "QuantumRouting/mecqkdnetwork.h"

#include "Support/split.h"
#include "Support/tostring.h"

#include <boost/graph/subgraph.hpp>
#include <fstream>
#include <glog/logging.h>
#include <sstream>
#include <stdexcept>

namespace uiiit {
namespace qr {

std::vector<MecQkdAlgo> allMecQkdAlgos() {
  static const std::vector<MecQkdAlgo> myAlgos({
      MecQkdAlgo::Random,
      MecQkdAlgo::Spf,
      MecQkdAlgo::BestFit,
      MecQkdAlgo::RandomFeas,
      MecQkdAlgo::SpfFeas,
      MecQkdAlgo::BestFitFeas,
  });
  return myAlgos;
}

std::string toString(const MecQkdAlgo aAlgo) {
  switch (aAlgo) {
    case MecQkdAlgo::Random:
      return "random";
    case MecQkdAlgo::Spf:
      return "spf";
    case MecQkdAlgo::BestFit:
      return "bestfit";
    case MecQkdAlgo::RandomFeas:
      return "randomfeas";
    case MecQkdAlgo::SpfFeas:
      return "spffeas";
    case MecQkdAlgo::BestFitFeas:
      return "bestfitfeas";
    default:; /* fall-through */
  }
  return "unknown";
}

MecQkdAlgo mecQkdAlgofromString(const std::string& aAlgo) {
  if (aAlgo == "random") {
    return MecQkdAlgo::Random;
  } else if (aAlgo == "bestfit") {
    return MecQkdAlgo::BestFit;
  } else if (aAlgo == "spf") {
    return MecQkdAlgo::Spf;
  } else if (aAlgo == "randomfeas") {
    return MecQkdAlgo::RandomFeas;
  } else if (aAlgo == "bestfitfeas") {
    return MecQkdAlgo::BestFitFeas;
  } else if (aAlgo == "spffeas") {
    return MecQkdAlgo::SpfFeas;
  }
  throw std::runtime_error(
      "invalid edge QKD algorithm: " + aAlgo + " (valid options are: " +
      ::toString(allMecQkdAlgos(),
                 ",",
                 [](const auto& aAlgo) { return toString(aAlgo); }) +
      ")");
}

std::string MecQkdWorkload::AppInfo::toString() const {
  std::stringstream ret;
  ret << "region " << theRegion << ", weight " << theWeight << ", load "
      << theLoad << ", rate " << theRate;
  return ret.str();
}

MecQkdWorkload::MecQkdWorkload(const std::vector<AppInfo>& aAppInfo,
                               support::RealRvInterface&   aRv)
    : theAppInfo(aAppInfo)
    , theRv(aRv)
    , theRegions()
    , theWeights(aAppInfo.size()) {
  if (aAppInfo.empty()) {
    throw std::runtime_error("invalid empty MecQkd workload");
  }
  for (std::size_t i = 0; i < aAppInfo.size(); i++) {
    theRegions.emplace(aAppInfo[i].theRegion);
    theWeights[i] = aAppInfo[i].theWeight;
  }

  if (VLOG_IS_ON(1)) {
    for (std::size_t i = 0; i < aAppInfo.size(); i++) {
      LOG(INFO) << aAppInfo[i].toString();
    }
  }
}

MecQkdWorkload MecQkdWorkload::fromCsvFile(const std::string&        aFilename,
                                           support::RealRvInterface& aRv) {
  std::vector<AppInfo> myAppInfo;

  std::ifstream myInfile(aFilename);
  if (not myInfile) {
    throw std::runtime_error("could not open file for reading: " + aFilename);
  }

  std::string myLine;
  std::size_t myLineNo = 0;
  while (myInfile) {
    ++myLineNo;
    std::getline(myInfile, myLine);
    if (myLine.empty() or myLine[0] == '#') {
      continue;
    }
    const auto myTokens = support::split<std::vector<std::string>>(myLine, ",");
    if (myTokens.size() != 4) {
      throw std::runtime_error("invalid input at file '" + aFilename +
                               "' line: " + std::to_string(myLineNo));
    }
    myAppInfo.emplace_back(AppInfo{std::stoull(myTokens[0]),
                                   std::stod(myTokens[1]),
                                   std::stod(myTokens[2]),
                                   std::stod(myTokens[3])});
  }

  return MecQkdWorkload(myAppInfo, aRv);
}

MecQkdWorkload::AppInfo MecQkdWorkload::operator()() {
  const auto res = support::sampleWeighted(theAppInfo, theWeights, 1, theRv);
  assert(res.size() == 1);
  return res[0];
}

MecQkdNetwork::MecQkdNetwork(
    const std::vector<std::pair<unsigned long, unsigned long>>& aEdges,
    support::RealRvInterface&                                   aWeightRv,
    const bool aMakeBidirectional)
    : CapacityNetwork(aEdges, aWeightRv, aMakeBidirectional) {
  // noop
}

MecQkdNetwork::MecQkdNetwork(const WeightVector& aEdgeWeights)
    : CapacityNetwork(aEdgeWeights) {
  // noop
}

void MecQkdNetwork::userNodes(const std::set<unsigned long>& aUserNodes) {
  const auto V = boost::num_vertices(theGraph);
  for (const auto& v : aUserNodes) {
    if (v >= V) {
      throw std::runtime_error("Invalid user node: " + std::to_string(v));
    }
  }

  theUserNodes = aUserNodes;
}

void MecQkdNetwork::edgeNodes(
    const std::map<unsigned long, double>& aEdgeProcessing) {
  const auto V = boost::num_vertices(theGraph);
  theEdgeNodes.clear();
  for (const auto& elem : aEdgeProcessing) {
    if (elem.second < 0) {
      throw std::runtime_error("Invalid negative processing capability: " +
                               std::to_string(elem.second));
    }
    if (elem.first >= V) {
      throw std::runtime_error("Invalid edge node: " +
                               std::to_string(elem.first));
    }
    theEdgeNodes.emplace(elem.first, EdgeNode{elem.second, 0.0});
  }
}

} // namespace qr
} // namespace uiiit
