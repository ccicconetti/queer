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

#include "QuantumRouting/mecqkdworkload.h"

#include "Support/split.h"

#include <algorithm>
#include <fstream>
#include <glog/logging.h>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace uiiit {
namespace qr {

std::string MecQkdWorkload::AppInfo::toString() const {
  std::stringstream ret;
  ret << "region " << theRegion << ", weight " << theWeight << ", load "
      << theLoad << ", rate " << theRate;
  return ret.str();
}

MecQkdWorkload::MecQkdWorkload(const std::vector<AppInfo>& aAppInfo,
                               support::RealRvInterface&   aRv,
                               const bool                  aWeighted)
    : theAppInfo(aAppInfo)
    , theRv(aRv)
    , theWeighted(aWeighted)
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
    for (const auto& myAppInfo : aAppInfo) {
      LOG(INFO) << myAppInfo.toString();
    }
  }
}

MecQkdWorkload MecQkdWorkload::fromCsvFile(const std::string&        aFilename,
                                           support::RealRvInterface& aRv,
                                           const bool aWeighted) {
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
    const auto myWeight = std::stod(myTokens[1]);
    if (myWeight > 0) {
      myAppInfo.emplace_back(AppInfo{std::stoull(myTokens[0]),
                                     myWeight,
                                     std::stod(myTokens[2]),
                                     std::stod(myTokens[3])});
    }
  }

  return {myAppInfo, aRv, aWeighted};
}

MecQkdWorkload::AppInfo MecQkdWorkload::operator()() {
  if (theWeighted) {
    const auto res = support::sampleWeighted(theAppInfo, theWeights, 1, theRv);
    assert(res.size() == 1);
    return res[0];
  } else {
    return support::choice(theAppInfo, theRv);
  }
}

} // namespace qr
} // namespace uiiit
