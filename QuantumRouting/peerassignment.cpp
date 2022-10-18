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

#include "QuantumRouting/peerassignment.h"

#include "Support/tostring.h"

#include <stdexcept>

namespace uiiit {
namespace qr {

std::vector<PeerAssignmentAlgo> allPeerAssignmentAlgos() {
  static const std::vector<PeerAssignmentAlgo> myAlgos({
      PeerAssignmentAlgo::Random,
      PeerAssignmentAlgo::ShortestPath,
      PeerAssignmentAlgo::Gap,
  });
  return myAlgos;
}

std::string toString(const PeerAssignmentAlgo aAlgo) {
  switch (aAlgo) {
    case PeerAssignmentAlgo::Random:
      return "random";
    case PeerAssignmentAlgo::ShortestPath:
      return "shortestpath";
    case PeerAssignmentAlgo::Gap:
      return "gap";
    default:; /* fall-through */
  }
  return "unknown";
}

PeerAssignmentAlgo peerAssignmentAlgofromString(const std::string& aAlgo) {
  if (aAlgo == "random") {
    return PeerAssignmentAlgo::Random;
  } else if (aAlgo == "shortestpath") {
    return PeerAssignmentAlgo::ShortestPath;
  } else if (aAlgo == "gap") {
    return PeerAssignmentAlgo::Gap;
  }
  throw std::runtime_error(
      "invalid peer assignment algorithm: " + aAlgo + " (valid options are: " +
      ::toString(allPeerAssignmentAlgos(),
                 ",",
                 [](const auto& aAlgo) { return toString(aAlgo); }) +
      ")");
}

} // namespace qr
} // namespace uiiit
