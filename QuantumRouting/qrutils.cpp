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

#include "QuantumRouting/qrutils.h"
#include "Support/random.h"

#include <cassert>
#include <cmath>
#include <memory>

namespace uiiit {
namespace qr {

double distance(const Coordinate& aLhs, const Coordinate& aRhs) {
  return std::sqrt(std::pow(std::get<0>(aLhs) - std::get<0>(aRhs), 2.0) +
                   std::pow(std::get<1>(aLhs) - std::get<1>(aRhs), 2.0) +
                   std::pow(std::get<2>(aLhs) - std::get<2>(aRhs), 2.0));
}

std::vector<std::pair<std::size_t, std::size_t>>
findLinks(const std::vector<Coordinate>& aItems,
          const double                   aThreshold,
          const double                   aProbability,
          const std::size_t              aSeed) {
  assert(aProbability >= 0 and aProbability <= 1);
  assert(aThreshold >= 0);

  const auto myRv =
      aProbability < 1 ?
          std::make_unique<support::UniformRv>(0, 1, aSeed, 0, 0) :
          nullptr;

  std::vector<std::pair<std::size_t, std::size_t>> ret;
  for (std::size_t i = 0; i < aItems.size(); i++) {
    for (std::size_t j = 0; j < i; j++) {
      if (distance(aItems[i], aItems[j]) < aThreshold and
          (myRv.get() == nullptr or (*myRv)() < aProbability)) {
        ret.push_back({i, j});
      }
    }
  }
  return ret;
}

double fidelitySwapping(const double      p1,
                        const double      p2,
                        const double      eta,
                        const std::size_t L,
                        const double      F) {
  assert(p1 > 0 and p1 <= 1);
  assert(p2 > 0 and p2 <= 1);
  assert(eta >= 0.5 and eta <= 1);
  assert(L > 0);
  assert(F >= 0 and F <= 1);

  return 1.0 / 4.0 +
         3.0 / 4.0 * std::pow(p1 * p1 * p2 * (4 * eta * eta - 1) / 3.0, L - 1) *
             std::pow((4.0 * F - 1.0) / 3.0, L);
}

} // namespace qr
} // namespace uiiit