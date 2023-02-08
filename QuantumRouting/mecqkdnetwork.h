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

#pragma once

#include "QuantumRouting/capacitynetwork.h"

namespace uiiit {
namespace qr {

enum class MecQkdAlgo {
  Random      = 0, //!< pick a random edge node
  Spf         = 1, //!< shortest path first
  BestFit     = 2, //!< best-fit allocation
  RandomFeas  = 3, //!< pick a random edge node among those feasible
  SpfFeas     = 4, //!< shortest path first, only on feasible edge nodes
  BestFitFeas = 5, //!< best-fit allocation, only along feasible paths

};

std::vector<MecQkdAlgo> allMecQkdAlgos();
std::string             toString(const MecQkdAlgo aAlgo);
MecQkdAlgo              mecQkdAlgofromString(const std::string& aAlgo);

/**
 * @brief A QKD network where some nodes also offer edge computing resources.
 */
class MecQkdNetwork final : public CapacityNetwork
{
 public:
  /**
   * @brief Create a network with given links and assign random weights
   *
   * The default measurement probability is 1.
   *
   * @param aEdges The edges of the network (src, dst).
   * @param aWeightRv The r.v. to draw the edge weights.
   * @param aMakeBidirectional If true then for each pair (A,B) two edges are
   * added A->B and B->A, with the same weight.
   */
  explicit MecQkdNetwork(const EdgeVector&         aEdges,
                         support::RealRvInterface& aWeightRv,
                         const bool                aMakeBidirectional);

  /**
   * @brief Create a network with given links and weights
   *
   * The default measurement probability is 1.
   *
   * @param aEdgeWeights The unidirectional edges and weights of the network
   * (src, dst, w).
   */
  explicit MecQkdNetwork(const WeightVector& aEdgeWeights);
};

} // namespace qr
} // namespace uiiit
