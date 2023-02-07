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

#include <cinttypes>
#include <istream>
#include <tuple>
#include <utility>
#include <vector>

namespace uiiit {
namespace qr {

//! Space coordinate (x, y, z)
using Coordinate = std::tuple<double, double, double>;

//! \return the Euclidean distance between two points.
double distance(const Coordinate& aLhs, const Coordinate& aRhs);

/**
 * @brief Detect items closer than a threshold distance, with a given
 * probability.
 *
 * Links are assumed to be reciprocal: if A and B are two items then either
 * there is no link with A and B or there is exactly one (A-B or B-A, but not
 * both).
 *
 * @param aItems the items to be considered
 * @param aThreshold the minimum distance to have a link between items
 * @param aProbability the probability that a link is created between two items,
 * provided that their Euclidean distance is within the threshold
 * @param aSeed the pseudo-random number generator seed to use, if aProbability
 * is < 1
 *
 * @return the list of pairs of items detected (returned through their indices)
 *
 * @pre aThreshold is non-negative
 * @pre aProbability is in [0, 1]
 */
std::vector<std::pair<unsigned long, unsigned long>>
findLinks(const std::vector<Coordinate>& aItems,
          const double                   aThreshold,
          const double                   aProbability = 1,
          const unsigned long            aSeed        = 0);

/**
 * @brief Return the links between vertices as read from a GraphML file.
 *
 * @param aGraphMl The stream containing the GraphML data.
 * @param aCoordinates The coordinates of the nodes read.
 *
 * @return the list of pairs of items detected (returned through their indices)
 */
std::vector<std::pair<unsigned long, unsigned long>>
findLinks(std::istream& aGraphMl, std::vector<Coordinate>& aCoordinates);

/**
 * @brief Detect if the bidirectional graph defined by the given edges is
 * connected.
 *
 * \return true if connected, false otherwise
 */
bool bigraphConnected(
    const std::vector<std::pair<unsigned long, unsigned long>>& aEdges);

/**
 * @brief Return the fidelity of a pair of entangled qubits through L
 * neighboring links, each performing entanglement swapping, without considering
 * decoherence and depolarization impairments, see:
 * https://arxiv.org/abs/quant-ph/9803056
 *
 * @param p1 the reliability on one-qubit operations
 * @param p2 the reliability of two-qubit operations
 * @param eta the probability of a wrong measurement
 * @param L the number of intermediate nodes, i.e., swaps
 * @param F the local entanglement fidelity
 * @return the resulting fidelity in [0,1]
 *
 * @pre p1, p2, and F are in (0,1]
 * @pre eta is in [0.5, 1]
 */
double fidelitySwapping(const double        p1,
                        const double        p2,
                        const double        eta,
                        const unsigned long L,
                        const double        F);

} // namespace qr
} // namespace uiiit