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
#include <tuple>

namespace uiiit {
namespace qr {

//! Space coordinated (x, y, z)
using Coordinate = std::tuple<double, double, double>;

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
 * @pre L is greater than or equal to 1
 */
double fidelitySwapping(const double      p1,
                        const double      p2,
                        const double      eta,
                        const std::size_t L,
                        const double      F);

} // namespace qr
} // namespace uiiit