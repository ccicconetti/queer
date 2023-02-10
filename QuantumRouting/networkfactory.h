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
#include "QuantumRouting/qrutils.h"

#include "QuantumRouting/poissonpointprocess.h"
#include "QuantumRouting/qrutils.h"
#include "Support/random.h"

#include <glog/logging.h>

#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace uiiit {
namespace qr {

/**
 * @brief Create a network from a Poisson Point Process.
 *
 * @tparam NETWORK The type of network to be create.
 * @param aEprRv The r.v. to draw the capacity of the edges.
 * @param aSeed The seed for random number generation.
 * @param aMu The average number of nodes.
 * @param aGridLength The grid length.
 * @param aThreshold The threshold to add an edge between two nodes.
 * @param aLinkProbability The probability that an edge is added.
 * @param aCoordinates The coordinateds of the nodes in a grid.
 * @return std::unique_ptr<CapacityNetwork> The network created.
 * @throw std::runtime_error if the network could not be generated.
 */
template <class NETWORK>
std::unique_ptr<NETWORK>
makeCapacityNetworkPpp(support::RealRvInterface& aEprRv,
                       const std::size_t         aSeed,
                       const double              aMu,
                       const double              aGridLength,
                       const double              aThreshold,
                       const double              aLinkProbability,
                       std::vector<Coordinate>&  aCoordinates) {
  const auto MANY_TRIES = 1000000u;

  auto myPppSeed = aSeed;
  for (unsigned myTry = 0; myTry < MANY_TRIES; myTry++) {
    auto myCoordinates =
        PoissonPointProcessGrid(aMu, myPppSeed, aGridLength, aGridLength)();
    const auto myEdges =
        findLinks(myCoordinates, aThreshold, aLinkProbability, aSeed);
    if (qr::bigraphConnected(myEdges)) {
      myCoordinates.swap(aCoordinates);
      return std::make_unique<NETWORK>(myEdges, aEprRv, true);

    } else {
      VLOG(1) << "graph with seed " << myPppSeed << " not connected, try again";
      myPppSeed += 1000000;
    }
  }

  throw std::runtime_error("Could not find a connected network after " +
                           std::to_string(MANY_TRIES) + " tries");
}

/**
 * @brief Create a network using the model proposed by Waxman in this paper
 *
 * B. M. Waxman, Routing of multipoint connections.
 * IEEE J. Select. Areas Commun. 6(9),(1988) 1617–1622.
 * https://doi.org/10.1109/49.12889
 *
 * @tparam NETWORK The type of network to be create.
 * @param aCapacityLambda The function to determine the capacity from the
 * actual distance between the nodes.
 * @param aSeed The seed for random number generation.
 * @param aNodes The number of nodes.
 * @param aL The maximum distance between two nodes.
 * @param aAlpha The larger this value, the higher the density of short edges
 * compared to longer ones.
 * @param aBeta The larger this value, the higher the edge density.
 * @param aCoordinates The coordinateds of the nodes in a grid.
 * @return std::unique_ptr<CapacityNetwork> The network created.
 * @throw std::range_error if aAlpha or aBeta are not in (0,1].
 * @throw std::runtime_error if the network could not be generated.
 */
template <class NETWORK>
std::unique_ptr<NETWORK> makeCapacityNetworkWaxman(
    const std::function<double(const double)>& aCapacityLambda,
    const std::size_t                          aSeed,
    const std::size_t                          aNodes,
    const double                               aL,
    const double                               aAlpha,
    const double                               aBeta,
    std::vector<Coordinate>&                   aCoordinates) {
  if (aAlpha <= 0 or aAlpha > 1) {
    throw std::range_error(
        "Value of alpha not in (0,1] in Waxman model network generation: " +
        std::to_string(aAlpha));
  }
  if (aBeta <= 0 or aBeta > 1) {
    throw std::range_error(
        "Value of beta not in (0,1] in Waxman model network generation: " +
        std::to_string(aBeta));
  }
  if (aL < 0) {
    throw std::range_error(
        "Negative value of L in Waxman model network generation: " +
        std::to_string(aL));
  }
  if (aNodes == 0) {
    throw std::range_error("Empty network");
  }
  const auto MANY_TRIES = 1000000u;

  const auto myProbLambda = [aAlpha, aBeta, aL](const double d) {
    return aBeta * std::exp(-d / (aAlpha * aL));
  };
  const auto myGridLength = aL / std::sqrt(2.0);

  for (unsigned myTry = 0; myTry < MANY_TRIES; myTry++) {

    // assign nodes to their coordinates
    support::UniformRv      myRvX(0, myGridLength, aSeed, 0, myTry);
    support::UniformRv      myRvY(0, myGridLength, aSeed, 1, myTry);
    std::vector<Coordinate> myCoordinates(aNodes, Coordinate{0, 0, 0});
    for (auto& myNode : myCoordinates) {
      std::get<0>(myNode) = myRvX();
      std::get<1>(myNode) = myRvY();
    }

    // create edges
    CapacityNetwork::WeightVector myEdges;
    support::UniformRv            myRxEdge(0, 1, aSeed, 2, myTry);
    for (std::size_t i = 0; i < aNodes; i++) {
      for (std::size_t j = (i + 1); j < aNodes; j++) {
        const auto myDistance = distance(myCoordinates[i], myCoordinates[j]);
        const auto myEdgeProb = myProbLambda(myDistance);
        assert(myEdgeProb >= 0 and myEdgeProb <= 1);
        if (myRxEdge() < myEdgeProb) {
          // add an edge between the two nodes in both directions
          const auto myCapacity = aCapacityLambda(myDistance);
          myEdges.push_back({i, j, myCapacity});
          myEdges.push_back({j, i, myCapacity});
        }
      }
    }

    std::vector<std::pair<unsigned long, unsigned long>> myEdgesUnweighted(
        myEdges.size());
    for (std::size_t i = 0; i < myEdges.size(); i++) {
      std::tie(myEdgesUnweighted[i].first, myEdgesUnweighted[i].second) = {
          std::get<0>(myEdges[i]), std::get<1>(myEdges[i])};
    }

    if (bigraphConnected(myEdgesUnweighted)) {
      myCoordinates.swap(aCoordinates);
      return std::make_unique<NETWORK>(myEdges);

    } else {
      VLOG(1) << "graph with seed " << myTry << " not connected, try again";
    }
  }

  throw std::runtime_error("Could not find a connected network after " +
                           std::to_string(MANY_TRIES) + " tries");
} // namespace qr

/**
 * @brief Create a network from a GraphML file.
 *
 * @tparam NETWORK The type of the network created.
 * @param aEprRv The r.v. to draw the capacity of the edges.
 * @param aGraphMl The input GraphML file.
 * @param aCoordinates The coordinateds of the nodes in a grid.
 * @return std::unique_ptr<NETWORK> The network created.
 * @throw std::runtime_error if the network in the file is not connected.
 */
template <class NETWORK>
std::unique_ptr<NETWORK>
makeCapacityNetworkGraphMl(support::RealRvInterface& aEprRv,
                           std::ifstream&            aGraphMl,
                           std::vector<Coordinate>&  aCoordinates) {
  const auto myEdges = findLinks(aGraphMl, aCoordinates);
  for (const auto& myEdge : myEdges) {
    VLOG(2) << '(' << myEdge.first << ',' << myEdge.second << ')';
  }
  if (bigraphConnected(myEdges)) {
    return std::make_unique<NETWORK>(myEdges, aEprRv, true);
  }

  throw std::runtime_error("The GraphML network is not fully connected");
}

} // namespace qr
} // namespace uiiit
