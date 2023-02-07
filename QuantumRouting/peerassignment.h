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

#include "QuantumRouting/esnetwork.h"
#include "Support/macros.h"
#include "Support/random.h"

#include <memory>
#include <string>
#include <vector>

namespace uiiit {
namespace qr {

//
// peer assignment algorithms
//

enum class PeerAssignmentAlgo : unsigned int {
  Random        = 0,
  ShortestPath  = 1,
  LoadBalancing = 2,
};

std::vector<PeerAssignmentAlgo> allPeerAssignmentAlgos();
std::string                     toString(const PeerAssignmentAlgo aAlgo);
PeerAssignmentAlgo peerAssignmentAlgofromString(const std::string& aAlgo);

//
// factory free function to create peer assignment instances
//

class PeerAssignment;
/**
 * @brief Construct a peer assignment instance bound to a given quantum network.
 *
 * @param aNetwork The reference quantum network for peer assignment.
 * @param aAlgo The assignment algorithm to be used.
 * @param aRv A random variable in [0,1] that might be used by the algorithm.
 * @param aCheckFunction The function that checks if a given path is valid.
 * @return std::unique_ptr<PeerAssignment> The peer assignment object.
 * @throw std::runtime_error if aAlgo is not known.
 */
std::unique_ptr<PeerAssignment>
makePeerAssignment(const EsNetwork&                   aNetwork,
                   const PeerAssignmentAlgo           aAlgo,
                   support::RealRvInterface&          aRv,
                   const EsNetwork::AppCheckFunction& aCheckFunction);

//
// peer assignment classes
//

class PeerAssignment
{
 public:
  struct AppDescriptor {
    AppDescriptor(const unsigned long aHost,
                  const double        aPriority,
                  const double        aFidelityThreshold) noexcept;

    const unsigned long theHost;     //!< the vertex that hosts the computation
    const double        thePriority; //!< priority weight
    const double        theFidelityThreshold; //!< fidelity threshold
  };

  virtual ~PeerAssignment() = default;

  /**
   * @brief Assign peers to the given apps. Pure virtual function.
   *
   * @param aApps The input applications.
   * @param aNumPeers The number of peers to be selected for each app.
   * @param aCandidatePeers The nodes that it is possible to select as peers.
   * @return std::vector<EsNetwork::AppDescriptor>
   */
  virtual std::vector<EsNetwork::AppDescriptor>
  assign(const std::vector<AppDescriptor>& aApps,
         const unsigned long               aNumPeers,
         const std::vector<unsigned long>& aCandidatePeers) = 0;

  //! @return The assignment algoritm.
  PeerAssignmentAlgo algo() const noexcept {
    return theAlgo;
  }

 protected:
  /**
   * @brief Construct a new Peer Assignment bound to a given quantum network.
   *
   * @param aNetwork The reference quantum network for peer assignment.
   * @param aAlgo The assignment algorithm to be used.
   */
  PeerAssignment(const EsNetwork& aNetwork, const PeerAssignmentAlgo aAlgo);

  void throwIfDuplicates(const std::vector<unsigned long>& aCandidatePeers);

 protected:
  const EsNetwork&         theNetwork;
  const PeerAssignmentAlgo theAlgo;
};

class PeerAssignmentRandom final : public PeerAssignment
{
 public:
  PeerAssignmentRandom(const EsNetwork&          aNetwork,
                       support::RealRvInterface& aRv);

  //! Assign peers choosing at random.
  std::vector<EsNetwork::AppDescriptor>
  assign(const std::vector<AppDescriptor>& aApps,
         const unsigned long               aNumPeers,
         const std::vector<unsigned long>& aCandidatePeers) override;

 private:
  support::RealRvInterface& theRv;
};

class PeerAssignmentShortestPath final : public PeerAssignment
{
 public:
  PeerAssignmentShortestPath(const EsNetwork&          aNetwork,
                             support::RealRvInterface& aRv);

  //! Pick peers from closest peers.
  std::vector<EsNetwork::AppDescriptor>
  assign(const std::vector<AppDescriptor>& aApps,
         const unsigned long               aNumPeers,
         const std::vector<unsigned long>& aCandidatePeers) override;

 private:
  support::RealRvInterface& theRv;
};

class PeerAssignmentLoadBalancing final : public PeerAssignment
{
 public:
  PeerAssignmentLoadBalancing(
      const EsNetwork&                   aNetwork,
      const EsNetwork::AppCheckFunction& aCheckFunction);

  //! Assign peers as the result of a generalized assignment problem.
  std::vector<EsNetwork::AppDescriptor>
  assign(const std::vector<AppDescriptor>& aApps,
         const unsigned long               aNumPeers,
         const std::vector<unsigned long>& aCandidatePeers) override;

 private:
  const EsNetwork::AppCheckFunction theCheckFunction;
};

} // namespace qr
} // namespace uiiit
