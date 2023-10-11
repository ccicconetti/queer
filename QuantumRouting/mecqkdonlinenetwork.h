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

#include <limits>

namespace uiiit {
namespace qr {

enum class MecQkdOnlineAlgo {
  NotInitialized  = 0,
  Policy014_k1    = 11, //!< all applications follow the same path
  Policy014_k3    = 13, //!< same as above, but keep k=3 alternatives
  Policy015       = 20, //!< each application is assigned own path
  Policy015_reuse = 21, //!< same as above, but try to re-use allocations
};

std::vector<MecQkdOnlineAlgo> allMecQkdOnlineAlgos();
std::string                   toString(const MecQkdOnlineAlgo aAlgo);
MecQkdOnlineAlgo mecQkdOnlineAlgofromString(const std::string& aAlgo);

/**
 * @brief A QKD network where some nodes also offer edge computing resources.
 *
 * Applications enter/leave the system dynamically.
 */
class MecQkdOnlineNetwork final : public CapacityNetwork
{

  // used within allocate()
  struct EdgeNode {
    // initialized from theEdgeProcessing
    unsigned long theId        = 0;
    double        theAvailable = 0;

    // working variables
    double                     theResidual = 0;
    double                     thePathSize = 0;
    std::vector<unsigned long> thePath; //!< the path from user to edge

    bool feasiblePath() const noexcept {
      return thePathSize > 0;
    }
    bool feasibleResidual() const noexcept {
      return theResidual >= 0;
    }
    bool feasible() const noexcept {
      return feasiblePath() and feasibleResidual();
    }
    std::string toString() const;
  };
  using Candidates = std::vector<EdgeNode>;

 public:
  struct Allocation {
    // input
    unsigned long theUserNode = 0; //!< the user node originating the request
    double        theRate     = 0; //!< the QKD rates requested
    double        theLoad     = 0; //!< the amount of processing load requested

    // output
    uint64_t theId = BLOCKED; //!< app id, if allocated
    Path     thePath;         //!< the current path, if allocated

    explicit Allocation(const unsigned long aUserNode,
                        const double        aRate,
                        const double        aLoad);
    std::string toString() const;
    bool        allocated() const noexcept {
      return theId != BLOCKED;
    }
    double totRate() const noexcept {
      return allocated() ? (static_cast<double>(thePath.size()) * theRate) :
                           0.0;
    }
    unsigned long edgeNode() const noexcept {
      assert(not allocated() or not thePath.empty());
      return allocated() ? thePath.rbegin()->m_target : 0;
    }
  };

  /**
   * @brief Create a network with given links and weights
   *
   * @param aEdgeWeights The unidirectional edges and weights of the network
   * (src, dst, w).
   */
  explicit MecQkdOnlineNetwork(const WeightVector& aEdgeWeights);

  /**
   * @brief Configure the network.
   *
   * @param aAlgorithm The algorithm to be used for the allocation.
   * @param aUserNodes The set of user node identifies.
   * @param aEdgeProcessing key: edge node identifier, value: processing.
   * @param aRv A r.v. in [0,1] to break ties.
   *
   * Clear any previous data. Can only be called before the arrival of apps.
   *
   * @throws std::runtime_error if trying to reconfigure an operational network.
   */
  void configure(const MecQkdOnlineAlgo&                     aAlgorithm,
                 std::unique_ptr<support::RealRvInterface>&& aRv,
                 const std::set<unsigned long>&              aUserNodes,
                 const std::map<unsigned long, double>&      aEdgeProcessing);

  //! @return the current availability on the edge nodes.
  const std::map<unsigned long, double>& edgeNodes() const {
    return theEdgeProcessing;
  }

  //! @return the total processing power of the edge nodes.
  double totProcessing() const noexcept {
    return theTotProcessing;
  }

  //! @return the totale net rate of currently active applications.
  double totNetRate() const noexcept {
    return theTotNetRate;
  }

  /**
   * @brief Allocate a new application, if possible.
   *
   * @param aApps The parameters of the application to be allocated, also
   * providing the structure for the output and an identifier of this
   * application to delete it later (set to BLOCKED if not allocated).
   */
  void add(Allocation& aApp);

  /**
   * @brief Remove the allocated application with the given identifier.
   */
  void del(const uint64_t aAppId);

  //! @return the number of links affected by signalling changes.
  std::size_t signalling() const noexcept {
    return theSignalling;
  }

  static constexpr uint64_t BLOCKED = std::numeric_limits<uint64_t>::max();

 private:
  static constexpr double EPSILON = 1e-5; // to break ties arbitrarily

  // data members initialized in configure()
  MecQkdOnlineAlgo                          theAlgorithm;
  std::unique_ptr<support::RealRvInterface> theRv;
  std::set<unsigned long>                   theUserNodes;
  std::map<unsigned long, double>           theEdgeProcessing;
  std::set<unsigned long>                   theEdgeNodes;
  double                                    theHighestProcessing = 0;

  // internal data structure
  uint64_t                       theNextId     = 0;
  std::size_t                    theSignalling = 0;
  std::map<uint64_t, Allocation> theActiveApps;
  double                         theTotProcessing = 0;
  double                         theTotNetRate    = 0;

  // used only with Policy014
  // key: user node
  // value:
  //   key: edge node
  //   value: vector of paths (at most K elements)
  // map of K paths from all user nodes to all edge nodes
  std::map<unsigned long, std::map<unsigned long, std::vector<Path>>> thePaths;

  // Used with Policy014 to compute at most aK paths from user to edge nodes.
  void computeAllUserEdgePaths(const std::size_t aK);

  // Try to allocate app with Policy014.
  void addPolicy014(Allocation& aApp);

  // Allocate a given app to a path.
  void allocate(Allocation& aApp, const Path& aPath);

  // Print the paths.
  std::string pathsToString();
};

} // namespace qr
} // namespace uiiit
