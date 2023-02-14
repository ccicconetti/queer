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
  Random       = 0, //!< pick a random edge node, among those feasible
  Spf          = 1, //!< shortest path first, among those feasible
  BestFit      = 2, //!< best-fit allocation, among those feasible
  RandomBlind  = 3, //!< pick a random edge node
  SpfBlind     = 4, //!< shortest path first
  BestFitBlind = 5, //!< best-fit allocation
  SpfStatic    = 6, //!< assign based on initial shortest-path
};

std::vector<MecQkdAlgo> allMecQkdAlgos();
std::string             toString(const MecQkdAlgo aAlgo);
MecQkdAlgo              mecQkdAlgofromString(const std::string& aAlgo);

class MecQkdWorkload final
{
 public:
  struct AppInfo {
    unsigned long theRegion = 0; //!< region identifier
    double        theWeight = 0; //!< weight to determine occurence probability
    double        theLoad   = 0; //!< load requested, in operations/second
    double        theRate   = 0; //!< QKD rate requested, in b/s

    std::string toString() const;
  };

  /**
   * @brief Construct a new MecQkdWorkload object from its app info data.
   *
   * @param aAppInfo The info to be copied into this data structure.
   * @param aRv The r.v. to draw randomly the apps.
   *
   * @throw std::runtime_error if no apps are passed.
   */
  MecQkdWorkload(const std::vector<AppInfo>& aAppInfo,
                 support::RealRvInterface&   aRv);

  /**
   * @brief Make a new MecQkdWorkload object from a csv file.
   *
   * @param aFilename The csv file where to read the data from.
   * @param aRv The r.v. to draw randomly the apps.
   * @return MecQkdWorkload
   */
  static MecQkdWorkload fromCsvFile(const std::string&        aFilename,
                                    support::RealRvInterface& aRv);

  /**
   * @brief Draw randomly an app.
   */
  AppInfo operator()();

  /**
   * @brief Return the regions in the app infos.
   *
   * @return std::vector<unsigned long>
   */
  const std::set<unsigned long>& regions() const {
    return theRegions;
  }

 private:
  const std::vector<AppInfo> theAppInfo;
  support::RealRvInterface&  theRv;
  std::set<unsigned long>    theRegions; // never changed after construction
  std::vector<double>        theWeights; // same
};

/**
 * @brief A QKD network where some nodes also offer edge computing resources.
 */
class MecQkdNetwork final : public CapacityNetwork
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
    bool          theAllocated = false; //!< true if the user has been allocated
    unsigned long theEdgeNode  = 0;     //!< the target edge user node assigned
    std::size_t   thePathLength = 0;    //!< the path length from user to edge

    explicit Allocation(const unsigned long aUserNode,
                        const double        aRate,
                        const double        aLoad);
    std::string toString() const;
    double      totRate() const noexcept {
      if (theAllocated) {
        return thePathLength * theRate;
      }
      return 0.0;
    }
  };

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

  /**
   * @brief Set the identifiers of user nodes.
   *
   * Clear any previous data.
   */
  void userNodes(const std::set<unsigned long>& aUserNodes);

  /**
   * @brief Set the edge nodes' identifiers and assign processing power values.
   *
   * Clear any previous data.
   *
   * @param aEdgeProcessing key: edge node identifier, value: processing
   */
  void edgeNodes(const std::map<unsigned long, double>& aEdgeProcessing);

  //! @return the current availability on the edge nodes.
  const std::map<unsigned long, double>& edgeNodes() const {
    return theEdgeProcessing;
  }

  /**
   * @brief Allocate the given user requests into the MEC/QKD resources.
   *
   * @param aApps The list of applications to be allocated, also providing the
   * structure for the output.
   * @param aAlgo The algorithm to be used.
   * @param aRv A r.v. in [0,1] to break ties.
   */
  void allocate(std::vector<Allocation>&  aApps,
                const MecQkdAlgo          aAlgo,
                support::RealRvInterface& aRv);

  //! @return the total processing power of the edge nodes.
  double totProcessing() const;

 private:
  std::set<unsigned long>         theUserNodes;
  std::map<unsigned long, double> theEdgeProcessing;
  std::set<unsigned long>         theEdgeNodes;

  /**
   * @brief Select the candidate edge node to be assigned.
   *
   * @param aCandidates The possible edge nodes.
   * @param aAlgo The algorithm used.
   * @param aRv A r.v. in [0,1] to break ties.
   * @return std::vector<EdgeNode>::const_iterator
   */
  static Candidates::iterator
  selectCandidate(std::vector<EdgeNode>&    aCandidates,
                  const MecQkdAlgo          aAlgo,
                  support::RealRvInterface& aRv);

  /**
   * @brief Allocate using MecQkdAlgo::SpfStatic.
   *
   * @param aApps The applications.
   */
  void allocateSpfStatic(std::vector<Allocation>& aApps);
};

} // namespace qr
} // namespace uiiit
