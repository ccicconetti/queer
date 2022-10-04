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

#include "QuantumRouting/network.h"
#include "QuantumRouting/qrutils.h"
#include "Support/random.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/detail/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/property_map.hpp>

#include <cinttypes>
#include <functional>
#include <list>
#include <map>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

namespace uiiit {
namespace qr {

enum class AppRouteAlgo {
  Random  = 0,
  BestFit = 1,
  Drr     = 2, // https://arxiv.org/abs/2203.05844v2
};

std::vector<AppRouteAlgo> allAppRouteAlgos();
std::string               toString(const AppRouteAlgo aAlgo);
AppRouteAlgo              appRouteAlgofromString(const std::string& aAlgo);

/**
 * @brief A quantum network where edges are characterized by their capacity
 * only, in terms of EPR pairs that can be generated by per second.
 *
 * Links are directional and, in principle, the capacity can be different
 * for the two directions.
 *
 * It is possible to route two types of resources:
 *
 * - flows: they are characterized by a source, a destination and a net EPR
 * rate; they represent metrology, sensing, and QKD applications that require a
 * constant rate of end-to-end entangled pairs
 *
 * - apps: they are characterized by a host node and a number of peers, as well
 * as a numeric priority; they represent elastic applications, e.g., for
 * distributed quantum computing
 */
class CapacityNetwork final : public Network
{
 public:
  using Graph =
      boost::adjacency_list<boost::listS,
                            boost::vecS,
                            boost::bidirectionalS,
                            boost::no_property,
                            boost::property<boost::edge_weight_t, double>,
                            boost::no_property,
                            boost::listS>;
  using VertexDescriptor = boost::graph_traits<Graph>::vertex_descriptor;
  using EdgeDescriptor   = boost::graph_traits<Graph>::edge_descriptor;

  struct FlowDescriptor {
    FlowDescriptor(const unsigned long aSrc,
                   const unsigned long aDst,
                   const double        aNetRate) noexcept;

    void movePathRateFrom(FlowDescriptor& aAnother);

    // input
    const unsigned long theSrc;     //!< the source vertex
    const unsigned long theDst;     //!< the destination vertex
    const double        theNetRate; //!< in EPR/s

    // output
    std::vector<unsigned long> thePath;      //!< hops not including src
    double                     theGrossRate; //!< in EPR/s
    std::size_t                theDijsktra;  //!< number of times called

    std::string toString() const;
  };

  struct AppDescriptor {
    using Path = std::list<EdgeDescriptor>;
    using Hops = std::vector<unsigned long>;

    AppDescriptor(const unsigned long               aHost,
                  const std::vector<unsigned long>& aPeers,
                  const double                      aPriority,
                  const double aFidelityThreshold) noexcept;

    // input
    const unsigned long theHost; //!< the vertex that hosts the computation
    const std::vector<unsigned long>
                 thePeers;             //!< the possible entanglement peers
    const double thePriority;          //! weight
    const double theFidelityThreshold; //! fidelity threhold

    // working variables
    std::map<unsigned long, std::list<Path>> theRemainingPaths;

    // output
    struct Output {
      explicit Output(const Path& aPath);

      Hops   theHops;
      double theNetRate   = 0; //!< in EPR-pairs/s
      double theGrossRate = 0; //!< in EPR-pairs/s
    };
    std::map<unsigned long, std::vector<Output>>
                theAllocated; //!< key: destination
    std::size_t theVisits;    //!< number of visits

    double      netRate() const;
    double      grossRate() const;
    std::string toString() const;

   private:
    double accumulate(
        const std::function<double(const AppDescriptor::Output&)>& aFn) const;
  };

  using FlowCheckFunction = std::function<bool(const FlowDescriptor&)>;
  using AppCheckFunction =
      std::function<bool(const AppDescriptor&, const AppDescriptor::Path&)>;

  // vector of (src, dst)
  using EdgeVector = std::vector<std::pair<unsigned long, unsigned long>>;
  // vector of (src, dst, weight)
  using WeightVector =
      std::vector<std::tuple<unsigned long, unsigned long, double>>;
  // map of host -> { peers }
  using ReachableNodes = std::map<unsigned long, std::set<unsigned long>>;

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
  explicit CapacityNetwork(const EdgeVector&         aEdges,
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
  explicit CapacityNetwork(const WeightVector& aEdgeWeights);

  /**
   * @brief Set the measurement probability.
   *
   * @param aMeasurementProbability the new measurement probability
   *
   * @throw std::runtime_error if the measurement probability is not in [0,1]
   */
  void measurementProbability(const double aMeasurementProbability);

  //! \return the measurement probability.
  double measurementProbability() const noexcept {
    return theMeasurementProbability;
  }

  //! \return the number of nodes.
  std::size_t numNodes() const;

  //! \return the number of edges.
  std::size_t numEdges() const;

  //! \return the min-max in-degree of the graph.
  std::pair<std::size_t, std::size_t> inDegree() const;

  //! \return the min-max out-degree of the graph.
  std::pair<std::size_t, std::size_t> outDegree() const;

  //! \return the total EPR capacity across all the edges.
  double totalCapacity() const;

  //! Save to a dot file.
  void toDot(const std::string& aFilename) const;

  //! \return the current weights, one per element in the return vector.
  WeightVector weights() const;

  /**
   * @brief For each node, find the reachable nodes with min/max distance.
   *
   * @param aMinHops The minimum distance, in hops.
   * @param aMaxHops The maximum distance, in hops.
   * @param aDiameter Return the network diameter, in hops.
   * @return std::map<unsigned long, std::set<unsigned long>>
   */
  ReachableNodes reachableNodes(const std::size_t aMinHops,
                                const std::size_t aMaxHops,
                                std::size_t&      aDiameter) const;

  /**
   * @brief Route the given flows in this network starting with current
   * capacities.
   *
   * The flows are routed one by one in the order in which they are passed.
   * The capacities are updated whenever a flow can be admitted, in which case
   * the corresponding descriptor is also updated with routing info.
   *
   * @param aFlows the flows to be routed (admitted flows are modified)
   * @param aCheckFunction the flow is considered feasible only if this
   * function returns true, otherwise it is inadmissible; the default is to
   * always accept the flow
   *
   * @throw std::runtime_error if aFlows contain an ill-formed request, in which
   * case we guarantee that the internal state is not changed
   */
  void route(
      std::vector<FlowDescriptor>& aFlows,
      const FlowCheckFunction&     aCheckFunction = [](const auto&) {
        return true;
      });

  /**
   * @brief Route the given elastic applications in the network.
   *
   * @param aApps the applications to be routed
   * @param aAlgo the algorithm to be used for resource allocation
   * @param aQuantum  the allocation quantum to be used
   * @param aK the maximum number of paths to be found for each app and peer
   * @param aCheckFunction the flow is considered feasible only if this
   * function returns true, otherwise it is inadmissible; the default is to
   * always accept the flow
   *
   * @throw std::runtime_error if aFlows contain an ill-formed request, in which
   * case we guarantee that the internal state is not changed
   */
  void route(
      std::vector<AppDescriptor>& aApps,
      const AppRouteAlgo          aAlgo,
      const double                aQuantum,
      const std::size_t           aK,
      const AppCheckFunction& aCheckFunction = [](const auto&, const auto&) {
        return true;
      });

  /**
   * @brief Add capacity on all the edges along a given path from a source node.
   *
   * @param aSrc the source node
   * @param aPath the path
   * @param aCapacity the capacity to be added
   *
   * @throw std::runtime_error if one of the edges in the path does not exist in
   * the graph
   */
  void addCapacityToPath(const VertexDescriptor               aSrc,
                         const std::vector<VertexDescriptor>& aPath,
                         const double                         aCapacity);
  /**
   * @brief Return the capacity for each node, defined as the sum of the
   * capacity of all its outgoing edges.
   *
   * @return std::vector<double> The capacity value for each node index.
   */
  std::vector<double> nodeCapacities() const;

  /**
   * @brief Save the nodes and vertices of the graph to two files that can be
   * read by Gnuplot.
   *
   * @param aFilename The base filename. If empty, then do nothing.
   * @param aCoordinates The coordinates to be printed for the nodes.
   */
  void toGnuplot(const std::string&             aFilename,
                 const std::vector<Coordinate>& aCoordinates) const;

 private:
  struct HopsFinder {
    HopsFinder(const std::vector<VertexDescriptor>& aPredecessors,
               const VertexDescriptor               aSource)
        : thePredecessors(aPredecessors)
        , theSource(aSource) {
      // noop
    }
    void operator()(std::vector<VertexDescriptor>& aHops,
                    const VertexDescriptor         aNext) {
      if (aNext == theSource) {
        return;
      }
      (*this)(aHops, thePredecessors[aNext]);
      aHops.push_back(aNext);
    }
    const std::vector<VertexDescriptor>& thePredecessors;
    const VertexDescriptor               theSource;
  };

  static bool checkCapacity(const VertexDescriptor               aSrc,
                            const std::vector<VertexDescriptor>& aPath,
                            const double                         aCapacity,
                            const Graph&                         aGraph);

  static void
  removeSmallestCapacityEdge(const VertexDescriptor               aSrc,
                             const std::vector<VertexDescriptor>& aPath,
                             Graph&                               aGraph);

  static void removeCapacityFromPath(const VertexDescriptor               aSrc,
                                     const std::vector<VertexDescriptor>& aPath,
                                     const double aCapacity,
                                     Graph&       aGraph);

  std::pair<std::size_t, std::size_t> minMaxVertexProp(
      const std::function<std::size_t(Graph::vertex_descriptor, const Graph&)>&
          aPropFunctor) const;

  //! \return the gross rate for a given path length, in num of edges.
  double toGrossRate(const double aNetRate, const std::size_t aNumEdges) const;

  //! \return the net rate for a given path length, in num of edges.
  double toNetRate(const double aGrossRate, const std::size_t aNumEdges) const;

  //! Resource allocation of apps using DRR.
  void routeDrr(
      std::vector<AppDescriptor>& aApps,
      const double                aQuantum,
      const std::size_t           aK,
      const AppCheckFunction& aCheckFunction = [](const auto&, const auto&) {
        return true;
      });

 private:
  Graph  theGraph;
  double theMeasurementProbability;
};

} // namespace qr
} // namespace uiiit
