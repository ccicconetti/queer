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

/**
 * @brief A network where edges are characterized by their capacity.
 *
 * Links are directional and, in principle, the capacity can be different
 * for the two directions.
 */
class CapacityNetwork : public Network
{
 public:
  FRIEND_TEST(TestCapacityNetwork, test_min_capacity_edges);

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
  using Path             = std::list<EdgeDescriptor>;

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
   * @param aEdgeWeights The unidirectional edges and weights of the network
   * (src, dst, w).
   */
  explicit CapacityNetwork(const WeightVector& aEdgeWeights);

  //! \return the number of nodes.
  std::size_t numNodes() const;

  //! \return the number of edges.
  std::size_t numEdges() const;

  //! \return the min-max in-degree of the graph.
  std::pair<std::size_t, std::size_t> inDegree() const;

  //! \return the min-max out-degree of the graph.
  std::pair<std::size_t, std::size_t> outDegree() const;

  //! \return the total capacity across all the edges.
  double totalCapacity() const;

  //! Save to a dot file.
  void toDot(const std::string& aFilename) const;

  //! \return the current weights, one per element in the return vector.
  WeightVector weights() const;

  /**
   * @brief Find the shortest path, in number of hops, from a source node to all
   * possible destinations passed, by following only the paths that satisfy a
   * minimum capacity requirement.
   *
   * @param aSource The source node.
   * @param aCapacity The capacity requirement.
   * @param aDestinations The set of destinations.
   * @return std::map<unsigned long, std::vector<unsigned long>> The constrained
   * shortest path found for each node. If there is no feasible path or if the
   * destination node is the same as the source, then the value in the map is
   * empty.
   */
  std::map<unsigned long, std::vector<unsigned long>>
  cspf(const unsigned long            aSource,
       const double                   aCapacity,
       const std::set<unsigned long>& aDestinations) const;

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
   * @brief Find the aNum nodes closest to a source node.
   *
   * @param aSrc Source node.
   * @param aNum Number of closest nodes to find
   * @param aRv A r.v. in [0,1] to break ties.
   * @param aWhiteList Only return nodes from this list if non-empty.
   * @return std::vector<unsigned long> The closest nodes found.
   * @post the size of the returned vector is smaller than or equal to aNum
   */
  std::vector<unsigned long>
  closestNodes(const unsigned long            aSrc,
               const unsigned long            aNum,
               support::RealRvInterface&      aRv,
               const std::set<unsigned long>& aWhiteList =
                   std::set<unsigned long>()) const;

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
   * @brief Find the minimum capacity along a path in the network.
   *
   * @param aSrc The source node.
   * @param aPath The path to the destination.
   * @return double The minimum capacity found.
   * @throw std::runtime_error if the path does not exist in the network.
   */
  double minCapacity(const VertexDescriptor               aSrc,
                     const std::vector<VertexDescriptor>& aPath) const;

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

 protected:
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

  /**
   * @brief Find the minimum capacity along a path in a graph.
   *
   * @param aPath The path.
   * @param aGraph The network.
   * @return double The minimum capacity found.
   */
  static double minCapacity(const Path& aPath, const Graph& aGraph);

  static void
  removeSmallestCapacityEdge(const VertexDescriptor               aSrc,
                             const std::vector<VertexDescriptor>& aPath,
                             Graph&                               aGraph);

  /**
   * @brief Remove the capacity from all the edges along a path.
   *
   * @param aSrc The source node.
   * @param aPath The path.
   * @param aCapacity The capacity to be subtracted.
   * @param aMinCapacity If specified,remove the edge with a smaller residual.
   * @param aGraph The graph.
   */
  static void removeCapacityFromPath(const VertexDescriptor               aSrc,
                                     const std::vector<VertexDescriptor>& aPath,
                                     const double                aCapacity,
                                     const std::optional<double> aMinCapacity,
                                     Graph&                      aGraph);

  std::pair<std::size_t, std::size_t> minMaxVertexProp(
      const std::function<std::size_t(Graph::vertex_descriptor, const Graph&)>&
          aPropFunctor) const;

 protected:
  Graph theGraph;
};

} // namespace qr
} // namespace uiiit
