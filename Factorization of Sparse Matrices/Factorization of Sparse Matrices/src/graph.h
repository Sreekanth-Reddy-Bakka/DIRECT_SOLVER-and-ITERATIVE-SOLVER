#ifndef _graph_
#define _graph_
#include <vector>
#include <set>
#include <list>
#include "dmatrix_denseCM.h"
#include "dmatrix_CCS.h"



class graph{
 public:
  typedef  std::set<int > neighbour;
  typedef   std::vector< neighbour > graph_storage;
  typedef   neighbour::const_iterator neighbour_iterator;
  /// add a node in the graph. return the node number
  int addNode();
  /// add n nodes in the graph, return the node number of the last added node
  int addNodes(int n);
  /// add an edge in the graph. nodes i and j are supposed already in the graph
  void addEdge(int i, int j);
  /// return all the direct neighbour of vertex i
  const neighbour & getNeighbour(int i) const;
  /// return iterator to the begining of the first neighbour of vertex i
  neighbour_iterator beginNeighbour( int i)const;
  neighbour_iterator endNeighbour( int i) const;
  /// clear the graph off all nodes and edges
  void clear();
  /// return the nuber of nodes in the graph
  int nbNodes() const;
  /// return the number of edges
  int nbEdges() const;
 private:
  graph_storage data;
};



/// Explore the graph starting at node start_node using breath first
//  return a vector containing all the visited nodes, in visiting order
std::vector<int > breathFirstSearch(const graph & g, int start_node);


/// Explore the graph starting at node start_node using depth first
//  return a vector containing all the visited nodes, in visiting order
std::vector<int > depthFirstSearch(const graph & g, int start_node);


/// return a new ordering of the vertices of the graph, such as the bandwith is smaller 
// using the cuthill McKee ordering scheme. You have to implement this function
std::vector<int > cuthillmckee(const graph &g);

 
/// compress the graph int  vector of int describing all the connectivity in a format similar to compress column storage
void compressGraph( const graph &graph, std::vector<int > &cpnt, std::vector<int> &lindex );


/// build the line graph of the input matrix( in dense storage). if |Aij| > eps, edgeij is added to the graph. 
graph buildLineGraph(const dmatrix_denseCM &A, double eps =0.);

/// build the column graph of the input matrix( in dense storage). if |Aij| > eps, edgeji is added to the graph. 
graph buildColumnGraph(const dmatrix_denseCM &A, double eps = 0.);

/// build the symetric graph of the input matrix( in dense storage). if |Aij| > eps, edge ij and edgeji is added to the graph. 
graph buildSymGraph(const dmatrix_denseCM &A, double eps = 0.);

/// build the line graph of the input matrix( in sparse storage). if Aij stored in A, edgeij added to the graph 
graph buildLineGraph(const dmatrix_CCS &A);

/// build the column graph of the input matrix( in sparse storage). if Aij stored in A, edgeji added to the graph 
graph buildColumnGraph(const dmatrix_CCS &A);

/// build the line graph of the input matrix( in sparse storage). if Aij stored in A, edgeij and edgeji added to the graph 
graph buildSymGraph(const dmatrix_CCS &A);


#endif
