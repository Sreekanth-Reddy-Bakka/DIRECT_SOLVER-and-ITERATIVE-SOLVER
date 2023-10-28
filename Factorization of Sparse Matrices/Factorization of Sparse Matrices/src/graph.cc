#include "graph.h"
#include <cassert>
#include <queue>
#include <math.h>
#include <random>

int graph::addNode(){
    data.push_back( std::set< int > ());
    return data.size()-1;
  };
int graph::addNodes(int n){
    int id;
    for( int i = 0; i < n ; ++i)   id = addNode();
    return id;
  };
  
void graph::addEdge(int i, int j){
    assert( i< data.size() && j < data.size());
    data[i].insert(j);
  }
const graph::neighbour & graph::getNeighbour(int i) const{
    assert (i < data.size());
    return data[i];
}
graph::neighbour_iterator graph::beginNeighbour( int i)const {
  return getNeighbour(i).begin();
}

graph::neighbour_iterator graph::endNeighbour( int i) const{
  return getNeighbour(i).end();
}

void graph::clear(){
  data.clear();
}

int graph::nbNodes()const{
  return data.size();
}

int graph::nbEdges()const{
  int ne = 0;
  for (int i = 0; i< nbNodes() ; ++i){
    ne += getNeighbour(i).size();
  } 
  return ne;
}

void df (int i, const graph &g, std::vector<int> & visited, std::vector<int > & order ){
  visited[i] = 1;
  order.push_back(i);
  std::list<int > neighbours;
  neighbours.insert(neighbours.begin(), g.beginNeighbour(i), g.endNeighbour(i));
  while(!neighbours.empty()){
    i = neighbours.front();
    if (!visited[i]) df(i, g, visited, order);
    neighbours.pop_front();
  }
}

std::vector<int > depthFirstSearch(const graph & g, int start_node){
  assert (start_node < g.nbNodes());
  std::vector<int > order;
  std::vector<int>  visited( g.nbNodes(), 0);
  df(start_node, g, visited, order);
  return order;
}

/// Explore the graph starting at node start_node.
//  return a vector containing all the visited nodes, in visiting order
std::vector<int > breathFirstSearch(const graph & g, int start_node){
  assert (start_node < g.nbNodes());
  std::vector<int > order;
  std::vector<int>  visited( g.nbNodes(), 0);
  std::queue<int > to_treat;
  to_treat.push(start_node);
  while(!to_treat.empty()){
    int i = to_treat.front();
    to_treat.pop();
    if (!visited[i]){
      visited[i] = 1;
      order.push_back(i);
      for(auto j : g.getNeighbour(i))  if(!visited[j]) to_treat.push(j);
    }
  }
  return order;
}

int give_minimum(const std::vector<int>& degree, const std::vector<int>& visited) 
{
	std::vector<int> not_visited;

	for (int i = 0; i < degree.size(); ++i) {
		if (visited[i] != 1) {
			not_visited.push_back(i);
		}
	}

	std::sort(not_visited.begin(), not_visited.end(), [&degree](int a, int b) {
		return degree[a] < degree[b];
	});
	if (!not_visited.empty()) {
		return not_visited[0]; 
	}

	return -1; 
}

std::vector<int > cuthillmckee(const graph &g)
{
	std::vector<int> order;
	std::queue<int> to_treat;

	int n = g.nbNodes();
	std::vector<int> degree(n, 0);
	std::vector<int> visited(n, 0);

	for (int i = 0; i < n; ++i) {
		degree[i] = g.getNeighbour(i).size();
	}

	while (order.size() != n) {
		int start_node = give_minimum(degree, visited);
		to_treat.push(start_node);

		while (!to_treat.empty()) {
			std::vector<int> unvisited_neighbors;
			int i = to_treat.front();
			to_treat.pop();

			if (!visited[i]) {
				visited[i] = 1;
				order.push_back(i);

				for (int neighbor : g.getNeighbour(i)) {
					if (!visited[neighbor]) {
						unvisited_neighbors.push_back(neighbor);
					}
				}

				// Sort in descending order of degree (maximum degree first).
				std::sort(unvisited_neighbors.begin(), unvisited_neighbors.end(), [&degree](int a, int b) {
					return degree[a] > degree[b];
				});

				for (int neighbor : unvisited_neighbors) {
					to_treat.push(neighbor);
				}
			}
        	}	
	}	
    return order;
}



/// return a new ordering of the vertex, such as the bandwith is smaller 
// using the cuthill McKee ordering scheme. You have to implement this function
// You can start by taking inspiration from the given breathFirstSearch implementation.
/* 
   You'll need some way to sort a vector in some assending order.
   example :
   We suggest to use algorithm from the standard c++ library.
   std::vector<int> a;
   //  .... a is filled with some values
   //  then a can be sorted like this :
   std::sort(a.begin(), a.end());
   // in the previous case a is sorted such as a[i] < a[i+1]
   // you can custumise the previous sort such as a is sorted such as 
   // f(a[i], a[i+1]) = true
   // to do that you must pass f to the sorting algorithm.
   // The easiest way to do that is to define a lambda function for f.
   // and pass it to the sort algorithm :
   std::sort(a.begin(), a.end(), [&weight](int a, int b){ return weight[a] < weight[b];} );
   // in the previous example, the vector weight is passed to the lambda function     and a is considered less than b if weight[a] < weight[b]
   

*/

   
void compressGraph( const graph &graph, std::vector<int > &cpnt, std::vector<int> &lindex ){
  int nv = graph.nbNodes();
  cpnt.resize(nv+1);
  cpnt[0] =0;
  int nnz =0;
  for (int i =0; i < nv; ++i){
    auto neibi = graph.getNeighbour(i);
    int ni = neibi.size();
    nnz += ni;
    cpnt[i+1] = cpnt[i] + neibi.size();
  }
  lindex.resize(nnz);
  int k =0;
  auto li = lindex.begin();
  for (int i =0; i < nv; ++i){
    auto neibi = graph.getNeighbour(i);
    int ni = neibi.size();
    std::copy(neibi.begin(), neibi.end(), li );
    li += ni; 
  }
}



graph buildLineGraph(const dmatrix_denseCM &A, double eps){
  int m = A.getNbLines();
  assert(m == A.getNbColumns());
  graph g;
  g.addNodes(m);
  for(int i = 0; i < m ; ++i){
    for(int j = 0; j < m ; ++j){
      if ( fabs(A(i,j)) > eps ) g.addEdge(i,j);
    }
  }
  return g;
}

graph buildColumnGraph(const dmatrix_denseCM &A, double eps){
  int m = A.getNbLines();
  assert(m == A.getNbColumns());
  graph g;
  g.addNodes(m);
  for(int i = 0; i < m ; ++i){
    for(int j = 0; j < m ; ++j){
      if (fabs (A(i,j)) > eps) g.addEdge(j,i);
    }
  }
  return g;
}

graph buildSymGraph(const dmatrix_denseCM &A){
  int m = A.getNbLines();
  assert(m == A.getNbColumns());
  graph g;
  g.addNodes(m);
  for(int i = 0; i < m ; ++i){
    for(int j = 0; j < m ; ++j){
      if (  fabs(A(i,j)) > 0. ) {
	g.addEdge(i,j);
	g.addEdge(j,i);
      }
    }
  }
  return g;
}


graph buildLineGraph(const dmatrix_CCS &A){
  int m = A.m;
  int n = A.n;
  graph gA;
  gA.addNodes(m);
  for (int j = 0; j < n; ++j){
    for(int k = A.columnptr[j]; k < A.columnptr[j+1]; ++k ){
      int i = A.lineindex[k];
      gA.addEdge(i,j);
    }
  }
  return gA;

}

graph buildColumnGraph(const dmatrix_CCS &A){
  int m = A.m;
  int n = A.n;
  graph gA;
  gA.addNodes(n);
  for (int j = 0; j < n; ++j){
    for(int k = A.columnptr[j]; k < A.columnptr[j+1]; ++k ){
      int i = A.lineindex[k];
      gA.addEdge(j,i);
    }
  }
  return gA;
}

graph buildSymGraph(const dmatrix_CCS &A){
  int m = A.m;
  int n = A.n;
  assert(m==n);
  graph gA;
  gA.addNodes(m);
  for (int j = 0; j < n; ++j){
    for(int k = A.columnptr[j]; k < A.columnptr[j+1]; ++k ){
      int i = A.lineindex[k];
      gA.addEdge(i,j);
      gA.addEdge(j,i);
    }
  }
  return gA;
}
