#ifndef _MESH_H
#define _MESH_H_

#include <vector>
#include <array>
#include <istream>
#include <string>
#include <iostream>
#include "dmatrix_denseCM.h"
#include "dmatrix_CCS.h"
#include "graph.h"

/// A class to store a mesh.
class mesh{
 public:
  /// a vector containing the nodes
  std::vector< std::array< double, 3> > nodes;
  /// a vector containing the node numbers of each triangle
  std::vector < std::array < int, 3 > > triangles; 
  /// a vector containing the node numbers of each tetrahedra
  std::vector < std::array < int, 4 > > tets; 
};

/// read a mesh from an input stream in .msh format
/*! (The mesh format as defined by gmsh, a pre and postprocessing meshing tool  )
    The mesh must contain only nodes, linear edges, triangles and tets. 
*/
std::istream& operator>>(std::istream& in, mesh& m);


/// This class is constructed from a mesh, and can return the associated mass matrix
/*! M = \int_mesh NiNj dv
  Where Ni is the linear shape function associated to node i 
*/
class fem{
 public: 
  typedef std::array<double,3> point;
  /// Constuctor from a mesh
  fem(const mesh &_m);
  /// Return the mass matrix in sparse Compressed Column Storage
  /*! if the mesh contain no tetrahedra, it return the mass matrix associated to all the triangles in the mesh.
      else it return the mass matrix associated to all the tetrahedra
   */
  dmatrix_CCS generateMassMatrix_CCS() const;
  /// Return the mass matrix in dense  column  major format.
  /*! if the mesh contain no tetrahedra, it return the mass matrix associated to all the triangles in the mesh.
    else it return the mass matrix associated to all the tetrahedra
  */
  dmatrix_denseCM generateMassMatrix_denseCM() const;
  
  double areaTriangle(const point &p0, const point& p1, const point &p2) const;
  double volumTet(const point &p0, const point& p1, const point &p2, const point & p3) const;
 private:
  const mesh &m;
  dmatrix_denseCM atri;
  dmatrix_denseCM atet;
  
};




#endif
