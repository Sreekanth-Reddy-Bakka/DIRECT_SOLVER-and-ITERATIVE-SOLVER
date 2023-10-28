#include "mesh.h"
#include <math.h>

std::istream& operator>>(std::istream& in, mesh& m){
  std::string tmp;
  in >> tmp;
  in >> tmp >> tmp >> tmp;
  in >> tmp;
  in >> tmp;
  int nnodes;
  in >> nnodes;
  m.nodes.reserve(nnodes);
  for (int i = 0; i < nnodes; ++i){
    //    int i;
    double x, y,z;
    in >> tmp >> x >> y >> z;
    m.nodes.push_back( std::array<double, 3>{x,y,z});
  }
  in >> tmp;
  in >> tmp;
  int nbelem;
  in >> nbelem;
  m.triangles.reserve(nbelem);
  for (int i =0; i < nbelem; ++i){
    int type, ntag,  i0, i1, i2, i3;
    in >> tmp >> type >> ntag;
    for (int tt = 0; tt <ntag ; ++tt) in >> tmp;
    switch (type){
    case 15: {in >> tmp; break;}
    case 1:  {in >> tmp >> tmp; break;}
    case 2:  { in >> i0 >> i1 >> i2; m.triangles.push_back(std::array<int, 3> {i0-1, i1-1,i2-1}); break; }
    case 4: {in >> i0>> i1 >> i2 >> i3; m.tets.push_back(std::array<int, 4> {i0-1, i1-1, i2-1, i3-1}); break; }
    default: {std::cout << "Untreated element type. " << __FILE__ << " : " << __LINE__ << std::endl; throw; break;}
    }
  }
  return in;
}

fem::fem(const mesh &_m):m(_m), atri(3,3,1.), atet(4,4,1.){
  atri(0,0) = 2.;
  atri(1,1) = 2.;
  atri(2,2) = 2.;
  atri/=24.;
  atet(0,0) = 2.;
  atet(1,1) = 2.;
  atet(2,2) = 2.;
  atet(3,3) = 2.;
  atet/=120.;
}
 

double fem::areaTriangle(const point &p0, const point& p1, const point &p2) const{
  point V01 ={{ p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] }};
  point V02 ={ p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2] };
  point N = { V01[1]*V02[2] - V01[2]*V02[1],
	      -V01[0]*V02[2] + V01[2]*V02[0],
	      V01[0]*V02[1] - V01[1]*V02[0]
  };
  return 0.5*sqrt( N[0]*N[0] + N[1]*N[1] + N[2]*N[2]);
}

double fem::volumTet(const point &p0, const point& p1, const point &p2, const point & p3) const{
  point V01 = {  p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] };
  point V02 = {  p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2] };
  point V03 = {  p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2] };
   
  point N = { V01[1]*V02[2] - V01[2]*V02[1],
	      -V01[0]*V02[2] + V01[2]*V02[0],
	      V01[0]*V02[1] - V01[1]*V02[0]
  };
  return fabs(N[0]*V03[0] + N[1]*V03[1]+ N[2]*V03[1])/6.;
}


dmatrix_denseCM fem::generateMassMatrix_denseCM() const{
    int nb_nodes = m.nodes.size();
    int M = nb_nodes;
    dmatrix_denseCM A(M,M,0.);
    if (m.tets.empty()){
      for(auto e_nodes : m.triangles){
	point X0 = m.nodes[e_nodes[0]];
	point X1 = m.nodes[e_nodes[1]];
	point X2 = m.nodes[e_nodes[2]];
	double s = areaTriangle(X0,X1,X2);
	// std::cout << " s " <<  s << " "; // << e_nodes[0] << " " <<e_nodes[1] << " "<< e_nodes[2] <<std::endl;
	for (int i = 0; i < 3; ++i)
	  for (int j = 0; j < 3; ++j)
	    A(e_nodes[i], e_nodes[j] ) += s*atri(i,j);
      }
    }
    else{ //3d case
      for(auto e_nodes : m.tets){
	point X0 = m.nodes[e_nodes[0]];
	point X1 = m.nodes[e_nodes[1]];
	point X2 = m.nodes[e_nodes[2]];
	point X3 = m.nodes[e_nodes[3]];
	
	double v = volumTet(X0,X1,X2,X3);
	for (int i = 0; i < 4; ++i)
	  for (int j = 0; j < 4; ++j)
	    A(e_nodes[i], e_nodes[j] ) += v*atet(i,j);
      }
    }
    return A;
}
 
dmatrix_CCS fem::generateMassMatrix_CCS() const{
    int nb_nodes = m.nodes.size();
    int M = nb_nodes;
    graph g;
    g.addNodes(nb_nodes);
    if (m.tets.empty()){
    for(auto e_nodes : m.triangles){
      for (auto i : e_nodes){
	for (auto j : e_nodes){
	  g.addEdge(i,j);
	  g.addEdge(j,i);
	}
      }
    }
    }
    else{ //3d case
     for(auto e_nodes : m.tets){
      for (auto i : e_nodes){
	for (auto j : e_nodes){
	  g.addEdge(i,j);
	  g.addEdge(j,i);
	}
      }
     } 
    }
    int NNZ = g.nbEdges();
    dmatrix_CCS A(M,M,NNZ);
    compressGraph(g, A.columnptr, A.lineindex);
    
    
    
    if (m.tets.empty()){
      for(auto e_nodes : m.triangles){
	point X0 = m.nodes[e_nodes[0]];
	point X1 = m.nodes[e_nodes[1]];
	point X2 = m.nodes[e_nodes[2]];
	double s = areaTriangle(X0,X1,X2);
    	for (int i = 0; i < 3; ++i)
	  for (int j = 0; j < 3; ++j)
	    A(e_nodes[i], e_nodes[j] ) += s*atri(i,j);
      }
    }

    else{ //3d case 
      for(auto e_nodes : m.tets){
	point X0 = m.nodes[e_nodes[0]];
	point X1 = m.nodes[e_nodes[1]];
	point X2 = m.nodes[e_nodes[2]];
	point X3 = m.nodes[e_nodes[3]];
	double v = volumTet(X0,X1,X2, X3);
      //    std::cout << s << std::endl;
	for (int i = 0; i < 4; ++i)
	  for (int j = 0; j < 4; ++j)
	    A(e_nodes[i], e_nodes[j] ) += v*atet(i,j);
      }
    }
    
    return A;
 }


 

