#include "mat2gnuplot.h"

#include "mesh.h"
#include "dmatrix_denseCM.h"
#include "dmatrix_CCS.h"
#include "dLLT_denseCM.h"

// this example can be used to plot the nonzero patern of matrices.
// it produces a file for each of the matrix tested here :
//      bcsstk14.plt_commands.txt
//      square4.plt_commands.txt
//      square4_LLT.plt_commands.txt
// you can latter read any of these files with  gnuplot to produce a nice jpg views of the nonzeroterm of your matrix  typing for example :
//      gnuplot   < bcsstk14.plt_commands.txt
int main(){
  //dmatrix_denseCM A;
  {
    dmatrix_CCS A;
    std::string filename("data/bcsstk14.mtx");
    std::ifstream in(filename.c_str());
    if(!in.is_open() ) std::cout << "Can't Open File " << filename  << std::endl;
    in >> A;
    spy_ge(A, "bcsstk14.plt" );
  }
  {
    //std::string inputfile ="data/cube.msh";
    //std::string inputfile ="data/cube2.msh";
    //std::string inputfile ="data/square.msh";
    //std::string inputfile ="data/square2.msh";
    //std::string inputfile ="data/square3.msh";
    std::string inputfile ="data/square4.msh";
    
    std::ifstream meshin(inputfile);
    if (!meshin.is_open()) {std::cout << "Can't open mesh file " << inputfile << std::endl;   throw;}
    std::cout << "Reading Mesh from file " << inputfile << std::endl; 
    mesh me;
    meshin >> me;
    std::cout << me.nodes.size() << " "<< me.triangles.size() << std::endl;
    fem fem1(me);
    dmatrix_denseCM A = fem1.generateMassMatrix_denseCM();
    dLLT_denseCM    ALLT (A);
    dmatrix_denseCM LpLT = ALLT.getL() + transpose(ALLT.getL());

    spy_ge(A, "square4.plt" );
    spy_ge(LpLT, "square4_LLT.plt" );
    /*  
    spy_ge(A, "square3.plt" );
    spy_ge(LpLT, "square3_LLT.plt" );
    spy_ge(A, "square4.plt" );
    spy_ge(LpLT, "square4_LLT.plt" );
    spy_ge(A, "cube.plt" );
    spy_ge(LpLT, "cube_LLT.plt" );
    spy_ge(A, "cube2.plt" );
    spy_ge(LpLT, "cube2_LLT.plt" );
     spy_ge(A, "square.plt" );
    spy_ge(LpLT, "square_LLT.plt" );
    */
  }
  
}
