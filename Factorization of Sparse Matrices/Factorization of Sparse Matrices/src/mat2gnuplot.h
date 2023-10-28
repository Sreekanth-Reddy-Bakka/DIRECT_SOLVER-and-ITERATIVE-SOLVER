#ifndef _mat2gnuplot_
#define _mat2gnuplot_

#include <string>
#include <fstream>
#include <iostream>
using std::string;
using std::ofstream;
using std::cout;

//****************************************************************************80
template < class MAT>
void spy_ge ( const MAT &A, string header )
  
//****************************************************************************80
//
//  Purpose:
//
//    SPY_GE plots a sparsity pattern for a  matrix of template type MAT
//     prerequist MAT class as member functions :
//      double operator (int i, int j ) const
//          return the i,j term of the matrix (0 based indexing)
//      int getNbLines() const :
//          return the number of line of the matrix
//      int getNbColumns() const :
//          return the number of columns in the matrix
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 Novembre 2018
//
//  Author:
//
//    John Burkardt
//    Modified by Nicolas Chevaugeon
//
//  Parameters:
//
//
//    Input, const MAT & the matrix.
//
//    Input, string HEADER, the name to be used for the
//    title of the plot, and as part of the names of the data, command
//    and plot files.
//
{
  string command_filename;
  ofstream command_unit;
  string data_filename;
  ofstream data_unit;
  int i;
  int j;
  int nz_num;
  string jpg_filename;
//
//  Create data file.
//
  data_filename = header + "_data.txt";
  data_unit.open ( data_filename.c_str ( ) );
  nz_num = 0;
  int  m = A.getNbLines();
  int  n = A.getNbColumns();
  
  for ( j = 0; j < n; ++j )
    for ( i = 0; i < m; ++i ){
      if ( A(i,j) != 0.0 ){
        data_unit << j << "  "
                  << i << "\n";
        nz_num +=  1;
      }
    }
  data_unit.close ( );
  cout << "\n";
  cout << "  Created sparsity data file '" << data_filename << "'\n";
//
//  Create command file.
//
  command_filename = header + "_commands.txt";
  command_unit.open ( command_filename.c_str ( ) );

  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "unset key\n";
  command_unit << "set term jpeg\n";

  jpg_filename = header + ".jpg";
  command_unit << "set output '" << jpg_filename << "'\n";
  command_unit << "set size ratio -1\n";
  command_unit << "set xlabel '<--- J --->'\n";
  command_unit << "set ylabel '<--- I --->'\n";
  command_unit << "set title '" << nz_num << " nonzeros for \"" 
               << header << "\"'\n";
  command_unit << "set timestamp\n";
  command_unit << "plot [y=0:" << n - 1 << "] [x="
               << m - 1 << ":0] '"
               << data_filename << "' with points pt 0\n";
 
  command_unit.close ( );
  cout << "  Created graphics command file '" << command_filename << "'\n";

  return;
}

#endif
