//-------------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//   Write a binary file and fill it with a constant value
//
// OUTPUT
//   output file -> 'build_model_homo.out.bin'
//
//-------------------------------------------------------------------------------------------------------

#include <fstream>
#include <iostream>

#include "type_def.h"

using namespace std ;
using namespace django ;

int main(int argc, char* argv[])
{

  cout << "\n=== build_model_homo started... ===\n" ; 

  //-------------------------------
  // read input parameters
  //-------------------------------

  Myint64 n_val ;
  Myint64 dim ;
  cin >> dim ;

  if (dim == ONE)

    {
      cout << " build 1D model\n" ;
      Myint64 nz ;
      cin >> nz ;
      n_val = nz ;

    }

  else if (dim == TWO)

    {
      cout << " build 2D model\n" ;
      Myint64 nz, nx ;
      cin >> nz >> nx ;
      n_val = nz * nx ;

    }

  else if (dim == THREE)

    {
      cout << " build 3D model\n" ;
      Myint64 nz, nx, ny ;
      cin >> nz >> nx >> ny ;
      n_val = nz * nx * ny ;

    }

  Myfloat val ; 
  cin >> val ;

  //-------------------------------
  // allocate val array
  //-------------------------------
  Myfloat32* val_array = new Myfloat32[n_val] ;
  
  //-------------------------------
  // fill the array
  //-------------------------------
  for (int ii=0; ii < n_val; ii++)
    {
      val_array[ii] = val ;
    }

  //-------------------------------
  // open output file
  //-------------------------------

  ofstream out_file("build_model_homo.out.bin", ios::binary | ios::trunc | ios::out) ;

  //-------------------------------
  // write output file
  //-------------------------------

  out_file.write((char*)val_array, n_val * sizeof(Myfloat32)) ;
 
  cout << " nb values in output file " << n_val << "\n" ;

  //-------------------------------
  // close output file
  //-------------------------------

  out_file.close();
 
  //-------------------------------
  // delete the array
  //-------------------------------
  delete[]val_array ;

  cout << "=== build_model_homo ended Ok ! ===\n\n" ; 

  return 0;
}
