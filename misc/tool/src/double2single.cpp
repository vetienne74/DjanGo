//-------------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//   convert a file from double precision to single precision
//
// OUTPUT
//   new file in single precision
//
//-------------------------------------------------------------------------------------------------------

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>

#include "type_def.h"

using namespace std ;
using namespace django ;

int main(int argc, char* argv[])
{

  cout << "\n=== double2single started... ===\n" ; 

  //-------------------------------
  // read input parameters
  //-------------------------------
  string file_in  = argv[1] ;
  
  //-------------------------------
  // read file
  //-------------------------------
  cout << "read input file " << file_in << "\n" ;
  ifstream in_file ;
  in_file.open(file_in.c_str(), ios::binary) ;
  assert(in_file.is_open());

  // get file size
  in_file.seekg(0,ios_base::end);
  Myint64 in_size = in_file.tellg();
  in_file.seekg(0,ios_base::beg);
  cout << "file size is :" << in_size << "\n" ;

  Myint64 nfloat = in_size / 8 ;
  cout << "assume #float (8 bytes):" << nfloat << "\n" ;
  
  // read file
  Myfloat64* val_array = new Myfloat64[nfloat] ;
  in_file.read((char*) val_array, in_size) ;
  in_file.close() ;
  
  // convert to double precision
  Myfloat32* val_array2 = new Myfloat32[nfloat] ;
  for (Myint64 ii=0; ii<nfloat; ii++)
    {
      val_array2[ii] = val_array[ii] ;
    }

  // write output file
  ofstream out_file("double2single.out.bin", ios::binary | ios::trunc | ios::out) ;
  out_file.write((char*) &val_array2[0], nfloat * sizeof(Myfloat32)) ;
  out_file.close();

  delete(val_array) ;
  
  cout << "=== double2single ended Ok ! ===\n\n" ; 

  return 0;
}
