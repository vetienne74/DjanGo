//-------------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//   Check file:
//   length should be > 0
//   should not contain NaN or Inf
//
// OUTPUT
//   display on the console
//
//-------------------------------------------------------------------------------------------------------

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>

#include "type_def.h"

using namespace std ;
using namespace django ;

int main(int argc, char* argv[])
{

  cout << "\n=== check file started... ===\n" ; 

  //-------------------------------
  // read input parameters
  //-------------------------------
  string file_in = argv[1] ;
  
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

  Myint64 nfloat = in_size / sizeof(Myfloat) ;
  if (sizeof(Myfloat) == 4)
    {
      cout << "assume #float (4 bytes):" << nfloat << "\n" ;
    }
  else
    {
      cout << "assume #double (8 bytes):" << nfloat << "\n" ;
    }
  
  // read file
  Myfloat* val_array = new Myfloat[nfloat] ;
  in_file.read((char*) val_array, in_size) ;

  Myfloat min_val, max_val ;
  Myfloat64 sum_val = 0 ;
  Myfloat n_inf, n_nan ;
 
  n_inf    = 0 ;
  n_nan    = 0 ;

  min_val = val_array[0] ;
  max_val = val_array[0] ;
  
  for (Myint64 ii=0; ii<nfloat; ii++)
    {
      min_val = min(val_array[ii], min_val) ;
      max_val = max(val_array[ii], max_val) ;
      if (isnan(val_array[ii])) n_nan++ ;
      if (isinf(val_array[ii])) n_inf++ ;
      sum_val += val_array[ii] ;	  
    }

  cout << "max val: " << max_val << "\n" ;
  cout << "average: " << sum_val / nfloat << "\n" ;
  cout << "min val: " << min_val << "\n" ;
  if (n_nan > 0) cout << "#nan: " << n_nan << "\n" ;
  if (n_inf > 0) cout << "#inf: " << n_inf << "\n" ;
  
  if ((in_size > 0) && (n_nan == 0) && (n_inf == 0))
    {
      cout << "FILE OK\n" ;
    }
  else
    {
      cout << "FILE KO\n" ;
    }
  
  in_file.close() ;

  delete(val_array) ;
  
  cout << "=== check file ended Ok ! ===\n\n" ; 

  return 0;
}
