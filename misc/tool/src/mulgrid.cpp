//-------------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//   multiply the content of the input binary file (assuming float32) by an input value
//
// OUTPUT
//   binary file
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

  cout << "\n=== mulgrid started... ===\n" ; 

  //-------------------------------
  // read input parameters
  //-------------------------------
  string file_in, file_out ;
  cin >> file_in ;
  cin >> file_out ; 
  Myfloat val ;
  cin >> val ;

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
  Myint64 nfloat = in_size / 4 ;
  cout << "assume #float:" << nfloat << "\n" ;
  
  // read file
  Myfloat* val_array = new Myfloat[nfloat] ;
  in_file.read((char*) val_array, in_size) ;

  Myfloat min_val, max_val ;
  Myfloat64 sum_val ;
  Myfloat n_inf, n_nan ;
 
  n_inf    = 0 ;
  n_nan    = 0 ;
  min_val = val_array[0] ;
  max_val = val_array[0] ;
  sum_val = 0 ;
  
  for (Myint64 ii=0; ii<nfloat; ii++)
    {
      min_val = min(val_array[ii], min_val) ;
      max_val = max(val_array[ii], max_val) ;
      if (isnan(val_array[ii])) n_nan++ ;
      if (isinf(val_array[ii])) n_inf++ ;
      sum_val += val_array[ii] ;	  
    }

  cout << "INITIAL values\n" ;
  cout << "max val: " << max_val << "\n" ;
  cout << "average: " << sum_val / nfloat << "\n" ;
  cout << "min val: " << min_val << "\n" ;
  if (n_nan > 0) cout << "#nan: " << n_nan << "\n" ;
  if (n_inf > 0) cout << "#inf: " << n_inf << "\n" ;
  
  if ((in_size > 0) && (n_nan == 0) && (n_inf == 0) && !((max_val == 0.) && (min_val == 0.)) )
    {
      cout << "FILE OK\n" ;
    }
  else
    {
      cout << "FILE KO\n" ;
    }
  
  in_file.close() ;

  // multiply
  for (Myint64 ii=0; ii<nfloat; ii++)
    {
      val_array[ii] *= val ;
    }

  n_inf    = 0 ;
  n_nan    = 0 ;
  min_val = val_array[0] ;
  max_val = val_array[0] ;
  sum_val = 0 ;

  for (Myint64 ii=0; ii<nfloat; ii++)
    {
      min_val = min(val_array[ii], min_val) ;
      max_val = max(val_array[ii], max_val) ;
      if (isnan(val_array[ii])) n_nan++ ;
      if (isinf(val_array[ii])) n_inf++ ;
      sum_val += val_array[ii] ;	  
    }

  cout << "NEW values\n" ;
  cout << "max val: " << max_val << "\n" ;
  cout << "average: " << sum_val / nfloat << "\n" ;
  cout << "min val: " << min_val << "\n" ;

  // write file
  cout << "write output file " << file_out << "\n" ;
  ofstream out_file(file_out.c_str(), ios::binary | ios::trunc | ios::out) ;
  out_file.write((char*) &val_array[0], nfloat*sizeof(Myfloat32)) ;
  out_file.close();

  delete(val_array) ;
  
  cout << "=== mulgrid ended Ok ! ===\n\n" ; 

  return 0;
}
