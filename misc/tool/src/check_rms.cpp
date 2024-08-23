//-------------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//   Compute NRMS between 2 files
//   files should be same length
//   normalisation of values by maximum is an option
//   check if NRMS is below the input threshold
//
// OUTPUT
//   display on the console
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

  cout << "\n=== check rms started... ===\n" ; 

  //-------------------------------
  // read input parameters
  //-------------------------------
  string file_in1, file_in2 ;
  cin >> file_in1 ;
  cin >> file_in2 ;

  Myfloat scal_val1 ;
  cin >> scal_val1 ;
  
  //-------------------------------
  // read file 1
  //-------------------------------
  
  cout << "read input file 1: " << file_in1 << "\n" ;
  ifstream in_file1 ;
  in_file1.open(file_in1.c_str(), ios::binary) ;
  assert(in_file1.is_open());

  // get file size 
  in_file1.seekg(0,ios_base::end);
  Myint64 in_size1 = in_file1.tellg();
  in_file1.seekg(0,ios_base::beg);
  cout << "file 1 size is :" << in_size1 << "\n" ;
  Myint64 nfloat = in_size1 / sizeof(Myfloat) ;
    if (sizeof(Myfloat) == 4)
    {
      cout << "assume #float (4 bytes):" << nfloat << "\n" ;
    }
  else
    {
      cout << "assume #double (8 bytes):" << nfloat << "\n" ;
    }

  cout << "apply factor: " << scal_val1 << "\n";
  
  // read file
  Myfloat* val_array1 = new Myfloat[nfloat] ;
  in_file1.read((char*) val_array1, in_size1) ;

  Myfloat min_val1, max_val1 ;
  Myfloat64 sum_val1 = 0 ;
  Myfloat n_inf, n_nan ;
 
  n_inf    = 0 ;
  n_nan    = 0 ;
  for (Myint64 ii=0; ii<nfloat; ii++)
    {
      val_array1[ii] *= scal_val1 ;
      if (ii == 0)
	{	  
	  min_val1 = val_array1[ii] ;
	  max_val1 = val_array1[ii] ;
	  if (isnan(val_array1[ii])) n_nan++ ;
	  if (isinf(val_array1[ii])) n_inf++ ;
	  sum_val1 += val_array1[ii] ;
	}
      if (val_array1[ii] > max_val1) max_val1 = val_array1[ii] ;
      if (val_array1[ii] < min_val1) min_val1 = val_array1[ii] ;
    }

  cout << "max val file 1: " << max_val1 << "\n" ;
  cout << "average file 1: " << sum_val1 / nfloat << "\n" ;
  cout << "min val: " << min_val1 << "\n" ;
  if (n_nan > 0) cout << "#nan: " << n_nan << "\n" ;
  if (n_inf > 0) cout << "#inf: " << n_inf << "\n" ;
  
  if ((in_size1 > 0) && (n_nan == 0) && (n_inf == 0) && !((max_val1 == 0.) && (min_val1 == 0.)) )
    {
      cout << "FILE OK\n" ;
    }
  else
    {
      cout << "FILE KO\n" ;
    }

  //-------------------------------
  // read file 2
  //-------------------------------
  
  cout << "\nread input file 2: " << file_in2 << "\n" ;
  ifstream in_file2 ;
  in_file2.open(file_in2.c_str(), ios::binary) ;
  assert(in_file2.is_open());

  // get file size
  in_file2.seekg(0,ios_base::end);
  Myint64 in_size2 = in_file2.tellg();
  in_file2.seekg(0,ios_base::beg);
  cout << "file 2 size is :" << in_size2 << "\n" ;

  if (in_size1 != in_size2)
   {
      cout << " error - files should have same size !\n" ;
      return 1;
    }
  
  // read file
  Myfloat* val_array2 = new Myfloat[nfloat] ;
  in_file2.read((char*) val_array2, in_size2) ;
  
  Myfloat min_val2, max_val2 ;
  Myfloat64 sum_val2 = 0 ;
 
  n_inf    = 0 ;
  n_nan    = 0 ;
  for (Myint64 ii=0; ii<nfloat; ii++)
    {
      if (ii == 0)
	{
	  min_val2 = val_array2[ii] ;
	  max_val2 = val_array2[ii] ;
	  if (isnan(val_array2[ii])) n_nan++ ;
	  if (isinf(val_array2[ii])) n_inf++ ;
	  sum_val2 += val_array2[ii] ;
	}
      if (val_array2[ii] > max_val2) max_val2 = val_array2[ii] ;
      if (val_array2[ii] < min_val2) min_val2 = val_array2[ii] ;
    }

  cout << "max val file 2: " << max_val2 << "\n" ;
  cout << "average file 2: " << sum_val2 / nfloat << "\n" ;
  cout << "min val: " << min_val2 << "\n" ;
  if (n_nan > 0) cout << "#nan: " << n_nan << "\n" ;
  if (n_inf > 0) cout << "#inf: " << n_inf << "\n" ;
  
  if ((in_size2 > 0) && (n_nan == 0) && (n_inf == 0) && !((max_val2 == 0.) && (min_val2 == 0.)) )
    {
      cout << "FILE OK\n" ;
    }
  else
    {
      cout << "FILE KO\n" ;
    }
  
  in_file1.close() ;
  in_file2.close() ;

  // output ratio
  cout << "\nratio max file1/file2: " << max_val1 / max_val2 << "\n" ;
  cout << "ratio min file1/file2: " << min_val1 / min_val2 << "\n" ;
  
  // normalise values
  Myint64 norm_rms ;
  cin >> norm_rms ;
  
  if (norm_rms == 0)
    {
      cout << "rms without normalisation\n" ;
    }
  else if (norm_rms == 1)
    {
     cout << "rms with normalisation\n" ;
     for (Myint64 ii=0; ii<nfloat; ii++)
       {
	 val_array1[ii] /= max_val1 ;
	 val_array2[ii] /= max_val2 ;
       }
    }
  else
    {
      cout << " error - rms normalisation should be 0 or 1 !\n" ;
      return 1;
    }

  // compute rms
 
  Myfloat tmp1 = 0. ;
  Myfloat tmp2 = 0. ;
  
  
  for (Myint64 ii=0; ii<nfloat; ii++)
    {
      tmp1 += pow((val_array1[ii] - val_array2[ii]), 2) ;
      tmp2 += pow(val_array2[ii], 2) ;
    }

  Myfloat rms = sqrt(tmp1 / tmp2) ;
  cout << "rms: " << rms << "\n" ;
  
  // check against threshold
  Myfloat rms_max ;
  cin >> rms_max ;

  cout << "Max RMS set to " << rms_max << "\n" ;
  if (rms <= rms_max)
    {
      cout << "*** RMS BELOW LIMIT ***\n" ;
    }
  else
    {
      cout << "*** RMS ABOVE LIMIT ***\n" ;
    }
  delete(val_array1) ;
  delete(val_array2) ;
  
  cout << "=== check rms ended Ok ! ===\n\n" ; 

  return 0;
}
