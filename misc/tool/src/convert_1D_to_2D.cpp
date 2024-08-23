//-------------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//   Read 1D profile and extend it to build a 2D model
//
// OUTPUT
//   output file 
//
//-------------------------------------------------------------------------------------------------------

#include <cassert>
#include <fstream>
#include <iostream>

#include "type_def.h"

using namespace std ;
using namespace django ;

template<class T> T* allocate_array(int mm) ;
template<class T> T** allocate_array(int mm, int nn) ;
template<class T> T*** allocate_array(int mm, int nn, int oo) ;

//-------------------------------------------------------------------------------------------------------

template<class T> T* allocate_array(int mm)
{
  T* A = new T[mm] ;
  assert(A != NULL) ;
  return(A) ;
} 

//-------------------------------------------------------------------------------------------------------

template<class T> T** allocate_array(int mm, int nn)
{ 
  T** A = new T*[mm] ;
  assert(A != NULL) ;
  T*  B = new T[mm*nn] ;
  assert(B != NULL) ;
  for (Myint ii=0; ii<mm; ii++)
    {
      A[ii] = B + ii * nn ;
    }
  return(A) ;
} 

//-------------------------------------------------------------------------------------------------------

template<class T> T*** allocate_array(int mm, int nn, int oo)
{ 
  T*** A = new T**[mm] ;
  assert(A != NULL) ;
  T**  B = new T*[mm*nn] ;
  assert(B != NULL) ;
  T* C = new T[mm*nn*oo] ; 
  assert(C != NULL) ;
  for (Myint ii=0; ii<mm; ii++)
    {
      for (Myint jj=0; jj<nn; jj++)
	{
	  B[nn*ii+jj] = C + (nn*ii+jj)*oo ;
	}
      A[ii] = B + nn*ii ;
    }
  return(A) ;
}

//-------------------------------------------------------------------------------------------------------

int main(int argc, char* argv[])
{

  cout << "\n=== convert_1D_to_2D started... ===\n" ; 

  //-------------------------------
  // read input parameters
  //-------------------------------
  string file_in, file_out ;
  cin >> file_in ;
  cin >> file_out ;

  Myint n1, n2 ;
 
  cin >> n1 >> n2 ;
  cout << "n1 " << n1 << "\n" ;
  cout << "n2 " << n2 << "\n" ;

  //-------------------------------
  // allocate arrays
  //--------------------------------
  Myfloat *array_in = allocate_array<Myfloat>(n1) ;
  if (array_in == NULL)
    {
      cout << " *** ERROR *** CAN NOT ALLOCATE MEMORY ! \n" ;
      return 0;
    }

  //-------------------------------
  // read inout file
  //-------------------------------
  cout << "read input file " << file_in << "\n" ;
  ifstream in_file ;
  in_file.open(file_in.c_str(), ios::binary) ;
  assert(in_file.is_open()); 
  Myint64 file_size = n1 * sizeof(Myfloat32);
  cout << "file size " << file_size << " bytes \n" ;
  in_file.read((char*) array_in, file_size) ;
  in_file.close() ;

  //-------------------------------
  // write output file
  //-------------------------------
  cout << "write output file " << file_out << "\n" ;
 
  ofstream out_file(file_out.c_str(), ios::binary | ios::trunc | ios::out) ;
  
  for (Myint64 i2 = 0; i2 < n2; i2++)
    {
      out_file.write((char*) &array_in[0], file_size) ;
    }
  out_file.close();

  //-------------------------------
  // delete the array
  //-------------------------------
  delete[]array_in ;

  cout << "=== convert_1D_to_2D ended Ok ! ===\n\n" ; 

  return 0;
}
