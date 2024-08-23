//-------------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//   Transpose a binary file 
//
// OUTPUT
//   output file 
//
//-------------------------------------------------------------------------------------------------------

#include <cassert>
#include <fstream>
#include <iostream>
#include <cstdlib>

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
  for (Myint64 ii=0; ii<mm; ii++)
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
  for (Myint64 ii=0; ii<mm; ii++)
    {
      for (Myint64 jj=0; jj<nn; jj++)
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

  cout << "\n=== transpose_2d started... ===\n" ; 

  //-------------------------------
  // read input parameters
  //-------------------------------
  string file_in = argv[1] ;
  Myint64 n1     = atoi(argv[2]) ;

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
  
  cout << "input n1 " << n1 << "\n" ;
  if (nfloat%n1 != 0)
    {
      cout << "*** error ***, file size is not a multiple of n1\n" ;
      return(1) ;
    }
  Myint64 n2 = nfloat / n1 ;
  cout << "====> n2 " << n2 << "\n" ;
  
  //-------------------------------
  // allocate val array
  //--------------------------------
  Myfloat **val_array = allocate_array<Myfloat>(n2, n1) ;
  if (val_array == NULL)
    {
      cout << " *** ERROR *** CAN NOT ALLOCATE MEMORY ! \n" ;
      return -1;
    }

  //-------------------------------
  // read file
  //-------------------------------
  cout << "read input file " << file_in << "\n" ;
  in_file.read((char*) &(val_array[0][0]), in_size) ;
  in_file.close() ;

  //-------------------------------
  // transpose
  //-------------------------------
  Myfloat **val_array2 = allocate_array<Myfloat>(n1, n2) ;
  if (val_array2 == NULL)
    {
      cout << " *** ERROR *** CAN NOT ALLOCATE MEMORY ! \n" ;
      return -1;
    }

  for (Myint64 i2=0; i2<n2; i2++)
    {
      for (Myint64 i1=0; i1<n1; i1++)
	{
	  val_array2[i1][i2] = val_array[i2][i1] ;
	}
    }
  
  //-------------------------------
  // write output file
  //-------------------------------
  //cout << "write output file " << file_out << "\n" ;
  ofstream out_file("transpose_2d.out", ios::binary | ios::trunc | ios::out) ;
  out_file.write((char*) &val_array2[0][0], nfloat * sizeof(Myfloat)) ;
  out_file.close();

  //-------------------------------
  // delete the array
  //-------------------------------
  delete[]val_array ;
  delete[]val_array2 ;
  
  cout << "=== transpose ended Ok ! ===\n\n" ; 

  return 0;
}
