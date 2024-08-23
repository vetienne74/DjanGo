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

  cout << "\n=== transpose started... ===\n" ; 

  //-------------------------------
  // read input parameters
  //-------------------------------
  string file_in, file_out ;
  cin >> file_in ;
  cin >> file_out ;

  Myint dim ;
  cin >> dim ;

  Myint n1, n2, n3 ;
  if (dim == THREE)

    {
      cout << "transpose 3D model\n" ;
      cin >> n1 >> n2 >> n3 ;
      cout << "n1 " << n1 << "\n" ;
      cout << "n2 " << n2 << "\n" ;
      cout << "n3 " << n3 << "\n" ;
    }

  else
    {
      cout << " *** ERROR *** NOT 3D ! \n" ; 
      return 0;
    }

  //-------------------------------
  // allocate val array
  //--------------------------------
  Myfloat ***val_array = allocate_array<Myfloat>(n3, n2, n1) ;
  if (val_array == NULL)
    {
      cout << " *** ERROR *** CAN NOT ALLOCATE MEMORY ! \n" ;
      return 0;
    }

  //-------------------------------
  // read file
  //-------------------------------
  cout << "read input file " << file_in << "\n" ;
  ifstream in_file ;
  in_file.open(file_in.c_str(), ios::binary) ;
  assert(in_file.is_open()); 
  Myint64 file_size = n1 * n2 * n3 * sizeof(Myfloat32);
  cout << "file size " << file_size << " bytes \n" ;
  in_file.read((char*) **val_array, file_size) ;
  in_file.close() ;

  //-------------------------------
  // transpose
  // and write output file
  //-------------------------------
  cout << "write output file " << file_out << "\n" ;
  ofstream out_file(file_out.c_str(), ios::binary | ios::trunc | ios::out) ;

  int option = 2 ;

  if (option == 1)
    {
      cout << "n1 " << n1 << "\n" ;
      cout << "n2 " << n3 << "\n" ;
      cout << "n3 " << n2 << "\n" ;
     
      for (Myint64 i2 = 0; i2 < n2; i2++)
	{
	  for (Myint64 i3 = 0; i3 < n3; i3++)
	    {
	      out_file.write((char*) &val_array[i3][i2][0], n1 * sizeof(Myfloat32)) ;
	    }
	}
    }
  else if (option == 2)
    {
      cout << "n1 " << n2 << "\n" ;
      cout << "n2 " << n1 << "\n" ;
      cout << "n3 " << n3 << "\n" ;     
  
      for (Myint64 i3 = 0; i3 < n3; i3++)
	{
	  for (Myint64 i1 = 0; i1 < n1; i1++)
	    {
	      for (Myint64 i2 = 0; i2 < n2; i2++)
		{
		  out_file.write((char*) &val_array[i3][i2][i1], sizeof(Myfloat32)) ;
		}
	    }
	}
    }

  out_file.close();

  //-------------------------------
  // delete the array
  //-------------------------------
  delete[]val_array ;

  cout << "=== transpose ended Ok ! ===\n\n" ; 

  return 0;
}
