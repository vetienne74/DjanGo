//-------------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//   Extract subgrid from a grid 
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

template<class T> T* allocate_array(Myint64 mm) ;
template<class T> T** allocate_array(Myint64 mm, Myint64 nn) ;
template<class T> T*** allocate_array(Myint64 mm, Myint64 nn, Myint64 oo) ;

//-------------------------------------------------------------------------------------------------------

template<class T> T* allocate_array(Myint64 mm)
{
  T* A = new T[mm] ;
  assert(A != NULL) ;
  return(A) ;
} 

//-------------------------------------------------------------------------------------------------------

template<class T> T** allocate_array(Myint64 mm, Myint64 nn)
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

template<class T> T*** allocate_array(Myint64 mm, Myint64 nn, Myint64 oo)
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

  cout << "\n=== subgrid started... ===\n" ; 

  //-------------------------------
  // read input parameters
  //-------------------------------
  string file_in, file_out ;
  cin >> file_in ;
  cin >> file_out ;  

  Myint64 n1, n2, n3 ;
  Myint64 d1, d2, d3 ;
  Myint64 o1, o2, o3 ;
        
  cin >> n1 >> n2 >> n3 ;
  cout << "n1 " << n1 << "\n" ;
  cout << "n2 " << n2 << "\n" ;
  cout << "n3 " << n3 << "\n" ;
  cin >> d1 >> d2 >> d3 ;
  cout << "d1 " << d1 << "\n" ;
  cout << "d2 " << d2 << "\n" ;
  cout << "d3 " << d3 << "\n" ;    
  cin >> o1 >> o2 >> o3 ;
  cout << "o1 " << o1 << "\n" ;
  cout << "o2 " << o2 << "\n" ;
  cout << "o3 " << o3 << "\n" ;    

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
  cout << "file size " << file_size << " byte \n" ;
  in_file.read((char*) **val_array, file_size) ;
  in_file.close() ;

  //-------------------------------
  // decimation  
  //-------------------------------
  cout << "write output file " << file_out << "\n" ;
  ofstream out_file(file_out.c_str(), ios::binary | ios::trunc | ios::out) ;
  
  for (Myint64 i3 = o3; i3 < n3; i3+=d3)
    {
      for (Myint64 i2 = o2; i2 < n2; i2+=d2)
	{
	  for (Myint64 i1 = o1; i1 < n1; i1+=d1)
	    {
	      out_file.write((char*) &val_array[i3][i2][i1], sizeof(Myfloat32)) ;
	    }
	}
      cout << "copied " << i3 << "/" << n3 << endl ;
    }
 
  Myint64 n11, n22, n33 ;
  n11 = 1 + ((n1-o1-1)/d1) ;
  n22 = 1 + ((n2-o2-1)/d2) ;
  n33 = 1 + ((n3-o3-1)/d3) ;

  cout << "new n1 " << n11 << "\n" ;
  cout << "new n2 " << n22 << "\n" ;
  cout << "new n3 " << n33 << "\n" ;

  out_file.close();

  //-------------------------------
  // delete the array
  //-------------------------------
  //delete[]val_array ;

  cout << "=== subgrid ended Ok ! ===\n\n" ; 

  return 0;
}
