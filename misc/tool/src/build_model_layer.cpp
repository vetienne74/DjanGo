//-------------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//   Build a layered model 
//
// OUTPUT
//   output file -> 'build_model_layer.out.bin'
//
//-------------------------------------------------------------------------------------------------------

#include <fstream>
#include <iostream>

#include "type_def.h"

using namespace std ;
using namespace django ;

int main(int argc, char* argv[])
{

  cout << "\n=== build_model_layer started... ===\n" ; 

  //-------------------------------
  // read input parameters
  //-------------------------------

  Myint64 nz, nx, ny ;
  Myint64 n_val ;
  Myint64 dim ;
  cin >> dim ;

  Myint64 nlayer ;
  cin >> nlayer ;

  Myint64* zlayer = new Myint64[nlayer] ;
  Myfloat32* vlayer = new Myfloat32[nlayer] ;

  nz = 0 ;
  nx = 0 ;
  ny = 0 ;

  if (dim == ONE)

    {
      cout << " build 1D model with no. layer " << nlayer << "\n" ;
      for (Myint64 ilayer = 0; ilayer < nlayer; ilayer++)
	{
	  cin >> zlayer[ilayer] ;
	  cin >> vlayer[ilayer] ;
	  nx = 1 ;
	  ny = 1 ;
	  nz += zlayer[ilayer] ;
	  cout << " layer " << ilayer + 1 << " with value " <<  vlayer[ilayer] << "\n";
	}
      cout << " grid size nz " << nz << "\n" ;
    }

  else if (dim == TWO)

    {
      cout << " build 2D model with no. layer " << nlayer << "\n" ;
      for (Myint64 ilayer = 0; ilayer < nlayer; ilayer++)
	{
	  cin >> zlayer[ilayer] >> nx ;
	  cin >> vlayer[ilayer] ;
	  ny = 1 ;
	  nz += zlayer[ilayer] ;
	  cout << " layer " << ilayer + 1 << " with value " <<  vlayer[ilayer] << "\n";
	}
      cout << " grid size nz " << nz << "\n" ;
      cout << " grid size nx " << nx << "\n" ;
    }

  else if (dim == THREE)

    {
      cout << " build 3D model with no. layer " << nlayer << "\n" ;
      for (Myint64 ilayer = 0; ilayer < nlayer; ilayer++)
	{
	  cin >> zlayer[ilayer] >> nx >> ny ;
	  cin >> vlayer[ilayer] ;
	  nz += zlayer[ilayer] ;
	  cout << " layer " << ilayer + 1 << " with value " <<  vlayer[ilayer] << "\n";
	}
      cout << " grid size nz " << nz << "\n" ;
      cout << " grid size nx " << nx << "\n" ;
      cout << " grid size ny " << ny << "\n" ;
    }

  n_val = nz * nx * ny ;

  //-------------------------------
  // allocate val array
  //-------------------------------
  Myfloat32* val_array = new Myfloat32[n_val] ;
  
  //-------------------------------
  // fill the array
  //-------------------------------
  Myint64 igrid = 0 ;

  for (Myint64 ix = 0; ix < nx; ix++)
    {
      for (Myint64 iy = 0; iy < ny; iy++)
	{
	  Myint64 izmin = 0 ;
	  for (Myint64 ilayer = 0; ilayer < nlayer; ilayer++)
	    {
	      Myint64 izmax = izmin + zlayer[ilayer] ; 
	      for (Myint64 iz = izmin; iz < izmax; iz++)
		{
		  val_array[igrid] = vlayer[ilayer] ;
		  igrid++ ;
		}
	      izmin = izmax ;
	    }	    
	}
    }

  //-------------------------------
  // open output file
  //-------------------------------

  ofstream out_file("build_model_layer.out.bin", ios::binary | ios::trunc | ios::out) ;

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

  cout << "=== build_model_layer ended Ok ! ===\n\n" ; 

  return 0;
}
