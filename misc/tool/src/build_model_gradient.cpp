//-------------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//   Build a layered gradient model 
//
// OUTPUT
//   output file -> 'build_model_gradient.out.bin'
//
//-------------------------------------------------------------------------------------------------------

#include <fstream>
#include <iostream>

#include "type_def.h"

using namespace std ;
using namespace django ;

int main(int argc, char* argv[])
{

  cout << "\n=== build_model_gradient started... ===\n" ; 

  //-------------------------------
  // read input parameters
  //-------------------------------

  Myint64 nz, nx, ny ;
  Myint64 n_val ;
  Myint64 dim ;
  cin >> dim ;

  Myint64 ngradient ;
  cin >> ngradient ;

  Myint64* zgradient = new Myint64[ngradient] ;
  Myfloat32* vgradient_min = new Myfloat32[ngradient] ;
  Myfloat32* vgradient_max = new Myfloat32[ngradient] ;

  nz = 0 ;
  nx = 0 ;
  ny = 0 ;

  if (dim == ONE)

    {
      cout << " build 1D model with no. gradient " << ngradient << "\n" ;
      for (Myint64 igradient = 0; igradient < ngradient; igradient++)
	{
	  cin >> zgradient[igradient] ;
	  cin >> vgradient_min[igradient] >> vgradient_max[igradient] ;
	  nx = 1 ;
	  ny = 1 ;
	  nz += zgradient[igradient] ;
	  cout << " gradient " << igradient + 1 << " with min " <<  vgradient_min[igradient] << " and max " << vgradient_max[igradient] << "\n";
	}
      cout << " grid size nz " << nz << "\n" ;
    }

  else if (dim == TWO)

    {
      cout << " build 2D model with no. gradient " << ngradient << "\n" ;
      for (Myint64 igradient = 0; igradient < ngradient; igradient++)
	{
	  cin >> zgradient[igradient] >> nx ;
	  cin >> vgradient_min[igradient] >> vgradient_max[igradient] ;
	  ny = 1 ;
	  nz += zgradient[igradient] ;
	  cout << " gradient " << igradient + 1 << " with min " <<  vgradient_min[igradient] << " and max " << vgradient_max[igradient] << "\n";
	}
      cout << " grid size nz " << nz << "\n" ;
      cout << " grid size nx " << nx << "\n" ;
    }

  else if (dim == THREE)

    {
      cout << " build 3D model with no. gradient " << ngradient << "\n" ;
      for (Myint64 igradient = 0; igradient < ngradient; igradient++)
	{
	  cin >> zgradient[igradient] >> nx >> ny ;
	  cin >> vgradient_min[igradient] >> vgradient_max[igradient] ;
	  nz += zgradient[igradient] ;
	  cout << " gradient " << igradient + 1 << " with min " <<  vgradient_min[igradient] << " and max " << vgradient_max[igradient] << "\n";
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
	  for (Myint64 igradient = 0; igradient < ngradient; igradient++)
	    {
	      Myint64 izmax = izmin + zgradient[igradient] ; 
	      for (Myint64 iz = izmin; iz < izmax; iz++)
		{
		  val_array[igrid] = vgradient_min[igradient] + (vgradient_max[igradient] - vgradient_min[igradient]) * (iz-izmin) / zgradient[igradient] ;
		  igrid++ ;
		}
	      izmin = izmax ;
	    }	    
	}
    }

  //-------------------------------
  // open output file
  //-------------------------------

  ofstream out_file("build_model_gradient.out.bin", ios::binary | ios::trunc | ios::out) ;

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

  cout << "=== build_model_gradient ended Ok ! ===\n\n" ; 

  return 0;
}
