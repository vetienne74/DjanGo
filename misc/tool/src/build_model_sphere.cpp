//-------------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//   Build a model with homogeneous background and spherical anomalies
//
// OUTPUT
//   output file -> 'build_model_sphere.out.bin'
//
//-------------------------------------------------------------------------------------------------------

#include <fstream>
#include <iostream>

#include "type_def.h"

using namespace std ;
using namespace django ;

int main(int argc, char* argv[])
{

  Myint64 nval ;
  Myfloat32* val_array ; 

  cout << "\n=== build_model_sphere started... ===\n" ; 

  //-------------------------------
  // read input parameters
  //-------------------------------

  Myint64 dim ;
  cin >> dim ;

  if (dim == TWO)

    {

      // model size
      Myint64 nz, nx ;
      cin >> nz >> nx ;
      nval = nz * nx ;

      // velocity background
      Myfloat32 val ; 
      cin >> val ;

      // nb sphere
      Myint64 nsz, nsx ;
      cin >> nsz >> nsx ;

      // origin of the sphere
      Myint64 osz, osx ;
      cin >> osz >> osx ;
      
      // interval between sphere
      Myint64 dsz, dsx ;
      cin >> dsz >> dsx ;

      // sphere radius
      Myint64 rs ;
      cin >> rs ;

      // sphere relative velocity
      Myfloat32 vals ; 
      cin >> vals ;

      //-------------------------------
      // allocate val array
      //-------------------------------
      val_array = new float[nval] ;
  
      //-------------------------------
      // build model
      //-------------------------------

      // background
      for (Myint64 iz=1; iz <= nz; iz++)
	{
	  for (Myint64 ix=1; ix <= nx; ix++)
	    {
	      Myint64 ii = (ix - 1) * nz + iz - 1 ;
	      val_array[ii] = val ;
	    }
	}

      // sphere
      const Myfloat32 rs2 = rs*rs ;
      Myint64 color = 1 ;

      for (Myint64 isz=1; isz <= nsz; isz++)
	{
	  for (Myint64 isx=1; isx <= nsx; isx++)
	    {
	      Myint64 csz = osz + (isz - 1) * dsz ;
	      Myint64 csx = osx + (isx - 1) * dsx ;

	      // cout << csz << " " << csx << "\n" ;

	      if (color == 1)
		{
		  color = -1 ;
		}
	      else
		{
		  color = 1 ;
		}

	      for (Myint64 iz=csz-rs; iz <= csz+rs; iz++)
		{
		  for (Myint64 ix=csx-rs; ix <= csx+rs; ix++)
		    {

		      Myint64 ii = (ix - 1) * nz + iz - 1 ;

		      if ((ii >= 1) && (ii<= nval))
			{
			  Myfloat dist2 = (iz - csz)*(iz - csz) + (ix - csx)*(ix - csx) ;
			  if (dist2 <= rs2)
			    {
			      val_array[ii] = val + color * vals;
			    }
			}
		    }
		}
	    }
	}
    }

  //-------------------------------
  // open output file
  //-------------------------------

  ofstream out_file("build_model_sphere.out.bin", ios::binary | ios::trunc | ios::out) ;

  //-------------------------------
  // write output file
  //-------------------------------

  out_file.write((char*)val_array, nval * sizeof(Myfloat32)) ;
 
  cout << " nb values in output file " << nval << "\n" ;

  //-------------------------------
  // close output file
  //-------------------------------

  out_file.close();
 
  //-------------------------------
  // delete the array
  //-------------------------------
  delete[]val_array ;

  cout << "=== build_model_sphere ended Ok ! ===\n\n" ; 

  return 0;
}
