//-------------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//   Write a ascii file that contains source and receiver positions
//
// INPUT  
//
// OUTPUT
//   output file -> 'acquisition.config'
//
//-------------------------------------------------------------------------------------------------------

#include <fstream>
#include <iostream>

#include "type_def.h"

using namespace std ;
using namespace django ;

int main(int argc, char* argv[])
{

  cout << "\n=== build_acqui_file_surround started... ===\n" ; 

  //-------------------------------
  // open output file
  //-------------------------------

  ofstream out_file("acquisition.config", ios::trunc | ios::out) ;

  //-------------------------------
  // read input parameters
  //-------------------------------

  Myint dim ;
  cin >> dim ;

  //-------------------------------------------------------------------------------------------------------

  Myint nsrc, nrec ;

  if (dim == ONE)

    {

      Myint nsrc_z ;
      Myfloat zsrc_or ;
      Myfloat delta_zsrc ;
      Myint nrec_z ;
      Myfloat zrec_or ;
      Myfloat delta_zrec ;

      cout << " create 1D acquisition file\n" ;

      cin >> nsrc_z ;
      cin >> zsrc_or ;
      cin >> delta_zsrc ;
      cin >> nrec_z ;
      cin >> zrec_or ;
      cin >> delta_zrec ;

      nsrc = nsrc_z ;
      nrec = nrec_z ;

      //-------------------------------
      // write output file
      //-------------------------------

      out_file << SEISMIC_FIXED << "\n" ;
      out_file << ONE << "\n" ;

      out_file << nsrc << "\n" ;
      Myint isrc = 0 ;
      for (Myint iz=0; iz < nsrc_z; iz++)
	{
	  out_file << ++isrc << " " << zsrc_or + iz * delta_zsrc << "\n" ;
	}

      out_file << nrec << "\n" ;
      Myint irec = 0 ;
      for (Myint iz=0; iz < nrec_z; iz++)
	{
	  out_file << ++irec << " " << zrec_or + iz * delta_zrec << "\n" ;
	}

    }

  //-------------------------------------------------------------------------------------------------------

  else if (dim == TWO)

    {

      Myint nsrc_z, nsrc_x ;
      Myfloat zsrc_or, xsrc_or ;
      Myfloat delta_zsrc, delta_xsrc ;
      Myint nrec_z, nrec_x ;
      Myfloat zrec_or, xrec_or ;
      Myfloat delta_zrec, delta_xrec ;

      cout << " create 2D acquisition file\n" ;

      cin >> nsrc_z >> nsrc_x ;
      cin >> zsrc_or >> xsrc_or ;
      cin >> delta_zsrc >>  delta_xsrc ;
      cin >> nrec_z >> nrec_x ;
      cin >> zrec_or >> xrec_or ;
      cin >> delta_zrec >> delta_xrec ;            

      //-------------------------------
      // write output file
      //-------------------------------

      out_file << SEISMIC_FIXED << "\n" ;
      out_file << TWO << "\n" ;

      // *** shots ***
      Myint n1, n2 ;
      if (nsrc_x == 1)
	{
	  n1 = 0 ;
	}
      else
	{
	  n1 = nsrc_x-2 ;
	}
      
      if (nsrc_z == 1)
	{
	  n2 = 0 ;
	}
      else
	{
	  n2 = nsrc_z-2 ;
	}
      
      nsrc = (nsrc_x * nsrc_z) - (n1 * n2) ;
      out_file << nsrc << "\n" ;
      Myint isrc = 0 ;      

      for (Myint ix=0; ix < nsrc_x; ix++)
	{
	  for (Myint iz=0; iz < nsrc_z; iz++)
	    {
	      if ((ix == 0) || (ix == nsrc_x-1) || (iz == 0) || (iz == nsrc_z-1))
		{		  
		  out_file << ++isrc << " " << zsrc_or + iz * delta_zsrc << " " << xsrc_or + ix * delta_xsrc << "\n" ;
		}	  	  
	    }
	}          

      if (isrc != nsrc)
	{
	  cout << "*** ERROR ON NO. SHOTS *** \n" ;
	}

      // *** receivers ***

      nrec = 2 * (nrec_z + nrec_x - 2) ;
      out_file << nrec << "\n" ;
      Myint irec = 0 ;      

      // up  
      for (Myint ix=0; ix < nrec_x; ix++)
	{
	  out_file << ++irec << " " << zrec_or + 0.           * delta_zrec << " " << xrec_or + ix * delta_xrec << "\n" ;	  
	}

      // down 
      for (Myint ix=0; ix < nrec_x; ix++)
	{	 
	  out_file << ++irec << " " << zrec_or + (nrec_z - 1) * delta_zrec << " " << xrec_or + ix * delta_xrec << "\n" ;
	}

      // left 
      for (Myint iz=1; iz < nrec_z-1; iz++)
	{
	  out_file << ++irec << " " << zrec_or + iz * delta_zrec << " " << xrec_or + 0.           * delta_xrec << "\n" ;	  
	}

      // right 
      for (Myint iz=1; iz < nrec_z-1; iz++)
	{	 
	  out_file << ++irec << " " << zrec_or + iz * delta_zrec << " " << xrec_or + (nrec_x - 1) * delta_xrec << "\n" ;
	}

      if (irec != nrec)
	{
	  cout << "*** ERROR ON NO. RECEIVERS *** \n" ;
	}
	
    }

  //-------------------------------------------------------------------------------------------------------

  else if (dim == THREE)

    {

      Myint nsrc_z, nsrc_x, nsrc_y ;
      Myfloat zsrc_or, xsrc_or, ysrc_or ;
      Myfloat delta_zsrc, delta_xsrc, delta_ysrc ;
      Myint nrec_z, nrec_x, nrec_y ;
      Myfloat zrec_or, xrec_or, yrec_or ;
      Myfloat delta_zrec, delta_xrec, delta_yrec ;

      cout << " create 3D acquisition file\n" ;

      cin >> nsrc_z >> nsrc_x >> nsrc_y ;
      cin >> zsrc_or >> xsrc_or >> ysrc_or ;
      cin >> delta_zsrc >>  delta_xsrc >> delta_ysrc ;

      cin >> nrec_z >> nrec_x >> nrec_y ;
      cin >> zrec_or >> xrec_or >> yrec_or ;
      cin >> delta_zrec >> delta_xrec >> delta_yrec ; 

      nsrc = nsrc_z * nsrc_x * nsrc_y ;
      nrec = nrec_z * nrec_x * nrec_y ;

      //-------------------------------
      // write output file
      //-------------------------------

      out_file << SEISMIC_FIXED << "\n" ;
      out_file << THREE << "\n" ;

      out_file << nsrc << "\n" ;
      Myint isrc = 0 ;
      for (Myint iz=0; iz < nsrc_z; iz++)
	{
	  for (Myint iy=0; iy < nsrc_y; iy++)
	    {
	      for (Myint ix=0; ix < nsrc_x; ix++)
		{
		  out_file << ++isrc << " " << zsrc_or + iz * delta_zsrc << " " << xsrc_or + ix * delta_xsrc << " " << ysrc_or + iy * delta_ysrc << "\n" ;
		}
	    }
	}

      out_file << nrec << "\n" ;
      Myint irec = 0 ;
      for (Myint iz=0; iz < nrec_z; iz++)
	{
	  for (Myint iy=0; iy < nrec_y; iy++)
	    {
	      for (Myint ix=0; ix < nrec_x; ix++)
		{
		  out_file << ++irec << " " << zrec_or + iz * delta_zrec << " " << xrec_or + ix * delta_xrec << " " << yrec_or + iy * delta_yrec << "\n" ;
		}
	    }

	}
    }

  cout << " nb shots:     " << nsrc << "\n" ;
  cout << " nb receivers: " << nrec << "\n" ;

  //-------------------------------
  // close output file
  //-------------------------------

  out_file.close();

  cout << "=== build_acqui_file_surround ended Ok ! ===\n\n" ; 

  return 0;
}
