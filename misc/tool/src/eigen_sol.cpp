//-------------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//   compute analytic solution for eigen mode 1d, 2d and 3d unity medium
//   l=1m, vp=1m/s, rho=1kg/m3
//
// INPUT
//   acquisition file
//
// OUTPUT
//   traces
//
//-------------------------------------------------------------------------------------------------------

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>

#include "constant.h"
#include "type_def.h"

using namespace std ;
using namespace django ;

int main(int argc, char* argv[])
{

  cout << "\n=== eigen_sol started... ===\n" ; 

  // single or double precision
  #ifdef _DOUBLE_PRECISION_
	  cout << " Computations in DOUBLE PRECISION\n" ;
  #else
		cout << " Computations in SINGLE PRECISION\n" ;
  #endif

  //-------------------------------
  // read input paramters
  //-------------------------------
  string acqui_file_name ;
  cin >> acqui_file_name ;
  cout << " acqui. file\t" << acqui_file_name << "\n" ;

  Myfloat tmin ;
  cin >> tmin ;
  cout << " tmin (s)\t" << tmin << "\n" ;

  Myfloat tmax ;
  cin >> tmax ;
  cout << " tmax (s)\t" << tmax << "\n" ; 

  Myfloat dt ;
  cin >> dt ;
  cout << " output dt (s)\t" << dt << "\n" ; 

  Myint nt = ceil((tmax-tmin)/dt + 1) ;
  cout << " output nt\t" << nt << "\n" ;
  
  //-------------------------------
  // read acquisition
  //-------------------------------
  cout << "\n read acquisition file \n" ;
   
  // open acquisition file
  ifstream acqui_file(acqui_file_name.c_str()) ;
  assert(acqui_file.is_open()) ;

  // acqui type
  Myint acqui_type ;
  acqui_file >> acqui_type ;
  if (acqui_type == 1)
    {
      cout << " acqui type \t FIXED\n" ;
    }
  else
    {
      cout << " error - acquisition file type is not supported\n" ;
      return -1 ;
    }

  // acqui dim
  Myint acqui_dim ;
  acqui_file >> acqui_dim ;
  if (acqui_dim == 1)
    {
      cout << " acqui dim \t 1D\n" ;
    }
  else if (acqui_dim == 2)
    {
      cout << " acqui dim \t 2D\n" ;
    }
  else
    {
      cout << " error - acquisition is not supported with dim=" << acqui_dim << "\n" ;
      return -1 ;
    }

  // get nb of sources
  Myint nsrc ;
  acqui_file >> nsrc ;

  if (nsrc > 0)
    {
      cout << " nb sources \t " << nsrc << "\n" ;
    }
  else
    {
      cout << " error - invalid nb of sources in acquisition file\n" ;
    }
  
  Myfloat *zsrc = new Myfloat[nsrc] ;
  Myfloat *xsrc = new Myfloat[nsrc] ;
  Myfloat *ysrc = new Myfloat[nsrc] ;
  
  // sources positions
  Myint itmp ;

  for (int isrc=0; isrc < nsrc; isrc++)
    {
      acqui_file >> itmp >> zsrc[isrc]  ;
      if (acqui_dim == 2) acqui_file >> xsrc[isrc] ;
      else if (acqui_dim == 3) acqui_file >> xsrc[isrc] >> ysrc[isrc] ;
    }

  // get nb of rec
  Myint nrec ;
  acqui_file >> nrec ;

  if (nrec > 0)
    {
      cout << " nb receivers \t " << nrec << "\n" ;
    }
  else
    {
      cout << " error - invalid nb of receivers in acquisition file\n" ;
    }
  
  Myfloat *zrec = new Myfloat[nrec] ;
  Myfloat *xrec = new Myfloat[nrec] ;
  Myfloat *yrec = new Myfloat[nrec] ;

  // rec positions
  for (int irec=0; irec < nrec; irec++)
    {
      acqui_file >> itmp >> zrec[irec]  ;
      if (acqui_dim == 2) acqui_file >> xrec[irec] ;
      else if (acqui_dim == 3) acqui_file >> xrec[irec] >> yrec[irec] ;
    }

  acqui_file.close() ; 
    
  //------------------------------------------------------
  // Compute analytic solution
  //------------------------------------------------------

  // open output file
  ofstream out_file("pr.time.rec.eigen_sol.out.bin", ios::binary | ios::trunc | ios::out) ;

  if (acqui_dim == 1)
    {
      for (Myint it = 0; it < nt; it++)
	{
	  for (Myint irec = 0; irec < nrec; irec++)
	    {	 
	      Myfloat time_pr = tmin + it*dt ;
	      Myfloat pr = -sin(M_PI*zrec[irec]) * sin(M_PI*time_pr) ;
	      out_file.write((char*) &pr, sizeof(Myfloat)) ;
	    }
	}
    }
  else if (acqui_dim == 2)
    {
      for (Myint it = 0; it < nt; it++)
	{
	  for (Myint irec = 0; irec < nrec; irec++)
	    {	
	      Myfloat time_pr = tmin + it*dt ;
	      Myfloat pr = -sqrt(2.0)*sin(M_PI*zrec[irec]) *sin(M_PI*xrec[irec]) * sin(sqrt(2.0)*M_PI*time_pr) ;
	      out_file.write((char*) &pr, sizeof(Myfloat)) ;
	    }
	}
    }

  out_file.close(); 
  
  cout << "=== eigen_sol ended Ok ! ===\n\n" ; 

  return 0;
}
