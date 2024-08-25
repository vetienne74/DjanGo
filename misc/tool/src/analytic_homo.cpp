//-------------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//   compute analytic solution in homogeneous 1D and 3D infinite medium
//   free surface option only for 1D
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

  cout << "\n=== analytic_homo started... ===\n" ; 

  //-------------------------------
  // read input paramters
  //-------------------------------
  string acqui_file_name ;
  cin >> acqui_file_name ;
  cout << " acqui. file\t" << acqui_file_name << "\n" ;

  Myfloat vp ;
  cin >> vp ;
  cout << " P-wave vel.\t" << vp << "\n" ;

  Myfloat src_freq ;
  cin >> src_freq ;
  cout << " Ricker freq\t" << src_freq << "\n" ;

  Myfloat dt ;
  cin >> dt ;
  cout << " ouput dt\t" << dt << "\n" ;

  Myint nt ;
  cin >> nt ;
  cout << " ouput nt\t" << nt << "\n" ;

  Myint free_surf ;
  cin >> free_surf ;
  if (free_surf == 1)
    {
      cout << " WITH free-surface\n" ;
    }
  else if (free_surf == 0)
    {
      cout << " WITHOUT free-surface\n" ;
    }
  else
    {
      cout << " error - invalid free-surface option\n" ;
      return -1 ;
    }
  
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
  else if (acqui_dim == 3)
    {
      cout << " acqui dim \t 3D\n" ;
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
      if (acqui_dim == 3) acqui_file >> xsrc[isrc] >> ysrc[isrc] ;
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
      if (acqui_dim == 3) acqui_file >> xrec[irec] >> yrec[irec] ;
    }

  acqui_file.close() ;

  //------------------------------------------------------
  // Ricker wavefield
  //------------------------------------------------------

  Myfloat *src_time_function = new Myfloat[nt+1] ;
  Myfloat amp       = DEFAULT_AMP_SRC ;

  // Ricker wavelet
  Myfloat da        = M_PI * src_freq ;
  Myfloat t0        = 1.5 * sqrt(6.) / da ;  
  Src_func src_func = RICKER ; 
    
  //------------------------------------------------------
  // Compute analytic solution
  //------------------------------------------------------

  // open output file
  ofstream out_file("pr.time.rec.analytic_homo.out.bin", ios::binary | ios::trunc | ios::out) ;
  
  for (Myint isrc = 0; isrc < nsrc; isrc++)
    {
      for (Myint irec = 0; irec < nrec; irec++)
	{

	  // direct wave
	  //---------------------------------------------------------------------------
	  {
	    Myfloat dist ;
	    Myfloat factor ;

	    if (acqui_dim == 1)
	      {
		dist = zrec[irec] - zsrc[isrc] ;
		factor = 0.5 / vp ;
	      }
	    else if (acqui_dim == 3)
	      {
		dist = sqrt(pow(zrec[irec] - zsrc[isrc], 2)
			    + pow(xrec[irec] - xsrc[isrc], 2)
			    + pow(yrec[irec] - ysrc[isrc], 2)) ;
		factor = 1. /( 4. * PI * dist *  vp * vp ) ; 
	      }	  
	
	    Myfloat delay = dist / vp ;
	    cout << " distance " << dist << " delay " << delay << "\n" ;

	    Myfloat max_amp   = 0. ;
	    Myfloat delay_src = 0. ;
	  
	    for (Myint it = 0; it <= nt; it ++)
	      {      

		Myfloat tt = it * dt ;
		Myfloat src_wavelet ;

		// Ricker wavelet
		Myfloat aa = M_PI * src_freq * (tt - t0 - abs(delay)) ;
		Myfloat a2 = pow(aa, 2.) ;     

		if (src_func == RICKER_PP)
		  {	  	  
		    src_wavelet = -(0.5 / (da*da)) * exp(-a2) ;   	  
		  }

		else if (src_func == RICKER_P)
		  {	  
		    src_wavelet = (aa/da) * exp(-a2) ;	 
		  }

		else if (src_func == RICKER)
		  {	  
		    src_wavelet = (1. - 2. * a2) * exp(-a2) ;	 
		  }

		else if (src_func == RICKER_D)
		  {	 
		    src_wavelet = -4. * aa * da * exp(-a2) -2. * aa * da * (1. -2. * a2) * exp(-a2)	; 
		  }

		else if (src_func == MONO_FREQ)
		  {
		    src_wavelet = sin(2. * M_PI * src_freq * tt) ;	 
		  }

		src_time_function[it] = amp * src_wavelet * factor  ;

		// retrieve time when source is max
		if (max_amp < abs(src_time_function[it])) 
		  {
		    max_amp = abs(src_time_function[it]) ;
		    delay_src = tt ;
		  }
	      }
	  } // direct wave

	  // reflection at free-surface
	  //---------------------------------------------------------------------------
	  if ((free_surf == 1) && (acqui_dim == 1)) {
	    Myfloat dist ;
	    Myfloat factor ;

	    if (acqui_dim == 1)
	      {
		dist = zrec[irec] + zsrc[isrc] ;
		factor = -1 * (0.5 / vp) ;
	      }
	    else if (acqui_dim == 3)
	      {
		dist = sqrt(pow(zrec[irec] - zsrc[isrc], 2)
			    + pow(xrec[irec] - xsrc[isrc], 2)
			    + pow(yrec[irec] - ysrc[isrc], 2)) ;
		factor = 1. /( 4. * PI * dist *  vp * vp ) ; 
	      }	  
	
	    Myfloat delay = dist / vp ;
	    cout << " distance " << dist << " delay " << delay << "\n" ;

	    Myfloat max_amp   = 0. ;
	    Myfloat delay_src = 0. ;
	  
	    for (Myint it = 0; it <= nt; it ++)
	      {      

		Myfloat tt = it * dt ;
		Myfloat src_wavelet ;

		// Ricker wavelet
		Myfloat aa = M_PI * src_freq * (tt - t0 - abs(delay)) ;
		Myfloat a2 = pow(aa, 2.) ;     

		if (src_func == RICKER_PP)
		  {	  	  
		    src_wavelet = -(0.5 / (da*da)) * exp(-a2) ;   	  
		  }

		else if (src_func == RICKER_P)
		  {	  
		    src_wavelet = (aa/da) * exp(-a2) ;	 
		  }

		else if (src_func == RICKER)
		  {	  
		    src_wavelet = (1. - 2. * a2) * exp(-a2) ;	 
		  }

		else if (src_func == RICKER_D)
		  {	 
		    src_wavelet = -4. * aa * da * exp(-a2) -2. * aa * da * (1. -2. * a2) * exp(-a2)	; 
		  }

		else if (src_func == MONO_FREQ)
		  {
		    src_wavelet = sin(2. * M_PI * src_freq * tt) ;	 
		  }

		src_time_function[it] += amp * src_wavelet * factor  ;

		// retrieve time when source is max
		if (max_amp < abs(src_time_function[it])) 
		  {
		    max_amp = abs(src_time_function[it]) ;
		    delay_src = tt ;
		  }
	      }
	  } // reflection at free surface
	  
	  out_file.write((char*)src_time_function, nt * sizeof(Myfloat)) ;
	}
    }

  out_file.close();

  // delete arrays
  delete[]xsrc ;
  delete[]ysrc ;
  delete[]zsrc ;
  delete[]xrec ;
  delete[]yrec ;
  delete[]zrec ;
  delete[]src_time_function ;
  
  cout << "=== analytic_homo ended Ok ! ===\n\n" ; 

  return 0;
}
