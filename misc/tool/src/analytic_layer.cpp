//-------------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//   compute analytic solution in 2-layer 1d infinite medium
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

  cout << "\n=== analytic_layer started... ===\n" ; 

  //-------------------------------
  // read input paramters
  //-------------------------------
  string acqui_file_name ;
  cin >> acqui_file_name ;
  cout << " acqui. file\t" << acqui_file_name << "\n" ;

  Myfloat vp1 ;
  cin >> vp1 ;
  cout << " P-wave vel1.\t" << vp1 << "\n" ;

  Myfloat vp2 ;
  cin >> vp2 ;
  cout << " P-wave vel2.\t" << vp2 << "\n" ;

  Myfloat interface_depth ;
  cin >> interface_depth ;
  cout << " Interf. depth\t" << interface_depth << "\n" ;
  
  Myfloat src_freq ;
  cin >> src_freq ;
  cout << " Ricker freq\t" << src_freq << "\n" ;

  Myfloat dt ;
  cin >> dt ;
  cout << " ouput dt\t" << dt << "\n" ;

  Myint nt ;
  cin >> nt ;
  cout << " ouput nt\t" << nt << "\n" ;
  
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
      return -1 ;
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
      return -1 ;
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

  Myfloat *src_time_function = new Myfloat[nt] ;
  Myfloat amp       = DEFAULT_AMP_SRC ;

  // Ricker wavelet
  Myfloat da        = M_PI * src_freq ;
  Myfloat t0        = 1.5 * sqrt(6.) / da ;  
  Src_func src_func = RICKER ; 
    
  //------------------------------------------------------
  // Compute analytic solution
  //------------------------------------------------------
  Myfloat factor ;
  Myfloat refl_coef = (vp2-vp1)/(vp2+vp1) ;
  Myfloat trans_coef = 1 + refl_coef ;
  cout << " refl_coef " << refl_coef << endl ;
  cout << " trans_coef " << trans_coef << endl ;
	      
  // open output file
  ofstream out_file("pr.time.rec.analytic_layer.out", ios::binary | ios::trunc | ios::out) ;
  
  for (Myint isrc = 0; isrc < nsrc; isrc++)
    {

      // check source above interface
      if (zsrc[isrc] > interface_depth)
	{
	  cout << " error - source below interface\n" ;
	  return -1 ;
	}

      // loop on receivers
      for (Myint irec = 0; irec < nrec; irec++)
	{
	  Myfloat dist = zrec[irec] - zsrc[isrc] ;
	  Myfloat delay ;	 
	  cout << " irec " << irec << "\n" ;
	  cout << " distance " << dist << "\n" ;
	  
	  // direct wave
	  if (zrec[irec] < interface_depth)
	    {
	      cout << " direct and reflected waves\n" ;
	      delay = dist / vp1 ;
	      factor = 0.5 / vp1 ;
	    }
	  else
	    // transmitted wave
	    {
	      Myfloat t1 = (interface_depth - zsrc[isrc]) / vp1 ;
	      Myfloat t2 = (zrec[irec] - interface_depth) / vp2 ;
	      delay = t1 + t2 ;
	      cout << " transmited waves\n" ;
	      factor = 0.5 / vp1 * trans_coef ;
	    }	  	
	  cout << " delay " << delay << endl ;
	  
	  // direct and transmitted wave
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
	    }

	  // reflected wave
	  if (zrec[irec] < interface_depth)
	    {
	      // delay
	      delay = ((interface_depth-zsrc[isrc])+(interface_depth-zrec[irec])) / vp1 ;
	      factor = 0.5 / vp1 * refl_coef ;
	      cout << " delay " << delay << endl ;	      
	      
	      // reflected wave wave
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
		}
	    }
	  
	  out_file.write((char*)src_time_function, nt * sizeof(Myfloat)) ;
	}
    }

  out_file.close();
  delete(xsrc) ;
  delete(ysrc) ;
  delete(zsrc) ;
  delete(xrec) ;
  delete(yrec) ;
  delete(zrec) ;
  delete(src_time_function) ;
  
  cout << "=== analytic_layer ended Ok ! ===\n\n" ; 

  return 0;
}
