//-------------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//   build source wavelet
//
// OUTPUT
//   one trace
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

  cout << "\n=== wavelet started... ===\n" ; 

  //-------------------------------
  // read input paramters
  //-------------------------------
  
  Myint src_tmp ;
  cin >> src_tmp ;

  Src_func src_func = (Src_func) (src_tmp) ;
  cout << " wavelet \t" ;
  if (src_func == RICKER_PP)
    {	  	  
      cout << "Gaussian \n" ;  	  
    }

  else if (src_func == RICKER_P)
    {	  
      cout << "derivative of Gaussian \n" ;  	 
    }

  else if (src_func == RICKER)
    {	  
      cout << "Ricker \n" ;  	
    }

  else if (src_func == RICKER_D)
    {	 
      cout << "derivative of Ricker \n" ;  
    }

  else if (src_func == MONO_FREQ)
    {
      cout << "mono chromatic \n" ;
    }
  
  Myfloat src_freq ;
  cin >> src_freq ;
  cout << " src. freq\t" << src_freq << "\n" ;

  Myfloat dt ;
  cin >> dt ;
  cout << " ouput dt\t" << dt << "\n" ;

  Myint nt ;
  cin >> nt ;
  cout << " ouput nt\t" << nt << "\n" ;    
  
  //------------------------------------------------------
  // build wavefield
  //------------------------------------------------------

  Myfloat *src_time_function = new Myfloat[nt+1] ;
  Myfloat amp       = DEFAULT_AMP_SRC ;

  // Ricker wavelet
  Myfloat da        = M_PI * src_freq ;
  Myfloat t0        = 1.5 * sqrt(6.) / da ;     
    
  //------------------------------------------------------
  // Compute analytic solution
  //------------------------------------------------------
  	 
  Myfloat max_amp   = 0. ;
  Myfloat delay_src = 0. ;
  Myfloat delay     = 0. ;
	  
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

      src_time_function[it] = amp * src_wavelet ;

      // retrieve time when source is max
      if (max_amp < abs(src_time_function[it])) 
	{
	  max_amp = abs(src_time_function[it]) ;
	  delay_src = tt ;
	}
    }

  cout << " Max time \t" << (nt-1)*dt << " s \n" ;
 
  // write binary file
  ofstream out_fileb("wavelet.out.bin", ios::binary | ios::trunc | ios::out) ;
  out_fileb.write((char*)src_time_function, nt * sizeof(Myfloat)) ;
  out_fileb.close();

  // write ascii file
  ofstream out_filea("wavelet.out.ascii", ios::trunc | ios::out) ;
  for (Myint it = 0; it <= nt; it ++)
    { 
      out_filea << it*dt << "\t" << src_time_function[it] << "\n" ;
    }
  out_filea.close();
  
  delete(src_time_function) ;
  
  cout << "=== wavelet ended Ok ! ===\n\n" ; 

  return 0;
}
