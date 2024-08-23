//-------------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//   compute local matrices used in FEM for 1D element (arbitrarily number of nodes)
//   - mass matrix
//   - stiffness matrix
//   - inverse of mass matrix
//
//   numerical integration is done with Gauss-Legendre quadrature or Gauss-Lobatto-Legendre 
//
//   GL points and matrix inverse are obtained with Numerical Recipes functions 
//
//-------------------------------------------------------------------------------------------------------

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>

#include "constant.h"
#include "type_def.h"

// numerical recipes
#include "nr3.h"
#include "gaussj.h"
#include "gamma.h"
#include "gauss_wgts.h"

#include "gll.h"

using namespace std ;
using namespace django ;

//-----------------------------------------------------------------------------------------
// compute mass matrix with quadrature rule such as Gauss-Legendre
//
// INPUT
// xgl (VecDoub) = array of GL coordinates
// wgl (VecDoub) = array of GL weights
// xnode (VecDoub) = array of node coordinates in element
//
// OUTPUT
// Mat_M (MatDoub) mass matrix
//-----------------------------------------------------------------------------------------
void compute_mass_matrix(VecDoub &xgl, VecDoub& wgl, VecDoub& xnode, MatDoub &Mat_M)
{   

  // # nodes in element
  Myint nnode = xnode.size() ;

  // # GL points
  Myint ngl   = xgl.size() ;

  // initialize mass matrix
  for (Myint ii=0; ii< nnode; ii++)
    {
      for (Myint jj=0; jj< nnode; jj++)
	{
	  Mat_M[ii][jj] = 0.0 ;
	}
    }
  
  // loop on element nodes 
  for (Myint in=0; in<nnode; in++)
    {
      for (Myint jn=0; jn<nnode; jn++)
	{
	  // loop on quadrature points
	  for (Myint ip=0; ip<ngl; ip++)
	    {
	      // compute denominator of Lagrange polynomials
	      Myfloat64 den_phi_in = 1.0 ;
	      for (Myint kn=0; kn<nnode; kn++)
		{
		  if (kn == in) continue ;			     
		  den_phi_in = den_phi_in * (xnode[in] - xnode[kn]) ;
		}  
	      Myfloat64 den_phi_jn = 1.0 ;
	      for (Myint kn=0; kn<nnode; kn++)
		{
		  if (kn == jn) continue ;			     
		  den_phi_jn = den_phi_jn * (xnode[jn] - xnode[kn]) ;
		}
            
	      // compute phi(in) at xgl(ip)
	      Myfloat64 phi_in = 1.0 ;
	      for (Myint kn=0; kn<nnode; kn++)
		{
		  if (kn == in) continue ;			     
		  phi_in = phi_in * (xgl[ip] - xnode[kn]) ;
		}        
	      phi_in = phi_in / den_phi_in ;
            
	      // compute phi(jn) at xgl(ip)
	      Myfloat64 phi_jn = 1.0 ;
	      for (Myint kn=0; kn<nnode; kn++)
		{
		  if (kn == jn) continue ;			     
		  phi_jn = phi_jn * (xgl[ip] - xnode[kn]) ;
		}         
	      phi_jn = phi_jn / den_phi_jn ;

	      Mat_M[in][jn] = Mat_M[in][jn] + wgl[ip] * phi_in * phi_jn ;
	    }
	}
    } 
} 

//-----------------------------------------------------------------------------------------
// compute stiffness matrix with quadrature rule such as Gauss-Legendre
//
// INPUT
// xgl (float*) = array of GL coordinates
// wgl (float*) = array of GL weights
// xnode (float*) = array of node coordinates in element
//
// OUTPUT
// Mat_K (MatDoub) stiffness matrix
//-----------------------------------------------------------------------------------------
void compute_stiffness_matrix(VecDoub &xgl, VecDoub& wgl, VecDoub& xnode, MatDoub &Mat_K)
{   

  // # nodes in element
  Myint nnode = xnode.size() ;

  // # GL points
  Myint ngl   = xgl.size() ;
  
  // initialize mass matrix
  Myfloat64 K_n[nnode][nnode] ;
  for (Myint ii=0; ii< nnode; ii++)
    {
      for (Myint jj=0; jj< nnode; jj++)
	{
	  Mat_K[ii][jj] = 0.0 ;
	}
    } 
  
  // loop on element nodes 
  for (Myint in=0; in<nnode; in++)
    {
      for (Myint jn=0; jn<nnode; jn++)
	{
	  // loop on quadrature points
	  for (Myint ip=0; ip<ngl; ip++)
	    {
	      // compute denominator of Lagrange polynomials
	      Myfloat64 den_phi_in = 1.0 ;
	      for (Myint kn=0; kn<nnode; kn++)
		{
		  if (kn == in) continue ;			     
		  den_phi_in = den_phi_in * (xnode[in] - xnode[kn]) ;
		}  
	      Myfloat64 den_phi_jn = 1.0 ;
	      for (Myint kn=0; kn<nnode; kn++)
		{
		  if (kn == jn) continue ;			     
		  den_phi_jn = den_phi_jn * (xnode[jn] - xnode[kn]) ;
		}
            
	      // compute phi(in) at xgl(ip)
	      Myfloat64 phi_in = 1.0 ;
	      Myfloat64 d_phi_in = 0.0 ;
	      for (Myint kn=0; kn<nnode; kn++)
		{
		  if (kn == in) continue ;
		  // derivative
		  d_phi_in = phi_in + d_phi_in * (xgl[ip] - xnode[kn]) ;

		  phi_in = phi_in * (xgl[ip] - xnode[kn]) ;
		}        
	      phi_in = phi_in / den_phi_in ;
	      d_phi_in = d_phi_in / den_phi_in ;
            
	      // compute phi(jn) at xgl(ip)
	      Myfloat64 phi_jn = 1.0 ;
	      for (Myint kn=0; kn<nnode; kn++)
		{
		  if (kn == jn) continue ;			     
		  phi_jn = phi_jn * (xgl[ip] - xnode[kn]) ;
		}         
	      phi_jn = phi_jn / den_phi_jn ;

	      Mat_K[in][jn] = Mat_K[in][jn] + wgl[ip] * d_phi_in * phi_jn ;
	    }
	}
    }
} 


//-----------------------------------------------------------------------------------------
// main program
//-----------------------------------------------------------------------------------------
int main(int argc, char* argv[])
{

  cout << "\n=== compute_element_matrices... ===\n" ; 

  // which quadrature ?
  cout << "\n Which quadrature? 1=GL, 2=GLL " ;
  Myint tmp ;
  cin >> tmp ;
  Node_integ node_integ ;
  if (tmp == NODE_INTEG_GL)
    {
      node_integ = NODE_INTEG_GL ;
    }
  else if (tmp == NODE_INTEG_GLL)
    {
      node_integ = NODE_INTEG_GLL ;
    }
  else
    {
      cout << "Error, invalid quadrature " << tmp << endl ;
      return -1 ;
    }      
        
  // read number of quadrature points
  Myint ng ;
  cout << "\n Number of quadrature points? " ;
  cin >> ng ;
    
  Doub x1 = -1.0 ;
  Doub x2 = +1.0 ;
  VecDoub xg(ng) ;
  VecDoub wg(ng) ;

  // compute quadrature points with nr3
  if (node_integ == NODE_INTEG_GL)
    {
      gauleg(x1, x2, xg, wg) ;
      cout << "GL obtained with nr3\n" ;    
    }

  if (node_integ == NODE_INTEG_GLL)
    {
      gll(xg, wg) ;
      cout << "GLL obtained with gll.h\n" ;
    }           

  for (Myint ii=0; ii<xg.size(); ii++)
    {    
      cout << "xg[" << ii << "] " << xg[ii] << " wg[" << ii << "] " << wg[ii] << "\n" ;
    }

   // which quadrature ?
  cout << "\n Node element? 1=EQD, 2=GLL " ;
  cin >> tmp ;
  Node_type node_type ;
  if (tmp == NODE_TYPE_EQD)
    {
      node_type = NODE_TYPE_EQD ;
    }
  else if (tmp == NODE_TYPE_GLL)
    {
      node_type = NODE_TYPE_GLL ;
    }
  else
    {
      cout << "Error, invalid node type " << tmp << endl ;
      return -1 ;
    }  
  
  // read number of element nodes
  Myint nnode ;
  cout << "\n Number of nodes in element? " ;
  cin >> nnode ;

  // compute coordinates of element nodes
  VecDoub xnode(nnode) ;
  if (node_type == NODE_TYPE_EQD)
    {      
      Doub hnode = (x2-x1) / (nnode-1) ;
      for (Myint ii=0; ii<nnode; ii++)
	{
	  xnode[ii] = x1 + ii*hnode ;      	  
	}
    }
  else if (node_type == NODE_TYPE_GLL)
    {      
      VecDoub dummy(nnode) ;
      gll(xnode, dummy) ;     
    }

  for (Myint ii=0; ii<nnode; ii++)
    {	   
      cout << "xnode[" << ii << "] " << xnode[ii] << "\n" ;
    }
  
  //return 0 ;

  // compute mass matrix
  MatDoub Mat_M(nnode, nnode) ;
  compute_mass_matrix(xg, wg, xnode, Mat_M) ;

  cout << "\n*** Matrix Mat_M *** \n" ;
  for (int ii=0; ii< Mat_M.nrows(); ii++)
    {
      for (int jj=0; jj< Mat_M.ncols(); jj++)
	{
	  cout << Mat_M[ii][jj] << " " ;
	}
      cout << endl ;
    }

  // write matrix
  ofstream out_file ;  
  out_file.open("M_inv.bin", ios::binary | ios::trunc) ;
  assert(out_file.is_open());
  for (int ii=0; ii< Mat_M.nrows(); ii++)
    {
      for (int jj=0; jj< Mat_M.ncols(); jj++)
	{
	  Myfloat tmp = Mat_M[ii][jj] ;
	  out_file.write((char*) &tmp, sizeof(tmp)) ;
	}
    }  
  out_file.close() ;
  
  // compute stiffness matrix
  MatDoub Mat_K(nnode, nnode) ;
  compute_stiffness_matrix(xg, wg, xnode, Mat_K) ;

  cout << "\n*** Matrix Mat_K *** \n" ;
  for (int ii=0; ii< Mat_K.nrows(); ii++)
    {
      for (int jj=0; jj< Mat_K.ncols(); jj++)
	{
	  cout << Mat_K[ii][jj] << " " ;
	}
      cout << endl ;
    }

  // compute invert of mass matrix
  gaussj(Mat_M) ;
  cout << "\n*** Invert Matrix Mat_M ***\n" ;
  for (int ii=0; ii< Mat_M.nrows(); ii++)
    {
      for (int jj=0; jj< Mat_M.ncols(); jj++)
	{
	  cout << Mat_M[ii][jj] << " " ;
	}
      cout << endl ;
    } 
  
  cout << "=== compute_element_matrices ended Ok ! ===\n\n" ; 

  return 0 ;
}
