//-------------------------------------------------------------------------------------------------------
//
// MODELLING WITH MGM IN 1D
//
// PARENT CLASS: Scheme
//   DERIVED CLASS: FEM
//     DERIVED CLASS: FEM_1D
//       DERIVED CLASS: FEM_1D_2nd
//
//-------------------------------------------------------------------------------------------------------

#include "FEM_1D_2nd.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>

#include "output_report.h"
#include "type_def.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------

  FEM_1D_2nd::FEM_1D_2nd(void)  : FEM_1D()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D_2nd::FEM_1D_2nd");

	eq_order = ORDER_2ND ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D_2nd::FEM_1D_2nd");
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_1D_2nd::initialize(Model* pModel, Myfloat fmax)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D_2nd::initialize");

	// call parent initialization
	Rtn_code rtn_code = FEM_1D::initialize(pModel, fmax) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D_2nd::initialize");
	return(RTN_CODE_OK) ;
} ;

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_1D_2nd::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D_2nd::finalize");

	// call parent finalize
	Rtn_code rtn_code = FEM_1D::finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// deallocate flux matrices
	for (Myint ineigh = 0; ineigh < 2; ineigh++)
	{
		for (Myint imat=0; imat <= MAX_POLY_ORDER; imat++)
		{
			if (pMat_Fz_i[ineigh][imat] != NULL) delete(pMat_Fz_i[ineigh][imat]) ;
			if (pMat_Fz_j[ineigh][imat] != NULL) delete(pMat_Fz_j[ineigh][imat]) ;
			if (pMat_Fz_c[ineigh][imat] != NULL) delete(pMat_Fz_c[ineigh][imat]) ;

			for (Myint imat2=0; imat2 <= MAX_POLY_ORDER; imat2++)
			{
				if (pMat_Gz_i[ineigh][imat][imat2] != NULL) delete(pMat_Gz_i[ineigh][imat][imat2]) ;
				if (pMat_Gz_j[ineigh][imat][imat2] != NULL) delete(pMat_Gz_j[ineigh][imat][imat2]) ;
				if (pMat_Gz_c[ineigh][imat][imat2] != NULL) delete(pMat_Gz_c[ineigh][imat][imat2]) ;
			}
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D_2nd::finalize");
	return(RTN_CODE_OK) ;
} ;

//-----------------------------------------------------------------------------------------
// compute flux matrix F with quadrature rule such as Gauss-Legendre
// in the reference element with xi in [-1, 1]
//
// (Fik)rj = Int_elem_face_ik (phi_ir x d(phi_ij)/dz ) dx
// with i = current element, k = neigbhour element
//      phi_ir = basis function in element i for node r
//      dx phi_ij = basis function in element i for node j 
//      r = line, j = column in matrix
//      Int_elem_face_ik = integral over surface shared between i and k
//
// INPUT
// xg (VecDoub) = array of quadrature coordinates
// wg (VecDoub) = array of quadrature weights
// xnode (VecDoub) = array of node coordinates in element
// neigh (Neigh_type) = face type (I_ZPREV, I_ZNEXT) for 1D
//
// OUTPUT
// Mat_F (MatDoub) flux matrix
//-----------------------------------------------------------------------------------------
Rtn_code FEM_1D_2nd::compute_flux_matrix_F(VecDoub &xg, VecDoub& wg, VecDoub& xnode, MatDoub &Mat_F, Neigh_type neigh)
{   
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D_2nd::compute_flux_matrix_F");

	// # nodes in element
	Myint nn = xnode.size() ;
	if (nn == 0)
	{
		print_error(" Error in FEM_1D_2nd::compute_flux_matrix_F, xnode.size() = 0") ;
		return(RTN_CODE_KO) ;
	}

	// check size of matrix
	if (nn != Mat_F.ncols())
	{
		print_error(" Error in FEM_1D_2nd::compute_flux_matrix_F, nn != Mat_F.ncols()") ;
		return(RTN_CODE_KO) ;
	}
	if (nn != Mat_F.nrows())
	{
		print_error(" Error in FEM_1D_2nd::compute_flux_matrix_F, nn != Mat_F.nrows()") ;
		return(RTN_CODE_KO) ;
	}

	// # GL points
	Myint ng = xg.size() ;
	if (ng == 0)
	{
		print_error(" Error in FEM_1D_2nd::compute_flux_matrix_F, xg.size() = 0") ;
		return(RTN_CODE_KO) ;
	}

	// initialize flux matrix
	for (Myint ii=0; ii< nn; ii++)
	{
		for (Myint jj=0; jj< nn; jj++)
		{
			Mat_F[ii][jj] = 0.0 ;
		}
	}

	// check nighbour type
	if ((neigh != I_ZPREV) && (neigh != I_ZNEXT))
	{
		print_error(" Error in FEM_1D_2nd::compute_flux_matrix_F, invalid neighbour type", neigh) ;
		return(RTN_CODE_KO) ;
	}

	// loop on element nodes i
	for (Myint in=0; in<nn; in++)
	{
		// loop on element nodes j
		for (Myint jn=0; jn<nn; jn++)
		{
			// loop on gradrature point

			// phi_in at xg
			Doub phi_in = 0.0 ;

			// phi_jn at xg
			Doub phi_jn = 0.0 ;

			for (Myint ig=0; ig<ng; ig++)
			{
				if ((neigh == I_ZPREV) && (ig == 0))
				{
					phi_in = lagran(xnode, in, xg[ig]) ;
					phi_jn = dlagran(xnode, jn, xg[ig]) ;
				}
				else if ((neigh == I_ZNEXT) && (ig == ng-1))
				{
					phi_in = lagran(xnode, in, xg[ig]) ;
					phi_jn = dlagran(xnode, jn, xg[ig]) ;
				}
				else
				{
					continue ;
				}

				// update mass matrix
				Mat_F[in][jn] = Mat_F[in][jn] + phi_in * phi_jn ;
			}
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D_2nd::compute_flux_matrix_F");
	return(RTN_CODE_OK) ;
}

//-----------------------------------------------------------------------------------------
// compute flux matrix G with quadrature rule such as Gauss-Legendre
// in the reference element with xi in [-1, 1]
//
// (Gik)rj = Int_elem_face_ik (phi_ir x d(phi_kj)/dz ) dxi
// with i = current element, k = neigbhour element
//      phi_ir = basis function in element i for node r
//      phi_kj = basis function in element k for node j 
//      r = line, j = column in matrix
//      Int_elem_face_ik = integral over surface shared between i and k
//
// INPUT
// xg (VecDoub) = array of quadrature coordinates
// wg (VecDoub) = array of quadrature weights
// xnode_i (VecDoub) = array of node coordinates in element i
// xnode_k (VecDoub) = array of node coordinates in element k
// neigh (Neigh_type) = face type (I_ZPREV, I_ZNEXT) for 1D
//
// OUTPUT
// Mat_G (MatDoub) flux matrix
//-----------------------------------------------------------------------------------------
Rtn_code FEM_1D_2nd::compute_flux_matrix_G(VecDoub &xg, VecDoub& wg,
		VecDoub& xnode_i, VecDoub& xnode_k, MatDoub &Mat_G, Neigh_type neigh)
{   
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D_2nd::compute_flux_matrix_G");

	// # nodes in element
	Myint ni = xnode_i.size() ;
	if (ni == 0)
	{
		print_error(" Error in FEM_1D_2nd::compute_flux_matrix_G, xnode_i.size() = 0") ;
		return(RTN_CODE_KO) ;
	}

	// check size of matrix
	if (ni != Mat_G.nrows())
	{
		print_error(" Error in FEM_1D_2nd::compute_flux_matrix_G, ni != Mat_F.nrows()") ;
		return(RTN_CODE_KO) ;
	}

	// # nodes in element
	Myint nk = xnode_k.size() ;
	if (nk == 0)
	{
		print_error(" Error in FEM_1D_2nd::compute_flux_matrix_G, xnode_k.size() = 0") ;
		return(RTN_CODE_KO) ;
	}

	// check size of matrix
	if (nk != Mat_G.ncols())
	{
		print_error(" Error in FEM_1D_2nd::compute_flux_matrix_G, nk != Mat_F.ncols()") ;
		return(RTN_CODE_KO) ;
	}

	// # GL points in element i
	Myint ngi = xg.size() ;
	if (ngi == 0)
	{
		print_error(" Error in FEM_1D_2nd::compute_flux_matrix_G, xg.size() = 0") ;
		return(RTN_CODE_KO) ;
	}

	// initialize flux matrix
	for (Myint ii = 0; ii < Mat_G.nrows(); ii++)
	{
		for (Myint jj = 0; jj < Mat_G.ncols(); jj++)
		{
			Mat_G[ii][jj] = 0.0 ;
		}
	}

	// check nighbour type
	if ((neigh != I_ZPREV) && (neigh != I_ZNEXT))
	{
		print_error(" Error in FEM_1D_2nd::compute_flux_matrix_G, invalid neighbour type", neigh) ;
		return(RTN_CODE_KO) ;
	}

	// loop on element i nodes r
	for (Myint ir=0; ir<ni; ir++)
	{
		// loop on element k nodes j
		for (Myint kj=0; kj<nk; kj++)
		{
			// loop on gradrature point
			for (Myint ig=0; ig<ngi; ig++)
			{
				// phi_in at xg[ig]
				Doub phi_in = 0.0 ;

				// phi_jn at xg[ig]
				Doub phi_jn = 0.0 ;

				if ((neigh == I_ZPREV) && (ig == 0))
				{
					phi_in = lagran(xnode_i, ir, xg[0]) ;
					phi_jn = dlagran(xnode_k, kj, xg[ngi-1]) ;
				}
				else if ((neigh == I_ZNEXT) && (ig == ngi-1))
				{
					phi_in = lagran(xnode_i, ir, xg[ngi-1]) ;
					phi_jn = dlagran(xnode_k, kj, xg[0]) ;
				}
				else
				{
					continue ;
				}

				// update mass matrix
				Mat_G[ir][kj] = Mat_G[ir][kj] + phi_in * phi_jn ;
			}
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D_2nd::compute_flux_matrix_G");
	return(RTN_CODE_OK) ;
}

//-----------------------------------------------------------------------------------------
Rtn_code FEM_1D_2nd::compute_flux_matrices(Myint imat, ofstream *out_file, Myint ni, VecDoub &xg, VecDoub& wg, VecDoub& xnode)
{   
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D_2nd::compute_flux_matrices");

	//----------------------
	// compute flux matrix F
	//----------------------
	for (Myint ineigh = 0; ineigh < 2; ineigh++)
	{
		MatDoub Mat_Fz(ni, ni) ;
		Rtn_code rtn_code = compute_flux_matrix_F(xg, wg, xnode, Mat_Fz, (Neigh_type) ineigh) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		*out_file << "\n*** Matrix F ineigh " << ineigh << " *** \n" ;
		for (int ii=0; ii< Mat_Fz.nrows(); ii++)
		{
			for (int jj=0; jj< Mat_Fz.ncols(); jj++)
			{
				*out_file << Mat_Fz[ii][jj] << " " ;
			}
			*out_file << endl ;
		}

		// convert full matrix to sparse matrix
		full2sparse(&pMat_Fz_i[ineigh][imat], &pMat_Fz_j[ineigh][imat], &pMat_Fz_c[ineigh][imat], Mat_Fz) ;

		// print sparse matrix
		Myint*   Mat_Fz_i = pMat_Fz_i[ineigh][imat]->pArray ;
		Myint*   Mat_Fz_j = pMat_Fz_j[ineigh][imat]->pArray ;
		Myfloat* Mat_Fz_c = pMat_Fz_c[ineigh][imat]->pArray ;
		*out_file << "\ni j c / sparse format\n" ;
		for (Myint ii=0; ii< pMat_Fz_i[ineigh][imat]->nz; ii++)
		{
			*out_file << Mat_Fz_i[ii] << " " << Mat_Fz_j[ii] << " " << Mat_Fz_c[ii] << endl ;
		}
	}

	//----------------------
	// compute flux matrix G
	//----------------------
	for (Myint ineigh = 0; ineigh < 2; ineigh++)
	{
		// loop on polynomial order of neighbour element
		for (Myint imat2=0; imat2 <= MAX_POLY_ORDER; imat2++)
		{
			Myint nk = imat2 + 1 ;
			MatDoub Mat_Gz(ni, nk) ;

			// quadrature with ng2 = max(ni,nk)
			Myint ng2 = max(ni,nk) ;
			Doub x1 = -1.0 ;
			Doub x2 = +1.0 ;
			VecDoub xg2(ng2) ;
			VecDoub wg2(ng2) ;

			//------------------------------------
			// compute quadrature points with nr3
			// --> case Gauss-Legendre
			//------------------------------------
			if (node_integ == NODE_INTEG_GL)
			{
				gauleg(x1, x2, xg2, wg2) ;
			}
			//------------------------------------
			// compute quadrature points with gll.h
			// --> case Gauss-Lobatto-Legendre
			//------------------------------------
			else if (node_integ == NODE_INTEG_GLL)
			{
				if (nk == 1)
				{
					gauleg(x1, x2, xg2, wg2) ;
				}
				else
				{
					gll(xg2, wg2) ;
				}
			}
			else
			{
				print_error(" Error in FEM_1D::compute_flux_matrices, numerical integration not supported", node_integ) ;
				return(RTN_CODE_KO) ;
			}

			//--------------------------------------
			// compute coordinates of element nodes
			//--------------------------------------
			VecDoub xnode_k(nk) ;

			//--------------------------------------
			// --> case equidistant node
			//--------------------------------------
			if (node_type == NODE_TYPE_EQD)
			{
				eqd(xnode_k) ;
			}
			//--------------------------------------
			// --> case GLL node
			//--------------------------------------
			else if (node_type == NODE_TYPE_GLL)
			{
				if (nk == 1)
				{
					xnode_k[0] = 0.0 ;
				}
				else
				{
					gll(xnode_k) ;
				}
			}
			else
			{
				print_error(" Error in FEM_1D::compute_flux_matrices, node type not supported", node_integ) ;
				return(RTN_CODE_KO) ;
			}

			Rtn_code rtn_code = compute_flux_matrix_G(xg2, wg2, xnode, xnode_k, Mat_Gz, (Neigh_type) ineigh) ;
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

			*out_file << "\n*** Matrix G ineigh " << ineigh << " *** \n" ;
			for (int ii=0; ii< Mat_Gz.nrows(); ii++)
			{
				for (int jj=0; jj< Mat_Gz.ncols(); jj++)
				{
					*out_file << Mat_Gz[ii][jj] << " " ;
				}
				*out_file << endl ;
			}

			// convert full matrix to sparse matrix
			full2sparse(&pMat_Gz_i[ineigh][imat][imat2], &pMat_Gz_j[ineigh][imat][imat2], &pMat_Gz_c[ineigh][imat][imat2], Mat_Gz) ;

			// print sparse matrix
			Myint*   Mat_Gz_i = pMat_Gz_i[ineigh][imat][imat2]->pArray ;
			Myint*   Mat_Gz_j = pMat_Gz_j[ineigh][imat][imat2]->pArray ;
			Myfloat* Mat_Gz_c = pMat_Gz_c[ineigh][imat][imat2]->pArray ;
			*out_file << "\ni j c / sparse format\n" ;
			for (Myint ii=0; ii< pMat_Gz_i[ineigh][imat][imat2]->nz; ii++)
			{
				*out_file << Mat_Gz_i[ii] << " " << Mat_Gz_j[ii] << " " << Mat_Gz_c[ii] << endl ;
			}
		}
	} // for (Myint ineigh = 0; ineigh < 2; ineigh++)

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D_2nd::compute_flux_matrices");
	return(RTN_CODE_OK) ;
}

} // namespace django
