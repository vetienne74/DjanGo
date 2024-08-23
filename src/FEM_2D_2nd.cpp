//-------------------------------------------------------------------------------------------------------
//
// MODELLING WITH MGM IN 2D
//
// PARENT CLASS: Scheme
//   DERIVED CLASS: FEM
//     DERIVED CLASS: FEM_2D
//       DERIVED CLASS: FEM_2D_2nd
//
//-------------------------------------------------------------------------------------------------------

#include "FEM_2D_2nd.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>

#include "output_report.h"
#include "type_def.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------

  FEM_2D_2nd::FEM_2D_2nd(void)  : FEM_2D()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_2nd::FEM_2D_2nd");

	eq_order = ORDER_2ND ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_2nd::FEM_2D_2nd");
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_2D_2nd::initialize(Model* pModel, Myfloat fmax)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_2nd::initialize");

	// call parent initialization
	Rtn_code rtn_code = FEM_2D::initialize(pModel, fmax) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_2nd::initialize");
	return(RTN_CODE_OK) ;
} ;

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_2D_2nd::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_2nd::finalize");

	// call parent finalize
	Rtn_code rtn_code = FEM_2D::finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// deallocate flux matrices
	for (Myint ineigh = 0; ineigh < 4; ineigh++)
	{
		for (Myint imat=0; imat <= MAX_POLY_ORDER; imat++)
		{
			if (pMat_Fz_i[ineigh][imat] != NULL) delete(pMat_Fz_i[ineigh][imat]) ;
			if (pMat_Fz_j[ineigh][imat] != NULL) delete(pMat_Fz_j[ineigh][imat]) ;
			if (pMat_Fz_c[ineigh][imat] != NULL) delete(pMat_Fz_c[ineigh][imat]) ;

			if (pMat_Fx_i[ineigh][imat] != NULL) delete(pMat_Fx_i[ineigh][imat]) ;
			if (pMat_Fx_j[ineigh][imat] != NULL) delete(pMat_Fx_j[ineigh][imat]) ;
			if (pMat_Fx_c[ineigh][imat] != NULL) delete(pMat_Fx_c[ineigh][imat]) ;

			for (Myint imat2=0; imat2 <= MAX_POLY_ORDER; imat2++)
			{
				if (pMat_Gz_i[ineigh][imat][imat2] != NULL) delete(pMat_Gz_i[ineigh][imat][imat2]) ;
				if (pMat_Gz_j[ineigh][imat][imat2] != NULL) delete(pMat_Gz_j[ineigh][imat][imat2]) ;
				if (pMat_Gz_c[ineigh][imat][imat2] != NULL) delete(pMat_Gz_c[ineigh][imat][imat2]) ;

				if (pMat_Gx_i[ineigh][imat][imat2] != NULL) delete(pMat_Gx_i[ineigh][imat][imat2]) ;
				if (pMat_Gx_j[ineigh][imat][imat2] != NULL) delete(pMat_Gx_j[ineigh][imat][imat2]) ;
				if (pMat_Gx_c[ineigh][imat][imat2] != NULL) delete(pMat_Gx_c[ineigh][imat][imat2]) ;
			}
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_2nd::finalize");
	return(RTN_CODE_OK) ;
} ;

//-----------------------------------------------------------------------------------------
// compute flux matrix F with quadrature rule such as Gauss-Legendre
// in the reference element with xi in [-1, 1]
//
// (Fx_i)rj = Int_elem_face_ik ( phi_ir x d(phi_ij)/dx ) dv
//
// (Fz_i)rj = Int_elem_face_ik ( phi_ir x d(phi_ij)/dz ) dv
//
// with i = current element, k = neigbhour element
//      phi_ir = basis function in element i for node r
//      phi_ij = basis function in element i for node j 
//      r = line, j = column in matrix
//      Int_elem_face_ik = integral over surface shared between i and k
//
// INPUT
// xg (VecDoub) = array of quadrature coordinates
// wg (VecDoub) = array of quadrature weights
// xnode (VecDoub) = array of node coordinates in element
// neigh (Neigh_type) = face type (I_ZPREV, I_ZNEXT, I_XPREV, I_XNEXT) for 2D
// axis = AXIS_Z or AXIS_X
//
// OUTPUT
// Mat_F (MatDoub) flux matrix
//-----------------------------------------------------------------------------------------
Rtn_code FEM_2D_2nd::compute_flux_matrix_F(VecDoub &xg, VecDoub& wg, VecDoub& xnode, MatDoub &Mat_F, Neigh_type neigh, Axis_type axis)
{   
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_2nd::compute_flux_matrix_F");

	// # nodes in 1D
	Myint nn = xnode.size() ;
	if (nn == 0)
	{
		print_error(" Error in FEM_2D_2nd::compute_flux_matrix_F, xnode.size() = 0") ;
		return(RTN_CODE_KO) ;
	}

	// check size of matrix
	Myint ni1= sqrt(Mat_F.ncols()) ;
	if (nn != ni1)
	{
		print_error(" Error in FEM_2D_2nd::compute_flux_matrix_F, nn != ni1") ;
		return(RTN_CODE_KO) ;
	}
	Myint ni2= sqrt(Mat_F.nrows()) ;
	if (nn != ni2)
	{
		print_error(" Error in FEM_2D_2nd::compute_flux_matrix_F, nn != ni2") ;
		return(RTN_CODE_KO) ;
	}

	// # GL points
	Myint ng = xg.size() ;
	if (ng == 0)
	{
		print_error(" Error in FEM_2D_2nd::compute_flux_matrix_F, xg.size() = 0") ;
		return(RTN_CODE_KO) ;
	}

	// initialize flux matrix
	for (Myint ii = 0; ii < Mat_F.nrows(); ii++)
	{
		for (Myint jj = 0; jj < Mat_F.ncols(); jj++)
		{
			Mat_F[ii][jj] = 0.0 ;
		}
	}

	Myint nxi   = nn ;
	Myint neta  = nn ;
	Myint netag = ng ;
	Myint nxig  = ng ;

	// check neighbour
	if ((neigh != I_ZPREV) && (neigh != I_ZNEXT) && (neigh != I_XPREV) && (neigh != I_XNEXT))
	{
		print_error(" Error in FEM_2D_2nd::compute_flux_matrix_F, invalid neighbour type", neigh) ;
		return(RTN_CODE_KO) ;
	}

	// loop on element nodes i
	Myint ii = -1 ;
	for (Myint ieta_i=0; ieta_i<neta; ieta_i++)
	{
		for (Myint ixi_i=0; ixi_i<nxi; ixi_i++)
		{
			ii++ ;

			// loop on element nodes j
			Myint jj = -1 ;
			for (Myint ieta_j=0; ieta_j<neta; ieta_j++)
			{
				for (Myint ixi_j=0; ixi_j<nxi; ixi_j++)
				{
					jj++ ;

					// loop on quadrature points
					for (Myint ieta_g=0; ieta_g<netag; ieta_g++)
					{
						for (Myint ixi_g=0; ixi_g<nxig; ixi_g++)
						{

							// compute phi_i at [xg(ixi_g), xg(ieta_g)]
							//=========================================

							// phi_ieta_i at xg[ieta_g]
							Doub phi_ieta_i = lagran(xnode, ieta_i, xg[ieta_g]) ;

							// phi_ixi_i at xg[ixi_g]
							Doub phi_ixi_i = lagran(xnode, ixi_i, xg[ixi_g]) ;

							// compute phi_i
							Doub phi_i = phi_ieta_i * phi_ixi_i ;

							// compute d(phi_j)/dxi at [xg(ixi_g), xg(ieta_g)]
							//================================================
							Doub dphi_j = 0.0 ;
							if (axis == X_AXIS)
							{
								// phi_ieta_i at xg[ieta_g]
								Doub phi_ieta_j = lagran(xnode, ieta_j, xg[ieta_g]) ;

								// d(phi_ixi_j)/dxi at xg[ixi_g]
								Doub dphi_ixi_j = dlagran(xnode, ixi_j, xg[ixi_g]) ;

								// compute phi_i
								dphi_j = phi_ieta_j * dphi_ixi_j ;
							}

							// compute d(phi_j)/deta at [xg(ixi_g), xg(ieta_g)]
							//================================================
							else if (axis == Z_AXIS)
							{
								// d(phi_ieta_j)deta at xg[ieta_g]
								Doub dphi_ieta_j = dlagran(xnode, ieta_j, xg[ieta_g]) ;

								// phi_ixi_j at xg[ixi_g]
								Doub phi_ixi_j = lagran(xnode, ixi_j, xg[ixi_g]) ;

								// compute phi_i
								dphi_j = dphi_ieta_j * phi_ixi_j ;
							}

							// update flux matrix
							//===================
							Myfloat weight = 0.0 ;
							if ((neigh == I_ZPREV) && (ieta_g == 0))
							{
								weight = wg[ixi_g] ;
							}
							else if ((neigh == I_ZNEXT) && (ieta_g == netag-1))
							{
								weight = wg[ixi_g] ;
							}
							else if ((neigh == I_XPREV) && (ixi_g == 0))
							{
								weight = wg[ieta_g] ;
							}
							else if ((neigh == I_XNEXT) && (ixi_g == nxig-1))
							{
								weight = wg[ieta_g] ;
							}
							Mat_F[ii][jj] = Mat_F[ii][jj] + weight * phi_i * dphi_j ;

						} // for (Myint ixi_g=0; ixi_g<nxig; ixi_g++)
					} // for (Myint ieta_g=0; ieta_g<netag; ieta_g++)

				} // for (Myint ixi_j=0; ixi_j<nxi; ixi_j++)
			} // for (Myint ieta_j=0; ieta_j<neta; ieta_j++)

		} // for (Myint ixi_i=0; ixi_i<nxi; ixi_i++)
	} // for (Myint ieta_i=0; ieta_i<neta; ieta_i++)

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_2nd::compute_flux_matrix_F");
	return(RTN_CODE_OK) ;
}

//-----------------------------------------------------------------------------------------
// compute flux matrix G with quadrature rule such as Gauss-Legendre
// in the reference element with xi in [-1, 1]
//
// (Gx_ik)rj = Int_elem_face_ik (phi_ir x d(phi_kj)/dx ) dv
//
// (Gz_ik)rj = Int_elem_face_ik (phi_ir x d(phi_kj)/dz ) dv
//
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
// neigh (Neigh_type) = face type (I_ZPREV, I_ZNEXT, I_XPREV, I_XNEXT) for 2D
//
// OUTPUT
// Mat_G (MatDoub) flux matrix
//-----------------------------------------------------------------------------------------
Rtn_code FEM_2D_2nd::compute_flux_matrix_G(VecDoub &xg, VecDoub& wg,
		VecDoub& xnode_i, VecDoub& xnode_k, MatDoub &Mat_G, Neigh_type neigh, Axis_type axis)
{   
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_2nd::compute_flux_matrix_G");

	// # nodes in 1D (element i)
	Myint ni = xnode_i.size() ;
	if (ni == 0)
	{
		print_error(" Error in FEM_2D_2nd::compute_flux_matrix_G, xnode_i.size() = 0") ;
		return(RTN_CODE_KO) ;
	}

	// check size of matrix
	// nrows = ni*ni
	Myint ni2= sqrt(Mat_G.nrows()) ;
	if (ni != ni2)
	{
		print_error(" Error in FEM_2D_2nd::compute_flux_matrix_G, ni != ni2") ;
		return(RTN_CODE_KO) ;
	}

	// # nodes in 1D (element k)
	Myint nk = xnode_k.size() ;
	if (nk == 0)
	{
		print_error(" Error in FEM_2D_2nd::compute_flux_matrix_G, xnode_k.size() = 0") ;
		return(RTN_CODE_KO) ;
	}

	// check size of matrix
	// ncols = nk*nk
	Myint nk2= sqrt(Mat_G.ncols()) ;
	if (nk != nk2)
	{
		print_error(" Error in FEM_2D_2nd::compute_flux_matrix_G, nk != nk2") ;
		return(RTN_CODE_KO) ;
	}

	// # GL points in element i
	Myint ng = xg.size() ;
	if (ng == 0)
	{
		print_error(" Error in FEM_2D_2nd::compute_flux_matrix_G, xg.size() = 0") ;
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

	Myint nxi_i  = ni ;
	Myint neta_i = ni ;
	Myint nxi_k  = nk ;
	Myint neta_k = nk ;
	Myint netag  = ng ;
	Myint nxig   = ng ;

	// check neighbour
	if ((neigh != I_ZPREV) && (neigh != I_ZNEXT) && (neigh != I_XPREV) && (neigh != I_XNEXT))
	{
		print_error(" Error in FEM_2D_2nd::compute_flux_matrix_G, invalid neighbour type", neigh) ;
		return(RTN_CODE_KO) ;
	}

	// loop on element i nodes r
	Myint ii = -1 ;
	for (Myint ieta_ir=0; ieta_ir<neta_i; ieta_ir++)
	{
		for (Myint ixi_ir=0; ixi_ir<nxi_i; ixi_ir++)
		{
			ii++ ;

			// loop on element k nodes j
			Myint jj = -1 ;
			for (Myint ieta_kj=0; ieta_kj<neta_k; ieta_kj++)
			{
				for (Myint ixi_kj=0; ixi_kj<nxi_k; ixi_kj++)
				{
					jj++ ;

					// loop on quadrature points
					for (Myint ieta_g=0; ieta_g<netag; ieta_g++)
					{
						for (Myint ixi_g=0; ixi_g<nxig; ixi_g++)
						{

							// update flux matrix
							//===================
							Doub weight = 0.0 ;
							Doub phi_ir = 0.0 ;
							Doub phi_kj = 0.0 ;

							if (axis == Z_AXIS)
							{
								if ((neigh == I_ZPREV) && (ieta_g == 0))
								{
									if (ieta_ir == 0)        phi_ir = lagran(xnode_i, ixi_ir, xg[ixi_g]) ;
									phi_kj = lagran(xnode_k, ixi_kj, xg[ixi_g]) * dlagran(xnode_k, ieta_kj, xg[neta_k-1]);
									weight = wg[ixi_g] ;
								}

								else if ((neigh == I_ZNEXT) && (ieta_g == netag-1))
								{
									if (ieta_ir == neta_i-1) phi_ir = lagran(xnode_i, ixi_ir, xg[ixi_g]) ;
									phi_kj = lagran(xnode_k, ixi_kj, xg[ixi_g]) * dlagran(xnode_k, ieta_kj, xg[0]);
									weight = wg[ixi_g] ;
								}

								else if ((neigh == I_XPREV) && (ixi_g == 0))
								{
									if (ixi_ir == 0)         phi_ir = lagran(xnode_i, ieta_ir, xg[ieta_g]) ;
									if (ixi_kj == nxi_k-1)   phi_kj = dlagran(xnode_k, ieta_kj, xg[ieta_g]) ;
									weight = wg[ieta_g] ;
								}
								else if ((neigh == I_XNEXT) && (ixi_g == nxig-1))
								{
									if (ixi_ir == nxi_i-1)   phi_ir = lagran(xnode_i, ieta_ir, xg[ieta_g]) ;
									if (ixi_kj == 0)         phi_kj = dlagran(xnode_k, ieta_kj, xg[ieta_g]) ;
									weight = wg[ieta_g] ;
								}
								Mat_G[ii][jj] = Mat_G[ii][jj] + weight * phi_ir * phi_kj ;
							}

							if (axis == X_AXIS)
							{
								if ((neigh == I_ZPREV) && (ieta_g == 0))
								{
									if (ieta_ir == 0)        phi_ir = lagran(xnode_i, ixi_ir, xg[ixi_g]) ;
									if (ieta_kj == neta_k-1) phi_kj = dlagran(xnode_k, ixi_kj, xg[ixi_g]) ;
									weight = wg[ixi_g] ;
								}

								else if ((neigh == I_ZNEXT) && (ieta_g == netag-1))
								{
									if (ieta_ir == neta_i-1) phi_ir = lagran(xnode_i, ixi_ir, xg[ixi_g]) ;
									if (ieta_kj == 0)        phi_kj = dlagran(xnode_k, ixi_kj, xg[ixi_g]) ;
									weight = wg[ixi_g] ;
								}

								else if ((neigh == I_XPREV) && (ixi_g == 0))
								{
									if (ixi_ir == 0)         phi_ir = lagran(xnode_i, ieta_ir, xg[ieta_g]) ;
									phi_kj = lagran(xnode_k, ieta_kj, xg[ieta_g]) * dlagran(xnode_k, ixi_kj, xg[nxi_k-1]) ;
									weight = wg[ieta_g] ;
								}
								else if ((neigh == I_XNEXT) && (ixi_g == nxig-1))
								{
									if (ixi_ir == nxi_i-1)   phi_ir = lagran(xnode_i, ieta_ir, xg[ieta_g]) ;
									phi_kj = lagran(xnode_k, ieta_kj, xg[ieta_g]) * dlagran(xnode_k, ixi_kj, xg[0]) ;
									weight = wg[ieta_g] ;
								}
								Mat_G[ii][jj] = Mat_G[ii][jj] + weight * phi_ir * phi_kj ;
							}


						} // for (Myint ixi_g=0; ixi_g<nxig; ixi_g++)
					} // for (Myint ieta_g=0; ieta_g<netag; ieta_g++)

				} // for (Myint ixi_kj=0; ixi_kj<nxi_k; ixi_kj++)
			} // for (Myint ieta_kj=0; ieta_kj<neta_k; ieta_kj++)

		} // for (Myint ixi_ir=0; ixi_ir<nxi_i; ixi_ir++)
	} // for (Myint ieta_ir=0; ieta_ir<neta_i; ieta_ir++)

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_2nd::compute_flux_matrix_G");
	return(RTN_CODE_OK) ;
}

//-----------------------------------------------------------------------------------------
Rtn_code FEM_2D_2nd::compute_flux_matrices(Myint imat, ofstream *out_file, Myint ni, VecDoub &xg, VecDoub& wg, VecDoub& xnode)
{   
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_2nd::compute_flux_matrices");

	//-----------------------
	// compute flux matrix Fz
	//-----------------------
	for (Myint ineigh = 0; ineigh < 4; ineigh++)
	{
		MatDoub Mat_Fz(ni*ni, ni*ni) ;
		Rtn_code rtn_code = compute_flux_matrix_F(xg, wg, xnode, Mat_Fz, (Neigh_type) ineigh, Z_AXIS) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		*out_file << "\n*** Matrix Fz ineigh " << ineigh << " *** \n" ;
		// for (int ii=0; ii< Mat_F.nrows(); ii++)
		//   {
		//     for (int jj=0; jj< Mat_F.ncols(); jj++)
		// 	{
		// 	  *out_file << Mat_F[ii][jj] << " " ;
		// 	}
		//     *out_file << endl ;
		//   }

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
	} // for (Myint ineigh = 0; ineigh < 4; ineigh++)

	//-----------------------
	// compute flux matrix Fx
	//-----------------------
	for (Myint ineigh = 0; ineigh < 4; ineigh++)
	{
		MatDoub Mat_Fx(ni*ni, ni*ni) ;
		Rtn_code rtn_code = compute_flux_matrix_F(xg, wg, xnode, Mat_Fx, (Neigh_type) ineigh, X_AXIS) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		*out_file << "\n*** Matrix Fx ineigh " << ineigh << " *** \n" ;
		// for (int ii=0; ii< Mat_F.nrows(); ii++)
		//   {
		//     for (int jj=0; jj< Mat_F.ncols(); jj++)
		// 	{
		// 	  *out_file << Mat_F[ii][jj] << " " ;
		// 	}
		//     *out_file << endl ;
		//   }

		// convert full matrix to sparse matrix
		full2sparse(&pMat_Fx_i[ineigh][imat], &pMat_Fx_j[ineigh][imat], &pMat_Fx_c[ineigh][imat], Mat_Fx) ;

		// print sparse matrix
		Myint*   Mat_Fx_i = pMat_Fx_i[ineigh][imat]->pArray ;
		Myint*   Mat_Fx_j = pMat_Fx_j[ineigh][imat]->pArray ;
		Myfloat* Mat_Fx_c = pMat_Fx_c[ineigh][imat]->pArray ;
		*out_file << "\ni j c / sparse format\n" ;
		for (Myint ii=0; ii< pMat_Fx_i[ineigh][imat]->nz; ii++)
		{
			*out_file << Mat_Fx_i[ii] << " " << Mat_Fx_j[ii] << " " << Mat_Fx_c[ii] << endl ;
		}
	} // for (Myint ineigh = 0; ineigh < 4; ineigh++)

	//-----------------------
	// compute flux matrix Gz
	//-----------------------
	for (Myint ineigh = 0; ineigh < 4; ineigh++)
	{
		// loop on polynomial order of neighbour element
		for (Myint imat2=0; imat2 <= MAX_POLY_ORDER; imat2++)
		{
			Myint nk = imat2 + 1 ;

			MatDoub Mat_Gz(ni*ni, nk*nk) ;
			MatDoub Mat_Gx(ni*ni, nk*nk) ;

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
				print_error(" Error in FEM_2D_2nd::compute_flux_matrices, numerical integration not supported", node_integ) ;
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
				print_error(" Error in FEM_2D_2nd::compute_flux_matrices, node type not supported", node_integ) ;
				return(RTN_CODE_KO) ;
			}

			// compute Gz
			Rtn_code rtn_code = compute_flux_matrix_G(xg2, wg2, xnode, xnode_k, Mat_Gz, (Neigh_type) ineigh, Z_AXIS) ;
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

			*out_file << "\n*** Matrix Gz ineigh " << ineigh << " *** \n" ;
			*out_file << "\n*** imat " << imat << " imat2 " << imat2 << " ***\n" ;
			// for (int ii=0; ii< Mat_G.nrows(); ii++)
			// 	{
			// 	  for (int jj=0; jj< Mat_G.ncols(); jj++)
			// 	    {
			// 	      *out_file << Mat_G[ii][jj] << " " ;
			// 	    }
			// 	  *out_file << endl ;
			// 	}

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

			// compute Gx
			rtn_code = compute_flux_matrix_G(xg2, wg2, xnode, xnode_k, Mat_Gx, (Neigh_type) ineigh, X_AXIS) ;
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

			*out_file << "\n*** Matrix Gx ineigh " << ineigh << " *** \n" ;
			*out_file << "\n*** imat " << imat << " imat2 " << imat2 << " ***\n" ;
			// for (int ii=0; ii< Mat_G.nrows(); ii++)
			// 	{
			// 	  for (int jj=0; jj< Mat_G.ncols(); jj++)
			// 	    {
			// 	      *out_file << Mat_G[ii][jj] << " " ;
			// 	    }
			// 	  *out_file << endl ;
			// 	}

			// convert full matrix to sparse matrix
			full2sparse(&pMat_Gx_i[ineigh][imat][imat2], &pMat_Gx_j[ineigh][imat][imat2], &pMat_Gx_c[ineigh][imat][imat2], Mat_Gx) ;

			// print sparse matrix
			Myint*   Mat_Gx_i = pMat_Gx_i[ineigh][imat][imat2]->pArray ;
			Myint*   Mat_Gx_j = pMat_Gx_j[ineigh][imat][imat2]->pArray ;
			Myfloat* Mat_Gx_c = pMat_Gx_c[ineigh][imat][imat2]->pArray ;
			*out_file << "\ni j c / sparse format\n" ;
			for (Myint ii=0; ii< pMat_Gx_i[ineigh][imat][imat2]->nz; ii++)
			{
				*out_file << Mat_Gx_i[ii] << " " << Mat_Gx_j[ii] << " " << Mat_Gx_c[ii] << endl ;
			}
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_2nd::compute_flux_matrices");
	return(RTN_CODE_OK) ;
}

} // namespace django
