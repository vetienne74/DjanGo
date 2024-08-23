//-------------------------------------------------------------------------------------------------------
//
// MODELLING WITH FEM IN 2D
//
// PARENT CLASS: Scheme
//   DERIVED CLASS: FEM
//     DERIVED CLASS: FEM_2D
//       DERIVED CLASS: FEM_2D_1st
//
//-------------------------------------------------------------------------------------------------------

#include "FEM_2D_1st.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>

#include "allocate_array.h"
#include "output_report.h"
#include "type_def.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------

  FEM_2D_1st::FEM_2D_1st(void)  : FEM_2D()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st::FEM_2D_1st");

	eq_order = ORDER_1ST ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st::FEM_2D_1st");
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_2D_1st::initialize(Model* pModel, Myfloat fmax)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st::initialize");

	// call parent initialization
	Rtn_code rtn_code = FEM_2D::initialize(pModel, fmax) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st::initialize");
	return(RTN_CODE_OK) ;
} ;

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_2D_1st::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st::finalize");

	// call parent finalize
	Rtn_code rtn_code = FEM_2D::finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// deallocate flux matrices
	for (Myint ineigh = 0; ineigh < 4; ineigh++)
	{
		for (Myint imat=0; imat <= MAX_POLY_ORDER; imat++)
		{
			if (pMat_F_i[ineigh][imat] != NULL) delete(pMat_F_i[ineigh][imat]) ;
			if (pMat_F_j[ineigh][imat] != NULL) delete(pMat_F_j[ineigh][imat]) ;
			if (pMat_F_c[ineigh][imat] != NULL) delete(pMat_F_c[ineigh][imat]) ;

			for (Myint imat2=0; imat2 <= MAX_POLY_ORDER; imat2++)
			{
				if (pMat_G_i[ineigh][imat][imat2] != NULL) delete(pMat_G_i[ineigh][imat][imat2]) ;
				if (pMat_G_j[ineigh][imat][imat2] != NULL) delete(pMat_G_j[ineigh][imat][imat2]) ;
				if (pMat_G_c[ineigh][imat][imat2] != NULL) delete(pMat_G_c[ineigh][imat][imat2]) ;
			}
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st::finalize");
	return(RTN_CODE_OK) ;
} ;

//-----------------------------------------------------------------------------------------
// compute flux matrix F with quadrature rule such as Gauss-Legendre
// in the reference element with xi in [-1, 1]
//
// (Fik)rj = Int_elem_face_ik (phi_ir x phi_ij) dxi
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
//
// OUTPUT
// Mat_F (MatDoub) flux matrix
//-----------------------------------------------------------------------------------------
Rtn_code FEM_2D_1st::compute_flux_matrix_F(VecDoub &xg, VecDoub& wg, VecDoub& xnode, MatDoub &Mat_F, Neigh_type neigh)
{   
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st::compute_flux_matrix_F");

	// # nodes in 1D
	Myint nn = xnode.size() ;
	if (nn == 0)
	{
		print_error(" Error in FEM_2D_1st::compute_flux_matrix_F, xnode.size() = 0") ;
		return(RTN_CODE_KO) ;
	}

	// check size of matrix
	Myint ni1= sqrt(Mat_F.ncols()) ;
	if (nn != ni1)
	{
		print_error(" Error in FEM_2D_1st::compute_flux_matrix_F, nn != ni1") ;
		return(RTN_CODE_KO) ;
	}
	Myint ni2= sqrt(Mat_F.nrows()) ;
	if (nn != ni2)
	{
		print_error(" Error in FEM_2D_1st::compute_flux_matrix_F, nn != ni2") ;
		return(RTN_CODE_KO) ;
	}

	// # GL points
	Myint ng = xg.size() ;
	if (ng == 0)
	{
		print_error(" Error in FEM_2D_1st::compute_flux_matrix_F, xg.size() = 0") ;
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
		print_error(" Error in FEM_2D_1st::compute_flux_matrix_F, invalid neighbour type", neigh) ;
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

							// compute phi_j at [xg(ixi_g), xg(ieta_g)]
							//=========================================

							// phi_ieta_j at xg[ieta_g]
							Doub phi_ieta_j = lagran(xnode, ieta_j, xg[ieta_g]) ;

							// phi_ixi_j at xg[ixi_g]
							Doub phi_ixi_j = lagran(xnode, ixi_j, xg[ixi_g]) ;

							// compute phi_j
							Doub phi_j = phi_ieta_j * phi_ixi_j ;

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
							Mat_F[ii][jj] = Mat_F[ii][jj] + weight * phi_i * phi_j ;

						} // for (Myint ixi_g=0; ixi_g<nxig; ixi_g++)
					} // for (Myint ieta_g=0; ieta_g<netag; ieta_g++)

				} // for (Myint ixi_j=0; ixi_j<nxi; ixi_j++)
			} // for (Myint ieta_j=0; ieta_j<neta; ieta_j++)

		} // for (Myint ixi_i=0; ixi_i<nxi; ixi_i++)
	} // for (Myint ieta_i=0; ieta_i<neta; ieta_i++)

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st::compute_flux_matrix_F");
	return(RTN_CODE_OK) ;
}

//-----------------------------------------------------------------------------------------
// compute flux matrix G with quadrature rule such as Gauss-Legendre
// in the reference element with xi in [-1, 1]
//
// (Gik)rj = Int_elem_face_ik (phi_ir x phi_kj) dxi
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
Rtn_code FEM_2D_1st::compute_flux_matrix_G(VecDoub &xg, VecDoub& wg,
		VecDoub& xnode_i, VecDoub& xnode_k, MatDoub &Mat_G, Neigh_type neigh)
{   
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st::compute_flux_matrix_G");

	// # nodes in 1D (element i)
	Myint ni = xnode_i.size() ;
	if (ni == 0)
	{
		print_error(" Error in FEM_2D_1st::compute_flux_matrix_G, xnode_i.size() = 0") ;
		return(RTN_CODE_KO) ;
	}

	// check size of matrix
	// nrows = ni*ni
	Myint ni2= sqrt(Mat_G.nrows()) ;
	if (ni != ni2)
	{
		print_error(" Error in FEM_2D_1st::compute_flux_matrix_G, ni != ni2") ;
		return(RTN_CODE_KO) ;
	}

	// # nodes in 1D (element k)
	Myint nk = xnode_k.size() ;
	if (nk == 0)
	{
		print_error(" Error in FEM_2D_1st::compute_flux_matrix_G, xnode_k.size() = 0") ;
		return(RTN_CODE_KO) ;
	}

	// check size of matrix
	// ncols = nk*nk
	Myint nk2= sqrt(Mat_G.ncols()) ;
	if (nk != nk2)
	{
		print_error(" Error in FEM_2D_1st::compute_flux_matrix_G, nk != nk2") ;
		return(RTN_CODE_KO) ;
	}

	// # GL points in element i
	Myint ng = xg.size() ;
	if (ng == 0)
	{
		print_error(" Error in FEM_2D_1st::compute_flux_matrix_G, xg.size() = 0") ;
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
		print_error(" Error in FEM_2D_1st::compute_flux_matrix_G, invalid neighbour type", neigh) ;
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
							if ((neigh == I_ZPREV) && (ieta_g == 0))
							{
								if (ieta_ir == 0)        phi_ir = lagran(xnode_i, ixi_ir, xg[ixi_g]) ;
								if (ieta_kj == neta_k-1) phi_kj = lagran(xnode_k, ixi_kj, xg[ixi_g]) ;
								weight = wg[ixi_g] ;
							}

							else if ((neigh == I_ZNEXT) && (ieta_g == netag-1))
							{
								if (ieta_ir == neta_i-1) phi_ir = lagran(xnode_i, ixi_ir, xg[ixi_g]) ;
								if (ieta_kj == 0)        phi_kj = lagran(xnode_k, ixi_kj, xg[ixi_g]) ;
								weight = wg[ixi_g] ;
							}

							else if ((neigh == I_XPREV) && (ixi_g == 0))
							{
								if (ixi_ir == 0)         phi_ir = lagran(xnode_i, ieta_ir, xg[ieta_g]) ;
								if (ixi_kj == nxi_k-1)   phi_kj = lagran(xnode_k, ieta_kj, xg[ieta_g]) ;
								weight = wg[ieta_g] ;
							}
							else if ((neigh == I_XNEXT) && (ixi_g == nxig-1))
							{
								if (ixi_ir == nxi_i-1)   phi_ir = lagran(xnode_i, ieta_ir, xg[ieta_g]) ;
								if (ixi_kj == 0)         phi_kj = lagran(xnode_k, ieta_kj, xg[ieta_g]) ;
								weight = wg[ieta_g] ;
							}
							Mat_G[ii][jj] = Mat_G[ii][jj] + weight * phi_ir * phi_kj ;

						} // for (Myint ixi_g=0; ixi_g<nxig; ixi_g++)
					} // for (Myint ieta_g=0; ieta_g<netag; ieta_g++)

				} // for (Myint ixi_kj=0; ixi_kj<nxi_k; ixi_kj++)
			} // for (Myint ieta_kj=0; ieta_kj<neta_k; ieta_kj++)

		} // for (Myint ixi_ir=0; ixi_ir<nxi_i; ixi_ir++)
	} // for (Myint ieta_ir=0; ieta_ir<neta_i; ieta_ir++)

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st::compute_flux_matrix_G");
	return(RTN_CODE_OK) ;
}

//-----------------------------------------------------------------------------------------
Rtn_code FEM_2D_1st::compute_flux_matrices(Myint imat, ofstream *out_file, Myint ni, VecDoub &xg, VecDoub& wg, VecDoub& xnode)
{   
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st::compute_flux_matrices");

	//----------------------
	// compute flux matrix F
	//----------------------
	for (Myint ineigh = 0; ineigh < 4; ineigh++)
	{
		MatDoub Mat_F(ni*ni, ni*ni) ;
		Rtn_code rtn_code = compute_flux_matrix_F(xg, wg, xnode, Mat_F, (Neigh_type) ineigh) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		*out_file << "\n*** Matrix F ineigh " << ineigh << " *** \n" ;
		// for (int ii=0; ii< Mat_F.nrows(); ii++)
		//   {
		//     for (int jj=0; jj< Mat_F.ncols(); jj++)
		// 	{
		// 	  *out_file << Mat_F[ii][jj] << " " ;
		// 	}
		//     *out_file << endl ;
		//   }

		// convert full matrix to sparse matrix
		full2sparse(&pMat_F_i[ineigh][imat], &pMat_F_j[ineigh][imat], &pMat_F_c[ineigh][imat], Mat_F) ;

		// print sparse matrix
		Myint*   Mat_F_i = pMat_F_i[ineigh][imat]->pArray ;
		Myint*   Mat_F_j = pMat_F_j[ineigh][imat]->pArray ;
		Myfloat* Mat_F_c = pMat_F_c[ineigh][imat]->pArray ;
		*out_file << "\ni j c / sparse format\n" ;
		for (Myint ii=0; ii< pMat_F_i[ineigh][imat]->nz; ii++)
		{
			*out_file << Mat_F_i[ii] << " " << Mat_F_j[ii] << " " << Mat_F_c[ii] << endl ;
		}
	} // for (Myint ineigh = 0; ineigh < 4; ineigh++)

	//----------------------
	// compute flux matrix G
	//----------------------
	for (Myint ineigh = 0; ineigh < 4; ineigh++)
	{
		// loop on polynomial order of neighbour element
		for (Myint imat2=0; imat2 <= MAX_POLY_ORDER; imat2++)
		{
			Myint nk = imat2 + 1 ;
			MatDoub Mat_G(ni*ni, nk*nk) ;

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
				print_error(" Error in FEM_2D_1st::compute_flux_matrices, numerical integration not supported", node_integ) ;
				return(RTN_CODE_KO) ;
			}

			//--------------------------------------
			// compute coordinates of element nodes
			//--------------------------------------
			VecDoub *xnode_k = &(xnode_ref_elem_1D[nk-1]) ;

			Rtn_code rtn_code = compute_flux_matrix_G(xg2, wg2, xnode, *xnode_k, Mat_G, (Neigh_type) ineigh) ;
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

			*out_file << "\n*** Matrix G ineigh " << ineigh << " *** \n" ;
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
			full2sparse(&pMat_G_i[ineigh][imat][imat2], &pMat_G_j[ineigh][imat][imat2], &pMat_G_c[ineigh][imat][imat2], Mat_G) ;

			// print sparse matrix
			Myint*   Mat_G_i = pMat_G_i[ineigh][imat][imat2]->pArray ;
			Myint*   Mat_G_j = pMat_G_j[ineigh][imat][imat2]->pArray ;
			Myfloat* Mat_G_c = pMat_G_c[ineigh][imat][imat2]->pArray ;
			*out_file << "\ni j c / sparse format\n" ;
			for (Myint ii=0; ii< pMat_G_i[ineigh][imat][imat2]->nz; ii++)
			{
				*out_file << Mat_G_i[ii] << " " << Mat_G_j[ii] << " " << Mat_G_c[ii] << endl ;
			}
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st::compute_flux_matrices");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_2D_1st::reset_scheme_from(FEM* schemeIn2)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st::reset_scheme_from");

	// call parent
	Rtn_code rtn_code = FEM::reset_scheme_from(schemeIn2) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	FEM_2D_1st* schemeIn = dynamic_cast<FEM_2D_1st*>(schemeIn2) ;
	if (schemeIn == NULL)
	{
		print_error("FEM_2D_1st::reset_scheme_from, not FEM_2D_1st") ;
		return(RTN_CODE_KO) ;
	}

	//--------------------------------------------------------------------------
	// free unused data
	//--------------------------------------------------------------------------

	// deallocate flux matrices in new scheme
	//---------------------------------------
	for (Myint ineigh = 0; ineigh < 4; ineigh++)
	{
		for (Myint imat=0; imat <= MAX_POLY_ORDER; imat++)
		{
			if (schemeIn->pMat_F_i[ineigh][imat] != NULL) delete(schemeIn->pMat_F_i[ineigh][imat]) ;
			if (schemeIn->pMat_F_j[ineigh][imat] != NULL) delete(schemeIn->pMat_F_j[ineigh][imat]) ;
			if (schemeIn->pMat_F_c[ineigh][imat] != NULL) delete(schemeIn->pMat_F_c[ineigh][imat]) ;

			for (Myint imat2=0; imat2 <= MAX_POLY_ORDER; imat2++)
			{
				if (schemeIn->pMat_G_i[ineigh][imat][imat2] != NULL) delete(schemeIn->pMat_G_i[ineigh][imat][imat2]) ;
				if (schemeIn->pMat_G_j[ineigh][imat][imat2] != NULL) delete(schemeIn->pMat_G_j[ineigh][imat][imat2]) ;
				if (schemeIn->pMat_G_c[ineigh][imat][imat2] != NULL) delete(schemeIn->pMat_G_c[ineigh][imat][imat2]) ;
			}
		}
	}

	// deallocate global node array in this scheme
	//--------------------------------------------
	deallocate_array<Node_2D_struct_type>(pNode, nnode) ;

	// copy tmin from initial to new scheme
	//-------------------------------------
	for (Myint ielem = 0; ielem < nelem; ielem++)
	{
		schemeIn->pElement[ielem].tmin = pElement[ielem].tmin ;
	}

	// deallocate local node array in this scheme
	//-------------------------------------------

	for (Myint ielem = 0; ielem < nelem; ielem++)
	{
		deallocate_array<Myint>(pElement[ielem].inode, pElement[ielem].nnode) ;
	}

	// deallocate element array in this scheme
	//----------------------------------------
	deallocate_array<Quad_struct_type>(pElement, nelem) ;

	// deallocate stiffness matrices in the new scheme
	//------------------------------------------------
	for (Myint imat = 0; imat <= MAX_POLY_ORDER; imat++)
	{
		// full matrix
		if (schemeIn->pMat_Dz[imat] != NULL) delete(schemeIn->pMat_Dz[imat]) ;
		if (schemeIn->pMat_Dx[imat] != NULL) delete(schemeIn->pMat_Dx[imat]) ;

		// sparse matrix
		if (schemeIn->pMat_Dz_i[imat] != NULL) delete(schemeIn->pMat_Dz_i[imat]) ;
		if (schemeIn->pMat_Dz_j[imat] != NULL) delete(schemeIn->pMat_Dz_j[imat]) ;
		if (schemeIn->pMat_Dz_c[imat] != NULL) delete(schemeIn->pMat_Dz_c[imat]) ;
		if (schemeIn->pMat_Dx_i[imat] != NULL) delete(schemeIn->pMat_Dx_i[imat]) ;
		if (schemeIn->pMat_Dx_j[imat] != NULL) delete(schemeIn->pMat_Dx_j[imat]) ;
		if (schemeIn->pMat_Dx_c[imat] != NULL) delete(schemeIn->pMat_Dx_c[imat]) ;
	}

	//--------------------------------------------------------------------------
	// reset scheme
	//--------------------------------------------------------------------------
	cout << "old nnode " << nnode << "new nnode " << schemeIn->nnode << "\n" ;
	nnode             = schemeIn->nnode ;
	pElement          = schemeIn->pElement ;
	pNode             = schemeIn->pNode ;
	fmax              = schemeIn->fmax ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st::reset_scheme_from");
	return(RTN_CODE_OK) ;
}

} // namespace django
