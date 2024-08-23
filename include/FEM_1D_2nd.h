#ifndef DJANGO_FEM_1D_2ND_H_
#define DJANGO_FEM_1D_2ND_H_

#include "FEM_1D.h"

#include "model.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class FEM_1D_2nd: public FEM_1D

{
public:

	FEM_1D_2nd(void) ;

	// initialize scheme
	virtual Rtn_code initialize(Model*, Myfloat fmax) ;

	// finalize scheme
	virtual Rtn_code finalize(void) ;

protected:

	// compute flux matrix
	virtual Rtn_code compute_flux_matrices(Myint imat, ofstream *out_file, Myint ni, VecDoub &xgl, VecDoub& wgl, VecDoub& xnode) ;

	// compute flux F matrix
	Rtn_code compute_flux_matrix_F(VecDoub &xgl, VecDoub& wgl, VecDoub& xnode, MatDoub &Mat_F, Neigh_type neigh) ;

	// compute flux G matrix
	Rtn_code compute_flux_matrix_G(VecDoub &xgl, VecDoub& wgl, VecDoub& xnode_i, VecDoub& xnode_k, MatDoub &Mat_G, Neigh_type neigh) ;

	// flux matrix Fz
	//===============
	// matrices are stored in sparse format (i, j, c) i=line, j=column, c=coefficient
	// pMat_Fz_i[ineigh][ii] is a pointer to the matrix (Grid_1D_float) to be used for:
	// - an element with polynomial order ii
	// - and a neighbour ineigh (I_ZPREV, I_ZNEXT)
	//
	// example:
	// Myint* const pMat_i = pMat_Fz_i[ineigh][ni-1]->pArray ;
	// Myint* const pMat_j = pMat_Fz_j[ineigh][ni-1]->pArray ;
	// Myfloat* const pMat_c = pMat_Fz_c[ineigh][ni-1]->pArray ;

	Grid_1D_int*   pMat_Fz_i[2][MAX_POLY_ORDER+1] ;
	Grid_1D_int*   pMat_Fz_j[2][MAX_POLY_ORDER+1] ;
	Grid_1D_float* pMat_Fz_c[2][MAX_POLY_ORDER+1] ;

	// flux matrix Gz
	//===============
	// matrices are stored in sparse format (i, j, c) i=line, j=column, c=coefficient
	// pMat_Gz_i[ineigh][ii][jj] is a pointer to the matrix (Grid_1D_float) to be used for:
	// - an element with polynomial order ii
	// - a neighbour element with polynomial order jj
	// - and a neighbour ineigh (I_ZPREV, I_ZNEXT)
	//
	// example:
	// Myint* const pMat_i = pMat_Gz_i[ineigh][ni-1][nk-1]->pArray ;
	// Myint* const pMat_j = pMat_Gz_j[ineigh][ni-1][nk-1]->pArray ;
	// Myfloat* const pMat_c = pMat_Gz_c[ineigh][ni-1][nk-1]->pArray ;

	Grid_1D_int*   pMat_Gz_i[2][MAX_POLY_ORDER+1][MAX_POLY_ORDER+1] ;
	Grid_1D_int*   pMat_Gz_j[2][MAX_POLY_ORDER+1][MAX_POLY_ORDER+1] ;
	Grid_1D_float* pMat_Gz_c[2][MAX_POLY_ORDER+1][MAX_POLY_ORDER+1] ;

} ;

} // namespace django

#endif
