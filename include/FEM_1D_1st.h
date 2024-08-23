#ifndef DJANGO_FEM_1D_1ST_H_
#define DJANGO_FEM_1D_1ST_H_

#include "FEM_1D.h"

#include "model.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class FEM_1D_1st: public FEM_1D

{
public:

	FEM_1D_1st(void) ;

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

	// reset scheme from another scheme
	// used for dynamic p-adaptivity
	virtual Rtn_code reset_scheme_from(FEM* schemeIn) ;

	// flux matrix F
	//==============
	// matrices are stored in sparse format (i, j, c) i=line, j=column, c=coefficient
	// pMat_F_i[ineigh][ii] is a pointer to the matrix (Grid_1D_float) to be used for:
	// - an element with polynomial order ii
	// - and a neighbour ineigh (I_ZPREV, I_ZNEXT)
	//
	// example:
	// Myint* const pMat_i = pMat_F_i[ineigh][ni-1]->pArray ;
	// Myint* const pMat_j = pMat_F_j[ineigh][ni-1]->pArray ;
	// Myfloat* const pMat_c = pMat_F_c[ineigh][ni-1]->pArray ;

	Grid_1D_int*   pMat_F_i[2][MAX_POLY_ORDER+1] ;
	Grid_1D_int*   pMat_F_j[2][MAX_POLY_ORDER+1] ;
	Grid_1D_float* pMat_F_c[2][MAX_POLY_ORDER+1] ;

	// flux matrix G
	//==============
	// matrices are stored in sparse format (i, j, c) i=line, j=column, c=coefficient
	// pMat_G_i[ineigh][ii][jj] is a pointer to the matrix (Grid_1D_float) to be used for:
	// - an element with polynomial order ii
	// - a neighbour element with polynomial order jj
	// - and a neighbour ineigh (I_ZPREV, I_ZNEXT)
	//
	// example:
	// Myint* const pMat_i = pMat_G_i[ineigh][ni-1][nk-1]->pArray ;
	// Myint* const pMat_j = pMat_G_j[ineigh][ni-1][nk-1]->pArray ;
	// Myfloat* const pMat_c = pMat_G_c[ineigh][ni-1][nk-1]->pArray ;

	Grid_1D_int*   pMat_G_i[2][MAX_POLY_ORDER+1][MAX_POLY_ORDER+1] ;
	Grid_1D_int*   pMat_G_j[2][MAX_POLY_ORDER+1][MAX_POLY_ORDER+1] ;
	Grid_1D_float* pMat_G_c[2][MAX_POLY_ORDER+1][MAX_POLY_ORDER+1] ;

} ;

} // namespace django

#endif
