#ifndef DJANGO_FEM_2D_H_
#define DJANGO_FEM_2D_H_

#include "FEM.h"

#include "acquisition.h"
#include "grid_1D_float.h"
#include "grid_1D_int.h"
#include "grid_2D_float.h"
#include "model.h"
#include "snapshot.h"
#include "type_def.h"

namespace django {

//=================
// Element quad
//=================

typedef struct Quad
{
	// number of nodes
	Myint nnode ;

	// number of nodes in one direction
	// -> used as index for element matrices
	Myint nnode_1D ;

	// node global indexes
	Myint *inode ;

	// element size (m)
	Myfloat size_x ;
	Myfloat size_z ;

	// min velocity (used for p-adaptivity)
	Myfloat vmin ;

	// max velocity (used for stable dt and dynamic front)
	Myfloat vmax ;

	// tmin (used for dynamic front)
	Myfloat tmin ;

	// region
	Region_type region ;

	// neighbour elements
	// I_ZPREV, I_ZNEXT, I_XPREV, I_XNEXT
	Myint neigh[4] ;

	// normal vector components
	// I_ZPREV, I_ZNEXT, I_XPREV, I_XNEXT
	Myfloat vnz[4] ;
	Myfloat vnx[4] ;

	// flux type
	// I_ZPREV, I_ZNEXT, I_XPREV, I_XNEXT
	Flux_type flux[4] ;

} Quad_struct_type ;

typedef struct Node_2D
{
	// coordinates
	Myfloat zcoord ;
	Myfloat xcoord ;

	// boundary type
	Boundary_type boundary ;

} Node_2D_struct_type ;

//------------------------------------------------------------------------------------

class FEM_2D: public FEM

{
public:

	FEM_2D(void) ;

	// initialize
	virtual Rtn_code initialize(Model*, Myfloat fmax) ;

	// finalize
	virtual Rtn_code finalize(void) ;

	// mesh statistics
	virtual Rtn_code mesh_info(void) ;

	// write mesh in VTK format
	virtual Rtn_code write_mesh_VTK(void) ;

	// write node in VTK format
	virtual Rtn_code write_node_VTK(void) ;

	// Write snapshot
	virtual Rtn_code write_snapshot(Variable*, Myint) ;

	// projet variable from one mesh to another
	Rtn_code project_variable(FEM_2D& schemeSrc, Myint varSrcId, Myint varDestId) ;

protected:

	// array of elements
	Quad_struct_type *pElement ;

	// array of nodes
	Node_2D_struct_type *pNode ;

	// mesh initialization
	virtual Rtn_code initialize_mesh(Model*, Myfloat fmax) ;

	// check time step
	virtual Myfloat compute_optimal_time_step(void) ;

	// locate src and rec in the mesh
	virtual Rtn_code locate_src_and_rec_in_mesh(Acquisition*) ;

	// compute tmin for dynamic front
	Rtn_code compute_tmin(Acquisition*) ;

	// locate snapshot pixel in the mesh
	virtual Rtn_code locate_pixel_in_mesh(Snapshot*) ;

	// write time solution at receivers
	virtual Rtn_code write_trace(Variable*, Myint) ;

	// compute error for eigen mode and write on disk
	Rtn_code compute_eigen_error(Myfloat*, Myint, Acquisition*) ;

	// interpolate value within element
	Myfloat interpolate_variable(Myfloat* pVal, Myint ielem, Myfloat eta_coord, Myfloat xi_coord) ;

	// compute element matrices
	virtual Rtn_code compute_element_matrices(void) ;

	// compute local mass matrix
	Rtn_code compute_mass_matrix(VecDoub &xgl, VecDoub& wgl, VecDoub& xnode, MatDoub &Mat_M) ;

	// compute global inverse mass matrix
	Rtn_code compute_global_inv_mass_matrix(void) ;

	// compute global inverse mass matrix
	Rtn_code compute_global_inv_mass_matrix( Model*, Var_type, Coef_type) ;

	// compute stiffness matrix
	Rtn_code compute_stiffness_matrix(VecDoub &xgl, VecDoub& wgl, VecDoub& xnode, MatDoub &Mat_M, Axis_type axis) ;

	// compute flux matrix
	virtual Rtn_code compute_flux_matrices(Myint imat, ofstream *out_file, Myint ni, VecDoub &xgl, VecDoub& wgl, VecDoub& xnode) = 0 ;

	// stiffness matrices
	//===================
	// pMat_Dz[ii] is a pointer to the matrix (Grid_2D_float) to be used
	// for an element with polynomial order ii
	// example:
	// Myfloat** const pMatK = pMat_Dz[ni-1]->pArray ;
	Grid_2D_float* pMat_Dz[MAX_POLY_ORDER+1] ;
	Grid_2D_float* pMat_Dx[MAX_POLY_ORDER+1] ;

	// matrices are stored in sparse format (i, j, c) i=line, j=column, c=coefficient
	// pMat_F_i[ineigh][ii] is a pointer to the matrix (Grid_1D_float) to be used for:
	// - an element with polynomial order ii
	//
	// example:
	// Myint* const pMat_i = pMat_Dz_i[ineigh][ni-1]->pArray ;
	// Myint* const pMat_j = pMat_Dz_j[ineigh][ni-1]->pArray ;
	// Myfloat* const pMat_c = pMat_Dz_c[ineigh][ni-1]->pArray ;

	Grid_1D_int*   pMat_Dz_i[MAX_POLY_ORDER+1] ;
	Grid_1D_int*   pMat_Dz_j[MAX_POLY_ORDER+1] ;
	Grid_1D_float* pMat_Dz_c[MAX_POLY_ORDER+1] ;

	Grid_1D_int*   pMat_Dx_i[MAX_POLY_ORDER+1] ;
	Grid_1D_int*   pMat_Dx_j[MAX_POLY_ORDER+1] ;
	Grid_1D_float* pMat_Dx_c[MAX_POLY_ORDER+1] ;

} ;

} // namespace django

#endif
