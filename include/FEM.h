#ifndef DJANGO_FEM_H_
#define DJANGO_FEM_H_

#include "scheme.h"

// numerical recipes
#include "nr3.h"

#include "grid_1D_float.h"
#include "grid_1D_int.h"
#include "grid_2D_float.h"
#include "model.h"
#include "snapshot.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------
// Constants definitions
//------------------------------------------------------------------------------------

const Myint MAX_POLY_ORDER = 9 ;
const Myfloat CFL_DGSE[MAX_POLY_ORDER+1] = {1.98, 0.97, 0.76, 0.77, 0.78, 0.8, 0.8, 0.8, 0.8, 0.8} ;
const Myfloat CFL_CGSE[MAX_POLY_ORDER+1] = {999, 1.84, 1.25, 1.27, 1.28, 1.36, 1.44, 1.44, 1.44, 1.44} ;

const Myint MAX_INV_MATRIX  = 10 ;

// optimal ratio for 20T
/* const Myfloat P0_RATIO = 1000.0 ; */
/* const Myfloat P1_RATIO = 100.0 ; */
/* const Myfloat P2_RATIO = 9.0 ; */
/* const Myfloat P3_RATIO = 3.3 ; */
/* const Myfloat P4_RATIO = 1.8 ; */
/* const Myfloat P5_RATIO = 1.2 ; */
/* const Myfloat P6_RATIO = 0.9 ; */

// optimal ratio for 100T
const Myfloat P0_RATIO = 1000.0 ;
const Myfloat P1_RATIO = 100.0 ;
const Myfloat P2_RATIO = 20.0 ;
const Myfloat P3_RATIO = 4.0 ;
const Myfloat P4_RATIO = 2.0 ;
const Myfloat P5_RATIO = 1.3 ;
const Myfloat P6_RATIO = 1.2 ;

//------------------------------------------------------------------------------------

class FEM: public Scheme

{
public:

	FEM(void) ;

	// initialize
	virtual Rtn_code initialize(void) ;

	// finalize
	virtual Rtn_code finalize(void) ;

	// display mesh info
	virtual Rtn_code mesh_info(void) ;

	// print FEM info
	virtual Rtn_code info(void) ;

	// write mesh in VTK format
	virtual Rtn_code write_mesh_VTK(void) = 0 ;

	// write node in VTK format
	virtual Rtn_code write_node_VTK(void) = 0 ;

protected:

	// total # elements
	Myint nelem ;

	// # elements in medium
	Myint nelem_med ;

	// # elements in layer
	Myint nelem_lay ;

	// # ghost elements
	Myint nelem_ghost ;

	// # acoustic iso elements
	Myint nelem_ac_iso ;

	// # acoustic lossy elements
	Myint nelem_ac_lossy ;

	// # elastic iso elements
	Myint nelem_el_iso ;

	// total # nodes
	Myint nnode ;

	// extreme mesh size
	Myfloat min_z_mesh, max_z_mesh ;
	Myfloat min_x_mesh, max_x_mesh ;

	// extreme distance between node
	Myfloat min_dist_node, max_dist_node ;

	// **** acquisition ****
	// number of nodes for the source excitation
	Myint  nnode_src ;
	// position of source in the grid (node index)
	Myint* inode_src ;
	// weight
	Myfloat* wnode_src ;

	// position of receivers in the grid (element and coordinates)
	Myint* ielem_rec ;
	Myfloat* eta_rec ;
	Myfloat* xi_rec ;

	// position of pixel of snapshot in the grid (element and coordinates)
	Myint* ielem_pixel ;
	Myfloat* eta_pixel ;
	Myfloat* xi_pixel ;

	// node type
	Node_type node_type ;

	// node distribution
	Node_distrib node_distrib ;

	// node integration
	Node_integ node_integ ;

	// flux type
	Flux_type flux_type ;

	// min. polynomial
	Myint pmin ;

	// max polynomial
	Myint pmax ;

	// no. element along axis in medium
	Myint nelem_x_med, nelem_z_med ;

	// no. element along axis in complete mesh
	Myint nelem_x, nelem_z ;

	// physical properties within elements
	Prop_type prop ;

	// source excitation
	Rtn_code source_excitation(Myfloat*, Myint, Wavefield_type, Data*, Myfloat) ;

	// mesh initialization
	virtual Rtn_code initialize_mesh(Model*, Myfloat fmax) = 0 ;

	// check time step
	virtual Myfloat compute_optimal_time_step(void) = 0 ;

	// locate src and rec in the mesh
	virtual Rtn_code locate_src_and_rec_in_mesh(Acquisition*) = 0 ;

	// locate snapshot pixel in the mesh
	virtual Rtn_code locate_pixel_in_mesh(Snapshot*) = 0 ;

	// write time solution at receivers
	virtual Rtn_code write_trace(Variable*, Myint) = 0 ;

	// free src and rec in the grid
	Rtn_code free_position_arrays(void) ;

	// compute element matrices
	virtual Rtn_code compute_element_matrices(void) = 0 ;

	// initialize node coordinates in ref element
	Rtn_code initialize_node_ref(void) ;

	// reset scheme from another scheme
	// used for dynamic p-adaptivity
	virtual Rtn_code reset_scheme_from(FEM* schemeIn) ;

	// compute Gauss Lobatto Legendre  quadrature points (weight and coordinates) in ref element
#include "gll.h"

	// compute equidistant point coordinates in ref element
#include "eqd.h"  

	// compute inverse of matrix
#include "gaussj.h"

	// compute Gauss Legendre quadrature points (weight and coordinates) in ref element
#include "gamma.h"
#include "gauss_wgts.h"

	// compute Lagrange polynomials and their derivative
#include "lagrange.h"

	// vector and matrix utility functions
#include "vecmat.h"

	// node coordinates in the 1D ref element
	VecDoub xnode_ref_elem_1D[MAX_POLY_ORDER+1] ;

	// inverse of mass matrix
	//========================
	// pMat_M_inv[ii] is a pointer to the matrix (Grid_2D_float) to be used
	// for an element with polynomial order ii
	// example:
	// Myfloat** const pMatM = pMat_M_inv[ni-1]->pArray ;
	Grid_2D_float* pMat_M_inv[MAX_POLY_ORDER+1] ;

	// global vector
	Grid_1D_float* pVec_k_glob ;

	// global mass matrix
	Grid_1D_float* pMat_M_inv_glob ;

	// global mass matrix embedding parameters
	Grid_1D_float* pMat_M_inv_glob_param[MAX_INV_MATRIX] ;

	// global vp vector
	Grid_1D_float* pVec_vp_glob ;

	// global vs vector
	Grid_1D_float* pVec_vs_glob ;

	// global rho vector
	Grid_1D_float* pVec_rho_glob ;

	// absorbing sponge coef
	Grid_1D_float *pSponge_coef_node ;

	// sigma_l (lossy term)
	Myfloat sigma_l ;

} ;

} // namespace django

#endif
