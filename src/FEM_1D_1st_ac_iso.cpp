//-------------------------------------------------------------------------------------------------------
//
// MODELLING WITH FEM IN 1D (CONTINUOUS/DISCONTINUOUS GALERKIN)
//
//  * ACOUSTIC ISOTROPIC MEDIA 
//  * 1ST ORDER WAVE EQUATION (PRESSURE / VELOCITY):
//
//            dvz/dt = coef2 dP/dz
//            dP /dt = coef1 dvz/dz
//
//            with
//
//            coef1 = kappa = (rho * vp * vp)
//            coef2 = 1/rho
//            where rho is CONSTANT
//
// FEM FORMULATION
// ===============
//
// pr = pr + dt * Mat_M_inv * vec_k 
// with vec_k += c1 { -Mat_Dz * vz(iel) + 0.5 * SUM_neigh [vnz * (Mat_F * vz(iel) + Mat_G * vz(iel_neigh)) ] }
//
// vz = vz + dt * Mat_M_inv * vec_k 
// with vec_k += c2 { -Mat_Dz * pr(iel) + 0.5 * SUM_neigh [vnz * (Mat_F * pr(iel) + Mat_G * pr(iel_neigh)) ] }
//
// PARENT CLASS: Scheme
//   DERIVED CLASS: FEM
//     DERIVED CLASS: FEM_1D
//       DERIVED CLASS: FEM_1D_1st
//         DERIVED CLASS: FEM_1D_1st_ac_iso
//
//-------------------------------------------------------------------------------------------------------

#include "FEM_1D_1st_ac_iso.h"

#include <cassert>
#include <fstream>
#include <iostream>

#include "mpi.h"

#include "grid_1D_float.h"
#include "grid_2D_complex.h"
#include "output_report.h"
#include "singleton.h"
#include "type_def.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------
FEM_1D_1st_ac_iso::FEM_1D_1st_ac_iso() : FEM_1D_1st()

{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D_1st_ac_iso::FEM_1D_1st_ac_iso");

	eq_type = ACOUSTIC ;
	CFL = 1.0 ;
	print_debug(ALL, MID_DEBUG, "CFL=", CFL) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D_1st_ac_iso::FEM_1D_1st_ac_iso");
}

//=======================================================================================================
//
// INITIALIZE MODELLING
//
//=======================================================================================================

Rtn_code FEM_1D_1st_ac_iso::initialize(Model* pModel, Myfloat fmax)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D_1st_ac_iso::initialize");

	// check space interpolation and time order
	Myint polynom_order = Singleton::Instance()->space_order ;
	if (polynom_order < 0)
	{
		print_error(" Can not initialize FEM_1D_1st_ac_iso, space_order < 0", polynom_order) ;
		return(RTN_CODE_KO) ;
	}
	Myint time_order = Singleton::Instance()->time_order ;
	if (time_order != 2)
	{
		print_error(" Can not initialize FEM_1D_1st_ac_iso, time_order != 2", time_order) ;
		return(RTN_CODE_KO) ;
	}

	// call parent initialiation
	Rtn_code rtn_code = FEM_1D_1st::initialize(pModel, fmax) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	//-------------------------------------------------------------------------------------------------------
	// unknown vectors
	//-------------------------------------------------------------------------------------------------------

	// velocity component (time domain)
	Variable* pVarVz = Singleton::Instance()->register_variable(VZ, "vz") ;
	pVarVz->allocate_grid(nnode, 0) ;

	// pressure component (time domain)
	Variable* pVarPr = Singleton::Instance()->register_variable(PR, "pr") ;
	pVarPr->allocate_grid(nnode, 0) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D_1st_ac_iso::initialize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_1D_1st_ac_iso::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D_1st_ac_iso::finalize");

	// delete components
	Singleton::Instance()->delete_variable(VZ) ;
	Singleton::Instance()->delete_variable(PR) ;

	// call parent finalize
	Rtn_code rtn_code = FEM_1D_1st::finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D_1st_ac_iso::finalize");
	return(RTN_CODE_OK) ;
}

//=======================================================================================================
//
// PERFORM FORWARD MODELLING FOR ONE SOURCE
//
//=======================================================================================================

Rtn_code FEM_1D_1st_ac_iso::solve_current_shot(Acquisition* pAcquisition, Data *pData, Wavefield_type wtype, Snapshot* pSnapshot)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D_1st_ac_iso::solve_current_shot");

	// get variables
	Variable* pr_var = Singleton::Instance()->get_variable(PR) ;
	Grid_1D_float *pr_grid = (Grid_1D_float*) pr_var->get_grid() ;
	Myfloat * const pr = pr_grid->pArray ;

	Variable* vz_var = Singleton::Instance()->get_variable(VZ) ;
	Grid_1D_float *vz_grid = (Grid_1D_float*) vz_var->get_grid() ;
	Myfloat * const vz = vz_grid->pArray ;

	// reset memory
	Rtn_code rtn_code = reset() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// locate src and rec in the grid
	rtn_code = locate_src_and_rec_in_mesh(pAcquisition) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// source scaling factor
	Myfloat src_factor = dt ;

	if (Singleton::Instance()->dryrun)
	{
		print_info(ALL, "\n*** DRY RUN - SKIP TIME LOOP ***") ;
	}
	else
	{

		//-------------------------------------------------------------------------------------------------------
		// start loop over time steps
		//-------------------------------------------------------------------------------------------------------

		for (Myint it = 0; it < nt; it ++)
		{
			//-------------------------------------------------------------------------------------------------------
			// computation of the pressure component
			//-------------------------------------------------------------------------------------------------------
			double t0 = MPI_Wtime() ;
			compute_pressure(pr, vz) ;
			time_in_kernel += MPI_Wtime() - t0 ;

			//-------------------------------------------------------------------------------------------------------
			// source excitation
			//-------------------------------------------------------------------------------------------------------
			if (src_type != NO_SRC_TYPE)
			{
				if (src_stype == EXPLOSIVE)
				{
					rtn_code = source_excitation(pr, it, wtype, pData, src_factor) ;
				}
				else if (src_stype == FORCE_Z)
				{
					rtn_code = source_excitation(vz, it, wtype, pData, src_factor) ;
				}
			}
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

			// compute energy
			//compute_energy(pr) ;

			// compute error for eigen mode
			if(Singleton::Instance()->pProgram->pModelling->get_case() == MODELLING_EIGEN)
			{
				rtn_code = compute_eigen_error(pr, it, pAcquisition) ;
				if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
			}

			//-------------------------------------------------------------------------------------------------------
			// computation of the velocity component
			//-------------------------------------------------------------------------------------------------------
			t0 = MPI_Wtime() ;
			compute_velocity(pr, vz) ;
			time_in_kernel += MPI_Wtime() - t0 ;

			//-------------------------------------------------------------------------------------------------------
			// write seismograms
			//-------------------------------------------------------------------------------------------------------
			rtn_code = write_trace(pr_var, it) ;
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
			rtn_code = write_trace(vz_var, it) ;
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

			//-------------------------------------------------------------------------------------------------------
			// write snapshot
			//-------------------------------------------------------------------------------------------------------
			rtn_code = write_snapshot(pr_var, it) ;
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
			rtn_code = write_snapshot(vz_var, it) ;
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

			// display remaining computation time
			display_remaining_computation_time(it, time_in_kernel) ;

		}

		//-------------------------------------------------------------------------------------------------------
		// end loop over time steps
		//-------------------------------------------------------------------------------------------------------

	}

	// free src and rec in the grid
	rtn_code = free_position_arrays() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D_1st_ac_iso::solve_current_shot");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_1D_1st_ac_iso::reset(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D_1st_ac_iso::reset");

	// call to parent class
	Rtn_code rtn_code = Scheme::reset() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// initialization for eigen mode test
	if(Singleton::Instance()->pProgram->pModelling->get_case() == MODELLING_EIGEN)
	{
		// initialization for eigen mode test
		init_eigen() ;
	}
	else
	{
		// reset grids
		rtn_code = Singleton::Instance()->get_variable(PR)->reset_grid() ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
		rtn_code = Singleton::Instance()->get_variable(VZ)->reset_grid() ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D_1st_ac_iso::reset");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_1D_1st_ac_iso::compute_pressure(Myfloat* pr, Myfloat* vz)
{
	print_debug(ALL, FULL_DEBUG, "IN FEM_1D_1st_ac_iso::compute_pressure");

	// update ncell
	ncell += nnode ;

	// initialize global vector
	//=========================

	Myfloat* const vec_k = pVec_k_glob->pArray ;
#pragma ivdep
	for (Myint inode = 0; inode < nnode; inode++)
	{
		vec_k[inode] = Myfloat(0.0) ;
	}

	// loop on elements
	//=================

	for (Myint iel = 0; iel < nelem; iel++)
	{
		// current element
		Myint const ni = pElement[iel].nnode ;
		Myint const mat_idx_i = pElement[iel].nnode_1D - 1 ;

		// compute -Mat_Dz * vz
		//=====================

		Myfloat** const Mat_Dz = pMat_Dz[mat_idx_i]->pArray ;
		for (Myint ii = 0; ii < ni; ii++) {
			Myint inode_ii = pElement[iel].inode[ii] ;
			for (Myint jj = 0; jj < ni; jj++) {
				vec_k[inode_ii] -= Mat_Dz[ii][jj] * vz[pElement[iel].inode[jj]] ;
			}
		}
		nb_op_kernel += ni * ni * 2 ;

		// add flux contributions
		// 0.5 * SUM_neigh [vnz * Mat_F * vz(iel) + vnz * Mat_G * vz(iel_neigh)]
		//======================================================================

		// loop on neighbours
		// I_ZPREV, I_ZNEXT
		for (Myint ineigh = 0; ineigh < 2; ineigh++)
		{
			if (pElement[iel].flux[ineigh] == NO_FLUX_TYPE) continue ;

			// neighbour element index
			Myint iel_neigh = pElement[iel].neigh[ineigh] ;

			if (iel_neigh != NO_NEIGH)
			{
				Myint const mat_idx_j = pElement[iel_neigh].nnode_1D - 1 ;
				Myfloat const vnz = pElement[iel].vnz[ineigh] ;

				// 0.5 * vnz * Mat_F * vz(iel)
				Myint    const n1      = pMat_F_i[ineigh][mat_idx_i]->nz ;
				Myint*   const Mat_F_i = pMat_F_i[ineigh][mat_idx_i]->pArray ;
				Myint*   const Mat_F_j = pMat_F_j[ineigh][mat_idx_i]->pArray ;
				Myfloat* const Mat_F_c = pMat_F_c[ineigh][mat_idx_i]->pArray ;
				for (Myint ii = 0; ii < n1; ii++) {
					vec_k[ pElement[iel].inode[ Mat_F_i[ii] ] ] += vnz * Myfloat(0.5) * Mat_F_c[ii] * vz[ pElement[iel].inode[ Mat_F_j[ii] ] ] ;
				}

				// 0.5 * vnz * Mat_G * vz(iel_neigh)
				Myint    const n2 = pMat_G_i[ineigh][mat_idx_i][mat_idx_j]->nz ;
				Myint*   const Mat_G_i = pMat_G_i[ineigh][mat_idx_i][mat_idx_j]->pArray ;
				Myint*   const Mat_G_j = pMat_G_j[ineigh][mat_idx_i][mat_idx_j]->pArray ;
				Myfloat* const Mat_G_c = pMat_G_c[ineigh][mat_idx_i][mat_idx_j]->pArray ;
				for (Myint ii = 0; ii < n2; ii++) {
					vec_k[ pElement[iel].inode[ Mat_G_i[ii] ] ] += vnz * Myfloat(0.5) * Mat_G_c[ii] * vz[ pElement[iel_neigh].inode[ Mat_G_j[ii] ] ] ;
				}

				nb_op_kernel += (n1+n2)*4 ;

			} // if (iel_neigh != NO_NEIGH)
		} // for (Myint ineigh = 0; ineigh < 2; ineigh++)

	} // for (Myint iel = 0; iel < nelem; iel++)

	// update pressure
	//================

	Myfloat* const Vec_vp_glob    = pVec_vp_glob->pArray ;
	Myfloat* const Mat_M_inv_glob = pMat_M_inv_glob->pArray ;
#pragma ivdep
	for (Myint inode = 0; inode < nnode; inode++)
	{
		pr[inode] += dt * Mat_M_inv_glob[inode] * vec_k[inode] * TEMP_RHO_CONST * Vec_vp_glob[inode] * Vec_vp_glob[inode] ;
	}
	nb_op_kernel += nnode * 6 ;

	// boundaries
	//===========
	Myfloat* sponge_coef_node ;
	if (is_there_boundary_type(SPG)) sponge_coef_node = pSponge_coef_node->pArray ;
#pragma ivdep
	for (Myint inode = 0; inode < nnode; inode++)
	{
		if (pNode[inode].boundary == FREESURF)
		{
			// free surface -> pressure is null
			pr[inode] = Myfloat(0.0) ;
		}
		else if (pNode[inode].boundary == SPG)
		{
			// sponge damping
			pr[inode] *= sponge_coef_node[inode] ;
		}
	}

	print_debug(ALL, FULL_DEBUG, "OUT FEM_1D_1st_ac_iso::compute_pressure");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_1D_1st_ac_iso::compute_velocity(Myfloat* pr, Myfloat* vz)
{
	print_debug(ALL, FULL_DEBUG, "IN FEM_1D_1st_ac_iso::compute_velocity");

	// c2 = 1/rho
	Myfloat c2 = 1.0 / TEMP_RHO_CONST ;

	// initialize global vector
	//=========================

	Myfloat* const vec_k = pVec_k_glob->pArray ;
#pragma ivdep
	for (Myint inode = 0; inode < nnode; inode++)
	{
		vec_k[inode] = Myfloat(0.0) ;
	}

	// loop on elements
	//=================

	for (Myint iel = 0; iel < nelem; iel++)
	{
		// current element
		Myint const ni = pElement[iel].nnode ;
		Myint const mat_idx_i = pElement[iel].nnode_1D - 1 ;

		// compute -Mat_Dz * pr
		//=====================

		Myfloat** const Mat_Dz = pMat_Dz[mat_idx_i]->pArray ;
		for (Myint ii = 0; ii < ni; ii++) {
			Myint inode_ii = pElement[iel].inode[ii] ;
			for (Myint jj = 0; jj < ni; jj++) {
				vec_k[inode_ii] -= Mat_Dz[ii][jj] * pr[pElement[iel].inode[jj]] ;
			}
		}
		nb_op_kernel += ni * ni * 2 ;

		// add flux contributions
		//=======================

		// loop on neighbours
		// I_ZPREV, I_ZNEXT
		for (Myint ineigh = 0; ineigh < 2; ineigh++)
		{
			if (pElement[iel].flux[ineigh] == NO_FLUX_TYPE) continue ;

			// neighbour element index
			Myint iel_neigh = pElement[iel].neigh[ineigh] ;

			if (iel_neigh != NO_NEIGH)
			{
				Myint const mat_idx_j = pElement[iel_neigh].nnode_1D - 1 ;
				Myfloat const vnz = pElement[iel].vnz[ineigh] ;

				// 0.5 * vnz * Mat_F * pr(iel)
				Myint    const n1      = pMat_F_i[ineigh][mat_idx_i]->nz ;
				Myint*   const Mat_F_i = pMat_F_i[ineigh][mat_idx_i]->pArray ;
				Myint*   const Mat_F_j = pMat_F_j[ineigh][mat_idx_i]->pArray ;
				Myfloat* const Mat_F_c = pMat_F_c[ineigh][mat_idx_i]->pArray ;
				for (Myint ii = 0; ii < n1; ii++) {
					vec_k[ pElement[iel].inode[ Mat_F_i[ii] ] ] += vnz * Myfloat(0.5) * Mat_F_c[ii] * pr[ pElement[iel].inode[ Mat_F_j[ii] ] ] ;
				}

				// 0.5 * vnz * Mat_G * pr(iel_neigh)
				Myint    const n2 = pMat_G_i[ineigh][mat_idx_i][mat_idx_j]->nz ;
				Myint*   const Mat_G_i = pMat_G_i[ineigh][mat_idx_i][mat_idx_j]->pArray ;
				Myint*   const Mat_G_j = pMat_G_j[ineigh][mat_idx_i][mat_idx_j]->pArray ;
				Myfloat* const Mat_G_c = pMat_G_c[ineigh][mat_idx_i][mat_idx_j]->pArray ;
				for (Myint ii = 0; ii < n2; ii++) {
					vec_k[ pElement[iel].inode[ Mat_G_i[ii] ] ] += vnz * Myfloat(0.5) * Mat_G_c[ii] * pr[ pElement[iel_neigh].inode[ Mat_G_j[ii] ] ] ;
				}

				nb_op_kernel += (n1+n2)*4 ;

			} // if (iel_neigh != NO_NEIGH)
		} // for (Myint ineigh = 0; ineigh < 2; ineigh++)
	} // for (Myint iel = 0; iel < nelem; iel++)

	// update velocity
	//================

	Myfloat* const Mat_M_inv_glob = pMat_M_inv_glob->pArray ;
#pragma ivdep
	for (Myint inode = 0; inode < nnode; inode++)
	{
		vz[inode] += dt * Mat_M_inv_glob[inode] * c2 * vec_k[inode] ;
	}
	nb_op_kernel += nnode * 4 ;

	// boundaries
	//===========
	Myfloat* sponge_coef_node ;
	if (is_there_boundary_type(SPG)) sponge_coef_node = pSponge_coef_node->pArray ;
	for (Myint inode = 0; inode < nnode; inode++)
	{
		if (pNode[inode].boundary == RIGID)
		{
			// rigid surface -> velocity is null
			vz[inode] = 0.0 ;
		}
		else if (pNode[inode].boundary == SPG)
		{
			// sponge damping
			vz[inode] *= sponge_coef_node[inode] ;
			nb_op_kernel++ ;
		}
	}

	print_debug(ALL, FULL_DEBUG, "OUT FEM_1D_1st_ac_iso::compute_velocity");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
void FEM_1D_1st_ac_iso::init_eigen()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D_1st_ac_iso::init_eigen");

	Myfloat time_vz = -dt/2.0 ;
	Myfloat time_pr = -dt ;

	Grid_1D_float *pr_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(PR))->get_grid() ;
	Myfloat * const pr = pr_grid->pArray ;

	Grid_1D_float *vz_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(VZ))->get_grid() ;
	Myfloat * const vz = vz_grid->pArray ;

	// get number mode number
	eigen_nmode = (Myint) Singleton::Instance()->pProgram->pModelling->get_param() ;

	// loop on elements
	for (Myint iel = 0; iel < nelem; iel++)
	{
		// loop on nodes
		for (Myint ii=0; ii< pElement[iel].nnode ; ii++)
		{
			Myint inode = pElement[iel].inode[ii] ;
			Myfloat znode = pNode[inode].zcoord ;
			pr[inode] = -sin(M_PI*znode * eigen_nmode) * sin(M_PI*time_pr * eigen_nmode) ;
			vz[inode] = cos(M_PI*znode * eigen_nmode) * cos(M_PI*time_vz * eigen_nmode) ;
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D_1st_ac_iso::init_eigen");
	return ;
}

} // namespace django