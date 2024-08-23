//-------------------------------------------------------------------------------------------------------
//
// MODELLING WITH FEM IN 2D (CONTINUOUS/DISCONTINUOUS GALERKIN)
//
//  * ACOUSTIC ISOTROPIC MEDIA 
//  * 1ST ORDER WAVE EQUATION (PRESSURE / VELOCITY):
//
//            coef2 dvx/dt = dP/dx
//            coef2 dvz/dt = dP/dz
//            coef1 dP /dt = (dvx/dx + dvz/dz) 
//
//            with
//
//            coef1 = 1/kappa = 1/(rho * vp * vp)
//            coef2 = rho
//
// FEM FORMULATION
// ==============
// pr = pr + dt * Mat_M_inv(coef1) * vec_k 
// with vec_k = c1 { -Mat_Dz * vz(iel) +  -Mat_Dx * vx(iel)
//                   + 0.5 * SUM_neigh [vnz * (Mat_F * vz(iel) + Mat_G * vz(iel_neigh))
//                                     +vnx * (Mat_F * vx(iel) + Mat_G * vx(iel_neigh)) ] }
//
// vz = vz + dt * Mat_M_inv(coef2) * vec_k 
// with vec_k = c2 * { -Mat_Dz * pr(iel)
//                   + 0.5 * SUM_neigh [vnz * (Mat_F * pr(iel) + Mat_G * pr(iel_neigh)) ] }
//
// vx = vx + dt * Mat_M_inv(coef2) * vec_k 
// with vec_k = c2 * { -Mat_Dx * pr(iel)
//                   + 0.5 * SUM_neigh [vnx * (Mat_F * pr(iel) + Mat_G * pr(iel_neigh)) ] }
//
// PARENT CLASS: Scheme
//   DERIVED CLASS: FEM
//     DERIVED CLASS: FEM_2D
//       DERIVED CLASS: FEM_2D_1s
//         DERIVED CLASS: FEM_2D_1st_ac_iso
//
//-------------------------------------------------------------------------------------------------------

#include "FEM_2D_1st_ac_iso.h"

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
FEM_2D_1st_ac_iso::FEM_2D_1st_ac_iso(void) : FEM_2D_1st()

{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st_ac_iso::FEM_2D_1st_ac_iso");

	eq_type = ACOUSTIC ;
	CFL = 1.0 ;
	print_debug(ALL, MID_DEBUG, "CFL=", CFL) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st_ac_iso::FEM_2D_1st_ac_iso");
}

//=======================================================================================================
//
// INITIALIZE MODELLING
//
//=======================================================================================================

Rtn_code FEM_2D_1st_ac_iso::initialize(Model* pModel, Myfloat fmax)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st_ac_iso::initialize");

	// check space interpolation and time order
	Myint polynom_order = Singleton::Instance()->space_order ;
	if (polynom_order < 0)
	{
		print_error(" Can not initialize FEM_2D_1st_ac_iso, space_order < 0", polynom_order) ;
		return(RTN_CODE_KO) ;
	}
	Myint time_order = Singleton::Instance()->time_order ;
	if (time_order != 2)
	{
		print_error(" Can not initialize FEM_2D_1st_ac_iso, time_order != 2", time_order) ;
		return(RTN_CODE_KO) ;
	}

	// call parent initialiation
	Rtn_code rtn_code = FEM_2D_1st::initialize(pModel, fmax) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// initialize global invert mass matrix
	//===============================================================

	// coef 1 = 1 / (rho * vp^2)
	rtn_code = compute_global_inv_mass_matrix(pModel, INV_RHOVP2, COEF1) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// coef 2 = rho
	rtn_code = compute_global_inv_mass_matrix(pModel, RHO, COEF2) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	//-------------------------------------------------------------------------------------------------------
	// unknown vectors
	//-------------------------------------------------------------------------------------------------------

	// velocity component (time domain)
	Variable* pVarVz = Singleton::Instance()->register_variable(VZ, "vz") ;
	pVarVz->allocate_grid(nnode, 0) ;
	Variable* pVarVx = Singleton::Instance()->register_variable(VX, "vx") ;
	pVarVx->allocate_grid(nnode, 0) ;

	// pressure component (time domain)
	Variable* pVarPr = Singleton::Instance()->register_variable(PR, "pr") ;
	pVarPr->allocate_grid(nnode, 0) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st_ac_iso::initialize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_2D_1st_ac_iso::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st_ac_iso::finalize");

	// delete components
	Singleton::Instance()->delete_variable(VZ) ;
	Singleton::Instance()->delete_variable(VX) ;
	Singleton::Instance()->delete_variable(PR) ;

	// call parent finalize
	Rtn_code rtn_code = FEM_2D_1st::finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st_ac_iso::finalize");
	return(RTN_CODE_OK) ;
}

//=======================================================================================================
//
// PERFORM FORWARD MODELLING FOR ONE SOURCE
//
//=======================================================================================================

Rtn_code FEM_2D_1st_ac_iso::solve_current_shot(Acquisition* pAcquisition, Data *pData, Wavefield_type wtype, Snapshot* pSnapshot)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st_ac_iso::solve_current_shot");

	// get variables
	Variable* pr_var = Singleton::Instance()->get_variable(PR) ;
	Grid_1D_float *pr_grid = (Grid_1D_float*) pr_var->get_grid() ;
	Myfloat * const pr = pr_grid->pArray ;

	Variable* vz_var = Singleton::Instance()->get_variable(VZ) ;
	Grid_1D_float *vz_grid = (Grid_1D_float*) vz_var->get_grid() ;
	Myfloat * const vz = vz_grid->pArray ;

	Variable* vx_var = Singleton::Instance()->get_variable(VX) ;
	Grid_1D_float *vx_grid = (Grid_1D_float*) vx_var->get_grid() ;
	Myfloat * const vx = vx_grid->pArray ;

	// reset memory
	Rtn_code rtn_code = reset() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// locate src and rec in the grid
	rtn_code = locate_src_and_rec_in_mesh(pAcquisition) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// compute tmin for dynamic front
	rtn_code = compute_tmin(pAcquisition) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// save mesh in VTK format
	write_mesh_VTK() ;

	// locate snapshot pixel in the mesh
	rtn_code = this->locate_pixel_in_mesh(pSnapshot) ;
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

		for (it = 0; it < nt; it ++)
		{
			//-------------------------------------------------------------------------------------------------------
			// computation of the pressure component
			//-------------------------------------------------------------------------------------------------------
			double t0 = MPI_Wtime() ;
			compute_pressure(pr, vz, vx) ;
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
			//compute_energy() ;

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
			compute_velocity(pr, vz, vx) ;
			time_in_kernel += MPI_Wtime() - t0 ;

			//-------------------------------------------------------------------------------------------------------
			// write seismograms
			//-------------------------------------------------------------------------------------------------------
			rtn_code = write_trace(pr_var, it) ;
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
			rtn_code = write_trace(vz_var, it) ;
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
			rtn_code = write_trace(vx_var, it) ;
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

			//-------------------------------------------------------------------------------------------------------
			// write snapshot
			//-------------------------------------------------------------------------------------------------------
			rtn_code = write_snapshot(pr_var, it) ;
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
			rtn_code = write_snapshot(vz_var, it) ;
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
			rtn_code = write_snapshot(vx_var, it) ;
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

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st_ac_iso::solve_current_shot");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_2D_1st_ac_iso::reset(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st_ac_iso::reset");

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
		rtn_code = Singleton::Instance()->get_variable(VX)->reset_grid() ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st_ac_iso::reset");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_2D_1st_ac_iso::compute_pressure(Myfloat* pr, Myfloat* vz, Myfloat* vx)
{
	print_debug(ALL, FULL_DEBUG, "IN FEM_2D_1st_ac_iso::compute_pressure");

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

	Myfloat current_time = it * dt ;

	for (Myint iel = 0; iel < nelem; iel++)
	{
		// dynamic front
		if (current_time <  pElement[iel].tmin) continue ;

		// current element
		const Myint ni = pElement[iel].nnode ;
		const Myint mat_idx_i = pElement[iel].nnode_1D - 1 ;

		// compute -Mat_Dz * vz
		//=====================

		const Myfloat cvol = pElement[iel].size_x / 2.0 ;

		// sparse matrix
		Myint    const n1      = pMat_Dz_i[mat_idx_i]->nz ;
		Myint*   const Mat_Dz_i = pMat_Dz_i[mat_idx_i]->pArray ;
		Myint*   const Mat_Dz_j = pMat_Dz_j[mat_idx_i]->pArray ;
		Myfloat* const Mat_Dz_c = pMat_Dz_c[mat_idx_i]->pArray ;
		for (Myint ii = 0; ii < n1; ii++) {
			vec_k[ pElement[iel].inode[Mat_Dz_i[ii]] ] -= cvol * Mat_Dz_c[ii] * vz[ pElement[iel].inode[Mat_Dz_j[ii]] ] ;
		}
		nb_op_kernel += n1 * 3 ;

		// compute -Mat_Dx * vx
		//=====================

		const Myfloat cvol2 = pElement[iel].size_z / 2.0 ;

		// sparse matrix
		Myint    const n2      = pMat_Dx_i[mat_idx_i]->nz ;
		Myint*   const Mat_Dx_i = pMat_Dx_i[mat_idx_i]->pArray ;
		Myint*   const Mat_Dx_j = pMat_Dx_j[mat_idx_i]->pArray ;
		Myfloat* const Mat_Dx_c = pMat_Dx_c[mat_idx_i]->pArray ;
		for (Myint ii = 0; ii < n2; ii++) {
			vec_k[ pElement[iel].inode[Mat_Dx_i[ii]] ] -= cvol2 * Mat_Dx_c[ii] * vx[ pElement[iel].inode[Mat_Dx_j[ii]] ] ;
		}
		nb_op_kernel += n2 * 3 ;

		// add flux contributions
		//=======================

		// loop on neighbours
		// I_ZPREV, I_ZNEXT, I_XPREV, I_XNEXT
		for (Myint ineigh = 0; ineigh < 4; ineigh++)
		{
			if (pElement[iel].flux[ineigh] == NO_FLUX_TYPE) continue ;

			// neighbour element index
			Myint iel_neigh = pElement[iel].neigh[ineigh] ;

			if (iel_neigh != NO_NEIGH)
			{
				Myint const mat_idx_j = pElement[iel_neigh].nnode_1D - 1 ;
				Myfloat const vnz = pElement[iel].vnz[ineigh] * Myfloat(0.5) ;
				Myfloat const vnx = pElement[iel].vnx[ineigh] * Myfloat(0.5) ;

				// 0.5 * vnz * Mat_F * vz(iel)
				// 0.5 * vnx * Mat_F * vx(iel)
				Myint    const n1      = pMat_F_i[ineigh][mat_idx_i]->nz ;
				Myint*   const Mat_F_i = pMat_F_i[ineigh][mat_idx_i]->pArray ;
				Myint*   const Mat_F_j = pMat_F_j[ineigh][mat_idx_i]->pArray ;
				Myfloat* const Mat_F_c = pMat_F_c[ineigh][mat_idx_i]->pArray ;
				for (Myint ii = 0; ii < n1; ii++) {
					vec_k[ pElement[iel].inode[ Mat_F_i[ii] ] ] += vnz * Myfloat(0.5) * Mat_F_c[ii] * vz[ pElement[iel].inode[ Mat_F_j[ii] ] ] ;
					vec_k[ pElement[iel].inode[ Mat_F_i[ii] ] ] += vnx * Myfloat(0.5) * Mat_F_c[ii] * vx[ pElement[iel].inode[ Mat_F_j[ii] ] ] ;
				}

				// 0.5 * vnz * Mat_G * vz(iel_neigh)
				// 0.5 * vnx * Mat_G * vx(iel_neigh)
				Myint    const n2 = pMat_G_i[ineigh][mat_idx_i][mat_idx_j]->nz ;
				Myint*   const Mat_G_i = pMat_G_i[ineigh][mat_idx_i][mat_idx_j]->pArray ;
				Myint*   const Mat_G_j = pMat_G_j[ineigh][mat_idx_i][mat_idx_j]->pArray ;
				Myfloat* const Mat_G_c = pMat_G_c[ineigh][mat_idx_i][mat_idx_j]->pArray ;
				for (Myint ii = 0; ii < n2; ii++) {
					vec_k[ pElement[iel].inode[ Mat_G_i[ii] ] ] += vnz * Myfloat(0.5) * Mat_G_c[ii] * vz[ pElement[iel_neigh].inode[ Mat_G_j[ii] ] ] ;
					vec_k[ pElement[iel].inode[ Mat_G_i[ii] ] ] += vnx * Myfloat(0.5) * Mat_G_c[ii] * vx[ pElement[iel_neigh].inode[ Mat_G_j[ii] ] ] ;
				}

				nb_op_kernel += (n1 + n2) * 8 ;

			} // if (iel_neigh != NO_NEIGH)
		} // for (Myint ineigh = 0; ineigh < 4; ineigh++)
	} // for (Myint iel = 0; iel < nelem; iel++)

	// update pressure
	//================

	Myfloat* const Vec_vp_glob    = pVec_vp_glob->pArray ;
	Myfloat* const Mat_M_inv_glob = pMat_M_inv_glob_param[COEF1]->pArray ;
#pragma ivdep
	for (Myint inode = 0; inode < nnode; inode++)
	{
		pr[inode] += dt * Mat_M_inv_glob[inode] * vec_k[inode] ;
	}
	nb_op_kernel += nnode * 3 ;

	//============
	// boundaries
	//============
	Myfloat* sponge_coef_node ;
	if (is_there_boundary_type(SPG))
	{
		if (pSponge_coef_node == NULL)
		{
			print_error(" FEM_2D_1st_ac_iso::compute_pressure, pSponge_coef_node == NULL") ;
			return(RTN_CODE_KO) ;
		}
		else
		{
			sponge_coef_node = pSponge_coef_node->pArray ;
#pragma ivdep
			for (Myint inode = 0; inode < nnode; inode++)
			{
				if (pNode[inode].boundary == SPG)
				{
					// sponge damping
					pr[inode] *= sponge_coef_node[inode] ;
					nb_op_kernel++ ;
				}
			}
		}
	}

#pragma ivdep
	for (Myint inode = 0; inode < nnode; inode++)
	{
		if (pNode[inode].boundary == FREESURF)
		{
			// free surface -> pressure is null
			pr[inode] = Myfloat(0.0) ;
		}
	}

	print_debug(ALL, FULL_DEBUG, "OUT FEM_2D_1st_ac_iso::compute_pressure");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_2D_1st_ac_iso::compute_velocity(Myfloat* pr, Myfloat* vz, Myfloat* vx)
{
	print_debug(ALL, FULL_DEBUG, "IN FEM_2D_1st_ac_iso::compute_velocity");

	// compute Vz
	//===========
	{
		// initialize global vector
		//=========================

		Myfloat* const vec_k = pVec_k_glob->pArray ;
#pragma ivdep
		for (Myint inode = 0; inode < nnode; inode++)
		{
			vec_k[inode] = Myfloat(0.0) ;
		}

		Myfloat current_time = it * dt ;

		// loop on elements
		//=================

		for (Myint iel = 0; iel < nelem; iel++)
		{
			// dynamic front
			if (current_time <  pElement[iel].tmin) continue ;

			// current element
			Myint const ni = pElement[iel].nnode ;
			Myint const mat_idx_i = pElement[iel].nnode_1D - 1 ;

			// compute -Mat_Dz * pr
			//=====================

			const Myfloat cvol = pElement[iel].size_x / 2.0 ;

			// sparse matrix
			Myint    const n1      = pMat_Dz_i[mat_idx_i]->nz ;
			Myint*   const Mat_Dz_i = pMat_Dz_i[mat_idx_i]->pArray ;
			Myint*   const Mat_Dz_j = pMat_Dz_j[mat_idx_i]->pArray ;
			Myfloat* const Mat_Dz_c = pMat_Dz_c[mat_idx_i]->pArray ;
			for (Myint ii = 0; ii < n1; ii++) {
				vec_k[ pElement[iel].inode[Mat_Dz_i[ii]] ] -= cvol * Mat_Dz_c[ii] * pr [ pElement[iel].inode[Mat_Dz_j[ii]] ] ;
			}
			nb_op_kernel += n1 * 3 ;

			// add flux contributions
			//=======================

			// loop on neighbours
			// I_ZPREV, I_ZNEXT, I_XPREV, I_XNEXT
			for (Myint ineigh = 0; ineigh < 4; ineigh++)
			{
				if (pElement[iel].flux[ineigh] == NO_FLUX_TYPE) continue ;

				// neighbour element index
				Myint iel_neigh = pElement[iel].neigh[ineigh] ;

				if (iel_neigh != NO_NEIGH)
				{
					Myint const mat_idx_j = pElement[iel_neigh].nnode_1D - 1 ;
					Myfloat const vnz = pElement[iel].vnz[ineigh] * Myfloat(0.5) ;

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

					nb_op_kernel += (n1 + n2) * 4 ;

				} // if (iel_neigh != NO_NEIGH)
			} // for (Myint ineigh = 0; ineigh < 4; ineigh++)

		} // for (Myint iel = 0; iel < nelem; iel++)

		// update velocity
		//================

		Myfloat* const Mat_M_inv_glob = pMat_M_inv_glob_param[COEF2]->pArray ;
#pragma ivdep
		for (Myint inode = 0; inode < nnode; inode++)
		{
			vz[inode] += dt * Mat_M_inv_glob[inode] * vec_k[inode] ;
		}
		nb_op_kernel += nnode * 3 ;
	}

	// compute Vx
	//===========
	{
		// initialize global vector
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

			// compute -Mat_Dx * pr
			//=====================

			const Myfloat cvol2 = pElement[iel].size_z / 2.0 ;

			// sparse matrix
			Myint    const n1      = pMat_Dx_i[mat_idx_i]->nz ;
			Myint*   const Mat_Dx_i = pMat_Dx_i[mat_idx_i]->pArray ;
			Myint*   const Mat_Dx_j = pMat_Dx_j[mat_idx_i]->pArray ;
			Myfloat* const Mat_Dx_c = pMat_Dx_c[mat_idx_i]->pArray ;
			for (Myint ii = 0; ii < n1; ii++) {
				vec_k[ pElement[iel].inode[Mat_Dx_i[ii]] ] -= cvol2 * Mat_Dx_c[ii] * pr [ pElement[iel].inode[Mat_Dx_j[ii]] ] ;
			}
			nb_op_kernel += n1 * 3 ;

			// add flux contributions
			//=======================

			// loop on neighbours
			// I_ZPREV, I_ZNEXT, I_XPREV, I_XNEXT
			for (Myint ineigh = 0; ineigh < 4; ineigh++)
			{
				if (pElement[iel].flux[ineigh] == NO_FLUX_TYPE) continue ;

				// neighbour element index
				Myint iel_neigh = pElement[iel].neigh[ineigh] ;

				if (iel_neigh != NO_NEIGH)
				{
					Myint const mat_idx_j = pElement[iel_neigh].nnode_1D - 1 ;
					Myfloat const vnx = pElement[iel].vnx[ineigh] * Myfloat(0.5) ;

					// 0.5 * vnx * Mat_F * pr(iel)
					Myint    const n1      = pMat_F_i[ineigh][mat_idx_i]->nz ;
					Myint*   const Mat_F_i = pMat_F_i[ineigh][mat_idx_i]->pArray ;
					Myint*   const Mat_F_j = pMat_F_j[ineigh][mat_idx_i]->pArray ;
					Myfloat* const Mat_F_c = pMat_F_c[ineigh][mat_idx_i]->pArray ;
					for (Myint ii = 0; ii < n1; ii++) {
						vec_k[ pElement[iel].inode[ Mat_F_i[ii] ] ] += vnx * Myfloat(0.5) * Mat_F_c[ii] * pr[ pElement[iel].inode[ Mat_F_j[ii] ] ] ;
					}

					// 0.5 * vnx * Mat_G * pr(iel_neigh)
					Myint    const n2 = pMat_G_i[ineigh][mat_idx_i][mat_idx_j]->nz ;
					Myint*   const Mat_G_i = pMat_G_i[ineigh][mat_idx_i][mat_idx_j]->pArray ;
					Myint*   const Mat_G_j = pMat_G_j[ineigh][mat_idx_i][mat_idx_j]->pArray ;
					Myfloat* const Mat_G_c = pMat_G_c[ineigh][mat_idx_i][mat_idx_j]->pArray ;
					for (Myint ii = 0; ii < n2; ii++) {
						vec_k[ pElement[iel].inode[ Mat_G_i[ii] ] ] += vnx * Myfloat(0.5) * Mat_G_c[ii] * pr[ pElement[iel_neigh].inode[ Mat_G_j[ii] ] ] ;
					}

					nb_op_kernel += (n1 + n2) * 4 ;

				} // if (iel_neigh != NO_NEIGH)
			} // for (Myint ineigh = 0; ineigh < 4; ineigh++)
		} // for (Myint iel = 0; iel < nelem; iel++)

		// update velocity
		// vx = vx + dt/rho * Mat_M_inv * vec_k
		Myfloat* const Mat_M_inv_glob = pMat_M_inv_glob_param[COEF2]->pArray ;
#pragma ivdep
		for (Myint inode = 0; inode < nnode; inode++)
		{
			vx[inode] += dt * Mat_M_inv_glob[inode] * vec_k[inode] ;
		}
		nb_op_kernel += nnode * 3 ;
	}

	//============
	// boundaries
	//============
	Myfloat* sponge_coef_node ;
	if (is_there_boundary_type(SPG))
	{
		if (pSponge_coef_node == NULL)
		{
			print_error(" FEM_2D_1st_ac_iso::compute_velocity, pSponge_coef_node == NULL") ;
			return(RTN_CODE_KO) ;
		}
		else
		{
			sponge_coef_node = pSponge_coef_node->pArray ;
#pragma ivdep
			for (Myint inode = 0; inode < nnode; inode++)
			{
				if (pNode[inode].boundary == SPG)
				{
					// sponge damping
					vz[inode] *= sponge_coef_node[inode] ;
					vx[inode] *= sponge_coef_node[inode] ;
					nb_op_kernel += 2 ;
				}
			}
		}
	}

	// #pragma ivdep
	//   for (Myint inode = 0; inode < nnode; inode++)
	//     {
	//       if (pNode[inode].boundary == RIGID)
	// 	{
	// 	  // rigid surface -> velocity is null
	// 	  vz[inode] = Myfloat(0.0) ;
	// 	  vx[inode] = Myfloat(0.0) ;
	// 	}
	//     }

	print_debug(ALL, FULL_DEBUG, "OUT FEM_2D_1st_ac_iso::compute_velocity");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
void FEM_2D_1st_ac_iso::init_eigen()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st_ac_iso::init_eigen");

	Myfloat time_vz = -dt/2.0 ;
	Myfloat time_pr = -dt ;

	Grid_1D_float *pr_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(PR))->get_grid() ;
	Myfloat * const pr = pr_grid->pArray ;

	Grid_1D_float *vz_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(VZ))->get_grid() ;
	Myfloat * const vz = vz_grid->pArray ;

	Grid_1D_float *vx_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(VX))->get_grid() ;
	Myfloat * const vx = vx_grid->pArray ;

	// get number mode number
	eigen_nmode = (Myint) Singleton::Instance()->pProgram->pModelling->get_param() ;

	Myfloat sq2 = sqrt(2.0) ;

	// loop on elements
	for (Myint iel = 0; iel < nelem; iel++)
	{
		// loop on nodes
		for (Myint ii=0; ii< pElement[iel].nnode ; ii++)
		{
			Myint inode = pElement[iel].inode[ii] ;
			Myfloat znode = pNode[inode].zcoord ;
			Myfloat xnode = pNode[inode].xcoord ;

			pr[inode] = -sq2 * sin(M_PI*xnode * eigen_nmode) * sin(M_PI*znode * eigen_nmode) * sin(sq2*M_PI*time_pr * eigen_nmode) ;
			vz[inode] = cos(M_PI*znode * eigen_nmode) * sin(M_PI*xnode * eigen_nmode) * cos(sq2*M_PI*time_vz * eigen_nmode) ;
			vx[inode] = cos(M_PI*xnode * eigen_nmode) * sin(M_PI*znode * eigen_nmode) * cos(sq2*M_PI*time_vz * eigen_nmode) ;
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st_ac_iso::init_eigen");
	return ;
}

} // namespace django
