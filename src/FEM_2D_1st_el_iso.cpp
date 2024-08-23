//-------------------------------------------------------------------------------------------------------
//
// MODELLING WITH FEM IN 2D (CONTINUOUS/DISCONTINUOUS GALERKIN)
//
//  * ELASTIC ISOTROPIC MEDIA 
//  * 1ST ORDER WAVE EQUATION (STRESS / VELOCITY):
//
//            coef1 dvx/dt      = d(tau+tau')/dx + d(tau'')/dz + coef1 Fx
//            coef1 dvz/dt      = d(tau'')/dx + d(tau-tau')/dz + coef1 Fz
//            coef2 d(tau)/dt   = dvx/dx + dvz/dz 
//            coef3 d(tau')/dt  = dvx/dx - dvz/dz
//            coef3 d(tau'')/dt = dvz/dx + dvx/dz
//
//            with
//
//            coef1 = rho
//            coef2 = 1/(rho * (vp^2-vs^2))
//            coef3 = 1/(rho * vs^2)
//
// FEM FORMULATION
// ===============
// vx = vx + dt * Mat_M_inv(coef1) * vec_k 
// with vec_k = -Mat_Dx * (tau(iel)+tau'(iel)) -Mat_Dz * (tau''(iel))
//              + 0.5 * SUM_neigh [vnx * (Mat_F * (tau(iel)+tau'(iel)) + Mat_G * (tau(iel_neigh)+tau'(iel_neigh)))
//              +vnz * (Mat_F * (tau''(iel)) + Mat_G * (tau''(iel_neigh))) ] 
//
// vz = vz + dt * Mat_M_inv(coef1) * vec_k 
// with vec_k = -Mat_Dx * (tau''(iel)) -Mat_Dz * (tau(iel)-tau'(iel))
//              + 0.5 * SUM_neigh [vnx * (Mat_F * (tau''(iel)) + Mat_G * (tau''(iel_neigh)))
//              +vnz * (Mat_F * (tau(iel)-tau'(iel)) + Mat_G * (tau(iel_neigh)-tau'(iel_neigh))) ]
//
// tau = tau + dt * Mat_M_inv(coef2) * vec_k 
// with vec_k = -Mat_Dx * (vx(iel)) -Mat_Dz * (vz(iel))
//              + 0.5 * SUM_neigh [vnx * (Mat_F * (vx(iel)) + Mat_G * (vx(iel_neigh)))
//              +vnz * (Mat_F * (vz(iel)) + Mat_G * (vz(iel_neigh))) ]
//
// tau' = tau' + dt * Mat_M_inv(coef3) * vec_k 
// with vec_k = -Mat_Dx * (vx(iel)) +Mat_Dz * (vz(iel))
//              + 0.5 * SUM_neigh [vnx * (Mat_F * (vx(iel)) + Mat_G * (vx(iel_neigh)))
//              -vnz * (Mat_F * (vz(iel)) + Mat_G * (vz(iel_neigh))) ]
//
// tau'' = tau'' + dt * Mat_M_inv(coef3) * vec_k 
// with vec_k = -Mat_Dx * (vz(iel)) -Mat_Dz * (vx(iel))
//              + 0.5 * SUM_neigh [vnx * (Mat_F * (vz(iel)) + Mat_G * (vz(iel_neigh)))
//              +vnz * (Mat_F * (vx(iel)) + Mat_G * (vx(iel_neigh))) ]
//
// PARENT CLASS: Scheme
//   DERIVED CLASS: FEM
//     DERIVED CLASS: FEM_2D
//       DERIVED CLASS: FEM_2D_1s
//         DERIVED CLASS: FEM_2D_1st_el_iso
//
//-------------------------------------------------------------------------------------------------------

#include "FEM_2D_1st_el_iso.h"

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
FEM_2D_1st_el_iso::FEM_2D_1st_el_iso(void) : FEM_2D_1st()

{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st_el_iso::FEM_2D_1st_el_iso");

	eq_type = ELASTIC ;
	CFL = 1.0 ;
	print_debug(ALL, MID_DEBUG, "CFL=", CFL) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st_el_iso::FEM_2D_1st_el_iso");
}

//=======================================================================================================
//
// INITIALIZE MODELLING
//
//=======================================================================================================

Rtn_code FEM_2D_1st_el_iso::initialize(Model* pModel, Myfloat fmax)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st_el_iso::initialize");

	// check space interpolation and time order
	Myint polynom_order = Singleton::Instance()->space_order ;
	if (polynom_order < 0)
	{
		print_error(" Can not initialize FEM_2D_1st_el_iso, space_order < 0", polynom_order) ;
		return(RTN_CODE_KO) ;
	}
	Myint time_order = Singleton::Instance()->time_order ;
	if (time_order != 2)
	{
		print_error(" Can not initialize FEM_2D_1st_el_iso, time_order != 2", time_order) ;
		return(RTN_CODE_KO) ;
	}

	// call parent initialiation
	Rtn_code rtn_code = FEM_2D_1st::initialize(pModel, fmax) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// initialize global invert mass matrix
	//===============================================================

	// coef 1 = rho
	rtn_code = compute_global_inv_mass_matrix(pModel, RHO, COEF1) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// coef 2 = 1 /(rho * (vp^2-vs^2))
	rtn_code = compute_global_inv_mass_matrix(pModel, INV_RHOVP2MINUSVS2, COEF2) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// coef 3 = 1 /(rho * vs^2))
	rtn_code = compute_global_inv_mass_matrix(pModel, INV_RHOVS2, COEF3) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	//-------------------------------------------------------------------------------------------------------
	// unknown vectors
	//-------------------------------------------------------------------------------------------------------

	// velocity component (time domain)
	Variable* pVarVz = Singleton::Instance()->register_variable(VZ, "vz") ;
	pVarVz->allocate_grid(nnode, 0) ;
	Variable* pVarVx = Singleton::Instance()->register_variable(VX, "vx") ;
	pVarVx->allocate_grid(nnode, 0) ;

	// stress component (time domain)
	Variable* pVarTau = Singleton::Instance()->register_variable(TAU, "pr") ;
	pVarTau->allocate_grid(nnode, 0) ;
	Variable* pVarTauP = Singleton::Instance()->register_variable(TAUP, "taup") ;
	pVarTauP->allocate_grid(nnode, 0) ;
	Variable* pVarTauPP = Singleton::Instance()->register_variable(TAUPP, "taupp") ;
	pVarTauPP->allocate_grid(nnode, 0) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st_el_iso::initialize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_2D_1st_el_iso::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st_el_iso::finalize");

	// delete components
	Singleton::Instance()->delete_variable(VZ) ;
	Singleton::Instance()->delete_variable(VX) ;
	Singleton::Instance()->delete_variable(TAU) ;
	Singleton::Instance()->delete_variable(TAUP) ;
	Singleton::Instance()->delete_variable(TAUPP) ;

	// call parent finalize
	Rtn_code rtn_code = FEM_2D_1st::finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st_el_iso::finalize");
	return(RTN_CODE_OK) ;
}

//=======================================================================================================
//
// PERFORM FORWARD MODELLING FOR ONE SOURCE
//
//=======================================================================================================

Rtn_code FEM_2D_1st_el_iso::solve_current_shot(Acquisition* pAcquisition, Data *pData, Wavefield_type wtype, Snapshot* pSnapshot)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st_el_iso::solve_current_shot");

	// get variables
	Variable* tau_var = Singleton::Instance()->get_variable(TAU) ;
	Grid_1D_float *tau_grid = (Grid_1D_float*) tau_var->get_grid() ;
	Myfloat * const tau = tau_grid->pArray ;

	Variable* tauP_var = Singleton::Instance()->get_variable(TAUP) ;
	Grid_1D_float *tauP_grid = (Grid_1D_float*) tauP_var->get_grid() ;
	Myfloat * const tauP = tauP_grid->pArray ;

	Variable* tauPP_var = Singleton::Instance()->get_variable(TAUPP) ;
	Grid_1D_float *tauPP_grid = (Grid_1D_float*) tauPP_var->get_grid() ;
	Myfloat * const tauPP = tauPP_grid->pArray ;

	Variable* vz_var = Singleton::Instance()->get_variable(VZ) ;
	Grid_1D_float *vz_grid = (Grid_1D_float*) vz_var->get_grid() ;
	Myfloat * const vz = vz_grid->pArray ;

	Variable* vx_var = Singleton::Instance()->get_variable(VX) ;
	Grid_1D_float *vx_grid = (Grid_1D_float*) vx_var->get_grid() ;
	Myfloat * const vx = vx_grid->pArray ;

	// reset memory
	Rtn_code rtn_code = reset() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// locate src and rec in the mesh
	rtn_code = locate_src_and_rec_in_mesh(pAcquisition) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

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

		for (Myint it = 0; it < nt; it ++)
		{
			//-------------------------------------------------------------------------------------------------------
			// computation of the velocity component at time = it * dt
			//-------------------------------------------------------------------------------------------------------
			double t0 = MPI_Wtime() ;
			compute_velocity(tau, tauP, tauPP, vz, vx) ;
			time_in_kernel += MPI_Wtime() - t0 ;

			//-------------------------------------------------------------------------------------------------------
			// source excitation on velocity component
			//-------------------------------------------------------------------------------------------------------
			if (src_type != NO_SRC_TYPE)
			{
				if (src_stype == FORCE_Z)
				{
					source_excitation(vz, it, wtype, pData, src_factor) ;
				}
			}

			//-------------------------------------------------------------------------------------------------------
			// computation of the stress component at time = [it+1/2] * dt
			//-------------------------------------------------------------------------------------------------------
			t0 = MPI_Wtime() ;
			compute_stress(tau, tauP, tauPP, vz, vx) ;
			time_in_kernel += MPI_Wtime() - t0 ;

			//-------------------------------------------------------------------------------------------------------
			// source excitation on the stress component
			//-------------------------------------------------------------------------------------------------------
			if (src_type != NO_SRC_TYPE)
			{
				if (src_stype == EXPLOSIVE)
				{
					rtn_code = source_excitation(tau, it, wtype, pData, src_factor) ;
				}
			}

			// compute error for eigen mode
			if(Singleton::Instance()->pProgram->pModelling->get_case() == MODELLING_EIGEN)
			{
				rtn_code = compute_eigen_error(tau, it, pAcquisition) ;
				if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
			}

			//-------------------------------------------------------------------------------------------------------
			// write seismograms
			//-------------------------------------------------------------------------------------------------------
			rtn_code = write_trace(vz_var, it) ;
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
			rtn_code = write_trace(vx_var, it) ;
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

			//-------------------------------------------------------------------------------------------------------
			// write snapshot
			//-------------------------------------------------------------------------------------------------------
			rtn_code = write_snapshot(vz_var, it) ;
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
			rtn_code = write_snapshot(vx_var, it) ;
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

			// compute_energy
			compute_energy(vx_var, vz_var, tau_var, tauP_var, tauPP_var, it) ;

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

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st_el_iso::solve_current_shot");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_2D_1st_el_iso::reset(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st_el_iso::reset");

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
		rtn_code = Singleton::Instance()->get_variable(TAU)->reset_grid() ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
		rtn_code = Singleton::Instance()->get_variable(TAUP)->reset_grid() ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
		rtn_code = Singleton::Instance()->get_variable(TAUPP)->reset_grid() ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
		rtn_code = Singleton::Instance()->get_variable(VZ)->reset_grid() ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
		rtn_code = Singleton::Instance()->get_variable(VX)->reset_grid() ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st_el_iso::reset");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_2D_1st_el_iso::compute_stress(Myfloat* tau, Myfloat* tauP, Myfloat* tauPP, Myfloat* vz, Myfloat* vx)
{
	print_debug(ALL, FULL_DEBUG, "IN FEM_2D_1st_el_iso::compute_stress");

	// free surface implementation
	bool flux_formulation = false ;
	//bool flux_formulation = true ;

	// get pointer to global arrays
	Myfloat* const Vec_rho_glob   = pVec_rho_glob->pArray ;
	Myfloat* const Vec_vp_glob    = pVec_vp_glob->pArray ;
	Myfloat* const Vec_vs_glob    = pVec_vs_glob->pArray ;
	Myfloat* const Mat_M_inv_glob = pMat_M_inv_glob->pArray ;
	Myfloat* const vec_k          = pVec_k_glob->pArray ;

	// update ncell
	ncell += nnode ;

	//=====================================================================================================
	// compute TAU
	//=====================================================================================================
	{
		// initialize global vector
		//=========================
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
			const Myint ni = pElement[iel].nnode ;
			const Myint mat_idx_i = pElement[iel].nnode_1D - 1 ;

			// compute -Mat_Dz * vz
			//=====================

			const Myfloat cvol = pElement[iel].size_x / 2.0 ;

			// full matrix
			// Myfloat** const Mat_Dz = pMat_Dz[mat_idx_i]->pArray ;
			// for (Myint ii = 0; ii < ni; ii++) {
			//   const Myint inode_ii = pElement[iel].inode[ii] ;
			//   for (Myint jj = 0; jj < ni; jj++) {
			//     vec_k[inode_ii] -= cvol * Mat_Dz[ii][jj] * vz[pElement[iel].inode[jj]] ;
			//   }
			// }
			// nb_op_kernel += ni * ni * 3 ;

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

			// Myfloat** const Mat_Dx = pMat_Dx[mat_idx_i]->pArray ;
			// for (Myint ii = 0; ii < ni; ii++) {
			//   const Myint inode_ii = pElement[iel].inode[ii] ;
			//   for (Myint jj = 0; jj < ni; jj++) {
			//     vec_k[inode_ii] -= cvol2 * Mat_Dx[ii][jj] * vx[pElement[iel].inode[jj]] ;
			//   }
			// }
			// nb_op_kernel += ni * ni * 3 ;

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

				// free surf
				else if (flux_formulation)
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
						vec_k[ pElement[iel].inode[ Mat_F_i[ii] ] ] += vnz * Mat_F_c[ii] * vz[ pElement[iel].inode[ Mat_F_j[ii] ] ] ;
						vec_k[ pElement[iel].inode[ Mat_F_i[ii] ] ] += vnx * Mat_F_c[ii] * vx[ pElement[iel].inode[ Mat_F_j[ii] ] ] ;
					}
					nb_op_kernel += n1 * 6 ;
				}

			} // for (Myint ineigh = 0; ineigh < 4; ineigh++)
		} // for (Myint iel = 0; iel < nelem; iel++)

		// update tau
		Myfloat* const Mat_M_inv_glob = pMat_M_inv_glob_param[COEF2]->pArray ;
#pragma ivdep
		for (Myint inode = 0; inode < nnode; inode++)
		{
			tau[inode] += dt * Mat_M_inv_glob[inode] * vec_k[inode] ;
		}
		nb_op_kernel += nnode * 3 ;
	}

	//=====================================================================================================
	// compute TAUP
	//=====================================================================================================
	{
		// initialize global vector
		//=========================
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
			const Myint ni = pElement[iel].nnode ;
			const Myint mat_idx_i = pElement[iel].nnode_1D - 1 ;

			// compute +Mat_Dz * vz
			//=====================

			const Myfloat cvol = pElement[iel].size_x / 2.0 ;

			// sparse matrix
			Myint    const n1      = pMat_Dz_i[mat_idx_i]->nz ;
			Myint*   const Mat_Dz_i = pMat_Dz_i[mat_idx_i]->pArray ;
			Myint*   const Mat_Dz_j = pMat_Dz_j[mat_idx_i]->pArray ;
			Myfloat* const Mat_Dz_c = pMat_Dz_c[mat_idx_i]->pArray ;
			for (Myint ii = 0; ii < n1; ii++) {
				vec_k[ pElement[iel].inode[Mat_Dz_i[ii]] ] += cvol * Mat_Dz_c[ii] * vz[ pElement[iel].inode[Mat_Dz_j[ii]] ] ;
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

					// -0.5 * vnz * Mat_F * vz(iel)
					// +0.5 * vnx * Mat_F * vx(iel)
					Myint    const n1      = pMat_F_i[ineigh][mat_idx_i]->nz ;
					Myint*   const Mat_F_i = pMat_F_i[ineigh][mat_idx_i]->pArray ;
					Myint*   const Mat_F_j = pMat_F_j[ineigh][mat_idx_i]->pArray ;
					Myfloat* const Mat_F_c = pMat_F_c[ineigh][mat_idx_i]->pArray ;
					for (Myint ii = 0; ii < n1; ii++) {
						vec_k[ pElement[iel].inode[ Mat_F_i[ii] ] ] -= vnz * Myfloat(0.5) * Mat_F_c[ii] * vz[ pElement[iel].inode[ Mat_F_j[ii] ] ] ;
						vec_k[ pElement[iel].inode[ Mat_F_i[ii] ] ] += vnx * Myfloat(0.5) * Mat_F_c[ii] * vx[ pElement[iel].inode[ Mat_F_j[ii] ] ] ;
					}

					// -0.5 * vnz * Mat_G * vz(iel_neigh)
					// +0.5 * vnx * Mat_G * vx(iel_neigh)
					Myint    const n2 = pMat_G_i[ineigh][mat_idx_i][mat_idx_j]->nz ;
					Myint*   const Mat_G_i = pMat_G_i[ineigh][mat_idx_i][mat_idx_j]->pArray ;
					Myint*   const Mat_G_j = pMat_G_j[ineigh][mat_idx_i][mat_idx_j]->pArray ;
					Myfloat* const Mat_G_c = pMat_G_c[ineigh][mat_idx_i][mat_idx_j]->pArray ;
					for (Myint ii = 0; ii < n2; ii++) {
						vec_k[ pElement[iel].inode[ Mat_G_i[ii] ] ] -= vnz * Myfloat(0.5) * Mat_G_c[ii] * vz[ pElement[iel_neigh].inode[ Mat_G_j[ii] ] ] ;
						vec_k[ pElement[iel].inode[ Mat_G_i[ii] ] ] += vnx * Myfloat(0.5) * Mat_G_c[ii] * vx[ pElement[iel_neigh].inode[ Mat_G_j[ii] ] ] ;
					}

					nb_op_kernel += (n1 + n2) * 8 ;

				} // if (iel_neigh != NO_NEIGH)

				// freesurf
				else if (flux_formulation)
				{
					Myint const mat_idx_j = pElement[iel_neigh].nnode_1D - 1 ;
					Myfloat const vnz = pElement[iel].vnz[ineigh] * Myfloat(0.5) ;
					Myfloat const vnx = pElement[iel].vnx[ineigh] * Myfloat(0.5) ;

					// -0.5 * vnz * Mat_F * vz(iel)
					// +0.5 * vnx * Mat_F * vx(iel)
					Myint    const n1      = pMat_F_i[ineigh][mat_idx_i]->nz ;
					Myint*   const Mat_F_i = pMat_F_i[ineigh][mat_idx_i]->pArray ;
					Myint*   const Mat_F_j = pMat_F_j[ineigh][mat_idx_i]->pArray ;
					Myfloat* const Mat_F_c = pMat_F_c[ineigh][mat_idx_i]->pArray ;
					for (Myint ii = 0; ii < n1; ii++) {
						vec_k[ pElement[iel].inode[ Mat_F_i[ii] ] ] -= vnz * Mat_F_c[ii] * vz[ pElement[iel].inode[ Mat_F_j[ii] ] ] ;
						vec_k[ pElement[iel].inode[ Mat_F_i[ii] ] ] += vnx * Mat_F_c[ii] * vx[ pElement[iel].inode[ Mat_F_j[ii] ] ] ;
					}
					nb_op_kernel += n1 * 6 ;
				}
			} // for (Myint ineigh = 0; ineigh < 4; ineigh++)
		} // for (Myint iel = 0; iel < nelem; iel++)

		// update taup
		Myfloat* const Mat_M_inv_glob = pMat_M_inv_glob_param[COEF3]->pArray ;
#pragma ivdep
		for (Myint inode = 0; inode < nnode; inode++)
		{
			tauP[inode] += dt * Mat_M_inv_glob[inode] * vec_k[inode] ;
		}
		nb_op_kernel += nnode * 3 ;
	}

	//=====================================================================================================
	// compute TAUPP
	//=====================================================================================================
	{
		// initialize global vector
		//=========================
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
			const Myint ni = pElement[iel].nnode ;
			const Myint mat_idx_i = pElement[iel].nnode_1D - 1 ;

			// compute -Mat_Dx * vz
			//=====================

			const Myfloat cvol = pElement[iel].size_z / 2.0 ;

			// sparse matrix
			Myint    const n2      = pMat_Dx_i[mat_idx_i]->nz ;
			Myint*   const Mat_Dx_i = pMat_Dx_i[mat_idx_i]->pArray ;
			Myint*   const Mat_Dx_j = pMat_Dx_j[mat_idx_i]->pArray ;
			Myfloat* const Mat_Dx_c = pMat_Dx_c[mat_idx_i]->pArray ;
			for (Myint ii = 0; ii < n2; ii++) {
				vec_k[ pElement[iel].inode[Mat_Dx_i[ii]] ] -= cvol * Mat_Dx_c[ii] * vz[ pElement[iel].inode[Mat_Dx_j[ii]] ] ;
			}
			nb_op_kernel += n2 * 3 ;

			// compute -Mat_Dz * vx
			//=====================

			const Myfloat cvol2 = pElement[iel].size_x / 2.0 ;

			// sparse matrix
			Myint    const n1      = pMat_Dz_i[mat_idx_i]->nz ;
			Myint*   const Mat_Dz_i = pMat_Dz_i[mat_idx_i]->pArray ;
			Myint*   const Mat_Dz_j = pMat_Dz_j[mat_idx_i]->pArray ;
			Myfloat* const Mat_Dz_c = pMat_Dz_c[mat_idx_i]->pArray ;
			for (Myint ii = 0; ii < n1; ii++) {
				vec_k[ pElement[iel].inode[Mat_Dz_i[ii]] ] -= cvol2 * Mat_Dz_c[ii] * vx[ pElement[iel].inode[Mat_Dz_j[ii]] ] ;
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
					Myfloat const vnx = pElement[iel].vnx[ineigh] * Myfloat(0.5) ;

					// 0.5 * vnz * Mat_F * vx(iel)
					// 0.5 * vnx * Mat_F * vz(iel)
					Myint    const n1      = pMat_F_i[ineigh][mat_idx_i]->nz ;
					Myint*   const Mat_F_i = pMat_F_i[ineigh][mat_idx_i]->pArray ;
					Myint*   const Mat_F_j = pMat_F_j[ineigh][mat_idx_i]->pArray ;
					Myfloat* const Mat_F_c = pMat_F_c[ineigh][mat_idx_i]->pArray ;
					for (Myint ii = 0; ii < n1; ii++) {
						vec_k[ pElement[iel].inode[ Mat_F_i[ii] ] ] += vnz * Myfloat(0.5) * Mat_F_c[ii] * vx[ pElement[iel].inode[ Mat_F_j[ii] ] ] ;
						vec_k[ pElement[iel].inode[ Mat_F_i[ii] ] ] += vnx * Myfloat(0.5) * Mat_F_c[ii] * vz[ pElement[iel].inode[ Mat_F_j[ii] ] ] ;
					}

					// 0.5 * vnz * Mat_G * vx(iel_neigh)
					// 0.5 * vnx * Mat_G * vz(iel_neigh)
					Myint    const n2 = pMat_G_i[ineigh][mat_idx_i][mat_idx_j]->nz ;
					Myint*   const Mat_G_i = pMat_G_i[ineigh][mat_idx_i][mat_idx_j]->pArray ;
					Myint*   const Mat_G_j = pMat_G_j[ineigh][mat_idx_i][mat_idx_j]->pArray ;
					Myfloat* const Mat_G_c = pMat_G_c[ineigh][mat_idx_i][mat_idx_j]->pArray ;
					for (Myint ii = 0; ii < n2; ii++) {
						vec_k[ pElement[iel].inode[ Mat_G_i[ii] ] ] += vnz * Myfloat(0.5) * Mat_G_c[ii] * vx[ pElement[iel_neigh].inode[ Mat_G_j[ii] ] ] ;
						vec_k[ pElement[iel].inode[ Mat_G_i[ii] ] ] += vnx * Myfloat(0.5) * Mat_G_c[ii] * vz[ pElement[iel_neigh].inode[ Mat_G_j[ii] ] ] ;
					}

					nb_op_kernel += (n1 + n2) * 8 ;

				} // if (iel_neigh != NO_NEIGH)

				else if (flux_formulation)
				{
					// free surf
					Myint const mat_idx_j = pElement[iel_neigh].nnode_1D - 1 ;
					Myfloat const vnz = pElement[iel].vnz[ineigh] * Myfloat(0.5) ;
					Myfloat const vnx = pElement[iel].vnx[ineigh] * Myfloat(0.5) ;

					// 0.5 * vnz * Mat_F * vx(iel)
					// 0.5 * vnx * Mat_F * vz(iel)
					Myint    const n1      = pMat_F_i[ineigh][mat_idx_i]->nz ;
					Myint*   const Mat_F_i = pMat_F_i[ineigh][mat_idx_i]->pArray ;
					Myint*   const Mat_F_j = pMat_F_j[ineigh][mat_idx_i]->pArray ;
					Myfloat* const Mat_F_c = pMat_F_c[ineigh][mat_idx_i]->pArray ;
					for (Myint ii = 0; ii < n1; ii++) {
						vec_k[ pElement[iel].inode[ Mat_F_i[ii] ] ] += vnz * Mat_F_c[ii] * vx[ pElement[iel].inode[ Mat_F_j[ii] ] ] ;
						vec_k[ pElement[iel].inode[ Mat_F_i[ii] ] ] += vnx * Mat_F_c[ii] * vz[ pElement[iel].inode[ Mat_F_j[ii] ] ] ;
					}
					nb_op_kernel += n1 * 6 ;
				}
			} // for (Myint ineigh = 0; ineigh < 4; ineigh++)
		} // for (Myint iel = 0; iel < nelem; iel++)

		// update tauPP
		Myfloat* const Mat_M_inv_glob = pMat_M_inv_glob_param[COEF3]->pArray ;
#pragma ivdep
		for (Myint inode = 0; inode < nnode; inode++)
		{
			tauPP[inode] += dt * Mat_M_inv_glob[inode] * vec_k[inode] ;
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
			print_error(" FEM_2D_1st_el_iso::compute_stress, pSponge_coef_node == NULL") ;
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
					tau[inode]   *= sponge_coef_node[inode] ;
					tauP[inode]  *= sponge_coef_node[inode] ;
					tauPP[inode] *= sponge_coef_node[inode] ;
				}
			}
			nb_op_kernel += nnode * 3 ;
		}
	}

	// only for CGM (flux conditions are unstable with CGM)
	if ((!flux_formulation) || (type == SCHEME_CGM))
	{
#pragma ivdep
		for (Myint inode = 0; inode < nnode; inode++)
		{
			if (pNode[inode].boundary == FREESURF)
			{
				// free surface -> stress is null
				tau[inode]   = Myfloat(0.0) ;
				tauP[inode]  = Myfloat(0.0) ;
				tauPP[inode] = Myfloat(0.0) ;
			}
		}
	}

	print_debug(ALL, FULL_DEBUG, "OUT FEM_2D_1st_el_iso::compute_stress");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_2D_1st_el_iso::compute_velocity(Myfloat* tau, Myfloat* tauP, Myfloat* tauPP, Myfloat* vz, Myfloat* vx)
{
	print_debug(ALL, FULL_DEBUG, "IN FEM_2D_1st_el_iso::compute_velocity");

	// get pointer to global arrays
	Myfloat* const vec_k          = pVec_k_glob->pArray ;

	//=====================================================================================================
	// compute VZ
	//=====================================================================================================
	{
		// initialize global vector
		//=========================
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

			// compute -Mat_Dz * (tau-tau')
			//=============================

			const Myfloat cvol = pElement[iel].size_x / 2.0 ;

			// sparse matrix
			Myint    const n1      = pMat_Dz_i[mat_idx_i]->nz ;
			Myint*   const Mat_Dz_i = pMat_Dz_i[mat_idx_i]->pArray ;
			Myint*   const Mat_Dz_j = pMat_Dz_j[mat_idx_i]->pArray ;
			Myfloat* const Mat_Dz_c = pMat_Dz_c[mat_idx_i]->pArray ;
			for (Myint ii = 0; ii < n1; ii++) {
				vec_k[ pElement[iel].inode[Mat_Dz_i[ii]] ] -= cvol * Mat_Dz_c[ii] * (tau[ pElement[iel].inode[Mat_Dz_j[ii]] ]-tauP[ pElement[iel].inode[Mat_Dz_j[ii]] ]) ;
			}
			nb_op_kernel += n1 * 4 ;

			// compute -Mat_Dx * (tau'')
			//=============================

			const Myfloat cvol2 = pElement[iel].size_z / 2.0 ;

			// sparse matrix
			Myint    const n2      = pMat_Dx_i[mat_idx_i]->nz ;
			Myint*   const Mat_Dx_i = pMat_Dx_i[mat_idx_i]->pArray ;
			Myint*   const Mat_Dx_j = pMat_Dx_j[mat_idx_i]->pArray ;
			Myfloat* const Mat_Dx_c = pMat_Dx_c[mat_idx_i]->pArray ;
			for (Myint ii = 0; ii < n2; ii++) {
				vec_k[ pElement[iel].inode[Mat_Dx_i[ii]] ] -= cvol2 * Mat_Dx_c[ii] * tauPP[ pElement[iel].inode[Mat_Dx_j[ii]] ] ;
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

					// 0.5 * vnz * Mat_F * [tau(iel)-tau'(iel)]
					// 0.5 * vnx * Mat_F * [tau''(iel)]
					Myint    const n1      = pMat_F_i[ineigh][mat_idx_i]->nz ;
					Myint*   const Mat_F_i = pMat_F_i[ineigh][mat_idx_i]->pArray ;
					Myint*   const Mat_F_j = pMat_F_j[ineigh][mat_idx_i]->pArray ;
					Myfloat* const Mat_F_c = pMat_F_c[ineigh][mat_idx_i]->pArray ;
					for (Myint ii = 0; ii < n1; ii++) {
						vec_k[ pElement[iel].inode[ Mat_F_i[ii] ] ] += vnz * Myfloat(0.5) * Mat_F_c[ii] * (tau[ pElement[iel].inode[ Mat_F_j[ii] ] ] - tauP[ pElement[iel].inode[ Mat_F_j[ii] ] ])  ;
						vec_k[ pElement[iel].inode[ Mat_F_i[ii] ] ] += vnx * Myfloat(0.5) * Mat_F_c[ii] * (tauPP[ pElement[iel].inode[ Mat_F_j[ii] ] ])  ;
					}

					// 0.5 * vnz * Mat_G * [tau(iel_neigh)-tau'(iel_neigh)]
					// 0.5 * vnx * Mat_G * [tau''(iel_neigh)]
					Myint    const n2 = pMat_G_i[ineigh][mat_idx_i][mat_idx_j]->nz ;
					Myint*   const Mat_G_i = pMat_G_i[ineigh][mat_idx_i][mat_idx_j]->pArray ;
					Myint*   const Mat_G_j = pMat_G_j[ineigh][mat_idx_i][mat_idx_j]->pArray ;
					Myfloat* const Mat_G_c = pMat_G_c[ineigh][mat_idx_i][mat_idx_j]->pArray ;
					for (Myint ii = 0; ii < n2; ii++) {
						vec_k[ pElement[iel].inode[ Mat_G_i[ii] ] ] += vnz * Myfloat(0.5) * Mat_G_c[ii] * (tau[ pElement[iel_neigh].inode[ Mat_G_j[ii] ] ] - tauP[ pElement[iel_neigh].inode[ Mat_G_j[ii] ] ]) ;
						vec_k[ pElement[iel].inode[ Mat_G_i[ii] ] ] += vnx * Myfloat(0.5) * Mat_G_c[ii] * (tauPP[ pElement[iel_neigh].inode[ Mat_G_j[ii] ] ]) ;
					}

					nb_op_kernel += (n1 + n2) * 10 ;

				} // if (iel_neigh != NO_NEIGH)
			} // for (Myint ineigh = 0; ineigh < 4; ineigh++)

		} // for (Myint iel = 0; iel < nelem; iel++)

		// update vz
		//==========

		Myfloat* const Mat_M_inv_glob = pMat_M_inv_glob_param[COEF1]->pArray ;
#pragma ivdep
		for (Myint inode = 0; inode < nnode; inode++)
		{
			vz[inode] += dt * Mat_M_inv_glob[inode] * vec_k[inode] ;
		}
		nb_op_kernel += nnode * 3 ;
	}

	//=====================================================================================================
	// compute VX
	//=====================================================================================================
	{
		// initialize global vector
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

			// compute -Mat_Dx * (tau(iel)+tau'(iel)
			//======================================

			const Myfloat cvol2 = pElement[iel].size_z / 2.0 ;

			// sparse matrix
			Myint    const n2      = pMat_Dx_i[mat_idx_i]->nz ;
			Myint*   const Mat_Dx_i = pMat_Dx_i[mat_idx_i]->pArray ;
			Myint*   const Mat_Dx_j = pMat_Dx_j[mat_idx_i]->pArray ;
			Myfloat* const Mat_Dx_c = pMat_Dx_c[mat_idx_i]->pArray ;
			for (Myint ii = 0; ii < n2; ii++) {
				vec_k[ pElement[iel].inode[Mat_Dx_i[ii]] ] -= cvol2 * Mat_Dx_c[ii] * (tau[ pElement[iel].inode[Mat_Dx_j[ii]] ]+tauP[ pElement[iel].inode[Mat_Dx_j[ii]] ])  ;
			}
			nb_op_kernel += n2 * 4 ;

			// compute -Mat_Dz * (tau''(iel))
			//===============================

			const Myfloat cvol = pElement[iel].size_x / 2.0 ;

			// sparse matrix
			Myint    const n1      = pMat_Dz_i[mat_idx_i]->nz ;
			Myint*   const Mat_Dz_i = pMat_Dz_i[mat_idx_i]->pArray ;
			Myint*   const Mat_Dz_j = pMat_Dz_j[mat_idx_i]->pArray ;
			Myfloat* const Mat_Dz_c = pMat_Dz_c[mat_idx_i]->pArray ;
			for (Myint ii = 0; ii < n1; ii++) {
				vec_k[ pElement[iel].inode[Mat_Dz_i[ii]] ] -= cvol * Mat_Dz_c[ii] * tauPP[ pElement[iel].inode[Mat_Dz_j[ii]] ] ;
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
					Myfloat const vnz = pElement[iel].vnz[ineigh] * Myfloat(0.5) ;

					// 0.5 * vnx * Mat_F * (tau(iel)+tau'(iel))
					// 0.5 * vnz * Mat_F * (tauPP(iel))
					Myint    const n1      = pMat_F_i[ineigh][mat_idx_i]->nz ;
					Myint*   const Mat_F_i = pMat_F_i[ineigh][mat_idx_i]->pArray ;
					Myint*   const Mat_F_j = pMat_F_j[ineigh][mat_idx_i]->pArray ;
					Myfloat* const Mat_F_c = pMat_F_c[ineigh][mat_idx_i]->pArray ;
					for (Myint ii = 0; ii < n1; ii++) {
						vec_k[ pElement[iel].inode[ Mat_F_i[ii] ] ] += vnx * Myfloat(0.5) * Mat_F_c[ii] * (tau[ pElement[iel].inode[ Mat_F_j[ii] ] ] + tauP[ pElement[iel].inode[ Mat_F_j[ii] ] ]) ;
						vec_k[ pElement[iel].inode[ Mat_F_i[ii] ] ] += vnz * Myfloat(0.5) * Mat_F_c[ii] * (tauPP[ pElement[iel].inode[ Mat_F_j[ii] ] ]) ;
					}

					// 0.5 * vnx * Mat_G * (tau(iel_neigh)+tau'(iel_neigh))
					// 0.5 * vnz * Mat_G * (tauPP(iel_neigh))
					Myint    const n2 = pMat_G_i[ineigh][mat_idx_i][mat_idx_j]->nz ;
					Myint*   const Mat_G_i = pMat_G_i[ineigh][mat_idx_i][mat_idx_j]->pArray ;
					Myint*   const Mat_G_j = pMat_G_j[ineigh][mat_idx_i][mat_idx_j]->pArray ;
					Myfloat* const Mat_G_c = pMat_G_c[ineigh][mat_idx_i][mat_idx_j]->pArray ;
					for (Myint ii = 0; ii < n2; ii++) {
						vec_k[ pElement[iel].inode[ Mat_G_i[ii] ] ] += vnx * Myfloat(0.5) * Mat_G_c[ii] * (tau[ pElement[iel_neigh].inode[ Mat_G_j[ii] ] ] + tauP[ pElement[iel_neigh].inode[ Mat_G_j[ii] ] ]);
						vec_k[ pElement[iel].inode[ Mat_G_i[ii] ] ] += vnz * Myfloat(0.5) * Mat_G_c[ii] * (tauPP[ pElement[iel_neigh].inode[ Mat_G_j[ii] ] ]);
					}

					nb_op_kernel += (n1 + n2) * 9 ;

				} // if (iel_neigh != NO_NEIGH)
			} // for (Myint ineigh = 0; ineigh < 4; ineigh++)
		} // for (Myint iel = 0; iel < nelem; iel++)

		// update vx

		Myfloat* const Mat_M_inv_glob = pMat_M_inv_glob_param[COEF1]->pArray ;
#pragma ivdep
		for (Myint inode = 0; inode < nnode; inode++)
		{
			vx[inode] += dt * Mat_M_inv_glob[inode] * vec_k[inode]  ;
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
			print_error(" FEM_2D_1st_el_iso::compute_velocity, pSponge_coef_node == NULL") ;
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
				}
			}
			nb_op_kernel += nnode * 2 ;
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

	print_debug(ALL, FULL_DEBUG, "OUT FEM_2D_1st_el_iso::compute_velocity");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
void FEM_2D_1st_el_iso::init_eigen()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st_el_iso::init_eigen");

	Myfloat time_vz = -dt/2.0 ;
	Myfloat time_tau = -dt ;

	Grid_1D_float *tau_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(TAU))->get_grid() ;
	Myfloat * const tau = tau_grid->pArray ;

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

			tau[inode] = -sq2 * sin(M_PI*xnode * eigen_nmode) * sin(M_PI*znode * eigen_nmode) * sin(sq2*M_PI*time_tau * eigen_nmode) ;
			vz[inode] = cos(M_PI*znode * eigen_nmode) * sin(M_PI*xnode * eigen_nmode) * cos(sq2*M_PI*time_vz * eigen_nmode) ;
			vx[inode] = cos(M_PI*xnode * eigen_nmode) * sin(M_PI*znode * eigen_nmode) * cos(sq2*M_PI*time_vz * eigen_nmode) ;
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st_el_iso::init_eigen");
	return ;
} 

//-------------------------------------------------------------------------------------------------------

void FEM_2D_1st_el_iso::compute_energy(Variable* vx_var, Variable* vz_var, Variable* tau_var,
		Variable* tauP_var, Variable* tauPP_var, Myint it)
{
	print_debug(ALL, FULL_DEBUG, "IN FEM_2D_1st_el_iso::compute_energy");

	if (compute_energy_flag)
	{
		if (dt_out != 0)
		{
			// check if ouput is required at the current time step
			Myint decim = round(dt_out / dt) ;

			if (it%decim == 0)
			{
				// temp variable
				Myfloat rho = 1 ;
				Myfloat vp  = 4000.0 ;
				Myfloat vs  = 2310.0 ;
				Myfloat vp2 = vp*vp ;
				Myfloat vs2 = vs*vs ;

				// define coef.
				Myfloat coef_a = vp2 / (rho*vs2*(3.0*vp2-vs2)) ;
				Myfloat coef_b = (vp2-vs2) / (rho*vs2*(6.0*vp2-2.0*vs2)) ;
				Myfloat coef_c = 1.0/(2.0*rho*vs2) ;

				// retrieve pointers ro grids
				Grid_1D_float *vx_grid = (Grid_1D_float*) vx_var->get_grid() ;
				Myfloat* const vx = vx_grid->pArray ;
				Grid_1D_float *vz_grid = (Grid_1D_float*) vz_var->get_grid() ;
				Myfloat* const vz = vz_grid->pArray ;
				Grid_1D_float *tau_grid = (Grid_1D_float*) tau_var->get_grid() ;
				Myfloat* const tau = tau_grid->pArray ;
				Grid_1D_float *tauP_grid = (Grid_1D_float*) tauP_var->get_grid() ;
				Myfloat* const tauP = tauP_grid->pArray ;
				Grid_1D_float *tauPP_grid = (Grid_1D_float*) tauPP_var->get_grid() ;
				Myfloat* const tauPP = tauPP_grid->pArray ;

				// compute energy
				energy_kin = 0.0 ;
				energy_pot = 0.0 ;

				// loop on the receivers
				Myfloat vxr, vzr, taur, tauPr, tauPPr ;
				for (Myint irec=0; irec<nrec; irec++)
				{
					if (ielem_rec[irec] != NOT_FOUND)
					{
						vxr    = interpolate_variable(vx, ielem_rec[irec], eta_rec[irec], xi_rec[irec]) ;
						vzr    = interpolate_variable(vz, ielem_rec[irec], eta_rec[irec], xi_rec[irec]) ;
						taur   = interpolate_variable(tau, ielem_rec[irec], eta_rec[irec], xi_rec[irec]) ;
						tauPr  = interpolate_variable(tauP, ielem_rec[irec], eta_rec[irec], xi_rec[irec]) ;
						tauPPr = interpolate_variable(tauPP, ielem_rec[irec], eta_rec[irec], xi_rec[irec]) ;

						energy_kin += 0.5 * rho * (vxr*vxr + vzr*vzr) ;
						energy_pot += coef_a * (taur*taur + tauPr*tauPr)
					- coef_b * (taur*taur - tauPr*tauPr) + coef_c * tauPPr*tauPPr ;
					}
				}
				energy_tot = energy_kin + energy_pot ;

				// write energy
				write_energy() ;
			}
		}
	}

	print_debug(ALL, FULL_DEBUG, "OUT FEM_2D_1st_el_iso::compute_energy");
}

} // namespace django
