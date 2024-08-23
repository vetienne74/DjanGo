//-------------------------------------------------------------------------------------------------------
//
// MODELLING WITH FEM IN 2D (CONTINUOUS/DISCONTINUOUS GALERKIN)
//
//  * ACOUSTIC LOSSY MEDIA
//  * 1ST ORDER WAVE EQUATION (PRESSURE / VELOCITY):
//
//            [Adapted from Ref Bilbao 2009, p177 eq 7-28]
//            dvx/dt = [1 + sigma_l d/dt] coef2 dP/dx
//            dvz/dt = [1 + sigma_l d/dt] coef2 dP/dz
//            dP /dt = coef1 [dvx/dx + dvz/dz - sigma_0 P]
//
//            with
//
//            coef1 = 1/kappa = 1/(rho * vp * vp)
//            coef2 = rho
//
// FEM FORMULATION
// ==============
// pr = prp + dt * Mat_M_inv(coef1) * vec_k
// with vec_k = c1 { -Mat_Dz * vz(iel) +  -Mat_Dx * vx(iel)
//                   + 0.5 * SUM_neigh [vnz * (Mat_F * vz(iel) + Mat_G * vz(iel_neigh))
//                                     +vnx * (Mat_F * vx(iel) + Mat_G * vx(iel_neigh)) ] }
//
// for the expression of vz and vx below, pr should be replace with:
//    (1+sigma_l/dt)pr - (sigma_l/dt)*prp
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
//         DERIVED CLASS: FEM_2D_1st_ac_lossy
//
//-------------------------------------------------------------------------------------------------------

#include "FEM_2D_1st_ac_lossy.h"

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
FEM_2D_1st_ac_lossy::FEM_2D_1st_ac_lossy(void) : FEM_2D_1st()

{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st_ac_lossy::FEM_2D_1st_ac_lossy");

	eq_type = AC_LOSSY ;
	CFL = 1.0 ;
	print_debug(ALL, MID_DEBUG, "CFL=", CFL) ;
	varPrId  = 0 ;
	varPrpId = 0 ;
	varVzId  = 0 ;
	varVxId  = 0 ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st_ac_lossy::FEM_2D_1st_ac_lossy");
}

//=======================================================================================================
//
// INITIALIZE MODELLING
//
//=======================================================================================================

Rtn_code FEM_2D_1st_ac_lossy::initialize(Model* pModel, Myfloat fmaxIn)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st_ac_lossy::initialize");

	fmax = fmaxIn ;

	// check space interpolation and time order
	Myint polynom_order = Singleton::Instance()->space_order ;
	if (polynom_order < 0)
	{
		print_error(" Can not initialize FEM_2D_1st_ac_lossy, space_order < 0", polynom_order) ;
		return(RTN_CODE_KO) ;
	}
	Myint time_order = Singleton::Instance()->time_order ;
	if (time_order != 2)
	{
		print_error(" Can not initialize FEM_2D_1st_ac_lossy, time_order != 2", time_order) ;
		return(RTN_CODE_KO) ;
	}

	// call parent initialiation
	Rtn_code rtn_code = FEM_2D_1st::initialize(pModel, fmax) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// get loss model
	//---------------
	Variable* loss2_var = pModel->get_parameter(LOSS2) ;
	if (loss2_var == NULL)
	{
		print_error("IN FEM_2D_1st_ac_lossy::initialize --> LOSS2 model not found");
		return(RTN_CODE_KO) ;
	}

	// get loss grid
	Grid_2D_float* loss2_grid = dynamic_cast<Grid_2D_float*>(loss2_var->get_grid()) ;
	if (loss2_grid == NULL)
	{
		print_error("IN FEM_2D_1st_ac_lossy::initialize --> LOSS2 grid not initialized");
		return(RTN_CODE_KO) ;
	}
	sigma_l = loss2_grid->get_max() ;

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
	varVzId = pVarVz->get_id() ;
	Variable* pVarVx = Singleton::Instance()->register_variable(VX, "vx") ;
	pVarVx->allocate_grid(nnode, 0) ;
	varVxId = pVarVx->get_id() ;

	// pressure component (time domain)
	Variable* pVarPr = Singleton::Instance()->register_variable(PR, "pr") ;
	pVarPr->allocate_grid(nnode, 0) ;
	varPrId = pVarPr->get_id() ;

	// pressure component (time domain) previous time step
	Variable* pVarPrp = Singleton::Instance()->register_variable(PRP, "prp") ;
	pVarPrp->allocate_grid(nnode, 0) ;
	varPrpId = pVarPrp->get_id() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st_ac_lossy::initialize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_2D_1st_ac_lossy::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st_ac_lossy::finalize");

	// delete components
	Singleton::Instance()->delete_variable(VZ) ;
	Singleton::Instance()->delete_variable(VX) ;
	Singleton::Instance()->delete_variable(PR) ;
	Singleton::Instance()->delete_variable(PRP) ;

	// call parent finalize
	Rtn_code rtn_code = FEM_2D_1st::finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st_ac_lossy::finalize");
	return(RTN_CODE_OK) ;
}

//=======================================================================================================
//
// PERFORM FORWARD MODELLING FOR ONE SOURCE
//
//=======================================================================================================

Rtn_code FEM_2D_1st_ac_lossy::solve_current_shot(Acquisition* pAcquisition, Data *pData, Wavefield_type wtype, Snapshot* pSnapshot)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st_ac_lossy::solve_current_shot");

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
			// dynamic p-adaptivity
			//-------------------------------------------------------------------------------------------------------
			rtn_code = dynamic_adapt_freq(it, pAcquisition, pSnapshot) ;
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

			//-------------------------------------------------------------------------------------------------------
			// get variables
			//-------------------------------------------------------------------------------------------------------
			Variable* pr_var = Singleton::Instance()->get_variable(PR) ;
			Grid_1D_float *pr_grid = (Grid_1D_float*) pr_var->get_grid() ;
			Myfloat * const pr = pr_grid->pArray ;

			Variable* prp_var = Singleton::Instance()->get_variable(PRP) ;
			Grid_1D_float *prp_grid = (Grid_1D_float*) prp_var->get_grid() ;
			Myfloat * const prp = prp_grid->pArray ;

			Variable* vz_var = Singleton::Instance()->get_variable(VZ) ;
			Grid_1D_float *vz_grid = (Grid_1D_float*) vz_var->get_grid() ;
			Myfloat * const vz = vz_grid->pArray ;

			Variable* vx_var = Singleton::Instance()->get_variable(VX) ;
			Grid_1D_float *vx_grid = (Grid_1D_float*) vx_var->get_grid() ;
			Myfloat * const vx = vx_grid->pArray ;


			//-------------------------------------------------------------------------------------------------------
			// computation of the pressure component
			//-------------------------------------------------------------------------------------------------------
			double t0 = MPI_Wtime() ;
			compute_pressure(pr, prp, vz, vx) ;
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
			compute_energy() ;

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
			compute_velocity(pr, prp, vz, vx) ;
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

			//-------------------------------------------------------------------------------------------------------
			// swap variables
			//-------------------------------------------------------------------------------------------------------
			Singleton::Instance()->swap_variable(PR, PRP) ;

			// reset variables Id
			varPrId = Singleton::Instance()->get_variable(PR)->get_id() ;
			varPrpId = Singleton::Instance()->get_variable(PRP)->get_id() ;

			// display remaining computation time
			display_remaining_computation_time(it, time_in_kernel) ;

			// write stat for current time step
			write_stat_current_step() ;

		}

		//-------------------------------------------------------------------------------------------------------
		// end loop over time steps
		//-------------------------------------------------------------------------------------------------------

	}

	// free src and rec in the grid
	rtn_code = free_position_arrays() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st_ac_lossy::solve_current_shot");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_2D_1st_ac_lossy::dynamic_adapt_freq(Myint it, Acquisition* pAcquisition, Snapshot* pSnapshot)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st_ac_lossy::dynamic_adapt_freq");

	// check adaptivity type
	if (adaptType != ADAPT_FREQTIME)
	{
		print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st_ac_lossy::dynamic_adapt_freq");
		return(RTN_CODE_OK) ;
	}

	// check if last frequency has been reached
	if (adaptTimeIdx >= adaptFmax.size())
	{
		print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st_ac_lossy::dynamic_adapt_freq");
		return(RTN_CODE_OK) ;
	}

	// check if p-adaptivity is needed at current time
	Myfloat currentTime = it * dt ;
	if (abs(currentTime - adaptTmin.at(adaptTimeIdx)) <= (MY_EPSILON * currentTime))
	{
		// START P-ADAPTIVITY WITH A NEW MESH
		double t0 = MPI_Wtime() ;

		print_line5() ;
		print_info(MASTER, " = Dynamic p-adaptivity triggered") ;
		Myfloat fmax = adaptFmax.at(adaptTimeIdx) ;
		print_info(MASTER, " = At (s) / Freq. (Hz)", currentTime, fmax) ;

		//------------------------------------------------------
		// Step 1: create a new scheme with all initializations
		// includng buiding mesh with new frequency
		//------------------------------------------------------
		FEM_2D_1st_ac_lossy newScheme = FEM_2D_1st_ac_lossy() ;
		newScheme.info() ;
		newScheme.initialize(Singleton::Instance()->pProgram->pDomain->pModel, fmax) ;
		newScheme.mesh_info() ;

		// increment index on freq. vector
		adaptTimeIdx++ ;

		//------------------------------------------------------
		// Step 2: project wavefields
		// From old mesh to new mesh
		//------------------------------------------------------

		// project variables
		Rtn_code rtn_code = newScheme.project_variable(*this, varPrId, newScheme.varPrId) ;
		if (rtn_code != RTN_CODE_OK) return rtn_code ;

		rtn_code = newScheme.project_variable(*this, varPrpId, newScheme.varPrpId) ;
		if (rtn_code != RTN_CODE_OK) return rtn_code ;

		rtn_code = newScheme.project_variable(*this, varVzId, newScheme.varVzId) ;
		if (rtn_code != RTN_CODE_OK) return rtn_code ;

		rtn_code = newScheme.project_variable(*this, varVxId, newScheme.varVxId) ;
		if (rtn_code != RTN_CODE_OK) return rtn_code ;

		// delete old variables and set to new ones
		Singleton::Instance()->delete_variable(varPrId) ;
		varPrId = newScheme.varPrId ;

		Singleton::Instance()->delete_variable(varPrpId) ;
		varPrpId = newScheme.varPrpId ;

		Singleton::Instance()->delete_variable(varVzId) ;
		varVzId = newScheme.varVzId ;

		Singleton::Instance()->delete_variable(varVxId) ;
		varVxId = newScheme.varVxId ;

		//------------------------------------------------------
		// Step 3: reset this scheme from new scheme
		//------------------------------------------------------

		rtn_code = this->free_position_arrays() ;
		if (rtn_code != RTN_CODE_OK) return rtn_code ;

		rtn_code = this->reset_scheme_from(&newScheme) ;
		if (rtn_code != RTN_CODE_OK) return rtn_code ;

		rtn_code = this->locate_src_and_rec_in_mesh(pAcquisition) ;
		if (rtn_code != RTN_CODE_OK) return rtn_code ;

		rtn_code = this->locate_pixel_in_mesh(pSnapshot) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		// END P-ADAPTIVITY WITH A NEW MESH
		double time_in_function = MPI_Wtime() - t0 ;
		print_info(MASTER, " Time set new mesh (s)", (Myfloat) time_in_function) ;

		print_line5() ;

	}
	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D_1st_ac_lossy::dynamic_adapt_freq");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_2D_1st_ac_lossy::reset(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st_ac_lossy::reset");

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
		rtn_code = Singleton::Instance()->get_variable(PRP)->reset_grid() ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
		rtn_code = Singleton::Instance()->get_variable(VZ)->reset_grid() ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
		rtn_code = Singleton::Instance()->get_variable(VX)->reset_grid() ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st_ac_lossy::reset");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_2D_1st_ac_lossy::compute_pressure(Myfloat* pr, Myfloat* prp, Myfloat* vz, Myfloat* vx)
{
	print_debug(ALL, FULL_DEBUG, "IN FEM_2D_1st_ac_lossy::compute_pressure");

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

	Myfloat* const Mat_M_inv_glob = pMat_M_inv_glob_param[COEF1]->pArray ;
#pragma ivdep
	for (Myint inode = 0; inode < nnode; inode++)
	{
		pr[inode] = prp[inode] + dt * Mat_M_inv_glob[inode] * vec_k[inode] ;
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
			print_error(" FEM_2D_1st_ac_lossy::compute_pressure, pSponge_coef_node == NULL") ;
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

	print_debug(ALL, FULL_DEBUG, "OUT FEM_2D_1st_ac_lossy::compute_pressure");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_2D_1st_ac_lossy::compute_velocity(Myfloat* pr, Myfloat* prp, Myfloat* vz, Myfloat* vx)
{
	print_debug(ALL, FULL_DEBUG, "IN FEM_2D_1st_ac_lossy::compute_velocity");

	// lossy terms
	const Myfloat loss_term1 = Myfloat(1.0) + sigma_l / dt ;
	const Myfloat loss_term2 = - sigma_l / dt ;

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
				vec_k[ pElement[iel].inode[Mat_Dz_i[ii]] ] -= cvol * Mat_Dz_c[ii] *
						(loss_term1 * pr [ pElement[iel].inode[Mat_Dz_j[ii]] ]
										   + loss_term2 * prp [ pElement[iel].inode[Mat_Dz_j[ii]] ] ) ;
			}
			nb_op_kernel += n1 * 6 ;

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
						vec_k[ pElement[iel].inode[ Mat_F_i[ii] ] ] += vnz * Myfloat(0.5)
						* Mat_F_c[ii] * (loss_term1 * pr[ pElement[iel].inode[ Mat_F_j[ii] ] ]
														  + loss_term2 * prp[ pElement[iel].inode[ Mat_F_j[ii] ] ] );
					}

					// 0.5 * vnz * Mat_G * pr(iel_neigh)
					Myint    const n2 = pMat_G_i[ineigh][mat_idx_i][mat_idx_j]->nz ;
					Myint*   const Mat_G_i = pMat_G_i[ineigh][mat_idx_i][mat_idx_j]->pArray ;
					Myint*   const Mat_G_j = pMat_G_j[ineigh][mat_idx_i][mat_idx_j]->pArray ;
					Myfloat* const Mat_G_c = pMat_G_c[ineigh][mat_idx_i][mat_idx_j]->pArray ;
					for (Myint ii = 0; ii < n2; ii++) {
						vec_k[ pElement[iel].inode[ Mat_G_i[ii] ] ] += vnz * Myfloat(0.5)
						* Mat_G_c[ii] * (loss_term1 * pr[ pElement[iel_neigh].inode[ Mat_G_j[ii] ] ]
														  + loss_term2 * prp[ pElement[iel_neigh].inode[ Mat_G_j[ii] ] ] );
					}

					nb_op_kernel += (n1 + n2) * 7 ;

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
				vec_k[ pElement[iel].inode[Mat_Dx_i[ii]] ] -= cvol2 *
						Mat_Dx_c[ii] * (loss_term1 * pr [ pElement[iel].inode[Mat_Dx_j[ii]] ]
														  + loss_term2 * prp [ pElement[iel].inode[Mat_Dx_j[ii]] ] );
			}
			nb_op_kernel += n1 * 6 ;

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
						vec_k[ pElement[iel].inode[ Mat_F_i[ii] ] ] += vnx * Myfloat(0.5)
						* Mat_F_c[ii] * (loss_term1 * pr[ pElement[iel].inode[ Mat_F_j[ii] ] ]
														  + loss_term2 * prp[ pElement[iel].inode[ Mat_F_j[ii] ] ] );
					}

					// 0.5 * vnx * Mat_G * pr(iel_neigh)
					Myint    const n2 = pMat_G_i[ineigh][mat_idx_i][mat_idx_j]->nz ;
					Myint*   const Mat_G_i = pMat_G_i[ineigh][mat_idx_i][mat_idx_j]->pArray ;
					Myint*   const Mat_G_j = pMat_G_j[ineigh][mat_idx_i][mat_idx_j]->pArray ;
					Myfloat* const Mat_G_c = pMat_G_c[ineigh][mat_idx_i][mat_idx_j]->pArray ;
					for (Myint ii = 0; ii < n2; ii++) {
						vec_k[ pElement[iel].inode[ Mat_G_i[ii] ] ] += vnx * Myfloat(0.5)
						* Mat_G_c[ii] * (loss_term1 * pr[ pElement[iel_neigh].inode[ Mat_G_j[ii] ] ]
														  + loss_term2 * prp[ pElement[iel_neigh].inode[ Mat_G_j[ii] ] ] );
					}

					nb_op_kernel += (n1 + n2) * 7 ;

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
			print_error(" FEM_2D_1st_ac_lossy::compute_velocity, pSponge_coef_node == NULL") ;
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

	print_debug(ALL, FULL_DEBUG, "OUT FEM_2D_1st_ac_lossy::compute_velocity");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
void FEM_2D_1st_ac_lossy::init_eigen()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st_ac_lossy::init_eigen");

	Myfloat time_vz = -dt/2.0 ;
	Myfloat time_pr = -dt ;

	// init P a -dt (*** Grid PRP ***)
	Grid_1D_float *pr_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(PRP))->get_grid() ;
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

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st_ac_lossy::init_eigen");
	return ;
}

//-------------------------------------------------------------------------------------------------------
void FEM_2D_1st_ac_lossy::compute_energy()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D_1st_ac_lossy::compute_energy");

	if (compute_energy_flag)
	{
		// check if ouput is required at the current time step
		Myint decim = round(dt_out / dt) ;

		if (it%decim != 0)
		{
			print_debug(ALL, FULL_DEBUG, "OUT FEM_2D_1st_ac_lossy::compute_energy");
			return ;
		}

		// kinetic energy
		// = 0.5 * rho * (vx^2 + vz^2)
		//----------------------------

		Grid_1D_float *vz_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(VZ))->get_grid() ;
		Myfloat * const vz = vz_grid->pArray ;

		Grid_1D_float *vx_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(VX))->get_grid() ;
		Myfloat * const vx = vx_grid->pArray ;

		// normalize energy by np
		energy_kin = 0.0 ;
		Myint64 np = vz_grid->nz ;
		const Myfloat coef = 0.5 / np ;
		for (Myint ip = 0; ip < np; ip++)
		{
			energy_kin += coef * (vz[ip] * vz[ip] + vx[ip] * vx[ip]) ;
		}

		// potential energy
		// = 0.5 * P^2 / (rho * vp^2)
		//---------------------------

		Grid_1D_float *pr_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(PR))->get_grid() ;
		Myfloat * const pr = pr_grid->pArray ;

		Grid_1D_float *prp_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(PRP))->get_grid() ;
		Myfloat * const prp = prp_grid->pArray ;

		Myfloat* const Vec_vp_glob    = pVec_vp_glob->pArray ;

		// to get accurate energy, pr should be at same time than vx, vz
		// use average (pr + prp) / 2 in this purpose
		// normalize energy by np
		energy_pot = 0.0 ;
		for (Myint ip = 0; ip < np; ip++)
		{
			energy_pot += coef * Vec_vp_glob[ip] * Vec_vp_glob[ip] *
					(pr[ip] + prp[ip]) * (pr[ip] + prp[ip]) / 4 ;
		}

		// total energy
		energy_tot = energy_kin + energy_pot ;

		// write energy
		write_energy() ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D_1st_ac_lossy::compute_energy");
}

} // namespace django
