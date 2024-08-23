//-------------------------------------------------------------------------------------------------------
//
// MODELLING WITH FEM IN 1D (CONTINUOUS/DISCONTINUOUS GALERKIN)
//
// **********************************************************************************************
// *** CURRENT DISCONTINUOUS GALERKIN FORMULATION IS UNSTABLE ABOVE P1 FOR 2ND ORDER WAVE EQ. ***
// **********************************************************************************************
//
//  * ACOUSTIC ISOTROPIC MEDIA 
//  * 2ND ORDER WAVE EQUATION (PRESSURE):
//
//            d2P /dt2 = coef1 d2P/dz2 
//
//            with
//
//            coef1 = Vp^2
//
//
// FEM FORMULATION
// ===============
//
// prn = (dt^2 * Mat_M_inv * vec_k) - prn + 2 * prc
// with vec_k += coef1 * { -Mat_Dz * prc(iel) + 0.5 * SUM_neigh [vnz * (Mat_Fz * prc(iel) +  Mat_Gz * prc(iel_neigh)) ] }
//
// PARENT CLASS: Scheme
//   DERIVED CLASS: FEM
//     DERIVED CLASS: FEM_1D
//       DERIVED CLASS: FEM_1D_2nd
//         DERIVED CLASS: FEM_1D_2nd_ac_iso
//
//-------------------------------------------------------------------------------------------------------

#include "FEM_1D_2nd_ac_iso.h"

#include <cassert>
#include <fstream>
#include <iostream>

#include "mpi.h"

#include "type_def.h"
#include "grid_1D_float.h"
#include "grid_2D_complex.h"
#include "output_report.h"
#include "singleton.h"
#include "type_def.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------
FEM_1D_2nd_ac_iso::FEM_1D_2nd_ac_iso(void) : FEM_1D_2nd()

{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D_2nd_ac_iso::FEM_1D_2nd_ac_iso");

	eq_type = ACOUSTIC ;
	CFL = 1.0 ;
	print_debug(ALL, MID_DEBUG, "CFL=", CFL) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D_2nd_ac_iso::FEM_1D_2nd_ac_iso");
}

//=======================================================================================================
//
// INITIALIZE MODELLING
//
//=======================================================================================================

Rtn_code FEM_1D_2nd_ac_iso::initialize(Model* pModel, Myfloat fmax)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D_2nd_ac_iso::initialize");

	// check space interpolation and time order
	if (Singleton::Instance()->space_order < 0)
	{
		print_error(" Can not initialize FEM_1D_2nd_ac_iso, space_order < 0", Singleton::Instance()->space_order) ;
		return(RTN_CODE_KO) ;
	}
	if (Singleton::Instance()->time_order != 2)
	{
		print_error(" Can not initialize FEM_1D_2nd_ac_iso, time_order != 2", Singleton::Instance()->time_order) ;
		return(RTN_CODE_KO) ;
	}

	// call parent initialiation
	Rtn_code rtn_code = FEM_1D_2nd::initialize(pModel, fmax) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	//-------------------------------------------------------------------------------------------------------
	// unknown vectors
	//-------------------------------------------------------------------------------------------------------

	// pressure component (time domain)
	// current time step
	Variable* pVarPrc = Singleton::Instance()->register_variable(PRC, "prc") ;
	pVarPrc->allocate_grid(nnode, 0) ;

	// new time step
	Variable* pVarPrn = Singleton::Instance()->register_variable(PRN, "pr") ;
	pVarPrn->allocate_grid(nnode, 0) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D_2nd_ac_iso::initialize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_1D_2nd_ac_iso::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D_2nd_ac_iso::finalize");

	// delete components
	Singleton::Instance()->delete_variable(PRC) ;
	Singleton::Instance()->delete_variable(PRN) ;

	// call parent finalize
	Rtn_code rtn_code = FEM_1D_2nd::finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D_2nd_ac_iso::finalize");
	return(RTN_CODE_OK) ;
}

//=======================================================================================================
//
// PERFORM FORWARD MODELLING FOR ONE SOURCE
//
//=======================================================================================================

Rtn_code FEM_1D_2nd_ac_iso::solve_current_shot(Acquisition* pAcquisition, Data *pData, Wavefield_type wtype, Snapshot* pSnapshot)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D_2nd_ac_iso::solve_current_shot");

	// reset memory
	Rtn_code rtn_code = reset() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// locate src and rec in the grid
	rtn_code = locate_src_and_rec_in_mesh(pAcquisition) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// source scaling factor
	Myfloat src_factor = dt * dt ;

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
			// get pointers to variables
			Variable* prc_var = Singleton::Instance()->get_variable(PRC) ;
			Grid_1D_float *prc_grid = (Grid_1D_float*) prc_var->get_grid() ;
			Myfloat * const prc = prc_grid->pArray ;

			Variable* prn_var = Singleton::Instance()->get_variable(PRN) ;
			Grid_1D_float *prn_grid = (Grid_1D_float*) prn_var->get_grid() ;
			Myfloat * const prn = prn_grid->pArray ;

			//-------------------------------------------------------------------------------------------------------
			// computation of the pressure component
			//-------------------------------------------------------------------------------------------------------
			double t0 = MPI_Wtime() ;
			compute_pressure(prc, prn) ;
			time_in_kernel += MPI_Wtime() - t0 ;

			//-------------------------------------------------------------------------------------------------------
			// source excitation
			//-------------------------------------------------------------------------------------------------------
			if (src_type != NO_SRC_TYPE)
			{
				if (src_stype == EXPLOSIVE)
				{
					rtn_code = source_excitation(prn, it-1, wtype, pData, src_factor) ;
				}
			}
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

			// compute energy
			//compute_energy(prn) ;

			// compute error for eigen mode
			if(Singleton::Instance()->pProgram->pModelling->get_case() == MODELLING_EIGEN)
			{
				rtn_code = compute_eigen_error(prn, it, pAcquisition) ;
				if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
			}

			//-------------------------------------------------------------------------------------------------------
			// write seismograms
			//-------------------------------------------------------------------------------------------------------
			rtn_code = write_trace(prn_var, it) ;
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

			//-------------------------------------------------------------------------------------------------------
			// write snapshot
			//-------------------------------------------------------------------------------------------------------
			rtn_code = write_snapshot(prn_var, it) ;
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

			//-------------------------------------------------------------------------------------------------------
			// swap variables
			//-------------------------------------------------------------------------------------------------------
			Singleton::Instance()->swap_variable(PRC, PRN) ;

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

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D_2nd_ac_iso::solve_current_shot");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_1D_2nd_ac_iso::reset(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D_2nd_ac_iso::reset");

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
		rtn_code = Singleton::Instance()->get_variable(PRC)->reset_grid() ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
		rtn_code = Singleton::Instance()->get_variable(PRN)->reset_grid() ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D_2nd_ac_iso::reset");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM_1D_2nd_ac_iso::compute_pressure(Myfloat* prc, Myfloat* prn)
{
	print_debug(ALL, FULL_DEBUG, "IN FEM_1D_2nd_ac_iso::compute_pressure");

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
		Myfloat const c1 = pElement[iel].vmin * pElement[iel].vmin ;
		Myfloat const cvol = 2.0 / pElement[iel].size ;

		// compute -Mat_Dz * prc
		//======================

		Myfloat** const Mat_K = pMat_Dz[mat_idx_i]->pArray ;
		for (Myint ii = 0; ii < ni; ii++) {
			Myint inode_ii = pElement[iel].inode[ii] ;
			for (Myint jj = 0; jj < ni; jj++) {
				vec_k[inode_ii] -= c1 * cvol * Mat_K[ii][jj] * prc[pElement[iel].inode[jj]] ;
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

				// 0.5 * vnz * Mat_Fz * vz(iel)
				Myint    const n1      = pMat_Fz_i[ineigh][mat_idx_i]->nz ;
				Myint*   const Mat_Fz_i = pMat_Fz_i[ineigh][mat_idx_i]->pArray ;
				Myint*   const Mat_Fz_j = pMat_Fz_j[ineigh][mat_idx_i]->pArray ;
				Myfloat* const Mat_Fz_c = pMat_Fz_c[ineigh][mat_idx_i]->pArray ;
				Myfloat  const c2       = 2.0 / pElement[iel].size ;
				for (Myint ii = 0; ii < n1; ii++) {
					vec_k[ pElement[iel].inode[ Mat_Fz_i[ii] ] ] += c1 * vnz * c2
							* Myfloat(0.5) * Mat_Fz_c[ii] * prc[ pElement[iel].inode[ Mat_Fz_j[ii] ] ] ;
				}

				// 0.5 * vnz * Mat_Gz * vz(iel_neigh)
				Myint    const n2 = pMat_Gz_i[ineigh][mat_idx_i][mat_idx_j]->nz ;
				Myint*   const Mat_Gz_i = pMat_Gz_i[ineigh][mat_idx_i][mat_idx_j]->pArray ;
				Myint*   const Mat_Gz_j = pMat_Gz_j[ineigh][mat_idx_i][mat_idx_j]->pArray ;
				Myfloat* const Mat_Gz_c = pMat_Gz_c[ineigh][mat_idx_i][mat_idx_j]->pArray ;
				Myfloat  const c3       = 2.0 / pElement[iel_neigh].size ;
				for (Myint ii = 0; ii < n2; ii++) {
					vec_k[ pElement[iel].inode[ Mat_Gz_i[ii] ] ] += c1 * vnz * c3
							* Myfloat(0.5) * Mat_Gz_c[ii] * prc[ pElement[iel_neigh].inode[ Mat_Gz_j[ii] ] ] ;
				}

			} // if (iel_neigh != NO_NEIGH)
		} // for (Myint ineigh = 0; ineigh < 2; ineigh++)

	} // for (Myint iel = 0; iel < nelem; iel++)

	// update pressure
	//================

	Myfloat* const Mat_M_inv_glob = pMat_M_inv_glob->pArray ;
#pragma ivdep
	for (Myint inode = 0; inode < nnode; inode++)
	{
		prn[inode] = (dt * dt * Mat_M_inv_glob[inode] * vec_k[inode]) - prn[inode] + Myfloat(2.0) * prc[inode] ;
	}
	nb_op_kernel += nnode * 5 ;

	// boundaries
	//===========
	Myfloat* sponge_coef_node ;
	if (is_there_boundary_type(SPG))
	{
		if (pSponge_coef_node == NULL)
		{
			print_error(" FEM_1D_2nd_ac_iso::compute_pressure, pSponge_coef_node == NULL") ;
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
					prn[inode] *= sponge_coef_node[inode] ;
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
			prn[inode] = Myfloat(0.0) ;
		}
	}

	print_debug(ALL, FULL_DEBUG, "OUT FEM_1D_2nd_ac_iso::compute_pressure");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
void FEM_1D_2nd_ac_iso::init_eigen()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D_2nd_ac_iso::init_eigen");

	Myfloat time_prn = -2.0*dt ;
	Myfloat time_prc = -dt ;

	Grid_1D_float *prc_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(PRC))->get_grid() ;
	Myfloat * const prc = prc_grid->pArray ;

	Grid_1D_float *prn_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(PRN))->get_grid() ;
	Myfloat * const prn = prn_grid->pArray ;

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
			prn[inode] = -sin(M_PI*znode * eigen_nmode) * sin(M_PI*time_prn * eigen_nmode) ;
			prc[inode] = -sin(M_PI*znode * eigen_nmode) * sin(M_PI*time_prc * eigen_nmode) ;
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D_2nd_ac_iso::init_eigen");
	return ;
}

} // namespace django
