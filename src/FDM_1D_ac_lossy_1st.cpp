//------------------------------------------------------------------------------------------------------
//
// MODELLING WITH FDM IN 1D
//
//  * ACOUSTIC ISOTROPIC LOSSY MEDIA
//  * 1ST ORDER WAVE EQUATION (PRESSURE / VELOCITY):
//
//            [Adapted from Ref Bilbao 2009, p177 eq 7-28]
//            dvz/dt = [1 + sigma_l d/dt] coef2 dP/dz
//            dP /dt = coef1 [dvz/dz - sigma_0 P]
//
//            with
//
//            coef1 = kappa = (rho * vp * vp) at pr location 
//            coef2 = 1/rho
//            where rho is CONSTANT
//            sigma_0 = freq. independent lossy term
//            sigma_l = freq. dependent lossy term
//
// * FD DISCRETISATION
//
//            Dt. vz = [1 + sigma_l Dt-] coef2 Dz. P
//            Dt. P  = coef1 [Dz. vz - sigma_o At. P]
//
//            Where Dt., Dz. = 1st derivative centered FD
//                  Dt-      = 1st derivative backward FD
//                  At.      = Centered average
//
// PARENT CLASS: Scheme
//   DERIVED CLASS: FDM
//     DERIVED CLASS: FDM_1D
//       DERIVED CLASS: FDM_1D_ac_lossy
//         DERIVED CLASS: FDM_1D_ac_lossy_1st
//
//-------------------------------------------------------------------------------------------------------

#include "FDM_1D_ac_lossy_1st.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>

#include "mpi.h"

#include "acquisition.h"
#include "data.h"
#include "data_std.h"
#include "grid.h"
#include "grid_2D_complex.h"
#include "model.h"
#include "output_report.h"
#include "singleton.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------
FDM_1D_ac_lossy_1st::FDM_1D_ac_lossy_1st(void) : FDM_1D_ac_lossy()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D_ac_lossy_1st::FDM_1D_ac_lossy_1st");
	eq_order = ORDER_1ST ;
	coef2 = 0.0 ;
	coef1 = NULL ;
	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D_ac_lossy_1st::FDM_1D_ac_lossy_1st");
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FDM_1D_ac_lossy_1st::update_physical_coef(Myint iz, Myint ix, Myint iy, Myfloat perturb)  
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D_ac_lossy_1st::update_physical_coef");

	// retrieve vp and rho
	Myfloat rho = TEMP_RHO_CONST ;
	Myfloat vp = sqrt(coef1->pArray[iz]/(dt * rho)) ;

	// add perturbation
	vp += perturb ;

	// update physical parameter
	coef1->pArray[iz] = dt * pow(vp, 2) * rho ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D_ac_lossy_1st::update_physical_coef");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FDM_1D_ac_lossy_1st::update_physical_coef(Model* pModel)  
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D_ac_lossy_1st::update_physical_coef");

	Myfloat rho = TEMP_RHO_CONST ;

	// get vp model
	Grid_1D_float *vp_grid = (Grid_1D_float*) (pModel->get_parameter(VP))->get_grid() ;
	Myfloat* vp = vp_grid->pArray ;

	// coef1_i = rho_i * vp_i**2 * dt (located on the pressure nodes)

	// loop on the grid points in medium
#pragma omp parallel for   
#pragma ivdep
	for (Myint iz = 0; iz < nz; iz++)
	{
		// update coef in the FDM grid point
		coef1->pArray[iz+izBeg2] = dt * pow(vp[iz], 2) * rho ;
	}

	// coef2_i = dt / rho_i (located on the velocity nodes)

	coef2 = dt / rho ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D_ac_lossy_1st::update_physical_coef");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FDM_1D_ac_lossy_1st::initialize(Model* pModel, Myfloat fmax)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D_ac_lossy_1st::initialize");

	// call parent initialiation
	Rtn_code rtn_code = FDM_1D_ac_lossy::initialize(pModel, fmax) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	//-------------------------------------------------------------------------------------------------------
	// unknown vectors
	//-------------------------------------------------------------------------------------------------------

	// velocity component (time domain)
	Variable* pVarVz = Singleton::Instance()->register_variable(VZ, "vz") ;
	pVarVz->allocate_grid(izEnd, dz) ;

	// pressure component (time domain)
	Variable* pVarPr = Singleton::Instance()->register_variable(PR, "pr") ;
	pVarPr->allocate_grid(izEnd+1, dz) ;

	// pressure component (time domain) previous time step
	Variable* pVarPrp = Singleton::Instance()->register_variable(PRP, "prp") ;
	pVarPrp->allocate_grid(izEnd+1, dz) ;

	//-------------------------------------------------------------------------------------------------------
	// initialize eq. coef. in the medium
	//-------------------------------------------------------------------------------------------------------

	// check model type is GRID and subtype is REGULAR
	if (pModel->get_type() != GRID)
	{
		print_error("IN FDM_1D_ac_lossy_1st::initialize --> model type is not GRID");
		return(RTN_CODE_KO) ;
	}
	if (pModel->get_sub_type() != REGULAR)
	{
		print_error("IN FDM_1D_ac_lossy_1st::initialize --> model subtype is not REGULAR");
		return(RTN_CODE_KO) ;
	}

	Myfloat rho = TEMP_RHO_CONST ;

	// coef1_i = rho_i * vp_i**2 * dt (located on the pressure nodes)

	coef1 = new Grid_1D_float(izEnd+1, dz) ;

	// get vp model
	Variable* vp_var = pModel->get_parameter(VP) ;
	if (vp_var == NULL)
	{
		print_error("IN FDM_1D_ac_lossy_1st::initialize --> VP model not found");
		return(RTN_CODE_KO) ;
	}
	// get vp grid
	Grid_1D_float* vp_grid = dynamic_cast<Grid_1D_float*>(vp_var->get_grid()) ;
	if (vp_grid == NULL)
	{
		print_error("IN FDM_1D_ac_lossy_1st::initialize --> VP grid not initialized");
		return(RTN_CODE_KO) ;
	}
	Myfloat* vp = vp_grid->pArray ;

	// loop on the grid points in medium
#pragma omp parallel for 
#pragma ivdep
	for (Myint iz = 0; iz < nz; iz++)
	{
		// update coef in the FDM grid point
		coef1->pArray[iz+izBeg2] = dt * vp[iz] * vp[iz] * rho ;
	}

	// coef2_i = dt / rho_i (located on the velocity nodes)

	coef2 = dt / rho ;

	//-------------------------------------------------------------------------------------------------------
	// initialize eq. coef. in the surrounding boundary
	//-------------------------------------------------------------------------------------------------------

	//---------------------------------------------------------------------------------------------------
	// random layer
	//---------------------------------------------------------------------------------------------------

	// z- layer
	if (get_boundary_type(ZBEG) == RANDOM)
	{
		Myfloat ratio, vp_rand, vp_tmp ;
		vp_tmp = vp[0] ;

		for (Myint iz = izBeg1; iz < izBeg2; iz++)
		{
			ratio = (float) (izBeg2-iz)/nlayer_zBeg ;
			vp_rand = vp_tmp - (rand() % (int) (RAND_BOUND_COEF1 * vp_tmp)) * ratio ;
			coef1->pArray[iz] = dt * vp_rand * vp_rand * rho ;
		}
	}

	// z+ layer
	if (get_boundary_type(ZEND) == RANDOM)
	{
		Myfloat ratio, vp_rand, vp_tmp ;
		vp_tmp = vp[nz-1] ;

		for (Myint iz = izEnd2; iz < izEnd1; iz++)
		{
			ratio = (float) (iz-(izEnd2-1))/nlayer_zEnd ;
			vp_rand = vp_tmp - (rand() % (int) (RAND_BOUND_COEF1 * vp_tmp)) * ratio ;
			coef1->pArray[iz] = dt * vp_rand * vp_rand * rho ;
		}
	}

	//---------------------------------------------------------------------------------------------------
	// sponge
	//---------------------------------------------------------------------------------------------------

	if (is_there_boundary_type(SPG))
	{
		// allocate sponge arrays
		sponge_coef_pr = new Grid_1D_float(izEnd+1, dz) ;
		sponge_coef_vz = new Grid_1D_float(izEnd, dz) ;

		// loop on all grid points
		for (Myint iz = izBeg1; iz < izEnd1; iz++)
		{
			// default value
			sponge_coef_pr->pArray[iz] = 1. ;
			sponge_coef_vz->pArray[iz] = 1. ;

			Myint iz2 = iz - izBeg2 ;
			iz2 = max(iz2, 0) ;
			iz2 = min(iz2, nz-1) ;

			Myfloat vp_tmp = vp[iz2] ;
			coef1->pArray[iz] = dt * vp_tmp * vp_tmp * rho ;

			// default value
			sponge_coef_pr->pArray[iz] = 1. ;
			sponge_coef_vz->pArray[iz] = 1. ;

			// z- layer
			if ((iz >= izBeg1) && (iz < izBeg2) && (get_boundary_type(ZBEG) == SPG))
			{
				Myfloat sponge_coef = get_boundary_coef(ZBEG) ;
				sponge_coef_pr->pArray[iz] *= exp(-pow(Myfloat(izBeg2 - iz) * sponge_coef, 2)) ;
				sponge_coef_vz->pArray[iz] *= exp(-pow(Myfloat(izBeg2 - iz - 0.5) * sponge_coef, 2)) ;
			}

			// z+ layer
			if ((iz >= izEnd2) && (iz < izEnd1) && (get_boundary_type(ZEND) == SPG))
			{
				Myfloat sponge_coef = get_boundary_coef(ZEND) ;
				sponge_coef_pr->pArray[iz]  *= exp(-pow(Myfloat(iz - izEnd2 + 1) * sponge_coef, 2)) ;
				sponge_coef_vz->pArray[iz-1] = exp(-pow(Myfloat(iz - izEnd2 + 0.5) * sponge_coef, 2)) ;
			}
		}
	}

	//---------------------------------------------------------------------------------------------------
	// CPML
	//---------------------------------------------------------------------------------------------------

	// z- layer
	if (get_boundary_type(ZBEG) == PML)
	{

		// compute some coefficients
		Myfloat cpml_alpha_max = CPML_FREQ * PI ;
		Myfloat cpml_pow       = CPML_POW ;
		Myfloat cpml_vmax      = vp_grid->get_max() ;
		Myfloat cpml_rcoef     = get_boundary_coef(ZBEG) ;

		// allocate CPML arrays
		Myint nzl = nlayer_zBeg+2*lstencil ;
		apml_zBeg       = new Grid_1D_float(nzl, dz) ;
		bpml_zBeg       = new Grid_1D_float(nzl, dz) ;
		apml_half_zBeg  = new Grid_1D_float(nzl, dz) ;
		bpml_half_zBeg  = new Grid_1D_float(nzl, dz) ;
		mem_vz_zBeg     = new Grid_1D_float(nzl, dz) ;
		mem_pr_zBeg     = new Grid_1D_float(nzl, dz) ;

		// fill CPML arrays
		Myfloat vp_tmp = vp[0] ;
		Myfloat d0     = -(cpml_pow+1.)*cpml_vmax*log10(cpml_rcoef) / (2.*nlayer_zBeg*dz) ;
		Myint iz1      = izBeg1 ;
		Myint iz2      = izBeg2 ;
#pragma omp parallel for
#pragma ivdep
		for (Myint iz = iz1; iz < iz2; iz++)
		{
			Myint ipml = iz ;

			// pressure node
			coef1->pArray[iz] = dt * vp_tmp * vp_tmp * rho ;

			Myfloat dnorm = Myfloat(izBeg2 - iz)/Myfloat(nlayer_zBeg) ;
			if (dnorm < 0) {dnorm = 0. ;}
			Myfloat alpha = cpml_alpha_max * (1. - dnorm) ;
			Myfloat dd = d0 * pow(dnorm, cpml_pow) ;

			bpml_zBeg->pArray[ipml] = exp(-(dd+alpha)*dt) ;
			apml_zBeg->pArray[ipml] = dd * (bpml_zBeg->pArray[ipml]-1.) / (dd+alpha) ;

			// velocity node
			dnorm = Myfloat(izBeg2 - iz -0.5)/Myfloat(nlayer_zBeg) ;
			if (dnorm < 0) {dnorm = 0. ;}
			alpha = cpml_alpha_max * (1. - dnorm) ;
			dd = d0 * pow(dnorm, cpml_pow) ;

			bpml_half_zBeg->pArray[ipml] = exp(-(dd+alpha)*dt) ;
			apml_half_zBeg->pArray[ipml] = dd * (bpml_half_zBeg->pArray[ipml]-1.) / (dd+alpha) ;
		}
	}

	// z+ layer
	if (get_boundary_type(ZEND) == PML)
	{
		// compute some coefficients
		Myfloat cpml_alpha_max = CPML_FREQ * PI ;
		Myfloat cpml_pow       = CPML_POW ;
		Myfloat cpml_vmax      = vp_grid->get_max() ;
		Myfloat cpml_rcoef     = get_boundary_coef(ZEND) ;

		// allocate CPML arrays
		Myint nzl = nlayer_zEnd+2*lstencil ;
		apml_zEnd       = new Grid_1D_float(nzl, dz) ;
		bpml_zEnd       = new Grid_1D_float(nzl, dz) ;
		apml_half_zEnd  = new Grid_1D_float(nzl, dz) ;
		bpml_half_zEnd  = new Grid_1D_float(nzl, dz) ;
		mem_vz_zEnd     = new Grid_1D_float(nzl, dz) ;
		mem_pr_zEnd     = new Grid_1D_float(nzl, dz) ;

		// fill CPML arrays
		Myint iz1 = izEnd2 ;
		Myint iz2 = izEnd1 ;
		Myfloat vp_tmp = vp[nz-1] ;
		Myfloat d0     = -(cpml_pow+1.)*cpml_vmax*log10(cpml_rcoef) / (2.*nlayer_zEnd*dz) ;
		Myint const offset_z = nlayer_zBeg + nz ;
#pragma omp parallel for
#pragma ivdep
		for (Myint iz = iz1; iz < iz2; iz++)
		{
			Myint ipml = iz - offset_z ;

			// pressure node
			coef1->pArray[iz] = dt * vp_tmp * vp_tmp * rho ;

			Myfloat dnorm = Myfloat(iz - izEnd2 + 1)/Myfloat(nlayer_zEnd) ;
			if (dnorm < 0) {dnorm = 0. ;}
			Myfloat alpha = cpml_alpha_max * (1. - dnorm) ;
			Myfloat dd = d0 * pow(dnorm, cpml_pow) ;

			bpml_zEnd->pArray[ipml] = exp(-(dd+alpha)*dt) ;
			apml_zEnd->pArray[ipml] = dd * (bpml_zEnd->pArray[ipml]-1.) / (dd+alpha) ;

			// velocity node
			dnorm = Myfloat(iz - izEnd2 + 0.5)/Myfloat(nlayer_zEnd) ;
			if (dnorm < 0) {dnorm = 0. ;}
			alpha = cpml_alpha_max * (1. - dnorm) ;
			dd = d0 * pow(dnorm, cpml_pow) ;

			bpml_half_zEnd->pArray[ipml] = exp(-(dd+alpha)*dt) ;
			apml_half_zEnd->pArray[ipml] = dd * (bpml_half_zEnd->pArray[ipml]-1.) / (dd+alpha) ;
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D_ac_lossy_1st::initialize");
	return(RTN_CODE_OK) ;
} ;

//-------------------------------------------------------------------------------------------------------
Rtn_code FDM_1D_ac_lossy_1st::reset(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D_ac_lossy_1st::reset");

	// call to parent class
	Rtn_code rtn_code = Scheme::reset() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// reset grids
	rtn_code = Singleton::Instance()->get_variable(PR)->reset_grid() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
	rtn_code = Singleton::Instance()->get_variable(PRP)->reset_grid() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
	rtn_code = Singleton::Instance()->get_variable(VZ)->reset_grid() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// initialize memory variables
	if (get_boundary_type(ZBEG) == PML)
	{
		mem_vz_zBeg->reset() ;
		mem_pr_zBeg->reset() ;
	}
	if (get_boundary_type(ZEND) == PML)
	{
		mem_vz_zEnd->reset() ;
		mem_pr_zEnd->reset() ;
	}

	// initialization for eigen mode test
	if(Singleton::Instance()->pProgram->pModelling->get_case() == MODELLING_EIGEN)
	{
		// initialization for eigen mode test
		init_eigen() ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D_ac_lossy_1st::reset");
	return(RTN_CODE_OK) ;

} ;

//=======================================================================================================
//
// FINALIZE MODELLING 
//
//=======================================================================================================

Rtn_code FDM_1D_ac_lossy_1st::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D_ac_lossy_1st::finalize");

	//---------------------------------------------------------------------------------------------------
	// deallocate arrays
	//---------------------------------------------------------------------------------------------------

	// physical parameter
	delete(coef1) ;

	// delete components
	Singleton::Instance()->delete_variable(VZ) ;
	Singleton::Instance()->delete_variable(PR) ;
	Singleton::Instance()->delete_variable(PRP) ;

	// call parent finalize
	Rtn_code rtn_code = FDM_1D_ac_lossy::finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D_ac_lossy_1st::finalize");
	return(RTN_CODE_OK) ;

} ;

//=======================================================================================================
//
// PERFORM FORWARD MODELLING FOR ONE SOURCE
//
//=======================================================================================================


Rtn_code FDM_1D_ac_lossy_1st::solve_current_shot(Acquisition* pAcquisition, Data *pData, Wavefield_type wtype, Snapshot* pSnapshot)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D_ac_lossy_1st::solve_current_shot");

	// reset memory
	Rtn_code rtn_code = reset() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// locate src and rec in the grid
	rtn_code = locate_src_and_rec_in_grid(pAcquisition) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// source scaling factor
	Myfloat src_factor = dt / dz ;

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
			// get variables
			Variable* pr_var = Singleton::Instance()->get_variable(PR) ;
			Grid_1D_float *pr_grid = (Grid_1D_float*) pr_var->get_grid() ;
			Myfloat * const pr = pr_grid->pArray ;

			Variable* vz_var = Singleton::Instance()->get_variable(VZ) ;
			Grid_1D_float *vz_grid = (Grid_1D_float*) vz_var->get_grid() ;
			Myfloat * const vz = vz_grid->pArray ;

			//-------------------------------------------------------------------------------------------------------
			// computation of the pressure component
			//-------------------------------------------------------------------------------------------------------
			double t0 = MPI_Wtime() ;
			compute_pressure() ;
			time_in_kernel += MPI_Wtime() - t0 ;

			//-------------------------------------------------------------------------------------------------------
			// source excitation
			//-------------------------------------------------------------------------------------------------------
			if (src_stype == EXPLOSIVE)
			{
				rtn_code = source_excitation(pr, it, wtype, pData, src_factor) ;
			}
			else if (src_stype == FORCE_Z)
			{
				rtn_code = source_excitation(vz, it, wtype, pData, src_factor) ;
			}
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

			// compute energy
			compute_energy(pr) ;

			//-------------------------------------------------------------------------------------------------------
			// computation of the velocity component
			//-------------------------------------------------------------------------------------------------------
			t0 = MPI_Wtime() ;
			compute_velocity() ;
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

			//-------------------------------------------------------------------------------------------------------
			// swap variables
			//-------------------------------------------------------------------------------------------------------
			Singleton::Instance()->swap_variable(PR, PRP) ;

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

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D_ac_lossy_1st::solve_current_shot");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
void FDM_1D_ac_lossy_1st::init_eigen()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D_ac_lossy_1st::init_eigen");

	{
		// init vz a -dt/2
		Grid_1D_float *vz_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(VZ))->get_grid() ;
		Myfloat * const vz = vz_grid->pArray ;
		Myfloat time_vz = -dt/2.0 ;
		Myint izmin = izBeg2 ;
		Myint izmax = izEnd2-1 ;
		Myint nz = izEnd2 - izBeg2 ;
		for (Myint iz=izmin; iz<izmax; iz++)
		{
			Myfloat zz = float(0.5+iz-izmin)/float(nz-1) ;
			vz[iz] = cos(M_PI*zz) * cos(M_PI*time_vz) ;
		}
		// ghost initilization
		Myint iz0 = izBeg2 ;
		for (Myint ii = 1; ii <= lstencil; ii++)
		{
			vz[iz0-ii] = vz[iz0+ii-1] ;
		}
		iz0 = izEnd2-1 ;
		for (Myint ii = 1; ii <= lstencil; ii++)
		{
			vz[iz0+ii-1] = vz[iz0-ii] ;
		}
	}
	{
		// init P a -dt (*** Grid PRP ***)
		Grid_1D_float *pr_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(PRP))->get_grid() ;
		Myfloat * const pr = pr_grid->pArray ;
		Myfloat time_pr = -dt ;
		Myint izmin = izBeg2 ;
		Myint izmax = izEnd2 ;
		Myint nz = izEnd2 - izBeg2 ;
		for (Myint iz=izmin; iz<izmax; iz++)
		{
			Myfloat zz = float(iz-izmin)/float(nz-1) ;
			pr[iz] = -sin(M_PI*zz) * sin(M_PI*time_pr) ;
		}
		Myint iz0 = izBeg2 ;
		pr[iz0] = 0. ;
		for (Myint ii = 1; ii <= lstencil; ii++)
		{
			pr[iz0-ii] = -pr[iz0+ii] ;
		}
		iz0 = izEnd2-1 ;
		pr[iz0] = 0. ;
		for (Myint ii = 1; ii <= lstencil; ii++)
		{
			pr[iz0+ii] = -pr[iz0-ii] ;
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D_ac_lossy_1st::init_eigen");
	return ;
} ;

} // namespace django
