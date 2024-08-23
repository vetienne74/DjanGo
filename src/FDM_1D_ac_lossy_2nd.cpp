//------------------------------------------------------------------------------------------------------
//
// MODELLING WITH FDM IN 1D
//
//  * ACOUSTIC LOSSY
//  * 2ND ORDER WAVE EQUATION (PRESSURE):
//
//            [Ref Bilbao 2009, p177 eq 7-28]
//            d2P /dt2 = coef1 [d2P/dz2 - sigma_0 dP/dt + sigma_l d3P/dtz2]
//
//            with
//
//            coef1 = Vp^2 at pr location
//            sigma_0 = freq. independent lossy term
//            sigma_l = freq. dependent lossy term
//
// * FD DISCRETISATION
//
//            Dtt P = coef1 [Dzz - sigma_0 Dt. + sigma_l Dt- Dzz] P
//
//            Where Dtt, Dzz = 2nd derivative centered FD
//                  Dt.      = 1st derivative centered FD
//                  Dt-      = 1st derivative backward FD
//
//            Remark: if sigma_l Dt. Dzz P is used instead, the system becomes implicit
//
// PARENT CLASS: Scheme
//   DERIVED CLASS: FDM
//     DERIVED CLASS: FDM_1D
//       DERIVED CLASS: FDM_1D_ac_lossy
//         DERIVED CLASS: FDM_1D_ac_lossy_2nd
//
//-------------------------------------------------------------------------------------------------------

#include "FDM_1D_ac_lossy_2nd.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>

#include "mpi.h"

#include "acquisition.h"
#include "data.h"
#include "data_std.h"
#include "grid.h"
#include "grid_2D_complex.h"
#include "model.h"
#include "output_report.h"
#include "singleton.h"
#include "type_def.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------
FDM_1D_ac_lossy_2nd::FDM_1D_ac_lossy_2nd(void) : FDM_1D_ac_lossy()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D_ac_lossy_2nd::FDM_1D_ac_lossy_2nd");
	eq_order = ORDER_2ND ;
	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D_ac_lossy_1st::FDM_1D_ac_lossy_2nd");
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FDM_1D_ac_lossy_2nd::update_physical_coef(Myint iz, Myint ix, Myint iy, Myfloat perturb)  
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D_ac_lossy_2nd::update_physical_coef");

	// retrieve vp
	Myfloat vp = sqrt(coef1->pArray[iz])/ dt ;

	// add perturbation
	vp += perturb ;

	// update physical parameter
	coef1->pArray[iz] = vp * vp * dt * dt ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D_ac_lossy_2nd::update_physical_coef");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FDM_1D_ac_lossy_2nd::initialize(Model* pModel, Myfloat fmax)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D_ac_lossy_2nd::initialize");

	// Force along z not yet available
	if (src_stype == FORCE_Z)
	{
		print_error(" Source force along z not yet available for 2nd order wave eq.");
		return(RTN_CODE_KO) ;
	}

	// call parent initialiation
	Rtn_code rtn_code = FDM_1D_ac_lossy::initialize(pModel, fmax) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	//-------------------------------------------------------------------------------------------------------
	// unknown vectors
	//-------------------------------------------------------------------------------------------------------

	// pressure component (time domain)
	Variable* pVarPrc = Singleton::Instance()->register_variable(PRC, "prc") ;
	pVarPrc->allocate_grid(izEnd+1, dz) ;

	Variable* pVarPrn = Singleton::Instance()->register_variable(PRN, "pr") ;
	pVarPrn->allocate_grid(izEnd+1, dz) ;

	Variable* pVarPrp = Singleton::Instance()->register_variable(PRP, "prp") ;
	pVarPrp->allocate_grid(izEnd+1, dz) ;

	//-------------------------------------------------------------------------------------------------------
	// initialize eq. coef. in the medium
	//-------------------------------------------------------------------------------------------------------

	// check model type is GRID and subtype is REGULAR
	if (pModel->get_type() != GRID)
	{
		print_error("IN FDM_1D_ac_lossy_2nd::initialize --> model type is not GRID");
		return(RTN_CODE_KO) ;
	}
	if (pModel->get_sub_type() != REGULAR)
	{
		print_error("IN FDM_1D_ac_lossy_2nd::initialize --> model subtype is not REGULAR");
		return(RTN_CODE_KO) ;
	}

	Myfloat rho = TEMP_RHO_CONST ;

	// coef_i = vp^2 * dt^2
	coef1 = new Grid_1D_float(izEnd+1, dz) ;

	// get vp model
	Variable* vp_var = pModel->get_parameter(VP) ;
	if (vp_var == NULL)
	{
		print_error("IN FDM_1D_ac_lossy_2nd::initialize --> VP model not found");
		return(RTN_CODE_KO) ;
	}
	// get vp grid
	Grid_1D_float* vp_grid = dynamic_cast<Grid_1D_float*>(vp_var->get_grid()) ;
	if (vp_grid == NULL)
	{
		print_error("IN FDM_1D_ac_lossy_2nd::initialize --> VP grid not initialized");
		return(RTN_CODE_KO) ;
	}

	Myfloat* vp = vp_grid->pArray ;

	for (Myint iz = 0; iz < nz; iz++)
	{
		coef1->pArray[iz+izBeg2] = vp[iz] * vp[iz] * dt * dt ;
	}

	print_debug(ALL, MID_DEBUG, "end coef compute");

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
			coef1->pArray[iz] = vp_rand * vp_rand * dt * dt ;
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
			coef1->pArray[iz] = vp_rand * vp_rand * dt * dt ;
		}
	}

	//---------------------------------------------------------------------------------------------------
	// sponge
	//---------------------------------------------------------------------------------------------------

	// if (layer_type == SPONGE)
	//   {

	//     sponge_coef_pr_z = new Grid_1D_float(izEnd+1, dz) ;

	//     Myfloat vp_tmp = vp[0] ;
	//     for (Myint iz = izBeg; iz <= izEnd; iz++)
	// 	{
	// 	  if ((iz >= izBeg1) && (iz < izBeg2))
	// 	    {
	// 	      sponge_coef_pr_z->pArray[iz] = exp(-pow(Myfloat(izBeg2 - iz) * sponge_coef, 2)) ;
	// 	      coef1->pArray[iz] = vp_tmp * vp_tmp * dt * dt ;

	// 	      print_info(ALL, "sponge_coef_pr_z[iz]",sponge_coef_pr_z->pArray[iz]) ;
	// 	    }
	// 	  else
	// 	    {
	// 	      sponge_coef_pr_z->pArray[iz] = 1. ;
	// 	    }
	// 	}
	//   }

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
			coef1->pArray[iz] = vp_tmp * vp_tmp * dt * dt ;

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
			coef1->pArray[iz] = vp_tmp * vp_tmp * dt * dt ;

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

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D_ac_lossy_2nd::initialize");
	return(RTN_CODE_OK) ;

} ;

//=======================================================================================================
//
// RESET MODELLING 
//
//=======================================================================================================

Rtn_code FDM_1D_ac_lossy_2nd::reset(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D_ac_lossy_2nd::reset");

	// call to parent class
	Rtn_code rtn_code = Scheme::reset() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// reset grids
	rtn_code = Singleton::Instance()->get_variable(PRN)->reset_grid() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
	rtn_code = Singleton::Instance()->get_variable(PRC)->reset_grid() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
	rtn_code = Singleton::Instance()->get_variable(PRP)->reset_grid() ;
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

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D_ac_lossy_2nd::reset");
	return(RTN_CODE_OK) ;

} 

//=======================================================================================================
//
// FINALIZE MODELLING 
//
//=======================================================================================================

Rtn_code FDM_1D_ac_lossy_2nd::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D_ac_lossy_2nd::finalize");

	//---------------------------------------------------------------------------------------------------
	// deallocate arrays
	//---------------------------------------------------------------------------------------------------

	// physical parameter
	delete(coef1) ;

	// delete components
	Singleton::Instance()->delete_variable(PRC) ;
	Singleton::Instance()->delete_variable(PRN) ;
	Singleton::Instance()->delete_variable(PRP) ;

	// call parent initialiation
	Rtn_code rtn_code = FDM_1D_ac_lossy::finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D_ac_lossy_2nd::finalize");
	return(RTN_CODE_OK) ;

} ;

//=======================================================================================================
//
// PERFORM FORWARD MODELLING FOR ONE SOURCE
//
//=======================================================================================================


Rtn_code FDM_1D_ac_lossy_2nd::solve_current_shot(Acquisition* pAcquisition, Data *pData, Wavefield_type wtype, Snapshot* pSnapshot)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D_ac_lossy_2nd::solve_current_shot");

	// reset memory
	Rtn_code rtn_code = reset() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// locate src and rec in the grid
	rtn_code = locate_src_and_rec_in_grid(pAcquisition) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// source scaling factor
	Myfloat src_factor = dt * dt / dz ;

	if (Singleton::Instance()->dryrun)
	{
		print_info(ALL, "\n*** DRY RUN - SKIP TIME LOOP ***") ;
	}
	else
	{

		// keep original grid index
		izBeg1_initial = izBeg1 ;
		izBeg2_initial = izBeg2 ;

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

			// update string length
			update_string_length(it, prn, prc) ;

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
				rtn_code = source_excitation(prn, it-1, wtype, pData, src_factor) ;
			}
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

			// compute energy
			compute_energy(prn) ;

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
			// PRN -> PRC
			// PRC -> PRP
			// PRP -> PRN (will be overwritten)
			//-------------------------------------------------------------------------------------------------------
			Singleton::Instance()->swap_variable(PRC, PRN) ;
			Singleton::Instance()->swap_variable(PRN, PRP) ;

			// display remaining computation time
			display_remaining_computation_time(it, time_in_kernel) ;

		}

		//-------------------------------------------------------------------------------------------------------
		// end loop over time steps
		//-------------------------------------------------------------------------------------------------------

		// reset initial string length
		izBeg1 = izBeg1_initial ;
		izBeg2 = izBeg2_initial ;

	}

	// free src and rec in the grid
	rtn_code = free_position_arrays() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D_ac_lossy_2nd::solve_current_shot");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
void FDM_1D_ac_lossy_2nd::update_string_length(Myint it, Myfloat *prn, Myfloat *prc)
{  
	print_debug(ALL, FULL_DEBUG, "IN FDM_1D_ac_lossy_2nd::update_string_length");

	if (string_length_vs_time_step != NULL)
	{

		// compute new izBeg1
		izBeg1 = izBeg1_initial + string_length_vs_time_step[it] / dz ;
		if (izBeg1 < izBeg1_initial) izBeg1 = izBeg1_initial ;
		if (izBeg1 > (izEnd2 - 1)) izBeg1 = izEnd2 - 1 ;
		izBeg2 = izBeg1 ;

		// initialize pressure to zero between izBeg1_initial and izBeg1
		for (Myint iz = izBeg1_initial; iz <= izBeg1; iz++)
		{
			prn[iz] = 0.0 ;
			prc[iz] = 0.0 ;
		}
	}

	print_debug(ALL, FULL_DEBUG, "IN FDM_1D_ac_lossy_2nd::update_string_length");
	return ;
}

//-------------------------------------------------------------------------------------------------------
void FDM_1D_ac_lossy_2nd::init_eigen()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D_ac_lossy_2nd::init_eigen");

	// initialize PRP at t = -2 * dt
	Grid_1D_float *prn_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(PRP))->get_grid() ;
	Myfloat * const prn = prn_grid->pArray ;

	// get number mode number
	eigen_nmode = (Myint) Singleton::Instance()->pProgram->pModelling->get_param() ;

	Myfloat time_prn = -2.0*dt ;
	Myint izmin = izBeg2 ;
	Myint izmax = izEnd2 ;
	Myint nz = izEnd2 - izBeg2 ;
	for (Myint iz=izmin; iz<izmax; iz++)
	{
		Myfloat zz = float(iz-izmin)/float(nz-1) ;
		prn[iz] = -sin(M_PI*zz*eigen_nmode) * sin(M_PI*time_prn*eigen_nmode) ;
	}
	Myint iz0 = izBeg2 ;
	prn[iz0] = 0. ;
	for (Myint ii = 1; ii <= lstencil; ii++)
	{
		prn[iz0-ii] = -prn[iz0+ii] ;
	}
	iz0 = izEnd2-1 ;
	prn[iz0] = 0. ;
	for (Myint ii = 1; ii <= lstencil; ii++)
	{
		prn[iz0+ii] = -prn[iz0-ii] ;
	}

	// initialize PRC at t = -dt
	Grid_1D_float *prc_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(PRC))->get_grid() ;
	Myfloat * const prc = prc_grid->pArray ;
	Myfloat time_prc = -dt ;
	izmin = izBeg2 ;
	izmax = izEnd2 ;
	nz = izEnd2 - izBeg2 ;
	for (Myint iz=izmin; iz<izmax; iz++)
	{
		Myfloat zz = float(iz-izmin)/float(nz-1) ;
		prc[iz] = -sin(M_PI*zz*eigen_nmode) * sin(M_PI*time_prc*eigen_nmode) ;
	}
	iz0 = izBeg2 ;
	prc[iz0] = 0. ;
	for (Myint ii = 1; ii <= lstencil; ii++)
	{
		prc[iz0-ii] = -prc[iz0+ii] ;
	}
	iz0 = izEnd2-1 ;
	prc[iz0] = 0. ;
	for (Myint ii = 1; ii <= lstencil; ii++)
	{
		prc[iz0+ii] = -prc[iz0-ii] ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D_ac_lossy_2nd::init_eigen");
	return ;
} ;

} // namespace django

