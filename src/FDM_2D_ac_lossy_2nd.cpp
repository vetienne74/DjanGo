//------------------------------------------------------------------------------------------------------
//
// MODELLING WITH FDM IN 2D
//
//  * ACOUSTIC LOSSY
//  * 2ND ORDER WAVE EQUATION (PRESSURE):
//
//            [Ref Bilbao, p177 eq 7-28]
//            d2P /dt2 = coef1 [d2P/dx2 + d2P/dz2 - sigma_0 dP/dt
//                       + sigma_l d/dt(d2P/dx2 + d2P/dz2)]
//
//            with
//
//            coef1 = Vp^2 at pr location
//            sigma_0 = freq. independent lossy term
//            sigma_l = freq. dependent lossy term
//
// * FD DISCRETISATION
//
//            Dtt P = coef1 [(Dxx + Dzz) - sigma_0 Dt. + sigma_l Dt- (Dxx + Dzz)] P
//
//            Where Dtt, Dzz, Dxx = 2nd derivative centered FD
//                  Dt.           = 1st derivative centered FD
//                  Dt-           = 1st derivative backward FD
//
//            Remark: if 2 sigma_l Dt. Dzz P is used instead, the system becomes implicit
//
// PARENT CLASS: Scheme
//   DERIVED CLASS: FDM
//     DERIVED CLASS FDM_1D
//       DERIVED CLASS: FDM_2D
//         DERIVED CLASS: FDM_2D_ac_lossy
//           DERIVED CLASS: FDM_2D_ac_lossy_2nd
//
//-------------------------------------------------------------------------------------------------------

#include "FDM_2D_ac_lossy_2nd.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>

#include "mpi.h"

#include "acquisition.h"
#include "data.h"
#include "data_std.h"
#include "grid.h"
#include "grid_3D_complex.h"
#include "model.h"
#include "output_report.h"
#include "singleton.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------
FDM_2D_ac_lossy_2nd::FDM_2D_ac_lossy_2nd(void) : FDM_2D_ac_lossy()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D_ac_lossy_2nd::FDM_2D_ac_lossy_2nd");
	eq_order = ORDER_2ND ;
	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D_ac_lossy_2nd::FDM_2D_ac_lossy_2nd");
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FDM_2D_ac_lossy_2nd::update_physical_coef(Myint iz, Myint ix, Myint iy, Myfloat perturb)  
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D_ac_lossy_2nd::update_physical_coef");

	// retrieve vp
	Myfloat vp = sqrt(coef1->pArray[ix][iz])/ dt ;

	// add perturbation
	vp += perturb ;

	// update physical parameter
	coef1->pArray[ix][iz] = vp * vp * dt * dt ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D_ac_lossy_2nd::update_physical_coef");
	return(RTN_CODE_OK) ;
}

//=======================================================================================================
//
// INITIALIZE MODELLING
//
//=======================================================================================================

Rtn_code FDM_2D_ac_lossy_2nd::initialize(Model* pModel, Myfloat fmax)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D_ac_lossy_2nd::initialize");

	// Force along z not yet available
	if (src_stype == FORCE_Z)
	{
		print_error(" Source force along z not yet available for 2nd order wave eq.");
		return(RTN_CODE_KO) ;
	}

	// call parent initialiation
	Rtn_code rtn_code = FDM_2D_ac_lossy::initialize(pModel, fmax) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	//-------------------------------------------------------------------------------------------------------
	// unknown vectors
	//-------------------------------------------------------------------------------------------------------

	// pressure component (time domain)
	Variable* pVarPrc = Singleton::Instance()->register_variable(PRC, "prc") ;
	pVarPrc->allocate_grid(izEnd+1, ixEnd+1, dz, dx) ;

	Variable* pVarPrn = Singleton::Instance()->register_variable(PRN, "pr") ;
	pVarPrn->allocate_grid(izEnd+1, ixEnd+1, dz, dx) ;

	Variable* pVarPrp = Singleton::Instance()->register_variable(PRP, "prp") ;
	pVarPrp->allocate_grid(izEnd+1, ixEnd+1, dz, dx) ;

	//-------------------------------------------------------------------------------------------------------
	// initialize eq. coef. in the medium
	//-------------------------------------------------------------------------------------------------------

	// check model type is GRID and subtype is REGULAR
	if (pModel->get_type() != GRID)
	{
		print_error("IN FDM_2D_ac_lossy_2nd::initialize --> model type is not GRID");
		return(RTN_CODE_KO) ;
	}
	if (pModel->get_sub_type() != REGULAR)
	{
		print_error("IN FDM_2D_ac_lossy_2nd::initialize --> model subtype is not REGULAR");
		return(RTN_CODE_KO) ;
	}

	Myfloat rho = TEMP_RHO_CONST ;

	// coef_i = vp^2 * dt^2

	coef1 = new Grid_2D_float(izEnd+1, ixEnd+1, dz, dx) ;

	// get vp model
	Variable* vp_var = pModel->get_parameter(VP) ;
	if (vp_var == NULL)
	{
		print_error("IN FDM_2D_ac_lossy_2nd::initialize --> VP model not found");
		return(RTN_CODE_KO) ;
	}
	// get vp grid
	Grid_2D_float* vp_grid = dynamic_cast<Grid_2D_float*>(vp_var->get_grid()) ;
	if (vp_grid == NULL)
	{
		print_error("IN FDM_2D_ac_lossy_2nd::initialize --> VP grid not Grid_2D_float*");
		return(RTN_CODE_KO) ;
	}

	Myfloat** vp = vp_grid->pArray ;

	// loop on the grid points in medium
#pragma omp parallel for
	for (Myint ix = 0; ix < nx; ix++)
	{
#pragma ivdep
		for (Myint iz = 0; iz < nz; iz++)
		{
			// update coef in the FDM grid point
			coef1->pArray[ix+ixBeg2][iz+izBeg2] = vp[ix][iz] * vp[ix][iz] * dt * dt ;
		}
	}

	print_debug(ALL, MID_DEBUG, "end coef compute");

	//-------------------------------------------------------------------------------------------------------
	// initialize eq. coef. in the surrounding boundary
	//-------------------------------------------------------------------------------------------------------

	//---------------------------------------------------------------------------------------------------
	// sponge
	//---------------------------------------------------------------------------------------------------

	if (is_there_boundary_type(SPG2))
	{
		// loop on all grid points
		for (Myint ix = ixBeg1; ix < ixEnd1; ix++)
		{
			for (Myint iz = izBeg1; iz < izEnd1; iz++)
			{
				Myint ix2 = ix - izBeg2 ;
				ix2 = max(ix2, 0) ;
				ix2 = min(ix2, nx-1) ;

				Myint iz2 = iz - izBeg2 ;
				iz2 = max(iz2, 0) ;
				iz2 = min(iz2, nz-1) ;

				Myfloat vp_tmp = vp[ix2][iz2] ;
				coef1->pArray[ix][iz] = dt * dt * vp_tmp * vp_tmp * rho ;
			}
		}
	}

	//---------------------------------------------------------------------------------------------------
	// random layer
	//---------------------------------------------------------------------------------------------------

	if (get_boundary_type(ZBEG) == RANDOM)
	{
		Myfloat ratio, vp_rand, vp_tmp ;

		for (Myint ix = izBeg1; ix < izEnd1; ix++)
		{
			Myint ix2 = ix - izBeg2 ;
			ix2 = max(ix2, 0) ;
			ix2 = min(ix2, nx-1) ;

			vp_tmp = vp[ix2][0] ;
			for (Myint iz = izBeg1; iz < izBeg2; iz++)
			{
				ratio = (float) (izBeg2-iz)/nlayer_zBeg ;
				vp_rand = vp_tmp - (rand() % (int) (RAND_BOUND_COEF1 * vp_tmp)) * ratio ;
				coef1->pArray[ix][iz] = vp_rand * vp_rand * dt * dt ;
			}
		}
	}

	if (get_boundary_type(ZEND) == RANDOM)
	{
		Myfloat ratio, vp_rand, vp_tmp ;

		for (Myint ix = izBeg1; ix < izEnd1; ix++)
		{
			Myint ix2 = ix - izBeg2 ;
			ix2 = max(ix2, 0) ;
			ix2 = min(ix2, nx-1) ;

			vp_tmp = vp[ix2][nz-1] ;
			for (Myint iz = izEnd2; iz < izEnd1; iz++)
			{
				ratio = (float) (iz-(izEnd2-1))/nlayer_zEnd ;
				vp_rand = vp_tmp - (rand() % (int) (RAND_BOUND_COEF1 * vp_tmp)) * ratio ;
				coef1->pArray[ix][iz] = vp_rand * vp_rand * dt * dt ;
			}
		}
	}

	if (get_boundary_type(XBEG) == RANDOM)
	{
		Myfloat ratio, vp_rand, vp_tmp ;

		for (Myint iz = izBeg1; iz < izEnd1; iz++)
		{
			Myint iz2 = iz - izBeg2 ;
			iz2 = max(iz2, 0) ;
			iz2 = min(iz2, nz-1) ;

			vp_tmp = vp[0][iz2] ;
			for (Myint ix = ixBeg1; ix < izBeg2; ix++)
			{
				ratio = (float) (ixBeg2-ix)/nlayer_xBeg ;
				vp_rand = vp_tmp - (rand() % (int) (RAND_BOUND_COEF1 * vp_tmp)) * ratio ;
				coef1->pArray[ix][iz] = vp_rand * vp_rand * dt * dt ;
			}
		}
	}

	if (get_boundary_type(XEND) == RANDOM)
	{
		Myfloat ratio, vp_rand, vp_tmp ;

		for (Myint iz = izBeg1; iz < izEnd1; iz++)
		{
			Myint iz2 = iz - izBeg2 ;
			iz2 = max(iz2, 0) ;
			iz2 = min(iz2, nz-1) ;

			vp_tmp = vp[nx-1][iz2] ;
			for (Myint ix = ixEnd2; ix < ixEnd1; ix++)
			{
				ratio = (float) (ix-(ixEnd2-1))/nlayer_xEnd ;
				vp_rand = vp_tmp - (rand() % (int) (RAND_BOUND_COEF1 * vp_tmp)) * ratio ;
				coef1->pArray[ix][iz] = vp_rand * vp_rand * dt * dt ;
			}
		}
	}

	// // output grid
	// ofstream coef1_file ;
	// coef1_file.open("grid_coef1.out", ios::binary) ;
	// assert(coef1_file.is_open());

	// for (Myint ix=ixBeg; ix < ixEnd; ix++)
	//   {
	//     for (Myint iz=izBeg; iz < izEnd; iz++)
	// 	{
	// 	  coef1_file.write((char*) &(coef1->pArray[ix][iz]), sizeof(Myfloat)) ;
	// 	}
	//   }
	// coef1_file.close() ;

	//---------------------------------------------------------------------------------------------------
	// CPML
	//---------------------------------------------------------------------------------------------------

	if (get_boundary_type(ZBEG) == PML)
	{

		// compute some coefficients
		Myfloat cpml_alpha_max = CPML_FREQ * PI ;
		Myfloat cpml_pow       = CPML_POW ;
		Myfloat cpml_vmax      = (vp_grid)->get_max() ;
		Myfloat cpml_rcoef     = get_boundary_coef(ZBEG) ;

		// allocate CPML arrays
		Myint nzl = nlayer_zBeg+2*lstencil ;
		Myint nxl = ixEnd+1 ;
		apml_zBeg         = new Grid_2D_float(nzl, nxl, dz, dx) ;
		bpml_zBeg         = new Grid_2D_float(nzl, nxl, dz, dx) ;
		apml_half_zBeg    = new Grid_2D_float(nzl, nxl, dz, dx) ;
		bpml_half_zBeg    = new Grid_2D_float(nzl, nxl, dz, dx) ;
		mem_vz_zBeg       = new Grid_2D_float(nzl, nxl, dz, dx) ;
		mem_pr_zBeg       = new Grid_2D_float(nzl, nxl, dz, dx) ;

		Myfloat d0     = -(cpml_pow+1.)*cpml_vmax*log10(cpml_rcoef) / (2.*nlayer_zBeg*dz) ;

		// fill CPML arrays
		Myint iz1 = izBeg1 ;
		Myint iz2 = izBeg2 ;
		Myint ix1 = ixBeg1 ;
		Myint ix2 = ixEnd1 ;
#pragma omp parallel for
		for (Myint ix = ix1; ix < ix2; ix++)
		{
			Myint ipmlx = ix ;
#pragma ivdep       
			for (Myint iz = iz1; iz < iz2; iz++)
			{
				Myint ipmlz = iz ;

				// pressure node
				Myint ix2 = ix - ixBeg2 ;
				ix2 = max(ix2, 0) ;
				ix2 = min(ix2, nx-1) ;
				Myfloat vp_tmp = vp[ix2][0] ;
				coef1->pArray[ix][iz] = vp_tmp * vp_tmp * dt * dt ;

				Myfloat dnorm = Myfloat(izBeg2 - iz)/Myfloat(nlayer_zBeg) ;
				if (dnorm < 0) {dnorm = 0. ;}
				Myfloat alpha = cpml_alpha_max * (1. - dnorm) ;
				Myfloat dd = d0 * pow(dnorm, cpml_pow) ;

				bpml_zBeg->pArray[ipmlx][ipmlz] = exp(-(dd+alpha)*dt) ;
				apml_zBeg->pArray[ipmlx][ipmlz] = dd * (bpml_zBeg->pArray[ipmlx][ipmlz]-1.) / (dd+alpha) ;

				// velocity node
				dnorm = Myfloat(izBeg2 - iz -0.5)/Myfloat(nlayer_zBeg) ;
				if (dnorm < 0) {dnorm = 0. ;}
				alpha = cpml_alpha_max * (1. - dnorm) ;
				dd = d0 * pow(dnorm, cpml_pow) ;

				bpml_half_zBeg->pArray[ipmlx][ipmlz] = exp(-(dd+alpha)*dt) ;
				apml_half_zBeg->pArray[ipmlx][ipmlz] = dd * (bpml_half_zBeg->pArray[ipmlx][ipmlz]-1.) / (dd+alpha) ;
			}
		}
	}

	// z+ layer
	if (get_boundary_type(ZEND) == PML)
	{
		// compute some coefficients
		Myfloat cpml_alpha_max = CPML_FREQ * PI ;
		Myfloat cpml_pow       = CPML_POW ;
		Myfloat cpml_vmax      = (vp_grid)->get_max() ;
		Myfloat cpml_rcoef     = get_boundary_coef(ZEND) ;

		// allocate CPML arrays
		Myint nzl = nlayer_zEnd+2*lstencil ;
		Myint nxl = ixEnd+1 ;
		apml_zEnd         = new Grid_2D_float(nzl, nxl, dz, dx) ;
		bpml_zEnd         = new Grid_2D_float(nzl, nxl, dz, dx) ;
		apml_half_zEnd    = new Grid_2D_float(nzl, nxl, dz, dx) ;
		bpml_half_zEnd    = new Grid_2D_float(nzl, nxl, dz, dx) ;
		mem_vz_zEnd       = new Grid_2D_float(nzl, nxl, dz, dx) ;
		mem_pr_zEnd       = new Grid_2D_float(nzl, nxl, dz, dx) ;

		Myfloat d0     = -(cpml_pow+1.)*cpml_vmax*log10(cpml_rcoef) / (2.*nlayer_zEnd*dz) ;

		// fill CPML arrays
		Myint iz1 = izEnd2 ;
		Myint iz2 = izEnd1 ;
		Myint ix1 = ixBeg1 ;
		Myint ix2 = ixEnd1 ;
		Myint const offset_z = nlayer_zBeg + nz ;
#pragma omp parallel for
		for (Myint ix = ix1; ix < ix2; ix++)
		{
			Myint ipmlx = ix ;
#pragma ivdep       
			for (Myint iz = iz1; iz < iz2; iz++)
			{
				Myint ipmlz = iz - offset_z ;

				// pressure node
				Myint ix2 = ix - ixBeg2 ;
				ix2 = max(ix2, 0) ;
				ix2 = min(ix2, nx-1) ;
				Myfloat vp_tmp = vp[ix2][nz-1] ;
				coef1->pArray[ix][iz] = vp_tmp * vp_tmp * dt * dt ;

				Myfloat dnorm = Myfloat(iz - izEnd2 + 1)/Myfloat(nlayer_zEnd) ;
				if (dnorm < 0) {dnorm = 0. ;}
				Myfloat alpha = cpml_alpha_max * (1. - dnorm) ;
				Myfloat dd = d0 * pow(dnorm, cpml_pow) ;

				bpml_zEnd->pArray[ipmlx][ipmlz] = exp(-(dd+alpha)*dt) ;
				apml_zEnd->pArray[ipmlx][ipmlz] = dd * (bpml_zEnd->pArray[ipmlx][ipmlz]-1.) / (dd+alpha) ;

				// velocity node
				dnorm = Myfloat(iz - izEnd2 + 0.5)/Myfloat(nlayer_zEnd) ;
				if (dnorm < 0) {dnorm = 0. ;}
				alpha = cpml_alpha_max * (1. - dnorm) ;
				dd = d0 * pow(dnorm, cpml_pow) ;

				bpml_half_zEnd->pArray[ipmlx][ipmlz] = exp(-(dd+alpha)*dt) ;
				apml_half_zEnd->pArray[ipmlx][ipmlz] = dd * (bpml_half_zEnd->pArray[ipmlx][ipmlz]-1.) / (dd+alpha) ;
			}
		}
	}

	// x- layer
	if (get_boundary_type(XBEG) == PML)
	{
		// compute some coefficients
		Myfloat cpml_alpha_max = CPML_FREQ * PI ;
		Myfloat cpml_pow       = CPML_POW ;
		Myfloat cpml_vmax      = (vp_grid)->get_max() ;
		Myfloat cpml_rcoef     = get_boundary_coef(XBEG) ;

		// allocate CPML arrays
		Myint nzl = izEnd+1 ;
		Myint nxl = nlayer_xBeg+2*lstencil ;
		apml_xBeg         = new Grid_2D_float(nzl, nxl, dz, dx) ;
		bpml_xBeg         = new Grid_2D_float(nzl, nxl, dz, dx) ;
		apml_half_xBeg    = new Grid_2D_float(nzl, nxl, dz, dx) ;
		bpml_half_xBeg    = new Grid_2D_float(nzl, nxl, dz, dx) ;
		mem_vx_xBeg       = new Grid_2D_float(nzl, nxl, dz, dx) ;
		mem_pr_xBeg       = new Grid_2D_float(nzl, nxl, dz, dx) ;

		Myfloat d0     = -(cpml_pow+1.)*cpml_vmax*log10(cpml_rcoef) / (2.*nlayer_xBeg*dx) ;

		// fill CPML arrays
		Myint iz1 = izBeg1 ;
		Myint iz2 = izEnd1 ;
		Myint ix1 = ixBeg1 ;
		Myint ix2 = ixBeg2 ;
#pragma omp parallel for
		for (Myint ix = ix1; ix < ix2; ix++)
		{
			Myint ipmlx = ix ;
#pragma ivdep       
			for (Myint iz = iz1; iz < iz2; iz++)
			{
				Myint ipmlz = iz ;

				// pressure node
				Myint iz2 = iz - izBeg2 ;
				iz2 = max(iz2, 0) ;
				iz2 = min(iz2, nz-1) ;
				Myfloat vp_tmp = vp[0][iz2] ;
				coef1->pArray[ix][iz] = vp_tmp * vp_tmp * dt * dt ;

				Myfloat dnorm = Myfloat(ixBeg2 - ix)/Myfloat(nlayer_xBeg) ;
				if (dnorm < 0) {dnorm = 0. ;}
				Myfloat alpha = cpml_alpha_max * (1. - dnorm) ;
				Myfloat dd = d0 * pow(dnorm, cpml_pow) ;

				bpml_xBeg->pArray[ipmlx][ipmlz] = exp(-(dd+alpha)*dt) ;
				apml_xBeg->pArray[ipmlx][ipmlz] = dd * (bpml_xBeg->pArray[ipmlx][ipmlz]-1.) / (dd+alpha) ;

				// velocity node
				dnorm = Myfloat(ixBeg2 - ix -0.5)/Myfloat(nlayer_xBeg) ;
				if (dnorm < 0) {dnorm = 0. ;}
				alpha = cpml_alpha_max * (1. - dnorm) ;
				dd = d0 * pow(dnorm, cpml_pow) ;

				bpml_half_xBeg->pArray[ipmlx][ipmlz] = exp(-(dd+alpha)*dt) ;
				apml_half_xBeg->pArray[ipmlx][ipmlz] = dd * (bpml_half_xBeg->pArray[ipmlx][ipmlz]-1.) / (dd+alpha) ;
			}
		}
	}

	// x+ layer
	if (get_boundary_type(XEND) == PML)
	{
		// compute some coefficients
		Myfloat cpml_alpha_max = CPML_FREQ * PI ;
		Myfloat cpml_pow       = CPML_POW ;
		Myfloat cpml_vmax      = (vp_grid)->get_max() ;
		Myfloat cpml_rcoef     = get_boundary_coef(XEND) ;

		// allocate CPML arrays
		Myint nzl = izEnd+1 ;
		Myint nxl = nlayer_xEnd+2*lstencil ;
		apml_xEnd         = new Grid_2D_float(nzl, nxl, dz, dx) ;
		bpml_xEnd         = new Grid_2D_float(nzl, nxl, dz, dx) ;
		apml_half_xEnd    = new Grid_2D_float(nzl, nxl, dz, dx) ;
		bpml_half_xEnd    = new Grid_2D_float(nzl, nxl, dz, dx) ;
		mem_vx_xEnd       = new Grid_2D_float(nzl, nxl, dz, dx) ;
		mem_pr_xEnd       = new Grid_2D_float(nzl, nxl, dz, dx) ;

		Myfloat d0     = -(cpml_pow+1.)*cpml_vmax*log10(cpml_rcoef) / (2.*nlayer_xEnd*dx) ;

		// fill CPML arrays
		Myint iz1 = izBeg1 ;
		Myint iz2 = izEnd1 ;
		Myint ix1 = ixEnd2 ;
		Myint ix2 = ixEnd1 ;
		Myint const offset_x = nlayer_xBeg + nx ;
#pragma omp parallel for
		for (Myint ix = ix1; ix < ix2; ix++)
		{
			Myint ipmlx = ix - offset_x ;
#pragma ivdep       
			for (Myint iz = iz1; iz < iz2; iz++)
			{
				Myint ipmlz = iz ;

				// pressure node
				Myint iz2 = iz - izBeg2 ;
				iz2 = max(iz2, 0) ;
				iz2 = min(iz2, nz-1) ;
				Myfloat vp_tmp = vp[nx-1][iz2] ;
				coef1->pArray[ix][iz] = vp_tmp * vp_tmp * dt * dt ;

				Myfloat dnorm = Myfloat(ix - ixEnd2 + 1)/Myfloat(nlayer_xEnd) ;
				if (dnorm < 0) {dnorm = 0. ;}
				Myfloat alpha = cpml_alpha_max * (1. - dnorm) ;
				Myfloat dd = d0 * pow(dnorm, cpml_pow) ;

				bpml_xEnd->pArray[ipmlx][ipmlz] = exp(-(dd+alpha)*dt) ;
				apml_xEnd->pArray[ipmlx][ipmlz] = dd * (bpml_xEnd->pArray[ipmlx][ipmlz]-1.) / (dd+alpha) ;

				// velocity node
				dnorm = Myfloat(ix - ixEnd2 + 0.5)/Myfloat(nlayer_xEnd) ;
				if (dnorm < 0) {dnorm = 0. ;}
				alpha = cpml_alpha_max * (1. - dnorm) ;
				dd = d0 * pow(dnorm, cpml_pow) ;

				bpml_half_xEnd->pArray[ipmlx][ipmlz] = exp(-(dd+alpha)*dt) ;
				apml_half_xEnd->pArray[ipmlx][ipmlz] = dd * (bpml_half_xEnd->pArray[ipmlx][ipmlz]-1.) / (dd+alpha) ;
			}
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D_ac_lossy_2nd::initialize");
	return(RTN_CODE_OK) ;

} ;


//=======================================================================================================
//
// RESET MODELLING 
//
//=======================================================================================================

Rtn_code FDM_2D_ac_lossy_2nd::reset(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D_ac_lossy_2nd::reset");

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

	// initialize CPML memory variables
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
	if (get_boundary_type(XBEG) == PML)
	{
		mem_vx_xBeg->reset() ;
		mem_pr_xBeg->reset() ;
	}
	if (get_boundary_type(XEND) == PML)
	{
		mem_vx_xEnd->reset() ;
		mem_pr_xEnd->reset() ;
	}

	if (Singleton::Instance()->pProgram->pModelling->get_case() == MODELLING_EIGEN)
	{
		// initialization for eigen mode test
		init_eigen() ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D_ac_lossy_2nd::reset");
	return(RTN_CODE_OK) ;

} ;


//=======================================================================================================
//
// FINALIZE MODELLING 
//
//=======================================================================================================

Rtn_code FDM_2D_ac_lossy_2nd::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D_ac_lossy_2nd::finalize");

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
	Rtn_code rtn_code = FDM_2D_ac_lossy::finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D_ac_lossy_2nd::finalize");
	return(RTN_CODE_OK) ;

} ;


//=======================================================================================================
//
// PERFORM FORWARD MODELLING FOR ONE SOURCE
//
//=======================================================================================================


Rtn_code FDM_2D_ac_lossy_2nd::solve_current_shot(Acquisition* pAcquisition, Data *pData, Wavefield_type wtype, Snapshot* pSnapshot)
{

	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D_ac_lossy_2nd::solve_current_shot");

	// reset memory
	Rtn_code rtn_code = this->reset() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// locate src and rec in the grid
	rtn_code = this->locate_src_and_rec_in_grid(pAcquisition) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// locate snapshot pixel in the grid
	rtn_code = this->locate_pixel_in_grid(pSnapshot) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// source scaling factor
	Myfloat src_factor = dt * dt / (dx * dz) ;

	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

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
			Variable* prn_var = Singleton::Instance()->get_variable(PRN) ;
			Grid_2D_float *prn_grid = (Grid_2D_float*) prn_var->get_grid() ;
			Myfloat ** const prn = prn_grid->pArray ;

			//-------------------------------------------------------------------------------------------------------
			// computation of the pressure component
			//-------------------------------------------------------------------------------------------------------
			double t0 = MPI_Wtime() ;
			this->compute_pressure() ;
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

	}

	// free src and rec in the grid
	rtn_code = this->free_position_arrays() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;


	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D_ac_lossy_2nd::solve_current_shot");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
void FDM_2D_ac_lossy_2nd::init_eigen()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D_ac_lossy_2nd::init_eigen");

	{
		// initialize PRC at t = -dt

		Grid_2D_float *pr_grid = (Grid_2D_float*) (Singleton::Instance()->get_variable(PRC))->get_grid() ;
		Myfloat ** const pr = pr_grid->pArray ;

		Myfloat time_p = -dt ;
		Myint izmin = izBeg2 ;
		Myint izmax = izEnd2 ;
		Myint nz = izEnd2 - izBeg2 ;
		Myint ixmin = ixBeg2 ;
		Myint ixmax = ixEnd2 ;
		Myint nx = ixEnd2 - ixBeg2 ;

		for (Myint ix=ixmin; ix<ixmax; ix++)
		{
			for (Myint iz=izmin; iz<izmax; iz++)
			{
				Myfloat zz = float(iz-izmin)/float(nz-1) ;
				Myfloat xx = float(ix-ixmin)/float(nx-1) ;
				pr[ix][iz] = - sqrt(2.0) * sin(M_PI*xx) * sin(M_PI*zz) * sin(M_PI*sqrt(2.0)*time_p) ;
			}
		}

		{
			//--------------------------------
			// free surface and rigid surface
			//--------------------------------

			// loop on boundaries
			Myint ii_sign, pr_sign, x_axe, z_axe ;
			Myint ixmin, ixmax, izmin, izmax ;
			bool switch1, switch2 ;
			for (Myint ib=1; ib<=nBoundary; ib++)
			{
				switch1 = false ;
				switch2 = false ;

				//  set corect index
				if (pBoundary[ib-1]->get_edge() == ZBEG)
				{
					izmin   = izBeg2 ;
					izmax   = izBeg2 ;
					ixmin   = ixBeg ;
					ixmax   = ixEnd ;
					ii_sign = 1 ;
					x_axe   = 0 ;
					z_axe   = 1 ;
					switch1 = true ;
				}
				else if (pBoundary[ib-1]->get_edge() == ZEND)
				{
					izmin   = izEnd2-1 ;
					izmax   = izEnd2-1 ;
					ixmin   = ixBeg ;
					ixmax   = ixEnd ;
					ii_sign = -1 ;
					x_axe   = 0 ;
					z_axe   = 1 ;
					switch1 = true ;
				}
				else if (pBoundary[ib-1]->get_edge() == XBEG)
				{
					izmin   = izBeg ;
					izmax   = izEnd ;
					ixmin   = ixBeg2 ;
					ixmax   = ixBeg2 ;
					ii_sign = 1 ;
					x_axe   = 1 ;
					z_axe   = 0 ;
					switch1 = true ;
				}
				else if (pBoundary[ib-1]->get_edge() == XEND)
				{
					izmin   = izBeg ;
					izmax   = izEnd ;
					ixmin   = ixEnd2-1 ;
					ixmax   = ixEnd2-1 ;
					ii_sign = -1 ;
					x_axe   = 1 ;
					z_axe   = 0 ;
					switch1 = true ;
				}

				// free surface with image method
				// anti-symetry of pressure wavefield
				if (pBoundary[ib-1]->get_type() == FREESURF)
				{
					pr_sign = -1.0 ;
					switch2 = true ;
				}

				// free surface with image method
				// anti-symetry of pressure wavefield
				else if (pBoundary[ib-1]->get_type() == RIGID)
				{
					pr_sign = 1.0 ;
					switch2 = true ;
				}

				// apply boundary condition
				if (switch1 && switch2)
				{
#pragma omp parallel for
					for (Myint ix = ixmin; ix <= ixmax; ix++)
					{
						for (Myint iz = izmin; iz <= izmax; iz++)
						{
							if (pBoundary[ib-1]->get_type() == FREESURF) pr[ix][iz] = 0. ;
#pragma ivdep
							for (Myint ii = 1; ii <= lstencil; ii++)
							{
								pr[ix-x_axe*(ii_sign*ii)][iz-z_axe*(ii_sign*ii)]
														  = pr_sign * pr[ix+x_axe*(ii_sign*ii)][iz+z_axe*(ii_sign*ii)] ;
							}
						}
					}
				}
			}
		}
	}

	{
		// initialize PRP at t = -2 * dt

		Grid_2D_float *pr_grid = (Grid_2D_float*) (Singleton::Instance()->get_variable(PRP))->get_grid() ;
		Myfloat ** const pr = pr_grid->pArray ;

		Myfloat time_p = -2.0 * dt ;
		Myint izmin = izBeg2 ;
		Myint izmax = izEnd2 ;
		Myint nz = izEnd2 - izBeg2 ;
		Myint ixmin = ixBeg2 ;
		Myint ixmax = ixEnd2 ;
		Myint nx = ixEnd2 - ixBeg2 ;

		for (Myint ix=ixmin; ix<ixmax; ix++)
		{
			for (Myint iz=izmin; iz<izmax; iz++)
			{
				Myfloat zz = float(iz-izmin)/float(nz-1) ;
				Myfloat xx = float(ix-ixmin)/float(nx-1) ;
				pr[ix][iz] = - sqrt(2.0) * sin(M_PI*xx) * sin(M_PI*zz) * sin(M_PI*sqrt(2.0)*time_p) ;
			}
		}

		{
			//--------------------------------
			// free surface and rigid surface
			//--------------------------------

			// loop on boundaries
			Myint ii_sign, pr_sign, x_axe, z_axe ;
			Myint ixmin, ixmax, izmin, izmax ;
			bool switch1, switch2 ;
			for (Myint ib=1; ib<=nBoundary; ib++)
			{
				switch1 = false ;
				switch2 = false ;

				//  set corect index
				if (pBoundary[ib-1]->get_edge() == ZBEG)
				{
					izmin   = izBeg2 ;
					izmax   = izBeg2 ;
					ixmin   = ixBeg ;
					ixmax   = ixEnd ;
					ii_sign = 1 ;
					x_axe   = 0 ;
					z_axe   = 1 ;
					switch1 = true ;
				}
				else if (pBoundary[ib-1]->get_edge() == ZEND)
				{
					izmin   = izEnd2-1 ;
					izmax   = izEnd2-1 ;
					ixmin   = ixBeg ;
					ixmax   = ixEnd ;
					ii_sign = -1 ;
					x_axe   = 0 ;
					z_axe   = 1 ;
					switch1 = true ;
				}
				else if (pBoundary[ib-1]->get_edge() == XBEG)
				{
					izmin   = izBeg ;
					izmax   = izEnd ;
					ixmin   = ixBeg2 ;
					ixmax   = ixBeg2 ;
					ii_sign = 1 ;
					x_axe   = 1 ;
					z_axe   = 0 ;
					switch1 = true ;
				}
				else if (pBoundary[ib-1]->get_edge() == XEND)
				{
					izmin   = izBeg ;
					izmax   = izEnd ;
					ixmin   = ixEnd2-1 ;
					ixmax   = ixEnd2-1 ;
					ii_sign = -1 ;
					x_axe   = 1 ;
					z_axe   = 0 ;
					switch1 = true ;
				}

				// free surface with image method
				// anti-symetry of pressure wavefield
				if (pBoundary[ib-1]->get_type() == FREESURF)
				{
					pr_sign = -1.0 ;
					switch2 = true ;
				}

				// free surface with image method
				// anti-symetry of pressure wavefield
				else if (pBoundary[ib-1]->get_type() == RIGID)
				{
					pr_sign = 1.0 ;
					switch2 = true ;
				}

				// apply boundary condition
				if (switch1 && switch2)
				{
#pragma omp parallel for        
					for (Myint ix = ixmin; ix <= ixmax; ix++)
					{
						for (Myint iz = izmin; iz <= izmax; iz++)
						{
							if (pBoundary[ib-1]->get_type() == FREESURF) pr[ix][iz] = 0. ;
#pragma ivdep
							for (Myint ii = 1; ii <= lstencil; ii++)
							{
								pr[ix-x_axe*(ii_sign*ii)][iz-z_axe*(ii_sign*ii)]
														  = pr_sign * pr[ix+x_axe*(ii_sign*ii)][iz+z_axe*(ii_sign*ii)] ;
							}
						}
					}
				}
			}
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D_ac_lossy_2nd::init_eigen");
	return ;
} ;

} // namespace django
