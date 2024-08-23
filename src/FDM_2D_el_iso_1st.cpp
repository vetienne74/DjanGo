//-------------------------------------------------------------------------------------------------------
//
// MODELLING WITH FDM IN 2D
//
//  * ELASTIC ISOTROPIC MEDIA 
//  * 1ST ORDER WAVE EQUATION (STRESS / VELOCITY):
//
//            dvx/dt  = coef1 (dSxx/dx + dSxz/dz) 
//            dvz/dt  = coef2 (dSxz/dx + dSzz/dz)
//            dSxx/dt = coef3 dvx/dx + coef4 dvz/dz
//            dSzz/dt = coef3 dvz/dz + coef4 dvx/dx
//            dSxz/dt = coef5 (dvx/dz + dvz/dx)
//
//            with
//
//            coef1 = (1/rho) at vx location
//            coef2 = (1/rho) at vz location
//            coef3 = (lambda + 2mu) at sxx and szz location
//            coef4 = (lambda) at sxx and szz location
//            coef5 = (mu) at sxz location
//
// PARENT CLASS: Scheme
//   DERIVED CLASS: FDM
//     DERIVED CLASS FDM_1D
//       DERIVED CLASS: FDM_2D
//         DERIVED CLASS: FDM_2D_el_iso
//           DERIVED CLASS: FDM_2D_el_iso_1st
//
//-------------------------------------------------------------------------------------------------------

#include "FDM_2D_el_iso_1st.h"

#include <cassert>
#include <cmath>
#include <iostream>

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
FDM_2D_el_iso_1st::FDM_2D_el_iso_1st(void) : FDM_2D_el_iso()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D_el_iso_1st::FDM_2D_el_iso_1st");
	eq_order = ORDER_1ST ;
	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D_el_iso_1st::FDM_2D_el_iso_1st");
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FDM_2D_el_iso_1st::initialize(Model* pModel, Myfloat fmax)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D_el_iso_1st::initialize");

	// call parent initialiation
	Rtn_code rtn_code = FDM_2D_el_iso::initialize(pModel, fmax) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	//-------------------------------------------------------------------------------------------------------
	// unknown vectors
	//-------------------------------------------------------------------------------------------------------

	// velocity component (time domain)
	Variable* pVarVz = Singleton::Instance()->register_variable(VZ, "vz") ;
	pVarVz->allocate_grid(izEnd+1, ixEnd+1, dz, dx) ;
	Variable* pVarVx = Singleton::Instance()->register_variable(VX, "vx") ;
	pVarVx->allocate_grid(izEnd+1, ixEnd+1, dz, dx) ;

	// stress component (time domain)
	Variable* pVarSxx = Singleton::Instance()->register_variable(SXX, "sxx") ;
	pVarSxx->allocate_grid(izEnd+1, ixEnd+1, dz, dx) ;
	Variable* pVarSzz = Singleton::Instance()->register_variable(SZZ, "szz") ;
	pVarSzz->allocate_grid(izEnd+1, ixEnd+1, dz, dx) ;
	Variable* pVarSxz = Singleton::Instance()->register_variable(SXZ, "sxz") ;
	pVarSxz->allocate_grid(izEnd+1, ixEnd+1, dz, dx) ;

	//-------------------------------------------------------------------------------------------------------
	// initialize eq. coef. in the medium
	//-------------------------------------------------------------------------------------------------------

	// check model type is GRID and subtype is REGULAR
	if (pModel->get_type() != GRID)
	{
		print_error("IN FDM_2D_el_iso_1st::initialize --> model type is not GRID");
		return(RTN_CODE_KO) ;
	}
	if (pModel->get_sub_type() != REGULAR)
	{
		print_error("IN FDM_2D_el_iso_1st::initialize --> model subtype is not REGULAR");
		return(RTN_CODE_KO) ;
	}

	// get vp model, grid and pointer
	Variable* vp_var = pModel->get_parameter(VP) ;
	if (vp_var == NULL)
	{
		print_error("IN FDM_2D_el_iso_1st::initialize --> VP model not found");
		return(RTN_CODE_KO) ;
	}
	Grid_2D_float* vp_grid = dynamic_cast<Grid_2D_float*>(vp_var->get_grid()) ;
	if (vp_grid == NULL)
	{
		print_error("IN FDM_2D_el_iso_1st::initialize --> VP grid not initialized");
		return(RTN_CODE_KO) ;
	}
	Myfloat** vp = vp_grid->pArray ;

	// get vs model, grid and pointer
	Variable* vs_var = pModel->get_parameter(VS) ;
	if (vs_var == NULL)
	{
		print_error("IN FDM_2D_el_iso_1st::initialize --> VS model not found");
		return(RTN_CODE_KO) ;
	}
	Grid_2D_float* vs_grid = dynamic_cast<Grid_2D_float*>(vs_var->get_grid()) ;
	if (vs_grid == NULL)
	{
		print_error("IN FDM_2D_el_iso_1st::initialize --> VS grid not initialized");
		return(RTN_CODE_KO) ;
	}
	Myfloat** vs = vs_grid->pArray ;

	// get rho model, grid and pointer
	Variable* rho_var = pModel->get_parameter(RHO) ;
	if (rho_var == NULL)
	{
		print_error("IN FDM_2D_el_iso_1st::initialize --> RHO model not found");
		return(RTN_CODE_KO) ;
	}
	Grid_2D_float* rho_grid = dynamic_cast<Grid_2D_float*>(rho_var->get_grid()) ;
	if (rho_grid == NULL)
	{
		print_error("IN FDM_2D_el_iso_1st::initialize --> RHO grid not initialized");
		return(RTN_CODE_KO) ;
	}
	Myfloat** rho = rho_grid->pArray ;

	// coef1 = dt * (1/rho) at vx location
	// INTERPOLATION with ARITHMETIC average from main grid
	coef1 = new Grid_2D_float(izEnd+1, ixEnd+1, dz, dx) ;

	// loop on the grid points in medium
	for (Myint ix = 0; ix < nx-1; ix++)
	{
		for (Myint iz = 0; iz < nz; iz++)
		{
			coef1->pArray[ix+ixBeg2][iz+izBeg2] = dt / (0.5 * (rho[ix][iz] + rho[ix+1][iz])) ;
		}
	}

	// coef2 = dt * (1/rho) at vz location
	// INTERPOLATION with ARITHMETIC average from main grid
	coef2 = new Grid_2D_float(izEnd+1, ixEnd+1, dz, dx) ;

	// loop on the grid points in medium
	for (Myint ix = 0; ix < nx; ix++)
	{
		for (Myint iz = 0; iz < nz-1; iz++)
		{
			coef2->pArray[ix+ixBeg2][iz+izBeg2] = dt / (0.5 * (rho[ix][iz] + rho[ix][iz+1])) ;
		}
	}

	// coef3 = dt * (lambda + 2mu) at sxx and szz location
	// lambda +2mu = vp**2 * rho
	// NO INTERPOLATION -> located on the main grid
	coef3 = new Grid_2D_float(izEnd+1, ixEnd+1, dz, dx) ;

	// loop on the grid points in medium
	for (Myint ix = 0; ix < nx; ix++)
	{
		for (Myint iz = 0; iz < nz; iz++)
		{
			coef3->pArray[ix+ixBeg2][iz+izBeg2] = dt * vp[ix][iz]*vp[ix][iz] * rho[ix][iz] ;
		}
	}

	// coef4 = dt * (lambda) at sxx and szz location
	// lambda = rho * (vp**2 - 2vs**2)
	// NO INTERPOLATION -> located on the main grid
	coef4 = new Grid_2D_float(izEnd+1, ixEnd+1, dz, dx) ;

	// loop on the grid points in medium
	for (Myint ix = 0; ix < nx; ix++)
	{
		for (Myint iz = 0; iz < nz; iz++)
		{
			coef4->pArray[ix+ixBeg2][iz+izBeg2] = dt * rho[ix][iz] * (vp[ix][iz]*vp[ix][iz] - 2*vs[ix][iz]*vs[ix][iz]) ;
		}
	}

	// coef5 = dt * (mu) at sxz location
	// mu = vs**2 * rho
	// INTERPOLATION with HARMONIC average from main grid
	coef5 = new Grid_2D_float(izEnd+1, ixEnd+1, dz, dx) ;

	// loop on the grid points in medium
	for (Myint ix = 0; ix < nx-1; ix++)
	{
		for (Myint iz = 0; iz < nz-1; iz++)
		{
			Myfloat mu_h = 4. / (
					1./(vs[ix][iz]     * vs[ix][iz]     * rho[ix][iz])    +
					1./(vs[ix+1][iz]   * vs[ix+1][iz]   * rho[ix+1][iz])  +
					1./(vs[ix][iz+1]   * vs[ix][iz+1]   * rho[ix][iz+1])  +
					1./(vs[ix+1][iz+1] * vs[ix+1][iz+1] * rho[ix+1][iz+1]) ) ;
			coef5->pArray[ix+ixBeg2][iz+izBeg2] = dt * mu_h ;
		}
	}

	print_debug(ALL, MID_DEBUG, "end coef compute");

	//-------------------------------------------------------------------------------------------------------
	// initialize eq. coef. in the surrounding boundary
	//-------------------------------------------------------------------------------------------------------

	//---------------------------------------------------------------------------------------------------
	// random layer
	//---------------------------------------------------------------------------------------------------

	// $$$ RANDOM LAYER TO DO

	//---------------------------------------------------------------------------------------------------
	// sponge
	//---------------------------------------------------------------------------------------------------

	// $$$ SPONGE TO DO

	//---------------------------------------------------------------------------------------------------
	// CPML
	//---------------------------------------------------------------------------------------------------

	if (get_boundary_type(ZBEG) == PML)
	{
		// compute some coefficients
		Myfloat cpml_alpha_max = CPML_FREQ * PI ;
		Myfloat cpml_pow       = CPML_POW ;
		Myfloat cpml_vmax      = vp_grid->get_max() ;
		Myfloat cpml_rcoef     = get_boundary_coef(ZBEG) ;

		// allocate CPML arrays
		Myint nzl = nlayer_zBeg+2*lstencil ;
		Myint nxl = ixEnd+1 ;
		apml_zBeg         = new Grid_2D_float(nzl, nxl, dz, dx) ;
		bpml_zBeg         = new Grid_2D_float(nzl, nxl, dz, dx) ;
		apml_half_zBeg    = new Grid_2D_float(nzl, nxl, dz, dx) ;
		bpml_half_zBeg    = new Grid_2D_float(nzl, nxl, dz, dx) ;

		mem_sxz_zBeg      = new Grid_2D_float(nzl, nxl, dz, dx) ;
		mem_szz_zBeg      = new Grid_2D_float(nzl, nxl, dz, dx) ;
		mem_vx_zBeg       = new Grid_2D_float(nzl, nxl, dz, dx) ;
		mem_vz_zBeg       = new Grid_2D_float(nzl, nxl, dz, dx) ;

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

				// extend physical parameters
				coef1->pArray[ix][iz] = coef1->pArray[ix][izBeg2] ;
				coef2->pArray[ix][iz] = coef2->pArray[ix][izBeg2] ;
				coef3->pArray[ix][iz] = coef3->pArray[ix][izBeg2] ;
				coef4->pArray[ix][iz] = coef4->pArray[ix][izBeg2] ;
				coef5->pArray[ix][iz] = coef5->pArray[ix][izBeg2] ;

				// main grid
				Myfloat dnorm = Myfloat(izBeg2 - iz)/Myfloat(nlayer_zBeg) ;
				if (dnorm < 0) {dnorm = 0. ;}
				Myfloat alpha = cpml_alpha_max * (1. - dnorm) ;
				Myfloat dd = d0 * pow(dnorm, cpml_pow) ;

				bpml_zBeg->pArray[ipmlx][ipmlz] = exp(-(dd+alpha)*dt) ;
				apml_zBeg->pArray[ipmlx][ipmlz] = dd * (bpml_zBeg->pArray[ipmlx][ipmlz]-1.) / (dd+alpha) ;

				// staggered grid
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
		Myfloat cpml_vmax      = vp_grid->get_max() ;
		Myfloat cpml_rcoef     = get_boundary_coef(ZEND) ;

		// allocate CPML arrays
		Myint nzl = nlayer_zEnd+2*lstencil ;
		Myint nxl = ixEnd+1 ;
		apml_zEnd         = new Grid_2D_float(nzl, nxl, dz, dx) ;
		bpml_zEnd         = new Grid_2D_float(nzl, nxl, dz, dx) ;
		apml_half_zEnd    = new Grid_2D_float(nzl, nxl, dz, dx) ;
		bpml_half_zEnd    = new Grid_2D_float(nzl, nxl, dz, dx) ;

		mem_sxz_zEnd      = new Grid_2D_float(nzl, nxl, dz, dx) ;
		mem_szz_zEnd      = new Grid_2D_float(nzl, nxl, dz, dx) ;
		mem_vx_zEnd       = new Grid_2D_float(nzl, nxl, dz, dx) ;
		mem_vz_zEnd       = new Grid_2D_float(nzl, nxl, dz, dx) ;

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

				// extend physical parameters
				coef1->pArray[ix][iz]   = coef1->pArray[ix][izEnd2-1] ;
				coef2->pArray[ix][iz-1] = coef2->pArray[ix][izEnd2-2] ;
				coef3->pArray[ix][iz]   = coef3->pArray[ix][izEnd2-1] ;
				coef4->pArray[ix][iz]   = coef4->pArray[ix][izEnd2-1] ;
				coef5->pArray[ix][iz-1] = coef5->pArray[ix][izEnd2-2] ;

				// main grid
				Myfloat dnorm = Myfloat(iz - izEnd2 + 1)/Myfloat(nlayer_zEnd) ;
				if (dnorm < 0) {dnorm = 0. ;}
				Myfloat alpha = cpml_alpha_max * (1. - dnorm) ;
				Myfloat dd = d0 * pow(dnorm, cpml_pow) ;

				bpml_zEnd->pArray[ipmlx][ipmlz] = exp(-(dd+alpha)*dt) ;
				apml_zEnd->pArray[ipmlx][ipmlz] = dd * (bpml_zEnd->pArray[ipmlx][ipmlz]-1.) / (dd+alpha) ;

				// staggered grid
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
		Myfloat cpml_vmax      = vp_grid->get_max() ;
		Myfloat cpml_rcoef     = get_boundary_coef(XBEG) ;

		// allocate CPML arrays
		Myint nzl = izEnd+1 ;
		Myint nxl = nlayer_xBeg+2*lstencil ;
		apml_xBeg         = new Grid_2D_float(nzl, nxl, dz, dx) ;
		bpml_xBeg         = new Grid_2D_float(nzl, nxl, dz, dx) ;
		apml_half_xBeg    = new Grid_2D_float(nzl, nxl, dz, dx) ;
		bpml_half_xBeg    = new Grid_2D_float(nzl, nxl, dz, dx) ;

		mem_sxx_xBeg      = new Grid_2D_float(nzl, nxl, dz, dx) ;
		mem_sxz_xBeg      = new Grid_2D_float(nzl, nxl, dz, dx) ;
		mem_vx_xBeg       = new Grid_2D_float(nzl, nxl, dz, dx) ;
		mem_vz_xBeg       = new Grid_2D_float(nzl, nxl, dz, dx) ;

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

				// extend physical parameters
				coef1->pArray[ix][iz] = coef1->pArray[ixBeg2][iz] ;
				coef2->pArray[ix][iz] = coef2->pArray[ixBeg2][iz] ;
				coef3->pArray[ix][iz] = coef3->pArray[ixBeg2][iz] ;
				coef4->pArray[ix][iz] = coef4->pArray[ixBeg2][iz] ;
				coef5->pArray[ix][iz] = coef5->pArray[ixBeg2][iz] ;

				// main grid
				Myfloat dnorm = Myfloat(ixBeg2 - ix)/Myfloat(nlayer_xBeg) ;
				if (dnorm < 0) {dnorm = 0. ;}
				Myfloat alpha = cpml_alpha_max * (1. - dnorm) ;
				Myfloat dd = d0 * pow(dnorm, cpml_pow) ;

				bpml_xBeg->pArray[ipmlx][ipmlz] = exp(-(dd+alpha)*dt) ;
				apml_xBeg->pArray[ipmlx][ipmlz] = dd * (bpml_xBeg->pArray[ipmlx][ipmlz]-1.) / (dd+alpha) ;

				// staggered grid
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
		Myfloat cpml_vmax      = vp_grid->get_max() ;
		Myfloat cpml_rcoef     = get_boundary_coef(XEND) ;

		// allocate CPML arrays
		Myint nzl = izEnd+1 ;
		Myint nxl = nlayer_xEnd+2*lstencil ;
		apml_xEnd         = new Grid_2D_float(nzl, nxl, dz, dx) ;
		bpml_xEnd         = new Grid_2D_float(nzl, nxl, dz, dx) ;
		apml_half_xEnd    = new Grid_2D_float(nzl, nxl, dz, dx) ;
		bpml_half_xEnd    = new Grid_2D_float(nzl, nxl, dz, dx) ;

		mem_sxx_xEnd      = new Grid_2D_float(nzl, nxl, dz, dx) ;
		mem_sxz_xEnd      = new Grid_2D_float(nzl, nxl, dz, dx) ;
		mem_vx_xEnd       = new Grid_2D_float(nzl, nxl, dz, dx) ;
		mem_vz_xEnd       = new Grid_2D_float(nzl, nxl, dz, dx) ;

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

				// extend physical parameters
				coef1->pArray[ix-1][iz] = coef1->pArray[ixEnd2-2][iz] ;
				coef2->pArray[ix][iz]   = coef2->pArray[ixEnd2-1][iz] ;
				coef3->pArray[ix][iz]   = coef3->pArray[ixEnd2-1][iz] ;
				coef4->pArray[ix][iz]   = coef4->pArray[ixEnd2-1][iz] ;
				coef5->pArray[ix-1][iz] = coef5->pArray[ixEnd2-2][iz] ;

				// main grid
				Myfloat dnorm = Myfloat(ix - ixEnd2 + 1)/Myfloat(nlayer_xEnd) ;
				if (dnorm < 0) {dnorm = 0. ;}
				Myfloat alpha = cpml_alpha_max * (1. - dnorm) ;
				Myfloat dd = d0 * pow(dnorm, cpml_pow) ;

				bpml_xEnd->pArray[ipmlx][ipmlz] = exp(-(dd+alpha)*dt) ;
				apml_xEnd->pArray[ipmlx][ipmlz] = dd * (bpml_xEnd->pArray[ipmlx][ipmlz]-1.) / (dd+alpha) ;

				// staggered grid
				dnorm = Myfloat(ix - ixEnd2 + 0.5)/Myfloat(nlayer_xEnd) ;
				if (dnorm < 0) {dnorm = 0. ;}
				alpha = cpml_alpha_max * (1. - dnorm) ;
				dd = d0 * pow(dnorm, cpml_pow) ;

				bpml_half_xEnd->pArray[ipmlx][ipmlz] = exp(-(dd+alpha)*dt) ;
				apml_half_xEnd->pArray[ipmlx][ipmlz] = dd * (bpml_half_xEnd->pArray[ipmlx][ipmlz]-1.) / (dd+alpha) ;
			}
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D_el_iso_1st::initialize");
	return(RTN_CODE_OK) ;

} ;

//=======================================================================================================
//
// RESET MODELLING 
//
//=======================================================================================================

Rtn_code FDM_2D_el_iso_1st::reset(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D_el_iso_1st::reset");

	// call to parent class
	Rtn_code rtn_code = Scheme::reset() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// reset grids
	rtn_code = Singleton::Instance()->get_variable(VZ)->reset_grid() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
	rtn_code = Singleton::Instance()->get_variable(VX)->reset_grid() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
	rtn_code = Singleton::Instance()->get_variable(SXX)->reset_grid() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
	rtn_code = Singleton::Instance()->get_variable(SZZ)->reset_grid() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
	rtn_code = Singleton::Instance()->get_variable(SXZ)->reset_grid() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// initialize CPML memory variables
	if (get_boundary_type(ZBEG) == PML)
	{
		mem_sxz_zBeg->reset() ;
		mem_szz_zBeg->reset() ;
		mem_vx_zBeg->reset() ;
		mem_vz_zBeg->reset() ;
	}
	if (get_boundary_type(ZEND) == PML)
	{
		mem_sxz_zEnd->reset() ;
		mem_szz_zEnd->reset() ;
		mem_vx_zEnd->reset() ;
		mem_vz_zEnd->reset() ;
	}
	if (get_boundary_type(XBEG) == PML)
	{
		mem_sxx_xBeg->reset() ;
		mem_sxz_xBeg->reset() ;
		mem_vx_xBeg->reset() ;
		mem_vz_xBeg->reset() ;
	}
	if (get_boundary_type(XEND) == PML)
	{
		mem_sxx_xEnd->reset() ;
		mem_sxz_xEnd->reset() ;
		mem_vx_xEnd->reset() ;
		mem_vz_xEnd->reset() ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D_el_iso_1st::reset");
	return(RTN_CODE_OK) ;

} ;

//=======================================================================================================
//
// FINALIZE MODELLING 
//
//=======================================================================================================

Rtn_code FDM_2D_el_iso_1st::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D_el_iso_1st::finalize");

	//---------------------------------------------------------------------------------------------------
	// deallocate arrays
	//---------------------------------------------------------------------------------------------------

	// physical parameter
	delete(coef1) ;
	delete(coef2) ;
	delete(coef3) ;
	delete(coef4) ;
	delete(coef5) ;

	// delete components
	Singleton::Instance()->delete_variable(VZ) ;
	Singleton::Instance()->delete_variable(VX) ;
	Singleton::Instance()->delete_variable(SXX) ;
	Singleton::Instance()->delete_variable(SZZ) ;
	Singleton::Instance()->delete_variable(SXZ) ;

	// call parent finalize
	Rtn_code rtn_code = FDM_2D_el_iso::finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D_el_iso_1st::finalize");
	return(RTN_CODE_OK) ;

} ;

//=======================================================================================================
//
// PERFORM FORWARD MODELLING FOR ONE SOURCE
//
//=======================================================================================================


Rtn_code FDM_2D_el_iso_1st::solve_current_shot(Acquisition* pAcquisition, Data *pData, Wavefield_type wtype, Snapshot* pSnapshot)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D_el_iso_1st::solve_current_shot");

	// get variables
	Variable* sxx_var = Singleton::Instance()->get_variable(SXX) ;
	Grid_2D_float *sxx_grid = (Grid_2D_float*) sxx_var->get_grid() ;
	Myfloat** const sxx = sxx_grid->pArray ;

	Variable* szz_var = Singleton::Instance()->get_variable(SZZ) ;
	Grid_2D_float *szz_grid = (Grid_2D_float*) szz_var->get_grid() ;
	Myfloat** const szz = szz_grid->pArray ;

	Variable* sxz_var = Singleton::Instance()->get_variable(SXZ) ;
	Grid_2D_float *sxz_grid = (Grid_2D_float*) sxz_var->get_grid() ;
	Myfloat** const sxz = sxz_grid->pArray ;

	Variable* vz_var = Singleton::Instance()->get_variable(VZ) ;
	Grid_2D_float *vz_grid = (Grid_2D_float*) vz_var->get_grid() ;
	Myfloat** const vz = vz_grid->pArray ;

	Variable* vx_var = Singleton::Instance()->get_variable(VX) ;
	Grid_2D_float *vx_grid = (Grid_2D_float*) vx_var->get_grid() ;
	Myfloat** const vx = vx_grid->pArray ;

	// reset unknowns to zero
	Rtn_code rtn_code = this->reset() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// locate src and rec in the grid
	rtn_code = this->locate_src_and_rec_in_grid(pAcquisition) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// locate snapshot pixel in the grid
	rtn_code = this->locate_pixel_in_grid(pSnapshot) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// source scaling factor
	Myfloat src_factor = dt / (dx * dz) ;

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
			this->compute_velocity() ;
			time_in_kernel += MPI_Wtime() - t0 ;

			//-------------------------------------------------------------------------------------------------------
			// source excitation on velocity component
			//-------------------------------------------------------------------------------------------------------
			if (src_stype == FORCE_Z)
			{
				this->source_excitation(vz, it, wtype, pData, src_factor) ;
			}

			//-------------------------------------------------------------------------------------------------------
			// computation of the stress component 1t time = [it + 1/2] * dt
			//-------------------------------------------------------------------------------------------------------
			t0 = MPI_Wtime() ;
			this->compute_stress() ;
			time_in_kernel += MPI_Wtime() - t0 ;

			//-------------------------------------------------------------------------------------------------------
			// source excitation on stress component
			//-------------------------------------------------------------------------------------------------------
			if (src_stype == EXPLOSIVE)
			{
				this->source_excitation(sxx, it, wtype, pData, src_factor) ;
				this->source_excitation(szz, it, wtype, pData, src_factor) ;
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

			// compute energy
			compute_energy(vx_var, vz_var, sxx_var, szz_var, sxz_var, it) ;

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

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D_el_iso_1st::solve_current_shot");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------

void FDM_2D_el_iso_1st::compute_energy(Variable* vx_var, Variable* vz_var, Variable* sxx_var,
		Variable* szz_var, Variable* sxz_var, Myint it)
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_2D_el_iso_1st::compute_energy");

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
				Grid_2D_float *vx_grid = (Grid_2D_float*) vx_var->get_grid() ;
				Myfloat** const vx = vx_grid->pArray ;
				Grid_2D_float *vz_grid = (Grid_2D_float*) vz_var->get_grid() ;
				Myfloat** const vz = vz_grid->pArray ;
				Grid_2D_float *sxx_grid = (Grid_2D_float*) sxx_var->get_grid() ;
				Myfloat** const sxx = sxx_grid->pArray ;
				Grid_2D_float *szz_grid = (Grid_2D_float*) szz_var->get_grid() ;
				Myfloat** const szz = szz_grid->pArray ;
				Grid_2D_float *sxz_grid = (Grid_2D_float*) sxz_var->get_grid() ;
				Myfloat** const sxz = sxz_grid->pArray ;

				// compute energy
				// loop on the receivers
				energy_kin = 0.0 ;
				energy_pot = 0.0 ;
				for (Myint irec=0; irec<nrec; irec++)
				{
					if ((ix_rec[irec] != NOT_FOUND) && (iz_rec[irec] != NOT_FOUND))
					{
						energy_kin += 0.5 * rho * (pow(vx[ix_rec[irec]][iz_rec[irec]], 2.0) + pow(vz[ix_rec[irec]][iz_rec[irec]], 2.0)) ;
						energy_pot += 0.5 * coef_a * (pow(sxx[ix_rec[irec]][iz_rec[irec]], 2.0) + pow(szz[ix_rec[irec]][iz_rec[irec]], 2.0))
					- coef_b * sxx[ix_rec[irec]][iz_rec[irec]] * szz[ix_rec[irec]][iz_rec[irec]]
																				   + coef_c * pow(sxz[ix_rec[irec]][iz_rec[irec]], 2.0) ;
					}
				}
				energy_tot = energy_kin + energy_pot ;

				// write energy
				write_energy() ;
			}
		}
	}

	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D_el_iso_1st::compute_energy");
}

} // namespace django
