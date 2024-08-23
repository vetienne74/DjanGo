//-------------------------------------------------------------------------------------------------------
//
// MODELLING WITH FDM IN 2D
//
//  * ACOUSTIC ISOTROPIC MEDIA 
//  * 1ST ORDER WAVE EQUATION (PRESSURE / VELOCITY):
//
//            dvx/dt = coef2 dP/dx
//            dvz/dt = coef2 dP/dz
//            dP /dt = coef1 (dvx/dx + dvz/dz) 
//
//            with
//
//            coef1 = kappa = (rho * vp * vp) at pr location 
//            coef2 = 1/rho
//            where rho is CONSTANT
//
// PARENT CLASS: Scheme
//   DERIVED CLASS: FDM
//     DERIVED CLASS FDM_1D
//       DERIVED CLASS: FDM_2D
//         DERIVED CLASS: FDM_2D_ac_iso
//           DERIVED CLASS: FDM_2D_ac_iso_1st
//
//-------------------------------------------------------------------------------------------------------

#include "FDM_2D_ac_iso_1st.h"

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
FDM_2D_ac_iso_1st::FDM_2D_ac_iso_1st(void) : FDM_2D_ac_iso()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D_ac_iso_1st::FDM_2D_ac_iso_1st");
	eq_order = ORDER_1ST ;
	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D_ac_iso_1st::FDM_2D_ac_iso_1st");
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FDM_2D_ac_iso_1st::update_physical_coef(Myint iz, Myint ix, Myint iy, Myfloat perturb)  
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D_ac_iso_1st::update_physical_coef");

	// retrieve vp and rho
	Myfloat rho = TEMP_RHO_CONST ;
	Myfloat vp = sqrt(coef1->pArray[ix][iz]/(dt * rho)) ;

	// add perturbation
	vp += perturb ;

	// update physical parameter
	coef1->pArray[ix][iz] = dt * pow(vp, 2) * rho ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D_ac_iso_1st::update_physical_coef");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FDM_2D_ac_iso_1st::update_physical_coef(Model* pModel)  
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D_ac_iso_1st::update_physical_coef");

	Myfloat rho = TEMP_RHO_CONST ;

	// get vp model
	Grid_2D_float *vp_grid = (Grid_2D_float*) (pModel->get_parameter(VP))->get_grid() ;
	Myfloat** vp = vp_grid->pArray ;

	// coef1_i = rho_i * vp_i**2 * dt (located on the pressure nodes)

	// loop on the grid points in medium
#pragma omp parallel for 
	for (Myint ix = 0; ix < nx; ix++)
	{
#pragma ivdep
		for (Myint iz = 0; iz < nz; iz++)
		{
			// update coef in the FDM grid point
			coef1->pArray[ix+ixBeg2][iz+izBeg2] = dt * pow(vp[ix][iz], 2) * rho ;
		}
	}

	// coef2_i = dt / rho_i (located on the velocity nodes)

	coef2 = dt / rho ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D_ac_iso_1st::update_physical_coef");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FDM_2D_ac_iso_1st::initialize(Model* pModel, Myfloat fmax)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D_ac_iso_1st::initialize");

	// call parent initialiation
	Rtn_code rtn_code = FDM_2D_ac_iso::initialize(pModel, fmax) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	//-------------------------------------------------------------------------------------------------------
	// unknown vectors
	//-------------------------------------------------------------------------------------------------------

	// velocity component (time domain)
	Variable* pVarVz = Singleton::Instance()->register_variable(VZ, "vz") ;
	pVarVz->allocate_grid(izEnd+1, ixEnd+1, dz, dx) ;
	Variable* pVarVx = Singleton::Instance()->register_variable(VX, "vx") ;
	pVarVx->allocate_grid(izEnd+1, ixEnd+1, dz, dx) ;

	// pressure component (time domain)
	Variable* pVarPr = Singleton::Instance()->register_variable(PR, "pr") ;
	pVarPr->allocate_grid(izEnd+1, ixEnd+1, dz, dx) ;

	//-------------------------------------------------------------------------------------------------------
	// initialize eq. coef. in the medium
	//-------------------------------------------------------------------------------------------------------

	// check model type is GRID and subtype is REGULAR
	if (pModel->get_type() != GRID)
	{
		print_error("IN FDM_2D_ac_iso_1st::initialize --> model type is not GRID");
		return(RTN_CODE_KO) ;
	}
	if (pModel->get_sub_type() != REGULAR)
	{
		print_error("IN FDM_2D_ac_iso_1st::initialize --> model subtype is not REGULAR");
		return(RTN_CODE_KO) ;
	}

	Myfloat rho = TEMP_RHO_CONST ;

	// coef1_i = rho_i * vp_i**2 * dt (located on the pressure nodes)

	coef1 = new Grid_2D_float(izEnd+1, ixEnd+1, dz, dx) ;

	// get vp model, grid and pointer
	Variable* vp_var = pModel->get_parameter(VP) ;
	if (vp_var == NULL)
	{
		print_error("IN FDM_2D_ac_iso_1st::initialize --> VP model not found");
		return(RTN_CODE_KO) ;
	}
	Grid_2D_float* vp_grid = dynamic_cast<Grid_2D_float*>(vp_var->get_grid()) ;
	if (vp_grid == NULL)
	{
		print_error("IN FDM_2D_ac_iso_1st::initialize --> VP grid not initialized");
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
			coef1->pArray[ix+ixBeg2][iz+izBeg2] = dt * pow(vp[ix][iz], 2) * rho ;
		}
	}

	// coef2_i = dt / rho_i (located on the velocity nodes)

	coef2 = dt / rho ;

	print_debug(ALL, MID_DEBUG, "end coef compute");

	//-------------------------------------------------------------------------------------------------------
	// initialize eq. coef. in the surrounding boundary
	//-------------------------------------------------------------------------------------------------------

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
				coef1->pArray[ix][iz] = dt * vp_rand * vp_rand * rho ;
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
				coef1->pArray[ix][iz] = dt * vp_rand * vp_rand * rho ;
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
				coef1->pArray[ix][iz] = dt * vp_rand * vp_rand * rho ;
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
				coef1->pArray[ix][iz] = dt * vp_rand * vp_rand * rho ;
			}
		}
	}

	//---------------------------------------------------------------------------------------------------
	// sponge
	//---------------------------------------------------------------------------------------------------

	if (is_there_boundary_type(SPG))
	{

		// allocate sponge arrays
		sponge_coef_pr = new Grid_2D_float(izEnd+1, ixEnd+1, dz, dx) ;
		sponge_coef_vz = new Grid_2D_float(izEnd,   ixEnd,   dz, dx) ;
		sponge_coef_vx = new Grid_2D_float(izEnd,   ixEnd,   dz, dx) ;

		// loop on all grid points
		for (Myint ix = ixBeg1; ix < ixEnd1; ix++)
		{
			for (Myint iz = izBeg1; iz < izEnd1; iz++)
			{

				// default value
				sponge_coef_pr->pArray[ix][iz] = 1. ;
				sponge_coef_vz->pArray[ix][iz] = 1. ;
				sponge_coef_vx->pArray[ix][iz] = 1. ;

				Myint ix2 = ix - izBeg2 ;
				ix2 = max(ix2, 0) ;
				ix2 = min(ix2, nx-1) ;

				Myint iz2 = iz - izBeg2 ;
				iz2 = max(iz2, 0) ;
				iz2 = min(iz2, nz-1) ;

				Myfloat vp_tmp = vp[ix2][iz2] ;
				coef1->pArray[ix][iz] = dt * vp_tmp * vp_tmp * rho ;

				// z- layer
				if ((iz >= izBeg1) && (iz < izBeg2) && (get_boundary_type(ZBEG) == SPG))
				{
					Myfloat sponge_coef = get_boundary_coef(ZBEG) ;
					sponge_coef_pr->pArray[ix][iz] *= exp(-pow(Myfloat(izBeg2 - iz) * sponge_coef, 2)) ;
					sponge_coef_vz->pArray[ix][iz] *= exp(-pow(Myfloat(izBeg2 - iz - 0.5) * sponge_coef, 2)) ;
					sponge_coef_vx->pArray[ix][iz] *= exp(-pow(Myfloat(izBeg2 - iz) * sponge_coef, 2)) ;
				}

				// z+ layer
				if ((iz >= izEnd2) && (iz < izEnd1) && (get_boundary_type(ZEND) == SPG))
				{
					Myfloat sponge_coef = get_boundary_coef(ZEND) ;
					sponge_coef_pr->pArray[ix][iz] *= exp(-pow(Myfloat(iz - izEnd2 + 1) * sponge_coef, 2)) ;
					sponge_coef_vz->pArray[ix][iz-1] *= exp(-pow(Myfloat(iz - izEnd2 + 0.5) * sponge_coef, 2)) ;
					sponge_coef_vx->pArray[ix][iz] *= exp(-pow(Myfloat(iz - izEnd2 + 1) * sponge_coef, 2)) ;
				}

				// x- layer
				if ((ix >= ixBeg1) && (ix < ixBeg2) && (get_boundary_type(XBEG) == SPG))
				{
					Myfloat sponge_coef = get_boundary_coef(XBEG) ;
					sponge_coef_pr->pArray[ix][iz] *= exp(-pow(Myfloat(ixBeg2 - ix) * sponge_coef, 2)) ;
					sponge_coef_vz->pArray[ix][iz] *= exp(-pow(Myfloat(ixBeg2 - ix) * sponge_coef, 2)) ;
					sponge_coef_vx->pArray[ix][iz] *= exp(-pow(Myfloat(ixBeg2 - ix - 0.5) * sponge_coef, 2)) ;
				}

				// x+ layer
				if ((ix >= ixEnd2) && (ix < ixEnd1) && (get_boundary_type(XEND) == SPG))
				{
					Myfloat sponge_coef = get_boundary_coef(XEND) ;
					sponge_coef_pr->pArray[ix][iz] *= exp(-pow(Myfloat(ix - ixEnd2 + 1) * sponge_coef, 2)) ;
					sponge_coef_vz->pArray[ix][iz] *= exp(-pow(Myfloat(ix - ixEnd2 + 1) * sponge_coef, 2)) ;
					sponge_coef_vx->pArray[ix-1][iz] *= exp(-pow(Myfloat(ix - ixEnd2 + 0.5) * sponge_coef, 2)) ;
				}
			}
		}
	}

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
				coef1->pArray[ix][iz] = coef1->pArray[ix][izBeg2] ;

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
		Myfloat cpml_vmax      = vp_grid->get_max() ;
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
				coef1->pArray[ix][iz] = coef1->pArray[ix][izEnd2-1] ;

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
		Myfloat cpml_vmax      = vp_grid->get_max() ;
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
				coef1->pArray[ix][iz] = coef1->pArray[ixBeg2][iz] ;

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
		Myfloat cpml_vmax      = vp_grid->get_max() ;
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
				coef1->pArray[ix][iz] = coef1->pArray[ixEnd2-1][iz] ;

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

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D_ac_iso_1st::initialize");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
Rtn_code FDM_2D_ac_iso_1st::reset(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D_ac_iso_1st::reset");

	// call to parent class
	Rtn_code rtn_code = Scheme::reset() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// reset grids
	rtn_code = Singleton::Instance()->get_variable(PR)->reset_grid() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
	rtn_code = Singleton::Instance()->get_variable(VZ)->reset_grid() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
	rtn_code = Singleton::Instance()->get_variable(VX)->reset_grid() ;
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

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D_ac_iso_1st::reset");
	return(RTN_CODE_OK) ;

} ;

//=======================================================================================================
//
// FINALIZE MODELLING 
//
//=======================================================================================================

Rtn_code FDM_2D_ac_iso_1st::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D_ac_iso_1st::finalize");

	//---------------------------------------------------------------------------------------------------
	// deallocate arrays
	//---------------------------------------------------------------------------------------------------

	// physical parameter
	delete(coef1) ;

	// delete components
	Singleton::Instance()->delete_variable(VZ) ;
	Singleton::Instance()->delete_variable(VX) ;
	Singleton::Instance()->delete_variable(PR) ;

	// call parent finalize
	Rtn_code rtn_code = FDM_2D_ac_iso::finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D_ac_iso_1st::finalize");
	return(RTN_CODE_OK) ;

} ;

//=======================================================================================================
//
// PERFORM FORWARD MODELLING FOR ONE SOURCE
//
//=======================================================================================================


Rtn_code FDM_2D_ac_iso_1st::solve_current_shot(Acquisition* pAcquisition, Data *pData, Wavefield_type wtype, Snapshot* pSnapshot)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D_ac_iso_1st::solve_current_shot");

	// get variables
	Variable* pr_var = Singleton::Instance()->get_variable(PR) ;
	Grid_2D_float *pr_grid = (Grid_2D_float*) pr_var->get_grid() ;
	Myfloat** const pr = pr_grid->pArray ;

	Variable* vz_var = Singleton::Instance()->get_variable(VZ) ;
	Grid_2D_float *vz_grid = (Grid_2D_float*) vz_var->get_grid() ;
	Myfloat** const vz = vz_grid->pArray ;

	Variable* vx_var = Singleton::Instance()->get_variable(VX) ;
	Grid_2D_float *vx_grid = (Grid_2D_float*) vx_var->get_grid() ;
	Myfloat** const vx = vx_grid->pArray ;

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
				rtn_code = this->source_excitation(pr, it, wtype, pData, src_factor) ;
			}
			else if (src_stype == FORCE_Z)
			{
				rtn_code = this->source_excitation(vz, it, wtype, pData, src_factor) ;
			}
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

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
			this->compute_velocity() ;
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

			// compute energy
			this->compute_energy(pr) ;

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

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D_ac_iso_1st::solve_current_shot");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
void FDM_2D_ac_iso_1st::init_eigen()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D_ac_iso_1st::init_eigen");

	{
		Grid_2D_float *vz_grid = (Grid_2D_float*) (Singleton::Instance()->get_variable(VZ))->get_grid() ;
		Myfloat ** const vz = vz_grid->pArray ;
		Grid_2D_float *vx_grid = (Grid_2D_float*) (Singleton::Instance()->get_variable(VX))->get_grid() ;
		Myfloat ** const vx = vx_grid->pArray ;

		// get number mode number
		eigen_nmode = (Myint) Singleton::Instance()->pProgram->pModelling->get_param() ;

		Myfloat time_v = -dt/2.0 ;
		Myint izmin = izBeg2 ;
		Myint izmax = izEnd2-1 ;
		Myint nz = izEnd2 - izBeg2 ;
		Myint ixmin = ixBeg2 ;
		Myint ixmax = ixEnd2-1 ;
		Myint nx = ixEnd2 - ixBeg2 ;

		for (Myint ix=ixmin; ix<ixmax; ix++)
		{
			for (Myint iz=izmin; iz<izmax; iz++)
			{
				Myfloat zz = float(0.5+iz-izmin)/float(nz-1) ;
				Myfloat xx = float(0.5+ix-ixmin)/float(nx-1) ;
				vx[ix][iz] = cos(M_PI*xx* eigen_nmode) * sin(M_PI*zz* eigen_nmode) * cos(M_PI*sqrt(2.0)*time_v* eigen_nmode) ;
				vz[ix][iz] = cos(M_PI*zz* eigen_nmode) * sin(M_PI*xx* eigen_nmode) * cos(M_PI*sqrt(2.0)*time_v* eigen_nmode) ;
			}
		}

		{
			//--------------------------------
			// free surface and rigid surface
			//--------------------------------

			// loop on boundaries
			Myint ii_sign, vel_sign, x_axe, z_axe ;
			Myint ixmin, ixmax, izmin, izmax ;
			Myint x_offset1, x_offset2, z_offset1, z_offset2 ;
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
					ixmax   = ixEnd-1 ;
					ii_sign = 1 ;
					x_axe   = 0 ;
					z_axe   = 1 ;
					x_offset1 = 0 ;
					x_offset2 = 0 ;
					z_offset1 = 0 ;
					z_offset2 = -1 ;
					switch1 = true ;
				}
				else if (pBoundary[ib-1]->get_edge() == ZEND)
				{
					izmin   = izEnd2-1 ;
					izmax   = izEnd2-1 ;
					ixmin   = ixBeg ;
					ixmax   = ixEnd-1 ;
					ii_sign = -1 ;
					x_axe   = 0 ;
					z_axe   = 1 ;
					x_offset1 = 0 ;
					x_offset2 = 0 ;
					z_offset1 = -1 ;
					z_offset2 = 0 ;
					switch1 = true ;
				}
				else if (pBoundary[ib-1]->get_edge() == XBEG)
				{
					izmin   = izBeg ;
					izmax   = izEnd-1 ;
					ixmin   = ixBeg2 ;
					ixmax   = ixBeg2 ;
					ii_sign = 1 ;
					x_axe   = 1 ;
					z_axe   = 0 ;
					x_offset1 = 0 ;
					x_offset2 = -1 ;
					z_offset1 = 0 ;
					z_offset2 = 0 ;
					switch1 = true ;
				}
				else if (pBoundary[ib-1]->get_edge() == XEND)
				{
					izmin   = izBeg ;
					izmax   = izEnd-1 ;
					ixmin   = ixEnd2-1 ;
					ixmax   = ixEnd2-1 ;
					ii_sign = -1 ;
					x_axe   = 1 ;
					z_axe   = 0 ;
					x_offset1 = -1 ;
					x_offset2 = 0 ;
					z_offset1 = 0 ;
					z_offset2 = 0 ;
					switch1 = true ;
				}

				// free surface with image method
				// anti-symetry of pressure wavefield
				if (pBoundary[ib-1]->get_type() == FREESURF)
				{
					vel_sign = 1.0 ;
					switch2 = true ;
				}

				// free surface with image method
				// anti-symetry of pressure wavefield
				else if (pBoundary[ib-1]->get_type() == RIGID)
				{
					vel_sign = -1.0 ;
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
#pragma ivdep
							for (Myint ii = 1; ii <= lstencil; ii++)
							{
								vx[ix-x_axe*(ii_sign*ii-x_offset1)][iz-z_axe*(ii_sign*ii-z_offset1)]
																	= vel_sign * vx[ix+x_axe*(ii_sign*ii+x_offset2)][iz+z_axe*(ii_sign*ii+z_offset2)] ;

								vz[ix-x_axe*(ii_sign*ii-x_offset1)][iz-z_axe*(ii_sign*ii-z_offset1)]
																	= vel_sign * vz[ix+x_axe*(ii_sign*ii+x_offset2)][iz+z_axe*(ii_sign*ii+z_offset2)] ;
							}
						}
					}
				}
			}

		}

	}
	{
		Grid_2D_float *pr_grid = (Grid_2D_float*) (Singleton::Instance()->get_variable(PR))->get_grid() ;
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
				pr[ix][iz] = - sqrt(2.0) * sin(M_PI*xx* eigen_nmode) * sin(M_PI*zz* eigen_nmode) * sin(M_PI*sqrt(2.0)*time_p* eigen_nmode) ;
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

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D_ac_iso_1st::init_eigen");
	return ;
} ;

} // namespace django
