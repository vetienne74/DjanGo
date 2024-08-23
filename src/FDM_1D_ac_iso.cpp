//------------------------------------------------------------------------------------------------------
//
// MODELLING WITH FDM IN 1D
//  * ACOUSTIC ISOTROPIC MEDIA 
//
// PARENT CLASS: Scheme
//   DERIVED CLASS: FDM
//     DERIVED CLASS: FDM_1D
//       DERIVED CLASS: FDM_1D_ac_iso
//
//-------------------------------------------------------------------------------------------------------

#include "FDM_1D_ac_iso.h"

#include "output_report.h"
#include "type_def.h"

namespace django {

//-------------------------------------------------------------------------------------------------------
FDM_1D_ac_iso::FDM_1D_ac_iso(void) : FDM_1D()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D_ac_iso::FDM_1D_ac_iso");
	eq_type = ACOUSTIC ;
	sponge_coef_pr = NULL ;
	sponge_coef_vz = NULL ;
	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D_ac_iso::FDM_1D_ac_iso");
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FDM_1D_ac_iso::initialize(Model* pModel, Myfloat fmax)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D_ac_iso::initialize");

	// get vp model
	Variable* vp_var = pModel->get_parameter(VP) ;
	if (vp_var == NULL)
	{
		print_error("IN FDM_1D_ac_iso_1st::initialize --> VP model not found");
		return(RTN_CODE_KO) ;
	}
	// get vp grid
	Grid* vp_grid = vp_var->get_grid() ;
	if (vp_grid == NULL)
	{
		print_error("IN FDM_1D_ac_iso_1st::initialize --> VP grid not initialized");
		return(RTN_CODE_KO) ;
	}

	// call parent initialization
	Rtn_code rtn_code = FDM_1D::initialize(vp_grid, fmax) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D_ac_iso::initialize");
	return(RTN_CODE_OK) ;
} ;

//-------------------------------------------------------------------------------------------------------
Rtn_code FDM_1D_ac_iso::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D_ac_iso::finalize");

	// cpml
	if (get_boundary_type(ZBEG) == PML)
	{
		delete(apml_zBeg) ;
		delete(bpml_zBeg) ;
		delete(apml_half_zBeg) ;
		delete(bpml_half_zBeg) ;
		delete(mem_vz_zBeg) ;
		delete(mem_pr_zBeg) ;
	}
	if (get_boundary_type(ZEND) == PML)
	{
		delete(apml_zEnd) ;
		delete(bpml_zEnd) ;
		delete(apml_half_zEnd) ;
		delete(bpml_half_zEnd) ;
		delete(mem_vz_zEnd) ;
		delete(mem_pr_zEnd) ;
	}

	// sponge
	if (is_there_boundary_type(SPG))
	{
		delete(sponge_coef_pr) ;
		delete(sponge_coef_vz) ;
	}

	// delete boundary
	//delete_all_boundary() ;

	// call parent finalize
	Rtn_code rtn_code = FDM_1D::finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D_ac_iso::finalize");
	return(RTN_CODE_OK) ;
} ;

} // namespace django
