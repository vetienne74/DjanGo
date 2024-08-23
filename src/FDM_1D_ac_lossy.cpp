//------------------------------------------------------------------------------------------------------
//
// MODELLING WITH FDM IN 1D
//  * ACOUSTIC LOSSY
//
// PARENT CLASS: Scheme
//   DERIVED CLASS: FDM
//     DERIVED CLASS: FDM_1D
//       DERIVED CLASS: FDM_1D_ac_lossy
//
//-------------------------------------------------------------------------------------------------------

#include "FDM_1D_ac_lossy.h"

#include <iostream>

#include "output_report.h"
#include "type_def.h"

namespace django {

//-------------------------------------------------------------------------------------------------------

FDM_1D_ac_lossy::FDM_1D_ac_lossy(void) : FDM_1D()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D_ac_lossy::FDM_1D_ac_lossy");
	eq_type = AC_LOSSY ;
	string_length_vs_time_step = NULL ;
	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D_ac_lossy::FDM_1D_ac_lossy");
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FDM_1D_ac_lossy::set_string_length_vs_time_step(Myfloat* string_length_vs_time_step_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D_ac_lossy::set_string_length_vs_time_step");
	string_length_vs_time_step = string_length_vs_time_step_in ;
	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D_ac_lossy::set_string_length_vs_time_step");
	return(RTN_CODE_OK) ;
} ;

//-------------------------------------------------------------------------------------------------------
Rtn_code FDM_1D_ac_lossy::initialize(Model* pModel, Myfloat fmax)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D_ac_lossy::initialize");

	// get vp model
	Variable* vp_var = pModel->get_parameter(VP) ;
	if (vp_var == NULL)
	{
		print_error("IN FDM_1D_ac_lossy::initialize --> VP model not found");
		return(RTN_CODE_KO) ;
	}
	// get vp grid
	Grid* vp_grid = vp_var->get_grid() ;
	if (vp_grid == NULL)
	{
		print_error("IN FDM_1D_ac_lossy::initialize --> VP grid not initialized");
		return(RTN_CODE_KO) ;
	}

	// call parent initialization
	Rtn_code rtn_code = FDM_1D::initialize(vp_grid, fmax) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// loss1
	loss1 = new Grid_1D_float(izEnd+1, dz) ;

	// get loss model
	Variable* loss1_var = pModel->get_parameter(LOSS1) ;
	if (loss1_var == NULL)
	{
		print_error("IN FDM_1D_ac_lossy::initialize --> LOSS1 model not found");
		return(RTN_CODE_KO) ;
	}

	// get loss grid
	Grid_1D_float* loss1_grid = dynamic_cast<Grid_1D_float*>(loss1_var->get_grid()) ;
	if (loss1_grid == NULL)
	{
		print_error("IN FDM_1D_ac_lossy::initialize --> LOSS1 grid not initialized");
		return(RTN_CODE_KO) ;
	}

	// loss2
	loss2 = new Grid_1D_float(izEnd+1, dz) ;

	// get loss model
	Variable* loss2_var = pModel->get_parameter(LOSS2) ;
	if (loss2_var == NULL)
	{
		print_error("IN FDM_1D_ac_lossy::initialize --> LOSS2 model not found");
		return(RTN_CODE_KO) ;
	}

	// get loss grid
	Grid_1D_float* loss2_grid = dynamic_cast<Grid_1D_float*>(loss2_var->get_grid()) ;
	if (loss2_grid == NULL)
	{
		print_error("IN FDM_1D_ac_lossy::initialize --> LOSS2 grid not initialized");
		return(RTN_CODE_KO) ;
	}

	// reset grid to initialize correctly the halos
	loss1->reset() ;
	loss2->reset() ;

	// loop on the grid points in medium
	for (Myint iz = 0; iz < nz; iz++)
	{
		loss1->pArray[iz+izBeg2] = loss1_grid->pArray[iz] ;
		loss2->pArray[iz+izBeg2] = loss2_grid->pArray[iz] ;
	}


	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D_ac_lossy::initialize");
	return(RTN_CODE_OK) ;
} ;

//-------------------------------------------------------------------------------------------------------
Rtn_code FDM_1D_ac_lossy::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D_ac_lossy::finalize");

	// physical parameter
	delete(loss1) ;
	delete(loss2) ;


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

	// delete boundary
	//delete_all_boundary() ;

	// call parent finalize
	Rtn_code rtn_code = FDM_1D::finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D_ac_lossy::finalize");
	return(RTN_CODE_OK) ;
} ;

} // namespace django
