#ifndef DJANGO_FDM_2D_AC_ISO_H_
#define DJANGO_FDM_2D_AC_ISO_H_

#include "FDM_2D.h"

#include "model.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class FDM_2D_ac_iso: public FDM_2D
{

protected:

	// constructor
	FDM_2D_ac_iso(void) ;

	// initialize scheme
	virtual Rtn_code initialize(Model*, Myfloat fmax) ;

	// finalize scheme
	virtual Rtn_code finalize(void) ;

	// pml coeffient
	Grid_2D_float* apml_zBeg ;
	Grid_2D_float* apml_xBeg ;
	Grid_2D_float* bpml_zBeg ;
	Grid_2D_float* bpml_xBeg ;
	Grid_2D_float* apml_half_zBeg;
	Grid_2D_float* apml_half_xBeg;
	Grid_2D_float* bpml_half_zBeg ;
	Grid_2D_float* bpml_half_xBeg ;
	Grid_2D_float* apml_zEnd ;
	Grid_2D_float* apml_xEnd ;
	Grid_2D_float* bpml_zEnd ;
	Grid_2D_float* bpml_xEnd ;
	Grid_2D_float* apml_half_zEnd;
	Grid_2D_float* apml_half_xEnd;
	Grid_2D_float* bpml_half_zEnd ;
	Grid_2D_float* bpml_half_xEnd ;

	// pml memory variable
	Grid_2D_float* mem_pr_zBeg ;
	Grid_2D_float* mem_pr_xBeg ;
	Grid_2D_float* mem_vz_zBeg ;
	Grid_2D_float* mem_vx_xBeg ;
	Grid_2D_float* mem_pr_zEnd ;
	Grid_2D_float* mem_pr_xEnd ;
	Grid_2D_float* mem_vz_zEnd ;
	Grid_2D_float* mem_vx_xEnd ;

	// sponge coefficient
	Grid_2D_float* sponge_coef_pr ;
	Grid_2D_float* sponge_coef_vz ;
	Grid_2D_float* sponge_coef_vx ;

} ;

} // namespace django

#endif
