#ifndef DJANGO_FDM_2D_EL_ISO_H_
#define DJANGO_FDM_2D_EL_ISO_H_

#include "FDM_2D.h"

#include "model.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class FDM_2D_el_iso: public FDM_2D
{

protected:

	// constructor
	FDM_2D_el_iso(void) ;

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
	Grid_2D_float* mem_sxz_zBeg ;
	Grid_2D_float* mem_szz_zBeg ;
	Grid_2D_float* mem_vx_zBeg ;
	Grid_2D_float* mem_vz_zBeg ;
	Grid_2D_float* mem_sxz_zEnd ;
	Grid_2D_float* mem_szz_zEnd ;
	Grid_2D_float* mem_vx_zEnd ;
	Grid_2D_float* mem_vz_zEnd ;
	Grid_2D_float* mem_sxx_xBeg ;
	Grid_2D_float* mem_sxz_xBeg ;
	Grid_2D_float* mem_vx_xBeg ;
	Grid_2D_float* mem_vz_xBeg ;
	Grid_2D_float* mem_sxx_xEnd ;
	Grid_2D_float* mem_sxz_xEnd ;
	Grid_2D_float* mem_vx_xEnd ;
	Grid_2D_float* mem_vz_xEnd ;

} ;

} // namespace django

#endif
