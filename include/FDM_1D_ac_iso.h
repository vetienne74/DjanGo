#ifndef DJANGO_FDM_1D_AC_ISO_H_
#define DJANGO_FDM_1D_AC_ISO_H_

#include "FDM_1D.h"

#include "model.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class FDM_1D_ac_iso: public FDM_1D
{

public:

	// constructor
	FDM_1D_ac_iso(void) ;

	// initialize scheme
	virtual Rtn_code initialize(Model*, Myfloat fmax) ;

	// finalize scheme
	virtual Rtn_code finalize(void) ;

	// pml coeffient
	Grid_1D_float *apml_zBeg ;
	Grid_1D_float *bpml_zBeg ;
	Grid_1D_float *apml_half_zBeg;
	Grid_1D_float *bpml_half_zBeg ;
	Grid_1D_float *apml_zEnd ;
	Grid_1D_float *bpml_zEnd ;
	Grid_1D_float *apml_half_zEnd;
	Grid_1D_float *bpml_half_zEnd ;

	// pml memory variable
	Grid_1D_float *mem_pr_zBeg ;
	Grid_1D_float *mem_vz_zBeg ;
	Grid_1D_float *mem_pr_zEnd ;
	Grid_1D_float *mem_vz_zEnd ;

	// sponge coefficient
	Grid_1D_float *sponge_coef_pr ;
	Grid_1D_float *sponge_coef_vz ;

} ;

} // namespace django

#endif
