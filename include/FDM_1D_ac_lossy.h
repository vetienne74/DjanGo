#ifndef DJANGO_FDM_1D_AC_LOSSY_H_
#define DJANGO_FDM_1D_AC_LOSSY_H_

#include "FDM_1D.h"

#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class FDM_1D_ac_lossy: public FDM_1D
{

public:

	Rtn_code set_string_length_vs_time_step(Myfloat*) ;

protected:

	// constructor
	FDM_1D_ac_lossy(void) ;

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

	// length of string vs time step
	Myfloat* string_length_vs_time_step ;
	Myint izBeg1_initial ;
	Myint izBeg2_initial ;

	// loss terms
	Grid_1D_float *loss1 ;
	Grid_1D_float *loss2 ;

} ;

} // namespace django

#endif
