#ifndef DJANGO_FDM_1D_AC_ISO_2ND_H_
#define DJANGO_FDM_1D_AC_ISO_2ND_H_

#include "FDM_1D_ac_iso.h"

#include "acquisition.h"
#include "data.h"
#include "grid_1D_float.h"
#include "snapshot.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class FDM_1D_ac_iso_2nd: public FDM_1D_ac_iso
{
public:

	// constructor
	FDM_1D_ac_iso_2nd(void) ;

	// initialize modelling
	Rtn_code initialize(Model*, Myfloat fmax) ;

	// reset modelling
	Rtn_code reset(void) ;

	// finalize modelling
	Rtn_code finalize(void) ;

	// solve forward problem for the current shot
	Rtn_code solve_current_shot(Acquisition*, Data*, Wavefield_type, Snapshot*)  ;

	// update model
	virtual Rtn_code update_physical_coef(Myint iz, Myint ix, Myint iy, Myfloat perturb) ;
	//virtual Rtn_code update_physical_coef(Model* pModel) ;

protected:

	// 1st coef in eq. = dt^2 * vp^2 / dt^2
	Grid_1D_float *coef1 ;

	// sponge coefficient
	Grid_1D_float *sponge_coef_pr_z ;

	// compute pressure component
	virtual void compute_pressure(void) = 0 ;

	// initialization for eigen mode
	void init_eigen(void) ;

} ;

} // namespace django

#endif
