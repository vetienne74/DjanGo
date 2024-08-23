#ifndef DJANGO_FDM_1D_AC_LOSSY_1ST_H_
#define DJANGO_FDM_1D_AC_LOSSY_1ST_H_

#include "FDM_1D_ac_lossy.h"

#include "acquisition.h"
#include "data.h"
#include "grid.h"
#include "grid_1D_float.h"
#include "model.h"
#include "snapshot.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class FDM_1D_ac_lossy_1st: public FDM_1D_ac_lossy
{
public:

	// constructor
	FDM_1D_ac_lossy_1st(void) ;

	// initialize scheme
	Rtn_code initialize(Model*, Myfloat fmax) ;

	// reset scheme
	Rtn_code reset(void) ;

	// finalize scheme
	Rtn_code finalize(void) ;

	// solve forward problem for the current shot
	Rtn_code solve_current_shot(Acquisition*, Data*, Wavefield_type, Snapshot*)  ;

	// update model
	virtual Rtn_code update_physical_coef(Myint iz, Myint ix, Myint iy, Myfloat perturb) ;
	virtual Rtn_code update_physical_coef(Model* pModel) ;

protected:

	// 1st coef in eq. = dt * vp * vp* rho
	Grid_1D_float *coef1 ;

	// 2nd coef in eq. = dt / rho
	Myfloat coef2 ;

	// compute pressure component
	virtual void compute_pressure(void) = 0 ;

	// compute velocity component
	virtual void compute_velocity(void) = 0 ;

	// initialization for eigen mode
	void init_eigen(void) ;

} ;

} // namespace django

#endif
