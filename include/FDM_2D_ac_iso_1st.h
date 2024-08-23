#ifndef DJANGO_FDM_2D_AC_ISO_1ST_H_
#define DJANGO_FDM_2D_AC_ISO_1ST_H_

#include "FDM_2D_ac_iso.h"

#include "acquisition.h"
#include "data.h"
#include "model.h"
#include "grid.h"
#include "grid_2D_float.h"
#include "snapshot.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class FDM_2D_ac_iso_1st: public FDM_2D_ac_iso
{
public:

	// constructor
	FDM_2D_ac_iso_1st(void) ;

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
	virtual Rtn_code update_physical_coef(Model* pModel) ;

protected:

	// 1st coef in eq. = dt * vp * vp* rho
	Grid_2D_float* coef1 ;

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
