#ifndef DJANGO_FDM_2D_EL_ISO_1ST_H_
#define DJANGO_FDM_2D_EL_ISO_1ST_H_

#include "FDM_2D_el_iso.h"

#include "acquisition.h"
#include "data.h"
#include "grid.h"
#include "grid_2D_float.h"
#include "model.h"
#include "snapshot.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class FDM_2D_el_iso_1st: public FDM_2D_el_iso
{
public:

	// constructor
	FDM_2D_el_iso_1st(void) ;

	// initialize modelling
	Rtn_code initialize(Model*, Myfloat fmax) ;

	// reset modelling
	Rtn_code reset(void) ;

	// finalize modelling
	Rtn_code finalize(void) ;

	// solve forward problem for the current shot
	Rtn_code solve_current_shot(Acquisition*, Data*, Wavefield_type, Snapshot*)  ;

	// compute energy
	void compute_energy(Variable*, Variable*, Variable*, Variable*, Variable*, Myint) ;

protected:

	// eq. coef1 = dt * (1/rho) at vx location
	Grid_2D_float* coef1 ;

	// eq. coef2 = dt * (1/rho) at vz location
	Grid_2D_float* coef2 ;

	// eq. coef3 = dt * (lambda + 2mu) at sxx and szz location
	Grid_2D_float* coef3 ;

	// eq. coef4 = dt * lambda at sxx and szz location
	Grid_2D_float* coef4 ;

	// eq. coef5 = dt * mu at sxz location
	Grid_2D_float* coef5 ;

	// compute stress component
	virtual void compute_stress(void) = 0 ;

	// compute velocity component
	virtual void compute_velocity(void) = 0 ;

} ;

} // namespace django

#endif
