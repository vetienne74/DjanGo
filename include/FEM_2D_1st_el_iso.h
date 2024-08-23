#ifndef DJANGO_FEM_2D_1ST_EL_ISO_H_
#define DJANGO_FEM_2D_1ST_EL_ISO_H_

#include "FEM_2D_1st.h"

#include "acquisition.h"
#include "data.h"
#include "grid.h"
#include "model.h"
#include "snapshot.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class FEM_2D_1st_el_iso: public FEM_2D_1st

{
public:

	FEM_2D_1st_el_iso(void) ;

	// initialize scheme
	Rtn_code initialize(Model*, Myfloat fmax) ;

	// reset scheme
	Rtn_code reset(void) ;

	// finalize scheme
	Rtn_code finalize(void) ;

	// solve forward problem for the current shot
	Rtn_code solve_current_shot(Acquisition*, Data*, Wavefield_type, Snapshot*)  ;

	// compute stress
	Rtn_code compute_stress(Myfloat* tau, Myfloat* tauP, Myfloat* tauPP, Myfloat* vz, Myfloat* vx) ;

	// compute velocity
	Rtn_code compute_velocity(Myfloat* tau, Myfloat* tauP, Myfloat* tauPP, Myfloat* vz, Myfloat* vx) ;

	// initialization for eigen mode
	void init_eigen(void) ;

	// compute energy
	void compute_energy(Variable*, Variable*, Variable*, Variable*, Variable*, Myint) ;


protected:

} ;

} // namespace django

#endif
