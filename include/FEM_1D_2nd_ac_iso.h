#ifndef DJANGO_FEM_1D_2ND_AC_ISO_H_
#define DJANGO_FEM_1D_2ND_AC_ISO_H_

#include "FEM_1D_2nd.h"

#include "acquisition.h"
#include "data.h"
#include "grid.h"
#include "model.h"
#include "snapshot.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class FEM_1D_2nd_ac_iso: public FEM_1D_2nd

{
public:

	FEM_1D_2nd_ac_iso(void) ;

	// initialize scheme
	Rtn_code initialize(Model*, Myfloat fmax) ;

	// reset scheme
	Rtn_code reset(void) ;

	// finalize scheme
	Rtn_code finalize(void) ;

	// solve forward problem for the current shot
	Rtn_code solve_current_shot(Acquisition*, Data*, Wavefield_type, Snapshot*)  ;

	// compute pressure
	Rtn_code compute_pressure(Myfloat* prn, Myfloat* prc) ;

	// initialization for eigen mode
	void init_eigen(void) ;

protected:

} ;

} // namespace django

#endif
