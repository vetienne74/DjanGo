#ifndef DJANGO_FEM_2D_1ST_AC_LOSSY_H_
#define DJANGO_FEM_2D_1ST_AC_LOSSY_H_

#include "FEM_2D_1st.h"

#include "acquisition.h"
#include "data.h"
#include "grid.h"
#include "model.h"
#include "snapshot.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class FEM_2D_1st_ac_lossy: public FEM_2D_1st

{
public:

	FEM_2D_1st_ac_lossy(void) ;

	// initialize scheme
	Rtn_code initialize(Model*, Myfloat fmax) ;

	// reset scheme
	Rtn_code reset(void) ;

	// finalize scheme
	Rtn_code finalize(void) ;

	// solve forward problem for the current shot
	Rtn_code solve_current_shot(Acquisition*, Data*, Wavefield_type, Snapshot*)  ;

	// compute pressure
	Rtn_code compute_pressure(Myfloat* pr, Myfloat* prp, Myfloat* vz, Myfloat* vx) ;

	// compute velocity
	Rtn_code compute_velocity(Myfloat* pr, Myfloat* prp, Myfloat* vz, Myfloat* vx) ;

	// initialization for eigen mode
	void init_eigen(void) ;

	// dynamic p-adaptivity
	Rtn_code dynamic_adapt_freq(Myint it, Acquisition* pAcquisition, Snapshot* pSnapshot) ;

	// compute energy
	void compute_energy(void) ;

protected:

	// variable id, used for projection from one mesh to another
	Myint varPrId, varPrpId, varVzId, varVxId ;

} ;

} // namespace django

#endif
