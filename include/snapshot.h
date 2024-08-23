#ifndef DJANGO_SNAPSHOT_H_
#define DJANGO_SNAPSHOT_H_

#include "grid_1D_float.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Snapshot

{
public:

	// construtor
	Snapshot(void) ;

	// destructor
	~Snapshot(void) ;

	Rtn_code set_time_param(Myfloat tmin, Myfloat tmax, Myfloat dt) ;
	Rtn_code set_pixel_param(Myfloat xmin, Myfloat xmax, Myint nx, Myfloat zmin, Myfloat zmax, Myint nz) ;
	Rtn_code info(void) ;

	// tmin, tmax and dt
	Myfloat tmin, tmax, dt ;

	// pixel coordinates
	Grid_1D_float* x_coord ;
	Grid_1D_float* z_coord ;

protected:

} ;

} // namespace django

#endif
