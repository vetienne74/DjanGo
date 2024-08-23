#ifndef DJANGO_GRID_1D_H_
#define DJANGO_GRID_1D_H_

#include "grid.h"

#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Grid_1D: public Grid

{
public:

	// nb grid point along z
	Myint nz ;

	// grid point spacing along axis z
	Myfloat dz ;

	// get nb of grid points
	virtual Myint64 get_nb_grid_point(void) ;

	// print info
	virtual void info(void) ;

} ;

} // namespace django

#endif
