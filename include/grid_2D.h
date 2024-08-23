#ifndef DJANGO_GRID_2D_H_
#define DJANGO_GRID_2D_H_

#include "grid_1D.h"

#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Grid_2D: public Grid_1D

{
public:

	// nb grid point along x
	Myint nx ;

	// grid point spacing along axis x
	Myfloat dx ;

	// get nb of grid points
	virtual Myint64 get_nb_grid_point(void) ;

	// print info
	virtual void info(void) ;

} ;

} // namespace django

#endif
