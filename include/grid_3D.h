#ifndef DJANGO_GRID_3D_H_
#define DJANGO_GRID_3D_H_

#include "grid_2D.h"

#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Grid_3D: public Grid_2D

{
public:

	// nb grid point along y
	Myint ny ;

	// grid point spacing along axis y
	Myfloat dy ;

	// get nb of grid points
	virtual Myint64 get_nb_grid_point(void) ;

	// print info
	virtual void info(void) ;

} ;

} // namespace django

#endif
