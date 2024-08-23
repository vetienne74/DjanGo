#ifndef DJANGO_GRID_3D_FLOAT_H_
#define DJANGO_GRID_3D_FLOAT_H_

#include "grid_3D.h"

#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Grid_3D_float: public Grid_3D

{
public:

	// constructor
	Grid_3D_float(Myint, Myint, Myint, Myfloat, Myfloat, Myfloat) ;

	// destructor
	~Grid_3D_float(void) ;

	// grid
	Myfloat ***pArray ;

	// read from disk
	Rtn_code read_from_disk(string) ;

	// write on disk
	Rtn_code write_on_disk(string) ;

	// get minimum value
	Myfloat get_min(void) ;

	// get minimum value
	Myfloat get_max(void) ;

	// print info
	virtual void info(void) ;

} ;

} // namespace django

#endif
