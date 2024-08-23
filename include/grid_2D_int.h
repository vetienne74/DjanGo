#ifndef DJANGO_GRID_2D_INT_H_
#define DJANGO_GRID_2D_INT_H_

#include "grid_2D.h"

#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Grid_2D_int: public Grid_2D

{
public:

	// constructor
	Grid_2D_int(Myint, Myint, Myfloat, Myfloat) ;

	// destructor
	~Grid_2D_int(void) ;

	// grid
	Myint **pArray ;

	// read from disk
	Rtn_code read_from_disk(string) ;

	// write on disk
	Rtn_code write_on_disk(string) ;

	// get minimum value
	Myint get_min(void) ;

	// get minimum value
	Myint get_max(void) ;

	// print info
	virtual void info(void) ;

	// reset grid
	virtual Rtn_code reset(void) ;
	virtual Rtn_code reset(Myfloat) ;

} ;

} // namespace django

#endif
