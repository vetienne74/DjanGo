#ifndef DJANGO_GRID_1D_FLOAT_H_
#define DJANGO_GRID_1D_FLOAT_H_

#include "grid_1D.h"

#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Grid_1D_float: public Grid_1D

{
public:

	// constructor (nelem, incremenent)
	Grid_1D_float(Myint, Myfloat) ;

	// constructor (nelem, start value, end value)
	Grid_1D_float(Myint, Myfloat, Myfloat) ;

	// destructor
	~Grid_1D_float(void) ;

	// grid
	Myfloat *pArray ;

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

	// reset grid
	virtual Rtn_code reset(void) ;
	virtual Rtn_code reset(Myfloat) ;

} ;

} // namespace django

#endif
