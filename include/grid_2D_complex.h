#ifndef DJANGO_GRID_2D_COMPLEX_H_
#define DJANGO_GRID_2D_COMPLEX_H_

#include "grid_2D.h"

#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Grid_2D_complex: public Grid_2D

{
public:

	// constructor
	Grid_2D_complex(Myint, Myint, Myfloat, Myfloat) ;

	// grid
	Mycomplex **pArray ;

	// get minimum value
	Myfloat get_min(void) ;

	// get minimum value
	Myfloat get_max(void) ;

	// print info
	virtual void info(void) ;

	// reset grid
	virtual Rtn_code reset(void) ;
	virtual Rtn_code reset(Myfloat) ;

	// write on disk
	Rtn_code write_on_disk(string) ;

} ;

} // namespace django

#endif
