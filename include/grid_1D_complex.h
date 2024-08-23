#ifndef DJANGO_GRID_1D_COMPLEX_H_
#define DJANGO_GRID_1D_COMPLEX_H_

#include "grid_1D.h"

#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Grid_1D_complex: public Grid_1D

{
public:

	// constructor
	Grid_1D_complex(Myint, Myfloat) ;

	// grid
	Mycomplex *pArray ;

	// get minimum value
	Myfloat get_min(void) ;

	// get minimum value
	Myfloat get_max(void) ;

	// print info
	virtual void info(void) ;

	// write on disk
	Rtn_code write_on_disk(string) ;

} ;

} // namespace django

#endif
