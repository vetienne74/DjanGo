#ifndef DJANGO_GRID_3D_COMPLEX_H_
#define DJANGO_GRID_3D_COMPLEX_H_

#include "grid_3D.h"

#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Grid_3D_complex: public Grid_3D

{
public:

	// constructor
	Grid_3D_complex(Myint, Myint, Myint, Myfloat, Myfloat, Myfloat) ;

	// destructor
	~Grid_3D_complex(void) ;

	// grid
	Mycomplex ***pArray ;

	// get minimum value
	Myfloat get_min(void) ;

	// get minimum value
	Myfloat get_max(void) ;

	// print info
	virtual void info(void) ;

	// reset grid
	virtual Rtn_code reset(void) ;
	virtual Rtn_code reset(Myfloat val) ;

	// write on disk
	Rtn_code write_on_disk(string) ;

} ;

} // namespace django

#endif
