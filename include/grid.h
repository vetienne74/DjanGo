#ifndef DJANGO_GRID_H_
#define DJANGO_GRID_H_

#include <string>

#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Grid

{
public:

	// pointer
	void* pArray  ;

	// virtual destructor
	virtual ~Grid(void) = 0 ;

	// get nb of grid points
	virtual Myint64 get_nb_grid_point(void) = 0 ;

	// print info
	virtual void info(void) ;

	// reset grid
	virtual Rtn_code reset(void) = 0 ;
	virtual Rtn_code reset(Myfloat val) = 0 ;

	// write grid
	virtual Rtn_code write_on_disk(string filename) = 0 ;

} ;

} // namespace django

#endif
