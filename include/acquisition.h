#ifndef DJANGO_ACQUISITION_H_
#define DJANGO_ACQUISITION_H_

#include <string>

#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Acquisition
{

public:

	// constructor
	Acquisition(string, Space_dim) ;

	// destructor
	~Acquisition(void) ;

	// acquisition type
	Acqui_type acquisition_type ;

	// acquisition dimension
	Space_dim acquisition_dim ;

	// nb shots
	Myint nsrc ;

	// current shot id
	Myint current_src ;

	// current shot coordinates ;
	Myfloat xsrc ;
	Myfloat ysrc ;
	Myfloat zsrc ;

	// source function
	Myfloat* src_function ;

	// nb receivers
	Myint nrec ;

	// current receiver coordinates ;
	Myfloat *xrec ;
	Myfloat *yrec ;
	Myfloat *zrec ;

	// initialize acquisition
	Rtn_code initialize(void) ;

	// read src and rec positions for current shot
	virtual Rtn_code get_position_for_current_shot(void) ;

	// flag to know if it remains some shot to compute
	bool remain_shot_to_compute(void) ;

	// move to next shot
	Rtn_code move_to_next_shot(void) ;

	// move to first shot
	Rtn_code move_to_first_shot(void) ;

	// write rec in VTK format
	Rtn_code write_rec_VTK(void) ;

	// write src in VTK format
	Rtn_code write_src_VTK(void) ;

private:

	// acquisition file name
	string acquisition_file_name ;

} ;

}  // namespace django

#endif
