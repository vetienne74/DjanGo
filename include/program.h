#ifndef DJANGO_PROGRAM_H_
#define DJANGO_PROGRAM_H_

#include "domain.h"
#include "inversion.h"
#include "modelling.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Program

{
public:

	Program(Prog_type) ;
	virtual Rtn_code initialize(void) = 0 ;
	virtual Rtn_code run(void) = 0 ;
	virtual Rtn_code finalize(void) = 0 ;
	virtual Rtn_code info(void) ;

	Prog_type type ;
	Domain* pDomain ;
	Modelling* pModelling ;
	Inversion* pInversion ;

} ;

} // namespace django

#endif
