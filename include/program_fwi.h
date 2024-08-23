#ifndef DJANGO_PROGRAM_FWI_H_
#define DJANGO_PROGRAM_FWI_H_


#include "domain.h"
#include "program.h"
#include "scheme.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Program_Fwi : public Program

{
public:

	Program_Fwi(Prog_type) ;
	Rtn_code initialize(void) ;
	Rtn_code run(void) ;
	Rtn_code finalize(void) ;
	Rtn_code info(void) ;

} ;

} // namespace django

#endif
