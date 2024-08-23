#ifndef DJANGO_PROGRAM_MODELLING_H_
#define DJANGO_PROGRAM_MODELLING_H_

#include "acquisition.h"
#include "data.h"
#include "domain.h"
#include "model.h"
#include "program.h"
#include "scheme.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Program_Modelling : public Program

{
public:

	Program_Modelling(Prog_type) ;
	Rtn_code initialize(void) ;
	Rtn_code run(void) ;
	Rtn_code finalize(void) ;
	Rtn_code info(void) ;

protected:
	Model* pModel ;
	Data*  pData ;

} ;

} // namespace django

#endif
