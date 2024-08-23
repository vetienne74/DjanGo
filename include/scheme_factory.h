#ifndef DJANGO_SCHEME_FACTORY_H_
#define DJANGO_SCHEME_FACTORY_H_

#include <vector>

#include "boundary.h"
#include "constant.h"
#include "freq_group.h"
#include "model.h"
#include "scheme.h"
#include "type_def.h"
#include "variable.h"

namespace django {

//------------------------------------------------------------------------------------

class Scheme_Factory

{
public:

	// constructor
	Scheme_Factory(void) ;

	// create scheme
	Scheme* create(void) ;

	// display info
	Rtn_code info(void) ;

	// set model
	Rtn_code set_model(Model*) ;


protected:

	// model
	Model* pModel ;

} ;

} // namespace django

#endif

