#ifndef DJANGO_DOMAIN_H_
#define DJANGO_DOMAIN_H_

#include "constant.h"
#include "model.h"
#include "scheme.h"
#include "scheme_factory.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Domain

{
public:

	Domain(void) ;
	Rtn_code initialize(void) ;
	Rtn_code info(void) ;

	Space_dim dim ;
	Model* pModel ;
	Scheme* pScheme ;
	Scheme_Factory* pScheme_Factory ;

} ;

} // namespace django

#endif
