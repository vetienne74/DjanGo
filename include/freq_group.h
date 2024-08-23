#ifndef DJANGO_FREQ_GROUP_H_
#define  DJANGO_FREQ_GROUP_H_

//
// configuration parameters
//

#include "type_def.h"

namespace django {

class Freq_group
{  
public:

	// destructor
	~Freq_group() ;

	// constructor
	Freq_group(Myfloat, Myfloat, Myfloat) ;

	// nb frequency groups
	Myint nb_freq ;

	// list of frequencies
	Myfloat* pFreq_list ;
} ;

} // namespace django

#endif
