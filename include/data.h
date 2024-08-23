#ifndef DJANGO_DATA_H_
#define DJANGO_DATA_H_

#include "acquisition.h"
#include "freq_group.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Data
{
public:

	// read data
	virtual Rtn_code read(Acquisition*) = 0 ;

	// initialize data management
	virtual Rtn_code initialize(Acquisition*, Prog_type, Domain_type, Myint) = 0 ;

	// finalize data management
	virtual Rtn_code finalize(Acquisition*, Prog_type) = 0 ;

	// reset data management
	virtual Rtn_code reset(Acquisition*, Prog_type) = 0 ;

	// read data configuration file
	virtual Rtn_code read_config(void) = 0 ;

protected:

	// data domain (time or frequency)
	Domain_type domain_in, domain_out ;

	// no. time step
	Myint nt_in, nt_out;

	// time step
	Myfloat dt_in, dt_out ;

	// frequency group
	Freq_group* pFreq_group_in, *pFreq_group_out ;
	Myint nf_in, nf_out ;

} ;

} // namespace django

#endif
