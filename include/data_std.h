#ifndef DJANGO_DATA_SEISMIC_H_
#define DJANGO_DATA_SEISMIC_H_

#include <string>

#include "acquisition.h"
#include "data.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Data_std: public Data
{
public:

	// computed / observed / adjoint term at receivers (pressure / freq. domain)
	Mycomplex **pr_freq_sol_rec ;
	Mycomplex **pr_freq_obs_rec ;
	Mycomplex **pr_freq_adj_rec ;

	// computed / observed / adjoint term at receivers (pressure / time domain)
	Myfloat   **pr_time_sol_rec ;
	Myfloat   **pr_time_obs_rec ;
	Myfloat   **pr_time_adj_rec ;

	// read data
	virtual Rtn_code read(Acquisition*) ;

	// initialize data management
	virtual Rtn_code initialize(Acquisition*, Prog_type, Domain_type, Myint) ;

	// finalize data management
	virtual Rtn_code finalize(Acquisition*, Prog_type) ;

	// reset data management
	virtual Rtn_code reset(Acquisition*, Prog_type) ;

	// read data configuration file
	virtual Rtn_code read_config(void) ;

	// data file name (pressure)
	string pr_data_filename ;

} ;

} // namespace django

#endif
