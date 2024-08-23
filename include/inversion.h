#ifndef DJANGO_INVERSION_H_
#define DJANGO_INVERSION_H_

#include "acquisition.h"
#include "data.h"
#include "gradient.h"
#include "model.h"
#include "modelling.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Inversion

{
public:

	// constructor
	Inversion(Modelling*) ;

	// initialize inversion
	virtual Rtn_code initialize(void) ;

	// perform FWI (Full waveform inversion)
	virtual Rtn_code perform_fwi(Acquisition*, Data*, Model*) ;

	// setters
	void set_niter(Myint) ;
	void set_ntry(Myint) ;
	void set_init_try(Myfloat) ;
	void set_freq_inv(Myfloat, Myfloat, Myfloat) ;

private:

	// modelling
	Modelling* pModelling ;

	// gradient
	Gradient* pGradient ;

	// compute fcost frequency domain
	Rtn_code compute_fcost_freq_domain(Acquisition*, Data*, Model*, Myfloat*) ;

	// inversion configuration file
	ifstream config_file ;

	// read inversion parameters
	virtual Rtn_code read_config(void) ;

	// open configuration file
	void open_config_file(const char*) ;

	// open configuration file
	void close_config_file(void) ;

	// no. iterations
	Myint niter ;

	// no. try for line search
	Myint ntry ;

	// initial step length for line search
	Myfloat init_try ;

	// frequency group
	Freq_group* pFreq_group ;

} ;

} // namespace django

#endif
