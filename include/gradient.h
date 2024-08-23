#ifndef DJANGO_GRADIENT_H_
#define DJANGO_GRADIENT_H_

#include "acquisition.h"
#include "data.h"
#include "freq_group.h"
#include "grid.h"
#include "model.h"
#include "modelling.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Gradient

{
public:

	// constructors
	Gradient(Acquisition*, Model*, Freq_group*) ;

	// destructor
	~Gradient(void) ;

	// compute gradient of objective function
	Rtn_code compute(Acquisition* , Data*, Modelling*, Model*, Myfloat*) ;

	// initialize modelling
	Rtn_code initialize(void) ;

	// update gradient
	Rtn_code update(Grid*, Grid*, Model*, Modelling*) ;

	// gather gradient contributions
	Rtn_code gather(void) ;

	// store value
	Rtn_code store_value(Myint iz, Myint ix, Myint iy, Myfloat grad_val) ;

	// gradient
	Model* pGrad_model ;

	// preconditionner
	Model* pPrecond ;

	// gradient preconditionned
	Model* pGrad_precond ;

	// gradient domain
	Domain_type domain ;

	// evaluate cost function, src wavelet and build corresponding adjoint source
	Rtn_code evaluate_fcost_wavelet_and_adjoint_source(Acquisition*, Data*, Modelling*, Myfloat*) ;

	// evaluate cost function
	Rtn_code evaluate_fcost(Acquisition*, Data*, Modelling*, Myfloat*) ;

	// compute FWI gradient time domain
	Rtn_code compute_gradient_time_domain(Acquisition*, Data*, Model*, Modelling*) ;

	// compute FWI gradient frequency domain with adjoint state method
	Rtn_code compute_gradient_freq_domain_with_adjoint_state_method(Acquisition*, Data*, Model*, Modelling*, Myfloat*) ;

	// compute FWI gradient frequency domain with fd of cost function
	Rtn_code compute_gradient_freq_domain_with_fd_fcost(Acquisition*, Data*, Model*, Modelling*) ;

private:

	// gradient configuration file
	ifstream config_file ;

	// read gradient parameters
	virtual Rtn_code read_config(void) ;

	// open configuration file
	void open_config_file(const char*) ;

	// open configuration file
	void close_config_file(void) ;

	// gradient type
	Grad_type type ;

	// fcost type
	Fcost_type fcost_type ;

	// source estimation flag
	bool src_estim_flag ;

	// min and max offset for source estimation
	Myfloat src_estim_offset_min, src_estim_offset_max ;

	// min and max offset for gradient computation
	Myfloat grad_offset_min, grad_offset_max ;

	// coef for gradient preconditionning
	Myfloat coef_grad_precond ;

	// coef for source estimation
	Mycomplex **src_coef ;

	// frequency group
	Freq_group* pFreq_group ;

	// frequency group
	Acquisition* pAcquisition ;

} ;

} // namespace django

#endif
