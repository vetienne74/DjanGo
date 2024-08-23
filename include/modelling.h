#ifndef DJANGO_MODELLING_H_
#define DJANGO_MODELLING_H_

#include <fstream>

#include "acquisition.h"
#include "data.h"
#include "freq_group.h"
#include "grid.h"
#include "model.h"
#include "scheme.h"
#include "snapshot.h"
#include "type_def.h"
#include "variable.h"

namespace django {

//------------------------------------------------------------------------------------

class Modelling

{
public:

	// constructor
	Modelling(void) ;

	// initialize modelling
	Rtn_code initialize(Model*) ;

	// finalize modelling
	Rtn_code finalize(void) ;

	// solve forward problem for all shots
	Rtn_code solve_all_shot(Data*, Model*) ;

	// set frequency list to process
	Rtn_code set_freq_list_to_extract(Freq_group* pFreq_group) ;

	// get frequency groups
	Freq_group* get_freq_group(void) ;

	// ouput initial source wavelet
	Rtn_code output_src_wavelet(void) ;

	// ouput initial source wavelet
	Rtn_code output_estimated_src_wavelet(Mycomplex*) ;

	// set time parameters
	Rtn_code set_time_param(Myfloat tmax, Myfloat dt) ;

	// set snapshot parameters
	Rtn_code set_snapshot_time_param(Myfloat tmin, Myfloat tmax, Myfloat dt) ;
	Rtn_code set_snapshot_pixel_param(Myfloat xmin, Myfloat xmax, Myint nx, Myfloat zmin, Myfloat zmax, Myint nz) ;

	// set acquisition
	Rtn_code set_acquisition(string file_name) ;

	// set frequency parameters
	Rtn_code set_freq_param(Myfloat fmin, Myfloat fmax, Myfloat df) ;

	// set modelling dt
	Rtn_code set_modelling_dt(Myfloat) ;

	// set source param
	Rtn_code set_src_param(Myfloat src_freq_in, Src_func src_func_in, Src_type src_type_in,
			Src_stype src_stype_in, string src_file_in, Myfloat src_dt_in,
			Myfloat src_t0_in, Myfloat src_sigma_in) ;

	// get source time function
	Myfloat* get_src_time_function(void) ;

	// set energy param
	Rtn_code set_energy_param(bool) ;

	// set nt modelling
	void set_nt(Myint) ;

	// set dt modelling
	void set_dt(Myfloat) ;

	// set modelling case
	void set_case(Modelling_case_type modelling_case) ;

	// set modelling parameter
	void set_param(Myfloat param) ;

	// get modelling case
	Modelling_case_type get_case(void) ;

	// get modelling case
	Myfloat get_param(void) ;

	// display info
	Rtn_code info(void) ;

	// scheme
	Scheme* pScheme ;

	// acquisition
	Acquisition* pAcquisition ;

	// snaphot
	Snapshot* pSnapshot ;

protected:

	// time parameters
	Myfloat tmax ;
	Myfloat dt_out ;
	Myfloat dt ;

	// modelling case
	Modelling_case_type modelling_case ;

	// modelling parameter
	Myfloat modelling_param ;

	// snapshot parameters
	bool    output_snapshot ;
	Myfloat tmin_snapshot ;
	Myfloat tmax_snapshot ;
	Myfloat dt_snapshot ;

	// frequency parameters
	bool    output_frequency ;
	Myfloat fmin ;
	Myfloat fmax ;
	Myfloat df ;
	Myint   nf ;

	// nb of time step for modelling
	Myint nt ;

	// nb of time step for output data
	Myint nt_out ;

	// frequency group
	Freq_group* pFreq_group ;

	// Source dominant frequency
	Myfloat src_freq ;

	// Source time function
	Src_func src_func ;

	// Source type
	Src_type src_type ;

	// Source subtype
	Src_stype src_stype ;

	// Source file name
	string src_file ;

	// Source file dt
	Myfloat src_dt ;

	// Source t0
	Myfloat src_t0 ;

	// Source gaussian smoothing
	Myfloat src_sigma ;

	// perc_completion
	Myint perc_completion ;

	// compute energy
	bool compute_energy_flag ;

	// CPML reflexion coef.
	Myfloat cpml_rcoef ;

	// sponge coef.
	Myfloat sponge_coef ;

	// source time function
	Myfloat* src_time_function ;

	// free source time function
	Rtn_code free_src_time_function(void) ;

} ;

} // namespace django
#endif
