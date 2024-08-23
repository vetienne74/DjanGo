//-------------------------------------------------------------------------------------------------------
//
// PARENT CLASS FOR ALL FORWARD MODELLING METHODS
//
//-------------------------------------------------------------------------------------------------------

#include "modelling.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>

#include "mpi.h"

#include "acquisition.h"
#include "allocate_array.h"
#include "data.h"
#include "freq_group.h"
#include "grid.h"
#include "output_report.h"
#include "parse_xml.h"
#include "singleton.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------

Modelling::Modelling(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Modelling::Modeling");

	// initialize members
	modelling_case      = MODELLING_STD ;
	modelling_param     = 0 ;

	tmax                = 0.0 ;
	dt                  = 0.0 ;
	dt_out              = 0.0 ;
	dt_snapshot         = 0.0 ;

	src_freq            = 0.0 ;
	src_func            = NO_SRC_FUNC ;
	src_type            = NO_SRC_TYPE ;
	src_stype           = NO_SRC_STYPE ;
	src_file            = "UNSPECIFIED" ;
	src_dt              = 0.0 ;
	src_t0              = 0.0 ;
	src_sigma           = 0.0 ;
	src_time_function   = NULL ;

	compute_energy_flag = false ;
	perc_completion     = 0 ;

	output_frequency    = false ;
	fmin                = 0.0 ;
	fmax                = 0.0 ;
	df                  = 0.0 ;
	nf                  = 0 ;

	output_snapshot     = false ;
	tmin_snapshot       = 0.0 ;
	tmax_snapshot       = 0.0 ;
	dt_snapshot         = 0.0 ;

	pAcquisition        = NULL ;
	pScheme             = NULL ;
	pSnapshot           = NULL ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Modelling::Modeling");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Modelling::set_acquisition(string file_name)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Modelling::set_acquisition");

	if (Singleton::Instance()->pProgram == NULL)
	{
		print_error(" IN Modelling::set_acquisition, pProgram not initialized") ;
		return(RTN_CODE_KO) ;
	}
	else
	{
		if (Singleton::Instance()->pProgram->pDomain == NULL)
		{
			print_error(" IN Modelling::set_acquisition, pDomain not initialized") ;
			return(RTN_CODE_KO) ;
		}
		else
		{
			if (Singleton::Instance()->pProgram->pDomain->pModel == NULL)
			{
				print_error(" IN Modelling::set_acquisition, pModel not initialized") ;
				return(RTN_CODE_KO) ;
			}
			else
			{
				Space_dim dim = Singleton::Instance()->pProgram->pDomain->pModel->get_dim() ;
				pAcquisition = new Acquisition(file_name, dim) ;
			}
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Modelling::set_acquisition");
	return (RTN_CODE_OK) ;
}


//-------------------------------------------------------------------------------------------------------
Rtn_code Modelling::set_time_param(Myfloat tmax_in, Myfloat dt_out_in) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Modelling::set_time_param");

	tmax        = tmax_in ;
	dt_out      = dt_out_in ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Modelling::set_time_param");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code Modelling::set_snapshot_time_param(Myfloat tmin_in, Myfloat tmax_in, Myfloat dt_in) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Modelling::set_snapshot_time_param");

	output_snapshot = true ;
	pSnapshot = new Snapshot ;
	pSnapshot->set_time_param(tmin_in, tmax_in, dt_in) ;
	tmin_snapshot   = tmin_in ;
	tmax_snapshot   = tmax_in ;
	dt_snapshot     = dt_in ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Modelling::set_snapshot_time_param");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code Modelling::set_snapshot_pixel_param(Myfloat xmin, Myfloat xmax, Myint nx,
		Myfloat zmin, Myfloat zmax, Myint nz)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Modelling::set_snapshot_pixel_param");
	pSnapshot->set_pixel_param(xmin,xmax,nx,zmin,zmax,nz) ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Modelling::set_snapshot_pixel_param");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
void Modelling::set_nt(Myint nt_in) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Modelling::set_nt");
	nt = nt_in ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Modelling::set_nt");
}

//-------------------------------------------------------------------------------------------------------
void Modelling::set_case(Modelling_case_type modelling_case_in) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Modelling::set_case");
	modelling_case = modelling_case_in ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Modelling::set_case");
}

//-------------------------------------------------------------------------------------------------------
void Modelling::set_param(Myfloat param_in) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Modelling::set_param");
	modelling_param = param_in ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Modelling::set_param");
}

//-------------------------------------------------------------------------------------------------------
Modelling_case_type Modelling::get_case(void) 
{
	print_debug(ALL, FULL_DEBUG, "IN Modelling::get_case");
	print_debug(ALL, FULL_DEBUG, "OUT Modelling::get_case");
	return(modelling_case) ;
}

//-------------------------------------------------------------------------------------------------------
Myfloat Modelling::get_param(void) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Modelling::get_param");
	print_debug(ALL, LIGHT_DEBUG, "OUT Modelling::get_param");
	return(modelling_param) ;
}

//-------------------------------------------------------------------------------------------------------
void Modelling::set_dt(Myfloat dt_in) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Modelling::set_dt");
	dt = dt_in ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Modelling::set_dt");
}

//-------------------------------------------------------------------------------------------------------
Rtn_code Modelling::set_freq_param(Myfloat fmin_in, Myfloat fmax_in, Myfloat df_in) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Modelling::set_freq_param");

	fmin = fmin_in ;
	fmax = fmax_in ;
	df   = df_in ;

	// init freq list

	print_debug(ALL, LIGHT_DEBUG, "OUT Modelling::set_freq_param");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code Modelling::set_modelling_dt(Myfloat dt_in) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Modelling::set_modelling_dt");
	dt = dt_in ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Modelling::set_modelling_dt");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code Modelling::set_energy_param(bool compute_energy_flag_in) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Modelling::set_energy_param");
	compute_energy_flag = compute_energy_flag_in ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Modelling::set_energy_param");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Modelling::set_src_param(Myfloat src_freq_in, Src_func src_func_in, Src_type src_type_in,
		Src_stype src_stype_in, string src_file_in, Myfloat src_dt_in,
		Myfloat src_t0_in, Myfloat src_sigma_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Modelling::set_src_param");

	// source frequency
	src_freq = src_freq_in ;

	// Source function ;
	src_func = src_func_in ;

	// Source type ;
	src_type = src_type_in ;

	// Source subtype ;
	src_stype = src_stype_in ;

	// Source file name
	src_file = src_file_in ;

	// Source file dt
	src_dt = src_dt_in ;

	// Source t0
	src_t0 = src_t0_in ;

	// source gaussian smoothing
	src_sigma = src_sigma_in ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Modelling::set_src_param");
	return(RTN_CODE_OK) ;
}


//-------------------------------------------------------------------------------------------------------

Rtn_code Modelling::free_src_time_function()
{
	print_debug(ALL, LIGHT_DEBUG, "IN Modelling::free_src_time_function");
	deallocate_array<Myfloat>(src_time_function, nt+1) ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Modelling::free_src_time_function");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat* Modelling::get_src_time_function()
{
	print_debug(ALL, LIGHT_DEBUG, "IN Modelling::get_src_time_function");

	if (src_func != NO_SRC_FUNC)
	{

		src_time_function = allocate_array<Myfloat>(nt+1) ;
		if (src_time_function == NULL)
		{
			print_error(" src_time_function is NULL") ;
			return(NULL) ;
		}

		print_info(MASTER, "") ;
		print_info(MASTER, " source initialization") ;

		Myfloat amp       = DEFAULT_AMP_SRC ;
		Myfloat max_amp   = 0. ;
		Myfloat delay_src = 0. ;

		// src function read from disk
		if (src_func == SRC_FILE)
		{
			// estimate nt required
			Myint nt_required = (nt * dt) / src_dt ;
			print_debug(ALL, MID_DEBUG, "nt required", nt_required) ;

			// read source file
			ifstream in_file ;
			in_file.open(src_file.c_str(), ios::binary) ;
			assert(in_file.is_open());

			// get length of file
			in_file.seekg (0, in_file.end);
			Myint nt_in_file = in_file.tellg() / sizeof(Myfloat) ;
			print_info(MASTER, " Nt in src file:", nt_in_file) ;
			in_file.seekg (0, in_file.beg);

			// determine # time step to read
			Myint nt_read = min(nt_in_file, nt_required) ;

			// allocate array
			print_info(MASTER, " Nt required:\t", nt_read) ;
			Myfloat* val_array = new Myfloat[nt_read] ;
			in_file.read((char*) val_array, nt_read * sizeof(Myfloat)) ;
			in_file.close() ;

			// interpolate
			for (Myint it = 0; it <= nt; it ++)
			{
				Myfloat t1 = it * dt ;
				Myint it_src = t1 / src_dt ;
				//cout << " * it_src " << it_src ;
				if (it_src < nt_read-1)
				{
					// nearest point
					//src_time_function[it] = val_array[it_src] ;

					// linear interpolation
					Myfloat64 t2 = it_src * src_dt ;
					Myfloat64 dd = t1 - t2 ;
					Myfloat64 a2 = dd / src_dt ;
					Myfloat64 a1 = 1. - a2 ;
					src_time_function[it] = a1*val_array[it_src] + a2*val_array[it_src+1];
					//cout << " a1 " << a1 << " a2 " << a2 ;
					//cout << " val_array[it_src] " << val_array[it_src] << "\n" ;
				}
				else if (it_src == nt_read-1)
				{
					src_time_function[it] = val_array[it_src] ;
				}
				else
				{
					src_time_function[it] = 0.0 ;
				}
				//cout << " * src_time_function " << src_time_function[it] << "\n" ;
			}
			delete(val_array) ;
		}

		else

		{

			// Ricker wavelet
			Myfloat da        = M_PI * src_freq ;
			Myfloat t00       = 1.5 * sqrt(6.) / da ;

			for (Myint it = 0; it <= nt; it ++)
			{

				Myfloat tt = it * dt ;
				Myfloat src_wavelet = 0.0 ;

				// Ricker wavelet
				Myfloat aa = M_PI * src_freq * (tt - t00 - src_t0) ;
				Myfloat a2 = pow(aa, 2.) ;

				if (src_func == RICKER_PP)
				{
					src_wavelet = -(0.5 / (da*da)) * exp(-a2) ;
				}

				else if (src_func == RICKER_P)
				{
					src_wavelet = (aa/da) * exp(-a2) ;
				}

				else if (src_func == RICKER)
				{
					src_wavelet = (1. - 2. * a2) * exp(-a2) ;
				}

				else if (src_func == RICKER_D)
				{
					src_wavelet = -4. * aa * da * exp(-a2) -2. * aa * da * (1. -2. * a2) * exp(-a2)	;
				}

				else if (src_func == MONO_FREQ)
				{
					src_wavelet = sin(2. * M_PI * src_freq * tt) ;
				}

				src_time_function[it] = amp * src_wavelet ;
			}
		}

		for (Myint it = 0; it <= nt; it ++)
		{
			// retrieve time when source is max
			if (max_amp < abs(src_time_function[it]))
			{
				max_amp = abs(src_time_function[it]) ;
				delay_src = it ;
			}
		}

		print_info(MASTER, " Max. src amp.:\t", max_amp) ;
		print_info(MASTER, " Max. src time (s):", (delay_src-1)*dt) ;

		// output src wavelet
		this->output_src_wavelet() ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Modelling::get_src_time_function");
	return(src_time_function) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Modelling::output_src_wavelet()
{
	print_debug(ALL, LIGHT_DEBUG, "IN Modelling::output_src_wavelet");

	ofstream out_file ;

	out_file.open(SRC_WAVELET_OUT_FILE, ios::binary) ;
	assert(out_file.is_open());
	out_file.write((char*) src_time_function, nt * sizeof(Myfloat)) ;
	out_file.close() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Modelling::output_src_wavelet");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Modelling::output_estimated_src_wavelet(Mycomplex* src_coef)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Modelling::output_estimated_src_wavelet");

	// convert src time wavelet to freq domain
	Mycomplex src_freq[pFreq_group->nb_freq] ;
	for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)
	{
		// convert data with DFT
		src_freq[ifreq] = ZERO_CMPLX ;
		for (Myint it = 0; it < nt; it++)
		{
			Myfloat coef_tmp = -2. * PI *  pFreq_group->pFreq_list[ifreq] * it * dt ;
			Mycomplex coef (cos(coef_tmp)*dt, sin(coef_tmp)*dt) ;
			src_freq[ifreq] += src_time_function[it] * coef ;
		}
	}

	// apply src correction on wavelet on freq domain
	Mycomplex estim_src_freq[pFreq_group->nb_freq] ;
	for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)
	{
		estim_src_freq[ifreq] = src_freq[ifreq] * src_coef[ifreq] ;
	}

	// rebuild estimated wavelet in time domain
	Myfloat estim_src_time_function[nt] ;
	Myfloat coef2 = 2. / nt ;

	for (Myint it = 0; it < nt; it++)
	{
		// convert data with IDFT
		estim_src_time_function[it] = 0. ;
		for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)
		{
			Myfloat coef_tmp = 2. * PI * pFreq_group->pFreq_list[ifreq] * it * dt;
			Mycomplex coef (cos(coef_tmp), sin(coef_tmp)) ;
			estim_src_time_function[it] += (estim_src_freq[ifreq] * coef).real() * coef2 ;
		}
	}

	// write time wavelet on disk
	ofstream out_file ;

	out_file.open(ESTIM_SRC_WAVELET_OUT_FILE, ios::binary) ;
	assert(out_file.is_open());
	out_file.write((char*) estim_src_time_function, nt * sizeof(Myfloat32)) ;
	out_file.close() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Modelling::output_estimated_src_wavelet");

	return(RTN_CODE_KO) ;
}


//-------------------------------------------------------------------------------------------------------

Freq_group* Modelling::get_freq_group()
{
	return(pFreq_group) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Modelling::initialize(Model* pModel)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Modelling::initialize");

	// check objects
	if (Singleton::Instance()->pProgram == NULL)
	{
		print_error(" Modelling::initialize, pProgram is NULL") ;
		return(RTN_CODE_KO) ;
	}
	if (Singleton::Instance()->pProgram->pDomain == NULL)
	{
		print_error(" Modelling::initialize, pDomain is NULL") ;
		return(RTN_CODE_KO) ;
	}
	if (Singleton::Instance()->pProgram->pDomain->pScheme == NULL)
	{
		print_error(" Modelling::initialize, pScheme is NULL") ;
		return(RTN_CODE_KO) ;
	}
	if (Singleton::Instance()->pProgram->pDomain->pModel == NULL)
	{
		print_error(" Modelling::initialize, pModel is NULL") ;
		return(RTN_CODE_KO) ;
	}

	// initialize acquisition
	Rtn_code rtn_code = pAcquisition->initialize() ;
	if (rtn_code != RTN_CODE_OK) return (rtn_code) ;

	// initialize scheme
	pScheme = Singleton::Instance()->pProgram->pDomain->pScheme ;

	rtn_code = pScheme->set_energy_param(compute_energy_flag) ;
	if (rtn_code != RTN_CODE_OK) return (rtn_code) ;

	rtn_code = pScheme->set_tmax(tmax) ;
	if (rtn_code != RTN_CODE_OK) return (rtn_code) ;

	rtn_code = pScheme->set_dt_out(dt_out) ;
	if (rtn_code != RTN_CODE_OK) return (rtn_code) ;

	rtn_code = pScheme->set_src_param(src_type, src_stype, src_sigma) ;
	if (rtn_code != RTN_CODE_OK) return (rtn_code) ;

	if (output_snapshot) {
		rtn_code = pScheme->set_snapshot_param(tmin_snapshot, tmax_snapshot, dt_snapshot) ;
		if (rtn_code != RTN_CODE_OK) return (rtn_code) ;
	}

	rtn_code = pScheme->initialize(pModel, Singleton::Instance()->fmax0) ;
	if (rtn_code != RTN_CODE_OK) return (rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Modelling::initialize");
	return(rtn_code) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Modelling::finalize()
{
	print_debug(ALL, LIGHT_DEBUG, "IN Modelling::finalize");

	// check objects
	if (Singleton::Instance()->pProgram == NULL)
	{
		print_error(" Modelling::finalize, pProgram is NULL") ;
		return(RTN_CODE_KO) ;
	}
	if (Singleton::Instance()->pProgram->pDomain == NULL)
	{
		print_error(" Modelling::finalize, pDomain is NULL") ;
		return(RTN_CODE_KO) ;
	}
	if (Singleton::Instance()->pProgram->pDomain->pScheme == NULL)
	{
		print_error(" Modelling::finalize, pScheme is NULL") ;
		return(RTN_CODE_KO) ;
	}
	if (Singleton::Instance()->pProgram->pDomain->pModel == NULL)
	{
		print_error(" Modelling::finalize, pModel is NULL") ;
		return(RTN_CODE_KO) ;
	}

	// delete member objects
	delete(pAcquisition) ;
	delete(pSnapshot) ;

	// finalize scheme
	pScheme = Singleton::Instance()->pProgram->pDomain->pScheme ;
	Rtn_code rtn_code = pScheme->finalize() ;
	if (rtn_code != RTN_CODE_OK) return (rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Modelling::finalize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Modelling::set_freq_list_to_extract(Freq_group* pFreq_group2) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Modelling::set_freq_list_to_extract");
	pFreq_group = pFreq_group2 ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Modelling::set_freq_list_to_extract");
	return(RTN_CODE_OK) ;
}


//=======================================================================================================
//
// SOLVE ALL SHOT
//
//=======================================================================================================

Rtn_code Modelling::solve_all_shot(Data* pData, Model* pModel) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Modelling::solve_all_shot");

	Grid* pIncident = NULL ;

	// initialize data
	Rtn_code rtn_code = pData->initialize(pAcquisition, MODELLING_PROG, FREQ, nf) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_info(MASTER, "") ;
	print_line2() ;
	print_info(MASTER, " START MODELLING") ;
	print_info(MASTER, "") ;

	double t0 = MPI_Wtime() ;

	//-------------------------------------------------------------------------------------------------------
	// loop on shots until all have been computed
	//-------------------------------------------------------------------------------------------------------

	while ( pAcquisition->remain_shot_to_compute() )

	{
		print_info(ALL, " SOLVE FOR SOURCE:", pAcquisition->current_src) ;

		// retrieve src and rec positions the current shot
		rtn_code = pAcquisition->get_position_for_current_shot() ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		// reset data
		rtn_code = pData->reset(pAcquisition, MODELLING_PROG) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		//-------------------------------------------------------------------------------------------------------
		// solve forward problem for the current shot
		//-------------------------------------------------------------------------------------------------------
		double t1 = MPI_Wtime() ;
		rtn_code = pScheme->solve_current_shot(pAcquisition, pData, INCIDENT, pSnapshot) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
		double t2 = MPI_Wtime() ;
		print_info(ALL, "") ;
		print_info(ALL, " Elapse time:\t", (Myfloat) (t2-t1)) ;

		// display stat
		pScheme->display_stat() ;

		// write stat on disk
		pScheme->log_stat() ;

		// increment shot counter
		rtn_code = pAcquisition->move_to_next_shot() ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
	}

	// finalize data
	rtn_code = pData->finalize(pAcquisition, MODELLING_PROG) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_info(MASTER, "") ;
	print_info(MASTER, " END MODELLING") ;
	print_info(MASTER, "") ;

	// if (nproc_world > 1)
	//   {
	//     // retrieve min and max computation time
	//     double min_comp_time, max_comp_time ;
	//     MPI::COMM_WORLD.Reduce(&time_in_kernel, &min_comp_time, 1, MPI_DOUBLE, MPI_MIN, 0) ;
	//     MPI::COMM_WORLD.Reduce(&time_in_kernel, &max_comp_time, 1, MPI_DOUBLE, MPI_MAX, 0) ;

	//     print_info(MASTER, " MIN. COMP. TIME:", (float) min_comp_time) ;
	//     print_info(MASTER, " MAX. COMP. TIME:", (float) max_comp_time) ;
	//     print_info(MASTER, " % DIFF COMP. TIME:", (float) ((max_comp_time-min_comp_time)/max_comp_time*100.)) ;
	//   }

	double t3 = MPI_Wtime() ;
	double total_time = t3 - t0 ;
	if (total_time < 60.0)
	{
		print_info(MASTER, " TOTAL TIME (sec):", (float) (total_time)) ;
	}
	else if (total_time < 3600.0)
	{
		print_info(MASTER, " TOTAL TIME (min):", (float) (total_time/60.)) ;
	}
	else
	{
		print_info(MASTER, " TOTAL TIME (hrs):", (float) (total_time/3600.)) ;
	}
	print_line2() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Modelling::solve_all_shot");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
// print info
// perform some checking (return error code)
//
Rtn_code Modelling::info(void)
{
	print_debug(ALL, FULL_DEBUG, "IN Modelling::info");

	print_info(MASTER, "") ;
	print_line2() ;
	print_info(MASTER, " MODELLING PARAMETERS") ;
	print_info(MASTER, "") ;

	// modelling case
	switch(modelling_case)
	{
	case MODELLING_STD:
		print_info(MASTER, " Modelling type:", "STANDARD") ;
		break;
	case MODELLING_EIGEN:
		print_info(MASTER, " Modelling type:", "EIGEN MODE") ;
		if ((Myint)modelling_param <= 0) modelling_param = 1 ;
		print_info(MASTER, " Number of modes:", (Myint) modelling_param) ;
		break;
	default:
		print_error(" Modelling case unknown", modelling_case) ;
		return(RTN_CODE_KO) ;
	}

	// tmax
	if (tmax <= 0)
	{
		print_error(" modelling time (s) should be > 0", tmax) ;
		return(RTN_CODE_KO) ;
	}
	else
	{
		print_info(MASTER, " Modelling time (s)", tmax) ;
	}

	// output dt
	if (dt_out < 0.)
	{
		print_error(" Ouput dt (s) should be > 0") ;
		return(RTN_CODE_KO) ;
	}
	else if (dt_out > tmax)
	{
		print_error(" Ouput dt (s) should be <= tmax") ;
		return(RTN_CODE_KO) ;
	}
	else
	{      
		print_info(MASTER, " Output dt (s)\t", dt_out) ;
		print_info(MASTER, " -> output rate (/s)", Myfloat(1.)/dt_out) ;
	}

	// freq parameters
	if (output_frequency)
	{
		print_info(MASTER, " Frequency extraction", "ON") ;
		if (fmin < 0.)
		{
			print_error(" Freq. min. should be > 0") ;
			return(RTN_CODE_KO) ;
		}
		else
		{
			print_info(MASTER, " Freq. min. (Hz)", fmin) ;
		}
		if (fmax < 0.)
		{
			print_error(" Freq. max. should be > 0") ;
			return(RTN_CODE_KO) ;
		}
		else
		{
			print_info(MASTER, " Freq. max. (Hz)", fmax) ;
		}
		if (df < 0.)
		{
			print_error(" Df should be > 0") ;
			return(RTN_CODE_KO) ;
		}
		else
		{
			print_info(MASTER, " Freq. df (Hz)\t", df) ;
		}
	}
	else
	{
		print_info(MASTER, " Frequency extraction", "OFF") ;
	}

	// snapshot
	if (output_snapshot)
	{
		print_info(MASTER, " Snapshot capture ", "ON") ;
		if (pSnapshot->info() != RTN_CODE_OK) return(RTN_CODE_KO) ;
	}
	else
	{
		print_info(MASTER, " Snapshot capture", "OFF") ;
	}

	// energy
	if (compute_energy_flag)
	{
		print_info(MASTER, " Compute energy ", "ON") ;
	}
	else
	{
		print_info(MASTER, " Compute energy ", "OFF") ;
	}

	print_info(MASTER, " Output time data", "ON") ;

	// source
	print_info(MASTER, "") ;

	// src function
	switch(src_func)
	{
	case RICKER_PP:
		print_info(MASTER, " Source function:", "RICKER 2ND PRIMITIVE") ;
		break;
	case RICKER_P:
		print_info(MASTER, " Source function:", "RICKER 1ST PRIMITIVE") ;
		break;
	case RICKER:
		print_info(MASTER, " Source function:", "RICKER") ;
		break;
	case RICKER_D:
		print_info(MASTER, " Source function:", "RICKER 1ST DERIVATIVE") ;
		break;
	case MONO_FREQ:
		print_info(MASTER, " Source function:", "MONO FREQUENCY") ;
		break;
	case SRC_FILE:
		print_info(MASTER, " Source function:", "READ FROM FILE") ;
		break;
	default:
		print_warning(" NO SOURCE FUNCTION HAS BEEN DEFINED") ;
	}

	if (src_func != NO_SRC_FUNC)
	{

		if (src_func == SRC_FILE)
		{
			// src file
			print_info(MASTER, " Input src file:", src_file.c_str()) ;

			// src dt
			if (src_dt <= 0.)
			{
				print_error(" Source dt should be > 0") ;
				return(RTN_CODE_KO) ;
			}
			else
			{
				print_info(MASTER, " Input src dt (s):", src_dt) ;
			}
		}
		else
		{
			// src freq
			if (src_freq <= 0.)
			{
				print_error(" Source frequency should be > 0") ;
				return(RTN_CODE_KO) ;
			}
			else
			{
				print_info(MASTER, " Source freq. (Hz):", src_freq) ;
			}
		}

		// src type
		switch(src_type)
		{
		case SRC_POINT:
			print_info(MASTER, " Source type: \t", "POINT SOURCE") ;
			break ;
		case SRC_GAUSSIAN:
			print_info(MASTER, " Source type: \t", "GAUSSIAN SOURCE") ;
			print_info(MASTER, " Source sigma (m):", src_sigma) ;
			break ;
		case SRC_SINC:
			print_info(MASTER, " Source type: \t", "SINC SOURCE") ;
			print_info(MASTER, " Sinc radius (m):", src_sigma) ;
			print_info(MASTER, " Sinc nw:\t", 999) ;
			break ;
		default:
			print_error(" Source type unknown") ;
			return(RTN_CODE_KO) ;
		}

		// src subtype
		switch(src_stype)
		{
		case EXPLOSIVE:
			print_info(MASTER, " Source subtype:", "EXPLOSIVE") ;
			break ;
		case FORCE_Z:
			print_info(MASTER, " Source subtype:", "FORCE ALONG Z AXIS") ;
			break ;
		default:
			print_error(" Source subtype unknown") ;
			return(RTN_CODE_KO) ;
		}
	}

	print_debug(ALL, FULL_DEBUG, "OUT Modelling::info");
	return(RTN_CODE_OK) ;
}

} // namespace django
