//-------------------------------------------------------------------------------------------------------
//
// SEISMIC DATA MANAGEMENT
//
// PARENT CLASS: Data
//   DERIVED CLASS: Data_std
//
//-------------------------------------------------------------------------------------------------------

#include "data_std.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <string>

#include "allocate_array.h"
#include "output_report.h"
#include "singleton.h"
#include "type_def.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------

Rtn_code Data_std::read(Acquisition* pAcquisition)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Data_std::read");

	// open file
	ifstream data_file ;
	data_file.open(pr_data_filename.c_str(), ios::binary) ;
	assert(data_file.is_open());

	// input data in freq.
	if (domain_in == FREQ)

	{
		Myint64 pos_in_file = Myint64(pAcquisition->current_src - 1) * pAcquisition->nrec * 8 ;
		data_file.seekg (pos_in_file) ;

		print_info(MASTER, " Read freq. data /src:", pAcquisition->current_src) ;
		print_info(MASTER, " Position in file: ", pos_in_file) ;

		// output data in freq. domain
		if (domain_out == FREQ)
		{
			for (Myint irec=0; irec<pAcquisition->nrec; irec++)
			{
				data_file.read((char*) &(pr_freq_obs_rec[irec]), sizeof(complex<float>)) ;
			}
		}

		// output data in time domain
		else
		{
			print_error(" Input data freq / output time not allowed !") ;
			return(RTN_CODE_KO) ;
		}

	}

	// input data in time domain
	else if (domain_in == TIME)

	{
		Myint64 pos_in_file = Myint64(pAcquisition->current_src - 1) * pAcquisition->nrec * 4 * nt_in ;
		data_file.seekg (pos_in_file) ;

		print_info(MASTER, " Read time data /src:", pAcquisition->current_src) ;
		print_info(MASTER, " Position in file: ", pos_in_file) ;

		// loop over the receivers
		for (Myint irec=0; irec<pAcquisition->nrec; irec++)
		{

			// output data in time domain
			if (domain_out == TIME)
			{
				data_file.read((char*) &(pr_time_obs_rec[irec,0]), sizeof(Myfloat32) * nt_in) ;
			}

			// output data in freq domain
			else
			{

				// read data
				Myfloat pr_time_obs_rec_tmp[nt_in] ;
				data_file.read((char*) &pr_time_obs_rec_tmp, sizeof(Myfloat) * nt_in) ;

				// loop over frequencies
				for (Myint ifreq = 0; ifreq < pFreq_group_out->nb_freq; ifreq++)
				{
					// convert data with DFT
					pr_freq_obs_rec[irec][ifreq] = ZERO_CMPLX ;
					for (Myint it = 0; it < nt_in; it++)
					{
						Myfloat coef_tmp = -2. * PI *  pFreq_group_out->pFreq_list[ifreq] * it * dt_in ;
						//Myfloat coef_tmp = -2. * PI *  pFreq_group_out->pFreq_list[ifreq] * it * float(nt_in-1)/float(nt_in) * dt_in ;

						Mycomplex coef (cos(coef_tmp)*dt_in, sin(coef_tmp)*dt_in) ;
						pr_freq_obs_rec[irec][ifreq] += pr_time_obs_rec_tmp[it] * coef ;
					}
				}
			}
		}
	}

	// close file
	data_file.close() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Data_std::read");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Data_std::reset(Acquisition* pAcquisition, Prog_type prog)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Data_std::reset");

	switch(prog)
	{

	case MODELLING_PROG:
		//if (Singleton::Instance()->pModelling->get_output_frequency_flag() == YES)
	{	  
		for (Myint irec=0; irec<pAcquisition->nrec; irec++)
		{
			for (Myint ifreq = 0; ifreq < nf_out; ifreq++)
			{
				pr_freq_sol_rec[irec][ifreq] = ZERO_CMPLX  ;
			}
		}
	}
	break ;

	case FWI_PROG:

		if (domain_out == FREQ)
		{
			for (Myint irec=0; irec<pAcquisition->nrec; irec++)
			{
				for (Myint ifreq = 0; ifreq < nf_out; ifreq++)
				{
					pr_freq_sol_rec[irec][ifreq] = ZERO_CMPLX  ;
					pr_freq_obs_rec[irec][ifreq] = ZERO_CMPLX  ;
					pr_freq_adj_rec[irec][ifreq] = ZERO_CMPLX  ;
				}
			}
		}

		else if (domain_out == TIME)
		{
			// TO BE DONE
		}

		else
		{
			print_error(" Data_std::finalize output data format not supported for this program") ;
			return(RTN_CODE_KO) ;
		}

		break ;
	default:
		print_error(" Data_std::finalize not supported for this program") ;
		return(RTN_CODE_KO) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Data_std::reset");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Data_std::initialize(Acquisition* pAcquisition, Prog_type prog, Domain_type domain_out2, Myint nfreq)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Data_std::initialize");

	if (pAcquisition == NULL)
	{
		print_error(" In Data_std::initialize, pAcquisition is NULL") ;
		return(RTN_CODE_KO) ;
	}

	// for output data
	domain_out = domain_out2 ;
	nf_out     = nfreq ;

	ifstream data_file(DATA_CONFIG_IN_FILE) ;

	switch(prog)
	{

	case MODELLING_PROG:
		// //if (Singleton::Instance()->pModelling->get_output_frequency_flag() == YES)
		// if (nfreq > 0)
		// 	{
		// 	  pr_freq_sol_rec = allocate_array<Mycomplex>(pAcquisition->nrec, nfreq) ;
		// 	}
		break ;

	case FWI_PROG:
	{
		// read data config file
		Rtn_code rtn_code = read_config() ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		if (domain_out == FREQ)
		{
			if (nfreq <= 0)
			{
				print_error(" Data_std::initialize nb. freq <= 0") ;
				return(RTN_CODE_KO) ;
			}
			pr_freq_sol_rec = allocate_array<Mycomplex>(pAcquisition->nrec, nfreq) ;
			pr_freq_obs_rec = allocate_array<Mycomplex>(pAcquisition->nrec, nfreq) ;
			pr_freq_adj_rec = allocate_array<Mycomplex>(pAcquisition->nrec, nfreq) ;
		}
		else if (domain_out == TIME)
		{
			if (nt_out <= 0)
			{
				print_error(" Data_std::initialize nt <= 0") ;
				return(RTN_CODE_KO) ;
			}
			pr_time_sol_rec = allocate_array<Myfloat>(pAcquisition->nrec, nt_out) ;
			pr_time_obs_rec = allocate_array<Myfloat>(pAcquisition->nrec, nt_out) ;
			pr_time_adj_rec = allocate_array<Myfloat>(pAcquisition->nrec, nt_out) ;
		}
		else
		{
			print_error(" Data_std::initialize output data format not supported for this program") ;
			return(RTN_CODE_KO) ;
		}

		break ;
	}
	default:
		print_error(" Data_std::initialize not supported for this program") ;
		return(RTN_CODE_KO) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Data_std::initialize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Data_std::finalize(Acquisition* pAcquisition, Prog_type prog)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Data_std::finalize");

	switch(prog)
	{

	case MODELLING_PROG:
		// //if (Singleton::Instance()->pModelling->get_output_frequency_flag() == YES)
		// if (nfreq > 0)
		// 	{
		// 	  deallocate_array<Mycomplex>(pr_freq_sol_rec, pAcquisition->nrec, nfreq) ;
		// 	}
		break ;

	case FWI_PROG:

		if (domain_out == FREQ)
		{
			deallocate_array<Mycomplex>(pr_freq_sol_rec, pAcquisition->nrec, pFreq_group_out->nb_freq) ;
			deallocate_array<Mycomplex>(pr_freq_obs_rec, pAcquisition->nrec, pFreq_group_out->nb_freq) ;
			deallocate_array<Mycomplex>(pr_freq_adj_rec, pAcquisition->nrec, pFreq_group_out->nb_freq) ;
		}
		else if (domain_out == TIME)
		{
			deallocate_array<Myfloat>(pr_time_sol_rec, pAcquisition->nrec, nt_out) ;
			deallocate_array<Myfloat>(pr_time_obs_rec, pAcquisition->nrec, nt_out) ;
			deallocate_array<Myfloat>(pr_time_adj_rec, pAcquisition->nrec, nt_out) ;
		}
		else
		{
			print_error(" Data_std::finalize output data format not supported for this program") ;
			return(RTN_CODE_KO) ;
		}

		break ;
	default:
		print_error(" Data_std::finalize not supported for this program") ;
		return(RTN_CODE_KO) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Data_std::finalize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Data_std::read_config()
{
	print_debug(ALL, LIGHT_DEBUG, "IN Data_std::read_config");

	ifstream config_file(DATA_CONFIG_IN_FILE) ;
	assert(config_file.is_open());

	print_info(MASTER, "") ;
	print_line2() ;
	print_info(MASTER, " DATA PARAMETERS") ;
	print_info(MASTER, "") ;

	// get data domain
	Myint int_tmp ;
	config_file >> int_tmp ;

	switch(int_tmp)
	{
	case TIME:
		domain_in = TIME ;
		print_info(MASTER, " Data domain: \t\tTIME DOMAIN") ;
		break ;
	case FREQ:
		domain_in = FREQ ;
		print_info(MASTER, " Data domain: \t\tFREQUENCY DOMAIN") ;
		break ;
	default:
		print_error(" Invalid data domain") ;
		return(RTN_CODE_KO) ;
	}

	// get nt and dt
	if (domain_in == TIME)
	{
		config_file >> nt_in ;
		print_info(MASTER, " Data nt:\t" , nt_in) ;
		config_file >> dt_in ;
		print_info(MASTER, " Data dt:\t" , dt_in) ;
	}

	// get data filename
	config_file >> pr_data_filename ;
	print_info(MASTER, " Data file name:",pr_data_filename.c_str()) ;

	print_line2() ;
	print_info(MASTER, "") ;

	config_file.close() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Data_std::read_config");
	return(RTN_CODE_OK) ;
}

} // namespce django
