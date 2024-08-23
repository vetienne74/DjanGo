//-------------------------------------------------------------------------------------------------------
//
// CLASS PROGRAM_GUITAR
//
//-------------------------------------------------------------------------------------------------------

#include "program_guitar.h"

#include <algorithm>

#include "acquisition.h"
#include "allocate_array.h"
#include "data.h"
#include "data_std.h"
#include "FDM_1D_ac_lossy.h"
#include "output_report.h"
#include "singleton.h"
#include "type_def.h"

namespace django {

//-------------------------------------------------------------------------------------------------------

Program_Guitar::Program_Guitar(Prog_type type_in) : Program(type_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program_Guitar::Program_Guitar");

	//type = type_in ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Program_Guitar::Program_Guitar");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Program_Guitar::initialize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program_Guitar::initialize");

	// initialize model
	pModel = Singleton::Instance()->pProgram->pDomain->pModel ;

	// initialise data
	pData  = new Data_std() ;

	// initialize modelling
	Rtn_code rtn_code = pModelling->initialize(pDomain->pModel) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// get left hand time function
	rtn_code = get_left_hand_time_function() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// set string length in Scheme
	FDM_1D_ac_lossy *pFDM_1D_ac_lossy = dynamic_cast<FDM_1D_ac_lossy*>(pDomain->pScheme) ;
	if (pFDM_1D_ac_lossy != NULL)
	{
		rtn_code = pFDM_1D_ac_lossy->set_string_length_vs_time_step(left_hand_time_function) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
	}
	else
	{
		print_error("IN Program_Guitar::initialize --> scheme is not FDM_1D_ac_lossy");
		return(RTN_CODE_KO) ;
	}


	print_debug(ALL, LIGHT_DEBUG, "OUT Program_Guitar::initialize");
	return(RTN_CODE_OK) ;
}


//-------------------------------------------------------------------------------------------------------

Rtn_code Program_Guitar::run(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program_Guitar::run");

	print_info(MASTER, "") ;
	print_line2() ;
	print_info(MASTER, " RUN PROGRAM GUITAR") ;
	print_info(MASTER, "") ;

	Rtn_code rtn_code = pModelling->solve_all_shot(pData, pModel) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Program_Guitar::run");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Program_Guitar::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program_Guitar::finalize");

	Rtn_code rtn_code = free_left_hand_time_function() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	rtn_code = pModelling->finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	delete(pModel) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Program_Guitar::finalize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Program_Guitar::free_left_hand_time_function(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program_Guitar::free_left_hand_time_function");
	deallocate_array<Myfloat>(left_hand_time_function, pDomain->pScheme->get_nt()+1) ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Program_Guitar::free_left_hand_time_function");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Program_Guitar::set_left_hand(string left_hand_file_in, Myfloat left_hand_file_dt_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program_Guitar::set_left_hand");
	left_hand_file = left_hand_file_in ;
	left_hand_file_dt = left_hand_file_dt_in ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Program_Guitar::set_left_hand");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Program_Guitar::info(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program_Guitar::info");

	// general prog info
	Rtn_code rtn_code = Program::info() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_info(MASTER, "") ;
	print_line2() ;
	print_info(MASTER, " SPECIFIC PROGRAM GUITAR PARAMETERS") ;
	print_info(MASTER, "") ;

	print_info(MASTER, " Left hand file\t", left_hand_file) ;
	print_info(MASTER, " Left hand file dt", left_hand_file_dt) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Program_Guitar::info");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Program_Guitar::get_left_hand_time_function(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program_Guitar::get_left_hand_time_function");

	Myint nt = pDomain->pScheme->get_nt() ;
	Myfloat dt = pDomain->pScheme->get_dt() ;

	left_hand_time_function = allocate_array<Myfloat>(nt+1) ;
	if (left_hand_time_function == NULL)
	{
		print_error(" left_hand_time_function is NULL") ;
		return(RTN_CODE_KO) ;
	}

	print_info(MASTER, "") ;
	print_info(MASTER, " left hand initialization") ;

	// left hand  function read from disk

	// estimate nt required
	Myint nt_required = (nt * dt) / left_hand_file_dt ;
	print_debug(ALL, MID_DEBUG, "nt required", nt_required) ;

	// read source file
	ifstream in_file ;
	in_file.open(left_hand_file.c_str(), ios::binary) ;
	assert(in_file.is_open());

	// get length of file
	in_file.seekg (0, in_file.end);
	Myint nt_in_file = in_file.tellg() / 4 ;
	in_file.seekg (0, in_file.beg);

	// determine # time step to read
	Myint nt_read = min(nt_in_file, nt_required) ;

	// allocate array
	print_info(MASTER, " Nt in left hand file:", nt_read) ;
	Myfloat32* val_array = new Myfloat32[nt_read] ;
	in_file.read((char*) val_array, nt_read * sizeof(Myfloat32)) ;
	in_file.close() ;

	// interpolate
	for (Myint it = 0; it <= nt; it ++)
	{
		Myfloat t1 = it * dt ;
		Myint it_left_hand = t1 / left_hand_file_dt ;

		if (it_left_hand < nt_read-1)
		{
			// nearest point
			//left_hand_time_function[it] = val_array[it_left_hand] ;

			// linear interpolation
			Myfloat t2 = it_left_hand * left_hand_file_dt ;
			Myfloat dd = t1 - t2 ;
			Myfloat a2 = dd / left_hand_file_dt ;
			Myfloat a1 = 1. - a2 ;
			left_hand_time_function[it] = a1*val_array[it_left_hand] + a2*val_array[it_left_hand+1];
		}
		else if (it_left_hand == nt_read-1)
		{
			left_hand_time_function[it] = val_array[it_left_hand] ;
		}
		else
		{
			left_hand_time_function[it] = 0.0 ;
		}
	}
	delete(val_array) ;

	Myfloat max_amp ;
	for (Myint it = 0; it <= nt; it ++)
	{
		// retrieve time when source is max
		if (max_amp < abs(left_hand_time_function[it]))
		{
			max_amp = abs(left_hand_time_function[it]) ;
		}
	}

	print_info(MASTER, " Max. left hand:", max_amp) ;

	// output src wavelet
	//this->output_src_wavelet() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Program_Guitar::get_left_hand_time_function");
	return(RTN_CODE_OK) ;
}

} // namespace django
