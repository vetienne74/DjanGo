//-------------------------------------------------------------------------------------------------------
//
// MODELLING WITH FDM
//
// PARENT CLASS: Scheme
//   DERIVED CLASS: FDM
//
//-------------------------------------------------------------------------------------------------------

#include "FDM.h"

#include <cassert>
#include <fstream>
#include <iostream>

#include "output_report.h"
#include "type_def.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------

FDM::FDM(void) : Scheme()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM::Scheme");

	method      = SCHEME_FDM ;
	type        = SCHEME_STAGGERED ;

	lstencil    = 0 ;
	npoint_lay  = 0 ;
	npoint_med  = 0 ;

	space_order = 0 ;
	time_order  = 0 ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM::Scheme");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM::initialize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM::initialize(void)");

	// display grid points stat
	print_info(MASTER, "") ;
	print_info(MASTER, " # medium grid points", npoint_med) ;
	print_info(MASTER, " # layer grid points ", npoint_lay) ;
	print_info(MASTER, " Total grid points   ", npoint) ;
	print_info(MASTER, " Layer volume (%)", (Myfloat)((Myfloat)npoint_lay/npoint*100.)) ;
	print_line2() ;

	// call parent initialization
	Rtn_code rtn_code = Scheme::initialize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM::initialize(void)");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM::finalize");

	// call parent finialize
	Rtn_code rtn_code = Scheme::finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM::finalize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM::info(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM::info");

	Scheme::info() ;

	print_info(MASTER, "") ;

	// space_order
	if (space_order <= 0)
	{
		print_error(" FDM::info Invalid spatial order", space_order) ;
		return(RTN_CODE_KO) ;
	}
	else
	{
		print_info(MASTER, " Spatial order\t", space_order) ;
	}

	// time order
	if (time_order < 0)
	{
		print_error(" FDM::info Invalid time order", time_order) ;
		return(RTN_CODE_KO) ;
	}
	else
	{
		print_info(MASTER, " Time order\t", time_order) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM::info");
	return(RTN_CODE_OK) ;
}

} // namespace django
