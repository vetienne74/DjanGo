//-------------------------------------------------------------------------------------------------------
//
// CLASS PROGRAM
//
//-------------------------------------------------------------------------------------------------------

#include "program.h"

#include "modelling.h"
#include "output_report.h"
#include "type_def.h"

namespace django {

//-------------------------------------------------------------------------------------------------------

Program::Program(Prog_type type_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program::Program");

	type = type_in ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Program::Program");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Program::info(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program::info");

	print_info(MASTER, "") ;
	print_line2() ;
	print_info(MASTER, " GENERAL PROGRAM PARAMETERS") ;
	print_info(MASTER, "") ;

	// type
	switch(type)
	{
	case MODELLING_PROG:
		print_info(MASTER, " Program: \t\tMODELLING") ;
		break ;
	case FWI_PROG:
		print_info(MASTER, " Program: \t\tFULL WAVEFORM INVERSION") ;
		break ;
	case RTM_PROG:
		print_info(MASTER, " Program: \t\tREVERSE TIME MIGRATION") ;
		break ;
	case GUITAR_PROG:
		print_info(MASTER, " Program: \t\tGUITAR") ;
		break ;
	default:
		print_error(" Invalid program in config file", type) ;
		return(RTN_CODE_KO) ;
	}

	// domain
	if (pDomain == NULL)
	{
		print_error(" Domain non initialized") ;
		return(RTN_CODE_KO) ;
	}
	Rtn_code rtn = pDomain->info() ;
	if (rtn != RTN_CODE_OK) return(rtn) ;

	// modelling
	if (pModelling == NULL)
	{
		print_error(" Modelling non initialized") ;
		return(RTN_CODE_KO) ;
	}
	rtn = pModelling->info() ;
	if (rtn != RTN_CODE_OK) return(rtn) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Program::info");
	return(RTN_CODE_OK) ;
}
} // namespace django
