//-------------------------------------------------------------------------------------------------------
//
// CLASS PROGRAM_RTM
//
//-------------------------------------------------------------------------------------------------------

#include "program_rtm.h"

#include "output_report.h"
#include "type_def.h"

namespace django {

//-------------------------------------------------------------------------------------------------------

Program_Rtm::Program_Rtm(Prog_type type_in) : Program(type_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program_Rtm::Program_Rtm");

	//type = type_in ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Program_Rtm::Program_Rtm");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Program_Rtm::initialize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program_Rtm::initialize");
	// TO DO
	print_debug(ALL, LIGHT_DEBUG, "OUT Program_Rtm::initialize");
	return(RTN_CODE_OK) ;
}


//-------------------------------------------------------------------------------------------------------

Rtn_code Program_Rtm::run(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program_Rtm::run");

	print_info(MASTER, "") ;
	print_line2() ;
	print_info(MASTER, " RUN PROGRAM RTM") ;
	print_info(MASTER, "") ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Program_Rtm::run");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Program_Rtm::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program_Rtm::finalize");
	// TO DO
	print_debug(ALL, LIGHT_DEBUG, "OUT Program_Rtm::finalize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Program_Rtm::info(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program_Rtm::info");

	// general prog info
	Program::info() ;

	print_info(MASTER, "") ;
	print_line2() ;
	print_info(MASTER, " SPECIFIC PROGRAM RTM PARAMETERS") ;
	print_info(MASTER, "") ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Program_Rtm::info");
	return(RTN_CODE_OK) ;
}

} // namespace django
