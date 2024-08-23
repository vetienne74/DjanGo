//-------------------------------------------------------------------------------------------------------
//
// CLASS PROGRAM_FWI
//
//-------------------------------------------------------------------------------------------------------

#include "program_fwi.h"

#include "output_report.h"
#include "type_def.h"

namespace django {

//-------------------------------------------------------------------------------------------------------

Program_Fwi::Program_Fwi(Prog_type type_in) : Program(type_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program_Fwi::Program_Fwi");

	//type = type_in ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Program_Fwi::Program_Fwi");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Program_Fwi::initialize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program_Fwi::initialize");
	// TO DO
	print_debug(ALL, LIGHT_DEBUG, "OUT Program_Fwi::initialize");
	return(RTN_CODE_OK) ;
}


//-------------------------------------------------------------------------------------------------------

Rtn_code Program_Fwi::run(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program_Fwi::run");

	print_info(MASTER, "") ;
	print_line2() ;
	print_info(MASTER, " RUN PROGRAM FWI") ;
	print_info(MASTER, "") ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Program_Fwi::run");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Program_Fwi::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program_Fwi::finalize");
	// TO DO
	print_debug(ALL, LIGHT_DEBUG, "OUT Program_Fwi::finalize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Program_Fwi::info(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program_Fwi::info");

	// general prog info
	Program::info() ;

	print_info(MASTER, "") ;
	print_line2() ;
	print_info(MASTER, " SPECIFIC PROGRAM FWI PARAMETERS") ;
	print_info(MASTER, "") ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Program_Fwi::info");
	return(RTN_CODE_OK) ;
}

} // namespace django
