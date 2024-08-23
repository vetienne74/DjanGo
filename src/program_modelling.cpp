//-------------------------------------------------------------------------------------------------------
//
// CLASS PROGRAM_MODELLING
//
//-------------------------------------------------------------------------------------------------------

#include "program_modelling.h"

#include "acquisition.h"
#include "data.h"
#include "data_std.h"
#include "output_report.h"
#include "singleton.h"
#include "type_def.h"

namespace django {

//-------------------------------------------------------------------------------------------------------

Program_Modelling::Program_Modelling(Prog_type type_in) : Program(type_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program_Modelling::Program_Modelling");

	//type = type_in ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Program_Modelling::Program_Modelling");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Program_Modelling::initialize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program_Modelling::initialize");

	pModel = Singleton::Instance()->pProgram->pDomain->pModel ;

	pData  = new Data_std() ;

	Rtn_code rtn_code = pModelling->initialize(pDomain->pModel) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Program_Modelling::initialize");
	return(RTN_CODE_OK) ;
}


//-------------------------------------------------------------------------------------------------------

Rtn_code Program_Modelling::run(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program_Modelling::run");

	print_info(MASTER, "") ;
	print_line2() ;
	print_info(MASTER, " RUN PROGRAM MODELLING") ;
	print_info(MASTER, "") ;

	Rtn_code rtn_code = pModelling->solve_all_shot(pData, pModel) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Program_Modelling::run");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Program_Modelling::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program_Modelling::finalize");

	Rtn_code rtn_code = pModelling->finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	delete(pModel) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Program_Modelling::finalize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Program_Modelling::info(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Program_Modelling::info");

	// general prog info
	Rtn_code rtn_code = Program::info() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_info(MASTER, "") ;
	print_line2() ;
	print_info(MASTER, " SPECIFIC PROGRAM MODELLING PARAMETERS") ;
	print_info(MASTER, "") ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Program_Modelling::info");
	return(RTN_CODE_OK) ;
}

} // namespace django
