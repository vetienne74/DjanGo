//-------------------------------------------------------------------------------------------------------
//
// PARENT CLASS FOR ALL BOUNDARY
//
//-------------------------------------------------------------------------------------------------------

#include "boundary.h"

#include "output_report.h"
#include "type_def.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------

Boundary::Boundary(Edge_type edge_in, Boundary_type type_in, Myint width_in, Myfloat coef_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Boundary::Boundary") ;

	edge  = edge_in ;
	type  = type_in ;
	if ((type == FREESURF) || (type == RIGID))
	{
		width = 0 ;
		coef  = 0.0 ;
	}
	else
	{
		width = width_in ;
		coef  = coef_in ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Boundary::Boundary") ;
}

//-------------------------------------------------------------------------------------------------------

void Boundary::info()
{
	print_debug(ALL, LIGHT_DEBUG, "IN Boundary::info") ;

	switch(edge)
	{
	case ZBEG:
		print_info(MASTER, " Boundary edge: \tZBEG") ;
		break ;
	case ZEND:
		print_info(MASTER, " Boundary edge: \tZEND") ;
		break ;
	case XBEG:
		print_info(MASTER, " Boundary edge: \tXBEG") ;
		break ;
	case XEND:
		print_info(MASTER, " Boundary edge: \tXEND") ;
		break ;
	case YBEG:
		print_info(MASTER, " Boundary edge: \tYBEG") ;
		break ;
	case YEND:
		print_info(MASTER, " Boundary edge: \tYEND") ;
		break ;
	default:
		print_error(" Invalid boundary edge", edge) ;
	}

	switch(type)
	{
	case PML:
		print_info(MASTER, " Boundary type: \tCPML") ;
		break ;
	case RANDOM:
		print_info(MASTER, " Boundary type: \tRANDOM LAYER") ;
		break ;
	case SPG:
		print_info(MASTER, " Boundary type: \tSPONGE (CERJAN)") ;
		break ;
	case SPG2:
		print_info(MASTER, " Boundary type: \tSPONGE (ISRAELI)") ;
		break ;
	case FREESURF:
		print_info(MASTER, " Boundary type: \tFREE SURFACE") ;
		break ;
	case RIGID:
		print_info(MASTER, " Boundary type: \tRIGID SURFACE") ;
		break ;
	default:
		print_error(" Invalid boundary type", type) ;
	}

	if ((type == PML) || (type == RANDOM) || (type == SPG) || (type == SPG2))
	{
		print_info(MASTER, " Boundary width:", width) ;
		print_info(MASTER, " Boundary coef:\t", coef) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Boundary::info") ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Boundary::get_coef(void)
{
	return coef ;
}

//-------------------------------------------------------------------------------------------------------

Boundary_type Boundary::get_type(void)
{
	return type ;
}

//-------------------------------------------------------------------------------------------------------

Myint Boundary::get_width(void)
{
	return width ;
}

//-------------------------------------------------------------------------------------------------------

Edge_type Boundary::get_edge(void)
{
	return edge ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Boundary::reset()
{
	return(RTN_CODE_OK) ;
}

} // namespace django
