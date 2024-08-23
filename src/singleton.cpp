//-------------------------------------------------------------------------------------------------------
//
// READ BASIC CONFIGURATION PARAMETERS
// AND INSTANTIATE APPROPRIATE OBJECTS
// THROUGH THE OBJECT FACTORY 
//
//-------------------------------------------------------------------------------------------------------

#include "singleton.h"

#include <cassert>
#include <cstdio>
#include <fstream>
#include <string>

#include "omp.h"

#include "FEM.h"
#include "global.h"
#include "output_report.h"
#include "parse_xml.h"
#include "program.h"

using namespace std;

namespace django {

// cache blocking size
static const Myint DEFAULT_CBX = 25 ;
static const Myint DEFAULT_CBY = 15 ;
static const Myint DEFAULT_CBZ = 10000 ;

// io buffer size
static const Myint DEFAULT_IO_BUFFER = 4096 ;

//-------------------------------------------------------------------------------------------------------

// Global static pointer used to ensure a single instance of the class.
Singleton* Singleton::m_pInstance = NULL;  

//-------------------------------------------------------------------------------------------------------

// This function is called to create an instance of the class. 
// Calling the constructor publicly is not allowed. The constructor 
// is private and is only called by this Instance function.

Singleton* Singleton::Instance()
{
	if (!m_pInstance)   // Only allow one instance of class to be generated.
		m_pInstance = new Singleton;

	return m_pInstance;
}

//-------------------------------------------------------------------------------------------------------

// constructor

Singleton::Singleton(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Singleton::Singleton");
	pProgram     = NULL ;
	nVariable    = 0 ;

	for (Myint ivar = 0; ivar < MAX_VARIABLE; ivar++)
	{
		pVariable[ivar] = NULL ;
	}

	dryrun       = false ;
	iobufsize    = DEFAULT_IO_BUFFER ;
	cbx          = DEFAULT_CBX ;
	cby          = DEFAULT_CBY ;
	cbz          = DEFAULT_CBZ ;
	nsubx        = 1 ;
	nsuby        = 1 ;
	nsubz        = 1 ;
	xml_version  = 0 ;
	eq_order     = NO_EQ_ORDER ;
	eq_type      = NO_EQ_TYPE ;
	ratio_cfl    = 1.0 ;
	dt           = 0.0 ;
	time_order   = 0 ;
	space_order  = 0 ;
	fmax0        = 0.0 ;
	type         = NO_SCHEME_TYPE ;
	nelem_x_med  = 0 ;
	nelem_z_med  = 0 ;
	flux_type    = NO_FLUX_TYPE ;
	adaptType    = NO_ADAPT ;
	dim          = NO_DIM ;
	scheme_type  = NO_SCHEME_TYPE ;
	node_type    = NO_NODE_TYPE ;
	front_type   = FRONT_STATIC ;
	prop         = PROP_CONST ;
	pmin         = 0 ;
	pmax         = MAX_POLY_ORDER ;
	node_distrib = NO_NODE_DISTRIB ;
	pModel       = NULL ;
	method       = NO_SCHEME_METHOD ;
	node_integ   = NODE_INTEG_GLL ;

	nBoundary    = 0 ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Singleton::Singleton");
}

//-------------------------------------------------------------------------------------------------------

Variable* Singleton::get_variable(Var_type type)
{
	print_debug(ALL, FULL_DEBUG, "IN Singleton::get_variable");

	Variable* pV = NULL ;

	print_debug(ALL, FULL_DEBUG, "nb variable in memory", nVariable);
	print_debug(ALL, FULL_DEBUG, "search for type", (Myint) type);

	for (Myint iVariable = 0; iVariable < MAX_VARIABLE; iVariable++)
	{
		print_debug(ALL, FULL_DEBUG, "iVariable ", iVariable);
		if (pVariable[iVariable] == NULL) continue ;
		print_debug(ALL, FULL_DEBUG, "type", (Myint) pVariable[iVariable]->get_type());
		if (pVariable[iVariable]->get_type() == type)
		{
			print_debug(ALL, FULL_DEBUG, "found required type");
			pV = pVariable[iVariable] ;
			break ;
		}
	}

	print_debug(ALL, FULL_DEBUG, "OUT Singleton::get_variable");
	return pV ;
}

//-------------------------------------------------------------------------------------------------------

Variable* Singleton::get_variable(Myint varId)
{
	print_debug(ALL, FULL_DEBUG, "IN Singleton::get_variable");

	Variable* pV = NULL ;

	if (varId < MAX_VARIABLE) pV = pVariable[varId] ;

	print_debug(ALL, FULL_DEBUG, "OUT Singleton::get_variable");
	return pV ;
}


//-------------------------------------------------------------------------------------------------------

Variable* Singleton::register_variable(Var_type type_in, string string_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Singleton::register_variable");

	// loop on variable array and find first empty slot
	Myint var_id = -1 ;
	for (Myint ivar = 0; ivar < MAX_VARIABLE; ivar++)
	{
		if (pVariable[ivar] == NULL)
		{
			var_id = ivar ;
			pVariable[var_id] = new Variable(type_in, string_in) ;
			pVariable[var_id]->set_id(var_id) ;
			nVariable++;
			break ;
		}
	}

	if (var_id == -1)
	{
		print_error(" Number of max variables has been reached");
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Singleton::register_variable");
	return pVariable[var_id] ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Singleton::swap_variable(Var_type v1, Var_type v2)
{
	print_debug(ALL, FULL_DEBUG, "IN Singleton::swap_variable");

	Variable* var1 = get_variable(v1) ;
	Variable* var2 = get_variable(v2) ;
	//var1->set_type(v2) ;
	//var2->set_type(v1) ;
	Grid* pGrid1 = var1->get_grid() ;
	Grid* pGrid2 = var2->get_grid() ;
	var1->set_grid(pGrid2) ;
	var2->set_grid(pGrid1) ;

	print_debug(ALL, FULL_DEBUG, "OUT Singleton::swap_variable");
	return(RTN_CODE_OK);
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Singleton::delete_all_variable(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Singleton::delete_all_variable");
	for (Myint ivar=0; ivar < MAX_VARIABLE; ivar++)
	{
		delete(pVariable[ivar]) ;
		nVariable-- ;
	}

	if (nVariable != 0)
	{
		print_error(" Error in Singleton::delete_all_variable") ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Singleton::delete_all_variable");
	return(RTN_CODE_OK);
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Singleton::delete_variable(Var_type type)
{
	print_debug(ALL, MID_DEBUG, "IN Singleton::delete_variable");

	for (Myint iVariable = 0; iVariable < MAX_VARIABLE; iVariable++)
	{
		print_debug(ALL, MID_DEBUG, "iVariable ", iVariable);
		if (pVariable[iVariable] != NULL)
		{
			if (pVariable[iVariable]->get_type() == type)
			{
				print_debug(ALL, MID_DEBUG, "found required type");
				delete(pVariable[iVariable]) ;
				pVariable[iVariable] = NULL ;
			}
		}
	}

	print_debug(ALL, MID_DEBUG, "OUT Singleton::delete_variable");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Singleton::delete_variable(Myint var_id)
{
	print_debug(ALL, MID_DEBUG, "IN Singleton::delete_variable");

	if ((var_id >= 0) && (var_id < MAX_VARIABLE))
	{
		print_debug(ALL, MID_DEBUG, "delete var_id", var_id);
		delete(pVariable[var_id]) ;
		pVariable[var_id] = NULL ;
		nVariable-- ;
	}

	print_debug(ALL, MID_DEBUG, "OUT Singleton::delete_variable");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Singleton::print_all_variable(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Singleton::print_all_variable");
	print_info(MASTER, " Variables in memory:") ;
	for (Myint ivar=0; ivar < MAX_VARIABLE; ivar++)
	{
		if (pVariable[ivar] != NULL) pVariable[ivar]->info() ;
	}
	print_debug(ALL, LIGHT_DEBUG, "OUT Singleton::print_all_variable");
	return(RTN_CODE_OK);
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Singleton::initialize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Singleton::initialize");

	// check config file name
	if (xml_config_file == "")
	{
		print_error(" XML config file name has not been specified (option -xml)") ;
		return(RTN_CODE_KO) ;
	}

	//-----------------------------------------
	// parse XML configuration file with EXPAT
	//-----------------------------------------
	char buf[BUFSIZ];
	XML_Parser parser = XML_ParserCreate(NULL);
	int done;
	int depth = 0;

	// init parser
	XML_SetUserData(parser, &depth);
	XML_SetElementHandler(parser, django_startElement, django_endElement);

	// open xml config file
	ifstream config_file(xml_config_file.c_str()) ;
	assert(config_file.is_open());

	// parse file byte after byte up to the end
	do {
		int len = 1 ;
		config_file.read (buf,len);
		done = config_file.eof() ;
		if (XML_Parse(parser, buf, len, done) == XML_STATUS_ERROR) {
			fprintf(stderr,
					"%s at line %" XML_FMT_INT_MOD "u\n",
					XML_ErrorString(XML_GetErrorCode(parser)),
					XML_GetCurrentLineNumber(parser));
			print_error(" Error while parsing XML configuration file") ;
			return(RTN_CODE_KO) ;
		}
	} while ( !done ) ;

	// free parser
	XML_ParserFree(parser);
	config_file.close() ;

	// print config info
	Rtn_code rtn = this->info() ;
	if (rtn != RTN_CODE_OK) return(rtn) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Singleton::initialize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code Singleton::info(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Singleton::info");

	print_info(MASTER, "") ;
	print_line2() ;
	print_info(MASTER, " SYSTEM CONFIGURATION") ;
	print_info(MASTER, "") ;

	// number of MPI process
	print_info(MASTER, " No. of MPI process", nproc_world) ;
	print_info(MASTER, " No. subdomains X", Singleton().Instance()->nsubx) ;
	print_info(MASTER, " No. subdomains Y", Singleton().Instance()->nsuby) ;
	print_info(MASTER, " No. subdomains Z", Singleton().Instance()->nsubz) ;

	if (nproc_world != (Singleton().Instance()->nsubx * Singleton().Instance()->nsuby
			* Singleton().Instance()->nsubz))
	{
		print_error("No. MPI inconsistent with subdomain decomposition") ;
		return(RTN_CODE_KO) ;
	}

	// number of OpenMP threads / MPI process
#ifdef _OPENMP
	print_info(MASTER, " OpenMP threads / MPI", omp_get_max_threads()) ;
#else
	print_info(MASTER, " OpenMP disabled") ;
#endif  

	if (xml_version != CURRENT_VERSION)
	{
		string err_msg = " Invalid XML file (not django / version " + to_string(CURRENT_VERSION) + ")\n" ;
		print_error(&err_msg) ;
		return(RTN_CODE_KO) ;
	}

	//-------------------------------------------------------------------------------------------------------
	// display configuration paramaters and perform some checkings
	//-------------------------------------------------------------------------------------------------------

	// program
	if (pProgram == NULL)
	{
		print_error(" Program non initialized") ;
		return(RTN_CODE_KO) ;
	}
	Rtn_code rtn = pProgram->info() ;
	if (rtn != RTN_CODE_OK) return(rtn) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Singleton::info");
	return(RTN_CODE_OK) ;

}

//-------------------------------------------------------------------------------------------------------

Rtn_code Singleton::delete_all_boundary(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Singleton::delete_all_boundary");
	for (Myint ibound=1; ibound <= nBoundary; ibound++)
	{
		delete(pBoundary[ibound-1]) ;
	}
	nBoundary = 0 ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Singleton::delete_all_boundary");
	return(RTN_CODE_OK);
}

//-------------------------------------------------------------------------------------------------------
Boundary* Singleton::get_boundary(Edge_type edge_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Singleton::get_boundary");

	Boundary* boundary_tmp = NULL ;
	for (Myint ib=1; ib<=nBoundary; ib++)
	{
		if (pBoundary[ib-1]->get_edge() == edge_in) boundary_tmp = pBoundary[ib-1] ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Singleton::get_boundary");
	return boundary_tmp ;
}

//-------------------------------------------------------------------------------------------------------
Myint Singleton::get_boundary_width(Edge_type edge_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Singleton::get_boundary_width");

	Myint width_tmp = 0 ;
	for (Myint ib=1; ib<=nBoundary; ib++)
	{
		if (pBoundary[ib-1]->get_edge() == edge_in) width_tmp = pBoundary[ib-1]->get_width() ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Singleton::get_boundary_width");
	return width_tmp ;
}

//-------------------------------------------------------------------------------------------------------
Myfloat Singleton::get_boundary_coef(Edge_type edge_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Singleton::get_boundary_coef");

	Myfloat coef_tmp = 0 ;
	for (Myint ib=1; ib<=nBoundary; ib++)
	{
		if (pBoundary[ib-1]->get_edge() == edge_in) coef_tmp = pBoundary[ib-1]->get_coef() ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Singleton::get_boundary_coef");
	return coef_tmp ;
}

//-------------------------------------------------------------------------------------------------------
Boundary_type Singleton::get_boundary_type(Edge_type edge_in)
{
	print_debug(ALL, FULL_DEBUG, "IN Singleton::get_boundary_type");

	Boundary_type type_tmp = NO_BOUNDARY ;
	for (Myint ib=1; ib<=nBoundary; ib++)
	{
		if (pBoundary[ib-1]->get_edge() == edge_in) type_tmp = pBoundary[ib-1]->get_type() ;
	}

	print_debug(ALL, FULL_DEBUG, "OUT Singleton::get_boundary_type");
	return type_tmp ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Singleton::create_boundary(Edge_type edge, Boundary_type type, Myint width, Myfloat coef)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Singleton::create_boundary");

	if ((nBoundary < 0) || (nBoundary >= MAX_BOUNDARY))
	{
		print_error(" Invalid nBoundary", nBoundary) ;
		return(RTN_CODE_KO) ;
	}

	pBoundary[nBoundary] = new Boundary(edge, type, width, coef) ;

	nBoundary++ ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Singleton::create_boundary");
	return(RTN_CODE_OK) ;
}

} // namespace django


