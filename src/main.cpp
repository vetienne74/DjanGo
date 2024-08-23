//-------------------------------------------------------------------------------------------------------
//
// GENERIC C++ PLATFORM FOR MODELLING AND INVERSION
//
//-------------------------------------------------------------------------------------------------------

#include <cassert>
#include <ctime>
#include <fstream>
#include <iostream>

#include "mpi.h" 

#include "global.h"
#include "output_report.h"
#include "parse_argument.h"
#include "program.h"
#include "singleton.h"
#include "type_def.h"

using namespace std;

//-------------------------------------------------------------------------------------------------------
// global variables
//-------------------------------------------------------------------------------------------------------

// debug level switch
django::Debug_level django::debug = NO_DEBUG ;

// global number of MPI process
int django::nproc_world ;

// global rank of MPI process
int django::myid_world ;

// memory management
django::Myint64 django::current_mem = 0 ;
django::Myint64 django::max_mem     = 0 ;

//-------------------------------------------------------------------------------------------------------
// start of main program
//-------------------------------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
	// start MPI environment
	MPI_Init(&argc, &argv) ;
	MPI_Comm_size(MPI_COMM_WORLD, &django::nproc_world) ;
	MPI_Comm_rank(MPI_COMM_WORLD, &django::myid_world) ;

	//-------------------------------------------------------------------------------------------------------
	// print header of output report
	//-------------------------------------------------------------------------------------------------------
	django::Rtn_code rtn_code = django::print_header_of_output_report() ;
	if (rtn_code != django::RTN_CODE_OK) {
		MPI_Finalize() ;
		return(rtn_code) ;
	}

	//-------------------------------------------------------------------------------------------------------
	// parse command line arguments
	//-------------------------------------------------------------------------------------------------------
	rtn_code = django::parse_argument(argc, argv) ;
	if (rtn_code != django::RTN_CODE_OK) {
		MPI_Finalize() ;
		return(django::RTN_CODE_OK) ;
	}

	//-------------------------------------------------------------------------------------------------------
	// initialize Singleton
	//-------------------------------------------------------------------------------------------------------
	rtn_code = django::Singleton::Instance()->initialize() ;
	if (rtn_code != django::RTN_CODE_OK)
	{
		MPI_Finalize() ;
		return(django::RTN_CODE_OK) ;
	}

	//-------------------------------------------------------------------------------------------------------
	// initialize program
	//-------------------------------------------------------------------------------------------------------
	django::Program* pProgram = django::Singleton::Instance()->pProgram ;
	rtn_code = pProgram->initialize() ;
	if (rtn_code != django::RTN_CODE_OK)
	{
		MPI_Finalize() ;
		return(django::RTN_CODE_OK) ;
	}

	//-------------------------------------------------------------------------------------------------------
	// run program
	//-------------------------------------------------------------------------------------------------------
	rtn_code = pProgram->run() ;
	if (rtn_code != django::RTN_CODE_OK)
	{
		MPI_Finalize() ;
		return(django::RTN_CODE_OK) ;
	}

	//-------------------------------------------------------------------------------------------------------
	// finalize program
	//-------------------------------------------------------------------------------------------------------
	rtn_code = pProgram->finalize() ;
	if (rtn_code != django::RTN_CODE_OK)
	{
		MPI_Finalize() ;
		return(django::RTN_CODE_OK) ;
	}
	delete(pProgram) ;

	//-------------------------------------------------------------------------------------------------------
	// print end of output report
	//-------------------------------------------------------------------------------------------------------
	django::print_end_of_output_report() ;

	// finalize MPI environment
	MPI_Finalize() ;

	return(django::RTN_CODE_OK);
}
