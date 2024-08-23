//-------------------------------------------------------------------------------------------------------
//
// CLASS TO HANDLE VARIABLES
//
//-------------------------------------------------------------------------------------------------------

#include "variable.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>

#include "mpi.h"

#include "allocate_array.h"
#include "grid_1D_float.h"
#include "output_report.h"
#include "singleton.h"

namespace django {

//-------------------------------------------------------------------------------------------------------

Variable::Variable(const char* name_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Variable::Variable");

	name           = name_in ;
	type           = NO_VAR_TYPE ;
	pGrid          = NULL ;

	nb_byte_write  = 0.0 ;
	time_in_write  = 0.0 ;
	pIO_buffer     = NULL ;
	idx_buffer     = 0 ;
	max_buffer     = 0 ;

	self_init      = NO_SELF_INIT ;
	self_init_file = UNSPECIFIED ;
	self_init_const = 0.0 ;

	nrec = 0 ;
	id   = 0 ;

	print_debug(ALL, LIGHT_DEBUG, "** create variable", name.c_str());
	print_debug(ALL, LIGHT_DEBUG, "** type", (Myint) type);
	print_debug(ALL, LIGHT_DEBUG, "OUT Variable::Variable");
}

//-------------------------------------------------------------------------------------------------------

Variable::Variable(Var_type type_in, string name_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Variable::Variable");

	name           = name_in ;
	type           = type_in ;
	pGrid          = NULL ;

	nb_byte_write  = 0.0 ;
	time_in_write  = 0.0 ;
	pIO_buffer     = NULL ;
	idx_buffer     = 0 ;
	max_buffer     = 0 ;

	self_init      = NO_SELF_INIT ;
	self_init_file = UNSPECIFIED ;
	self_init_const = 0.0 ;

	nrec = 0 ;
	id   = 0 ;

	print_debug(ALL, LIGHT_DEBUG, "** create variable", name.c_str());
	print_debug(ALL, LIGHT_DEBUG, "** type", (Myint) type);
	print_debug(ALL, LIGHT_DEBUG, "OUT Variable::Variable");
}

//-------------------------------------------------------------------------------------------------------

Variable::~Variable(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Variable::~Variable");
	print_debug(ALL, LIGHT_DEBUG, "** delete variable", name.c_str());

	if (pIO_buffer != NULL)
	{
		flush_rec() ;
		delete(pIO_buffer) ;

		if (nb_byte_write > 0)
		{
			print_info(ALL, "") ;
			print_info(ALL, " IO stat for Variable", name) ;
			print_info(ALL, " Bytes written:\t", nb_byte_write) ;
			print_info(ALL, " Write time:\t", (Myfloat) time_in_write) ;
			Myfloat io_speed = nb_byte_write / (Myfloat) time_in_write ;
			if (io_speed < 1.e+3)
			{
				print_info(ALL, " Write speed (B/s):", io_speed) ;
			}
			else if (io_speed < 1.e+6)
			{
				print_info(ALL, " Write speed (KB/s):", io_speed / Myfloat(1.e+3)) ;
			}
			else
			{
				print_info(ALL, " Write speed (MB/s):", io_speed / Myfloat(1.e+6)) ;
			}
		}
	}

	if (pGrid != NULL) delete(pGrid) ;
	pGrid = NULL ;
	name  = "*** DELETED ***" ;
	type  = NO_VAR_TYPE ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Variable::~Variable");
}

//-------------------------------------------------------------------------------------------------------

Grid* Variable::get_grid(void)
{
	print_debug(ALL, FULL_DEBUG, "IN Variable::get_grid");
	print_debug(ALL, FULL_DEBUG, "OUT Variable::get_grid");
	return(pGrid) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Variable::allocate_grid(Myint nz, Myfloat dz)
{
	print_debug(ALL, FULL_DEBUG, "IN Variable::allocate_grid(Myint nz, Myfloat dz)");

	// grid already allocated
	if (pGrid != NULL)
	{
		print_error(" Error IN Variable::allocate_grid, grid already allocated") ;
		pGrid->info() ;
		return(RTN_CODE_KO) ;
	}

	pGrid = new Grid_1D_float(nz, dz) ;

	print_debug(ALL, FULL_DEBUG, "OUT Variable::allocate_grid(Myint nz, Myfloat dz)");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Variable::allocate_grid(Myint nz, Myint nx, Myfloat dz, Myfloat dx)
{
	print_debug(ALL, FULL_DEBUG, "IN Variable::allocate_grid(Myint nz, Myint nx, Myfloat dz, Myfloat dx)");

	// grid already allocated
	if (pGrid != NULL)
	{
		print_error(" Error IN Variable::allocate_grid, grid already allocated") ;
		pGrid->info() ;
		return(RTN_CODE_KO) ;
	}

	pGrid = new Grid_2D_float(nz, nx, dz, dx) ;

	print_debug(ALL, FULL_DEBUG, "OUT Variable::allocate_grid(Myint nz, Myint nx, Myfloat dz, Myfloat dx)");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Variable::initialize_grid(Myint nz, Myfloat dz)
{
	print_debug(ALL, FULL_DEBUG, "IN Variable::initialize_grid(Myint nz, Myfloat dz)");

	// initialize grid if needed
	if (pGrid == NULL)
	{
		// allocate grid
		this->allocate_grid(nz, dz) ;
		Grid_1D_float* vp_grid = (Grid_1D_float *) pGrid ;

		// initialize grid
		print_info(MASTER, " Initialize parameter", name) ;
		if (self_init == NO_SELF_INIT)
		{
			print_error(" Error IN Variable::get_grid, variable can not self-initialized") ;
			return(RTN_CODE_KO) ;
		}
		else if (self_init == INIT_FROM_CONST)
		{
			print_info(MASTER, " With constant value", self_init_const) ;
			vp_grid->reset(self_init_const) ;
		}
		else if (self_init == INIT_FROM_FILE)
		{
			print_info(MASTER, " From file\t", self_init_file) ;
			vp_grid->read_from_disk(self_init_file) ;
			print_info(MASTER, " Min. value \t", vp_grid->get_min()) ;
			print_info(MASTER, " Max. value \t", vp_grid->get_max()) ;
		}
		print_info(MASTER, "") ;
	}

	print_debug(ALL, FULL_DEBUG, "OUT Variable::initialize_grid(Myint nz, Myfloat dz)");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Variable::initialize_grid(Myint nz, Myint nx, Myfloat dz, Myfloat dx)
{
	print_debug(ALL, FULL_DEBUG, "IN Variable::initialize_grid(Myint nz, Myint nx, Myfloat dz, Myfloat dx)");

	// initialize grid if needed
	if (pGrid == NULL)
	{
		// allocate grid
		this->allocate_grid(nz, nx, dz, dx) ;
		Grid_2D_float* vp_grid = (Grid_2D_float *) pGrid ;

		// initialize grid
		print_info(MASTER, " Initialize parameter", name) ;
		if (self_init == NO_SELF_INIT)
		{
			print_error(" Error IN Variable::get_grid, variable can not self-initialized") ;
			return(RTN_CODE_KO) ;
		}
		else if (self_init == INIT_FROM_CONST)
		{
			print_info(MASTER, " With constant value", self_init_const) ;
			vp_grid->reset(self_init_const) ;
		}
		else if (self_init == INIT_FROM_FILE)
		{
			print_info(MASTER, " From file\t", self_init_file) ;
			vp_grid->read_from_disk(self_init_file) ;
		}
		print_info(MASTER, " Min. (m/s) \t", vp_grid->get_min()) ;
		print_info(MASTER, " Max. (m/s) \t", vp_grid->get_max()) ;
		print_info(MASTER, "") ;
	}

	print_debug(ALL, FULL_DEBUG, "OUT Variable::initialize_grid(Myint nz, Myint nx, Myfloat dz, Myfloat dx)");
	return(RTN_CODE_OK) ;
}


//-------------------------------------------------------------------------------------------------------

Var_type Variable::get_type(void)
{
	print_debug(ALL, FULL_DEBUG, "IN Variable::get_type");
	print_debug(ALL, FULL_DEBUG, "type", type);
	print_debug(ALL, FULL_DEBUG, "OUT Variable::get_type");
	return(type) ;
}

//-------------------------------------------------------------------------------------------------------

Myint Variable::get_id(void)
{
	print_debug(ALL, FULL_DEBUG, "IN Variable::get_id");
	print_debug(ALL, FULL_DEBUG, "id", id);
	print_debug(ALL, FULL_DEBUG, "OUT Variable::get_id");
	return(id) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Variable::set_type(Var_type type_in)
{
	print_debug(ALL, FULL_DEBUG, "IN Variable::set_type");
	print_debug(ALL, FULL_DEBUG, "old type", type);
	type = type_in ;
	print_debug(ALL, FULL_DEBUG, "new type", type);
	print_debug(ALL, FULL_DEBUG, "OUT Variable::set_type");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Variable::set_id(Myint id_in)
{
	print_debug(ALL, FULL_DEBUG, "IN Variable::set_id");
	id = id_in ;
	print_debug(ALL, FULL_DEBUG, "OUT Variable::set_id");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Variable::set_self_init_mode(Self_init_type self_init_in)
{
	print_debug(ALL, FULL_DEBUG, "IN Variable::set_self_init_mode");
	self_init = self_init_in ;
	print_debug(ALL, FULL_DEBUG, "OUT Variable::set_self_init_mode");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Variable::set_self_init_file(string self_init_file_in)
{
	print_debug(ALL, FULL_DEBUG, "IN Variable::set_self_init_file");
	self_init_file = self_init_file_in ;
	print_debug(ALL, FULL_DEBUG, "OUT Variable::set_self_init_file");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Variable::set_self_init_const(Myfloat self_init_const_in)
{
	print_debug(ALL, FULL_DEBUG, "IN Variable::set_self_init_const");
	self_init_const = self_init_const_in ;
	print_debug(ALL, FULL_DEBUG, "OUT Variable::set_self_init_const");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Variable::set_grid(Grid* pGrid_in)
{
	print_debug(ALL, FULL_DEBUG, "IN Variable::set_grid");
	pGrid = pGrid_in ;
	print_debug(ALL, FULL_DEBUG, "OUT Variable::set_grid");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Variable::reset_grid()
{
	print_debug(ALL, FULL_DEBUG, "IN Variable::reset_grid");

	Rtn_code rtn_code = pGrid->reset() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, FULL_DEBUG, "OUT Variable::reset_grid");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

string Variable::get_name(void)
{
	print_debug(ALL, FULL_DEBUG, "IN Variable::get_name");
	print_debug(ALL, FULL_DEBUG, "type", name);
	print_debug(ALL, FULL_DEBUG, "OUT Variable::get_name");
	return(name) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Variable::read(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Variable::read");

	// read grid
	//pGrid->read() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Variable::read");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Variable::write(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Variable::write");

	// write grid
	//pGrid->write() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Variable::write");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Variable::write(string filename)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Variable::write");

	// write grid
	Rtn_code rtn_code = pGrid->write_on_disk(filename) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Variable::write");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

void Variable::info(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Variable::info");

	// type
	switch(type)
	{
	case NO_VAR_TYPE:
		print_info(MASTER, " Type\t\t", "UNSPECIFIED") ;
		break ;
	case VP:
		print_info(MASTER, " Type\t\t", "P-wave velocity (m/s)") ;
		break ;
	case VS:
		print_info(MASTER, " Type\t\t", "S-wave velocity (m/s)") ;
		break ;
	case RHO:
		print_info(MASTER, " Type\t\t", "Density (kg/m3)") ;
		break ;
	case LOSS1:
		print_info(MASTER, " Type\t\t", "Loss1 term (???)") ;
		break ;
	case LOSS2:
		print_info(MASTER, " Type\t\t", "Loss2 term (???)") ;
		break ;
	case VX:
		print_info(MASTER, " Type\t\t", "Particle velocity along x (m/s)") ;
		break ;
	case VY:
		print_info(MASTER, " Type\t\t", "Particle velocity along y (m/s)") ;
		break ;
	case VZ:
		print_info(MASTER, " Type\t\t", "Particle velocity along z (m/s)") ;
		break ;
	case PR:
		print_info(MASTER, " Type\t\t", "Pressure (Pa)") ;
		break ;
	case SXX:
		print_info(MASTER, " Type\t\t", "Stress sigma_xx (Pa)") ;
		break ;
	case SZZ:
		print_info(MASTER, " Type\t\t", "Stress sigma_zz (Pa)") ;
		break ;
	case SXZ:
		print_info(MASTER, " Type\t\t", "Stress sigma_xz (Pa)") ;
		break ;
	default:
		print_info(MASTER, " Type\t\t", "INVALID") ;
		break ;
	}

	// name
	print_info(MASTER, " Name\t\t", name) ;

	// grid
	if (pGrid == NULL)
	{
		print_info(MASTER, " Grid not initialized") ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Variable::info");
}

//-------------------------------------------------------------------------------------------------------
// 1d grid
// write at recceivers locations
// use buffer mechanism
Rtn_code Variable::write_rec(Myint nrec_in, Myint* iz_rec)
{
	print_debug(ALL, FULL_DEBUG, "IN Variable::write_rec(Myint nrec_in, Myint* iz_rec)");

	// initialize buffer if not done already
	if (pIO_buffer == NULL)
	{
		nrec = nrec_in ;
		max_buffer = max(1, Singleton::Instance()->iobufsize / (int)sizeof(Myfloat) / nrec) ;
		print_debug(ALL, LIGHT_DEBUG, " Max IO buffer size", Singleton::Instance()->iobufsize) ;
		print_debug(ALL, LIGHT_DEBUG, " -> Nt in IO buffer", max_buffer) ;
		pIO_buffer = new Grid_2D_float (nrec, max_buffer, 0, 0) ;
		print_debug(ALL, LIGHT_DEBUG," -> IO buffer size", max_buffer * nrec * (int)sizeof(Myfloat)) ;

		idx_buffer = 0 ;
	}

	// write buffer if full
	if (idx_buffer == max_buffer)
	{
		flush_rec() ;
	}

	// add data to buffer
	Grid_1D_float* grid_u = dynamic_cast<Grid_1D_float*>(pGrid) ;
	if (grid_u == NULL)
	{
		print_error(" Error IN Variable::write_rec -> pGrid is not Grid_1D_float") ;
		return(RTN_CODE_KO) ;
	}
	Myfloat* time_u = grid_u->pArray ;

	for (Myint irec=0; irec<nrec; irec++)
	{
		if (iz_rec[irec] == NOT_FOUND)
		{
			pIO_buffer->pArray[idx_buffer][irec] = 0.0 ;
		}
		else
		{
			pIO_buffer->pArray[idx_buffer][irec] = time_u[iz_rec[irec]] ;
		}
	}

	// increment index
	idx_buffer++ ;

	print_debug(ALL, FULL_DEBUG, "OUT Variable::write_rec(Myint nrec_in, Myint* iz_rec)");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
// 2d grid
// write at recceivers locations
// use buffer mechanism
Rtn_code Variable::write_rec(Myint nrec_in, Myint* iz_rec, Myint* ix_rec)
{
	print_debug(ALL, FULL_DEBUG, "IN Variable::write_rec(Myint nrec_in, Myint* iz_rec. Myint* ix_rec)");

	// initialize buffer if not done already
	if (pIO_buffer == NULL)
	{
		nrec = nrec_in ;
		max_buffer = max(1, Singleton::Instance()->iobufsize / (int)sizeof(Myfloat) / nrec) ;
		print_debug(ALL, LIGHT_DEBUG, " Max IO buffer size", Singleton::Instance()->iobufsize) ;
		print_debug(ALL, LIGHT_DEBUG, " -> Nt in IO buffer", max_buffer) ;
		pIO_buffer = new Grid_2D_float (nrec, max_buffer, 0, 0) ;
		print_debug(ALL, LIGHT_DEBUG," -> IO buffer size", max_buffer * nrec * (int)sizeof(Myfloat)) ;

		idx_buffer = 0 ;
	}

	// write buffer if full
	if (idx_buffer == max_buffer)
	{
		flush_rec() ;
	}

	// add data to buffer
	Grid_2D_float* grid_u = dynamic_cast<Grid_2D_float*>(pGrid) ;
	if (grid_u == NULL)
	{
		print_error(" Error IN Variable::write_rec -> pGrid is not Grid_2D_float") ;
		return(RTN_CODE_KO) ;
	}
	Myfloat** time_u = grid_u->pArray ;
	for (Myint irec=0; irec<nrec; irec++)
	{
		if ((iz_rec[irec] == NOT_FOUND) || (ix_rec[irec] == NOT_FOUND))
		{
			pIO_buffer->pArray[idx_buffer][irec] = 0.0 ;
		}
		else
		{
			pIO_buffer->pArray[idx_buffer][irec] = time_u[ix_rec[irec]][iz_rec[irec]] ;
		}
	}

	// increment index
	idx_buffer++ ;

	print_debug(ALL, FULL_DEBUG, "OUT Variable::write_rec(Myint nrec_in, Myint* iz_rec. Myint* ix_rec)");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Variable::flush_rec(void)
{
	print_debug(ALL, FULL_DEBUG, "IN Variable::flush_rec");

	// time measure for IO
	double t0 = MPI_Wtime() ;

	// check buffer has been allocated
	if (pIO_buffer != NULL)
	{
		// open, write and close
		if (idx_buffer != 0)
		{
			string file_name = name + TIME_REC_OUT_FILE ;
			ofstream pFile ;
			pFile.open(file_name.c_str(), ios::binary | ios::app | ios::out) ;
			assert(pFile.is_open());
			Myint byte_size = nrec * idx_buffer * sizeof(Myfloat) ;
			pFile.write((char*)&(pIO_buffer->pArray[0][0]), byte_size) ;
			pFile.close() ;

			print_debug(ALL, FULL_DEBUG, "write bytes", byte_size) ;

			// increment IO counter
			nb_byte_write += byte_size ;
			time_in_write += MPI_Wtime() - t0 ;

			// reinit buffer index
			idx_buffer = 0 ;
		}
	}

	print_debug(ALL, FULL_DEBUG, "OUT Variable::flush_rec");
	return(RTN_CODE_OK) ;
}

} // namespace django
