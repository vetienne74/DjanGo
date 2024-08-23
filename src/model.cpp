//-------------------------------------------------------------------------------------------------------
//
// PARENT CLASS FOR ALL MODELS
//
//-------------------------------------------------------------------------------------------------------

#include "model.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <fstream>
#include <string>

#include "allocate_array.h"
#include "grid_1D_float.h"
#include "output_report.h"
#include "parse_xml.h"
#include "singleton.h"
#include "type_def.h"
#include "variable.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------

Model::Model()
{
	print_debug(ALL, LIGHT_DEBUG, "IN Model::Model");
	dim      = NO_DIM ;
	type     = NO_MODEL_TYPE ;
	sub_type = NO_MODEL_SUBTYPE ;
	nx       = 0 ;
	ny       = 0 ;
	nz       = 0 ;
	dx       = 0.0 ;
	dy       = 0.0 ;
	dz       = 0.0 ;
	dxrandom = 0.0 ;
	dyrandom = 0.0 ;
	dzrandom = 0.0 ;
	xcoord   = NULL ;
	ycoord   = NULL ;
	zcoord   = NULL ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Model::Model");
}

//-------------------------------------------------------------------------------------------------------

Model::Model(Space_dim dim_in, Model_type type_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Model::Model");

	dim      = dim_in ;
	type     = type_in ;
	print_debug(ALL, LIGHT_DEBUG, "dim", dim);
	print_debug(ALL, LIGHT_DEBUG, "type", type);
	sub_type = NO_MODEL_SUBTYPE ;
	nx       = 0 ;
	ny       = 0 ;
	nz       = 0 ;
	dx       = 0.0 ;
	dy       = 0.0 ;
	dz       = 0.0 ;
	dxrandom = 0.0 ;
	dyrandom = 0.0 ;
	dzrandom = 0.0 ;
	xcoord   = NULL ;
	ycoord   = NULL ;
	zcoord   = NULL ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Model::Model");
}

//-------------------------------------------------------------------------------------------------------

Model::Model(const Model& pModel1, const Model& pModel2, Myfloat alpha)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Model::Model");
	nx       = 0 ;
	ny       = 0 ;
	nz       = 0 ;
	dx       = 0.0 ;
	dy       = 0.0 ;
	dz       = 0.0 ;
	dxrandom = 0.0 ;
	dyrandom = 0.0 ;
	dzrandom = 0.0 ;
	xcoord   = NULL ;
	ycoord   = NULL ;
	zcoord   = NULL ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Model::Model");
}

//-------------------------------------------------------------------------------------------------------

Model::~Model()
{
	print_debug(ALL, LIGHT_DEBUG, "IN ~Model::Model");
	Singleton::Instance()->delete_variable(VP) ;
	Singleton::Instance()->delete_variable(VS) ;
	Singleton::Instance()->delete_variable(RHO) ;
	Singleton::Instance()->delete_variable(LOSS1) ;
	Singleton::Instance()->delete_variable(LOSS2) ;

	if (xcoord != NULL) delete(xcoord) ;
	if (ycoord != NULL) delete(ycoord) ;
	if (zcoord != NULL) delete(zcoord) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT ~Model::Model");
}

//-------------------------------------------------------------------------------------------------------

Variable* Model::get_parameter(Var_type type)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Model::get_parameter");

	// retrieve variable
	Variable* pVariable = Singleton::Instance()->get_variable(type) ;

	// variable not found in config
	if (pVariable == NULL)
	{
		print_error(" No appropriate parameter found in config") ;
		return(NULL) ;
	}

	// allocate grid if null
	else
	{
		if (pVariable->get_grid() == NULL) {
			Rtn_code rtn_code ;
			if (dim == ONE)
			{
				rtn_code = pVariable->initialize_grid(nz, dz) ;
			}
			else if (dim == TWO)
			{
				rtn_code = pVariable->initialize_grid(nz, nx, dz, dx) ;
			}
			if (rtn_code != RTN_CODE_OK)
			{
				print_error(" Error IN Model::get_parameter, variable can not be initialized") ;
				return(NULL) ;
			}
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Model::get_parameter");
	return pVariable ;
}

//-------------------------------------------------------------------------------------------------------

Variable* Model::register_parameter(Var_type type, string file_name)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Model::register_parameter");

	// register parameter as a variable in Singleton
	Variable* pVar = Singleton::Instance()->register_variable(type, file_name) ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Model::register_parameter");

	return(pVar) ;
}

//-------------------------------------------------------------------------------------------------------

Space_dim Model::get_dim(void) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Model::get_dim");
	print_debug(ALL, LIGHT_DEBUG, "OUT Model::get_dim");
	return dim ;
}

//-------------------------------------------------------------------------------------------------------
Myfloat Model::get_dx(void) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Model::get_dx");
	print_debug(ALL, LIGHT_DEBUG, "OUT Model::get_dx");
	return dx ;
}
Myfloat Model::get_dy(void) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Model::get_dy");
	print_debug(ALL, LIGHT_DEBUG, "OUT Model::get_dy");
	return dy ;
}
Myfloat Model::get_dz(void) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Model::get_dz");
	print_debug(ALL, LIGHT_DEBUG, "OUT Model::get_dz");
	return dz ;
}

//-------------------------------------------------------------------------------------------------------
Myint Model::get_nx(void) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Model::get_nx");
	print_debug(ALL, LIGHT_DEBUG, "OUT Model::get_nx");
	return nx ;
}
Myint Model::get_ny(void) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Model::get_ny");
	print_debug(ALL, LIGHT_DEBUG, "OUT Model::get_ny");
	return ny ;
}
Myint Model::get_nz(void) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Model::get_nz");
	print_debug(ALL, LIGHT_DEBUG, "OUT Model::get_nz");
	return nz ;
}

//-------------------------------------------------------------------------------------------------------
Grid_1D_float* Model::get_xcoord(void) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Model::get_xcoord");
	print_debug(ALL, LIGHT_DEBUG, "OUT Model::get_xcoord");
	return xcoord ;
}
Grid_1D_float* Model::get_ycoord(void) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Model::get_ycoord");
	print_debug(ALL, LIGHT_DEBUG, "OUT Model::get_ycoord");
	return ycoord ;
}
Grid_1D_float* Model::get_zcoord(void) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Model::get_zcoord");
	print_debug(ALL, LIGHT_DEBUG, "OUT Model::get_zcoord");
	return zcoord ;
}

//-------------------------------------------------------------------------------------------------------
Model_type Model::get_type(void) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Model::get_type");
	return type ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Model::get_type");
}

//-------------------------------------------------------------------------------------------------------
Model_sub_type Model::get_sub_type(void) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Model::get_sub_type");
	return sub_type ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Model::get_sub_type");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Model::set_size(Myint nx_in, Myint ny_in, Myint nz_in) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Model::set_size");

	nx = nx_in ;
	ny = ny_in ;
	nz = nz_in ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Model::set_size");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Model::set_sampling(Myfloat dx_in, Myfloat dy_in, Myfloat dz_in) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Model::set_sampling");

	dx = dx_in ;
	dy = dy_in ;
	dz = dz_in ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Model::set_sampling");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Model::set_sampling_random(Myfloat dxrandom_in, Myfloat dyrandom_in, Myfloat dzrandom_in) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Model::set_sampling_random");

	dxrandom = dxrandom_in ;
	dyrandom = dyrandom_in ;
	dzrandom = dzrandom_in ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Model::set_sampling_random");
	return(RTN_CODE_OK) ;
}


//-------------------------------------------------------------------------------------------------------

Rtn_code Model::initialize()
{
	print_debug(ALL, LIGHT_DEBUG, "IN Model::initialize");

	sub_type = REGULAR ;

	// initialize zcoord
	//------------------

	if (dim >= ONE)
	{
		zcoord = new Grid_1D_float(nz, dz) ;

		if (dzrandom == 0.0)
		{
			for (Myint ii=0; ii<nz; ii++)
			{
				zcoord->pArray[ii] = ii*dz ;
			}
		}
		else
		{
			sub_type = UNREGULAR ;
			Myfloat coord_tmp = 0 ;

			// distribute random coordinates
			for (Myint ii=0; ii<nz; ii++)
			{
				zcoord->pArray[ii] = coord_tmp ;
				Myfloat aa = Myfloat(rand()) / RAND_MAX ;
				Myint   bb = rand() % 2 ;
				if (bb == 0) aa *= -1 ;
				Myfloat hh = dz * (1. + dzrandom * aa) ;
				coord_tmp += hh ;
			}

			// reajust if needed
			Myfloat64 alpha = (nz-1)*dz / zcoord->pArray[nz-1] ;
			for (Myint ii=0; ii<nz; ii++)
			{
				zcoord->pArray[ii] *= alpha ;
			}
		}
	}

	// initialize xcoord
	//------------------

	if (dim >= TWO)
	{
		xcoord = new Grid_1D_float(nx, dx) ;

		if (dxrandom == 0.0)
		{
			for (Myint ii=0; ii<nx; ii++)
			{
				xcoord->pArray[ii] = ii*dx ;
			}
		}
		else
		{
			sub_type = UNREGULAR ;
			Myfloat coord_tmp = 0 ;

			// distribute random coordinates
			for (Myint ii=0; ii<nx; ii++)
			{
				xcoord->pArray[ii] = coord_tmp ;
				Myfloat aa = Myfloat(rand()) / RAND_MAX ;
				Myint   bb = rand() % 2 ;
				if (bb == 0) aa *= -1 ;
				Myfloat hh = dx * (1. + dxrandom * aa) ;
				coord_tmp += hh ;
			}

			// reajust if needed
			Myfloat alpha = (nx-1)*dx / xcoord->pArray[nx-1] ;
			for (Myint ii=0; ii<nx; ii++)
			{
				xcoord->pArray[ii] *= alpha ;
			}
		}
	}

	// initialize ycoord
	//------------------

	if (dim >= THREE)
	{
		ycoord = new Grid_1D_float(ny, dy) ;

		if (dyrandom == 0.0)
		{
			for (Myint ii=0; ii<ny; ii++)
			{
				ycoord->pArray[ii] = ii*dy ;
			}
		}
		else
		{
			sub_type = UNREGULAR ;
			Myfloat coord_tmp = 0 ;

			// distribute random coordinates
			for (Myint ii=0; ii<ny; ii++)
			{
				ycoord->pArray[ii] = coord_tmp ;
				Myfloat aa = Myfloat(rand()) / RAND_MAX ;
				Myint   bb = rand() % 2 ;
				if (bb == 0) aa *= -1 ;
				Myfloat hh = dy * (1. + dyrandom * aa) ;
				coord_tmp += hh ;
			}

			// reajust if needed
			Myfloat alpha = (ny-1)*dy / ycoord->pArray[ny-1] ;
			for (Myint ii=0; ii<ny; ii++)
			{
				ycoord->pArray[ii] *= alpha ;
			}
		}
	}

	// ouput grid point in ascii file
	ofstream out_file ;
	out_file.open(GRID_POINT_OUT_FILE, ios::trunc | ios::out) ;

	if (dim == ONE)
	{
		for (Myint iz=0; iz<nz; iz++)
		{
			out_file << zcoord->pArray[iz] << "\n";
		}
	}
	else if (dim == TWO)
	{
		for (Myint ix=0; ix<nx; ix++)
		{
			for (Myint iz=0; iz<nz; iz++)
			{
				out_file << zcoord->pArray[iz] << " " << xcoord->pArray[ix] << "\n";
			}
		}
	}
	else if (dim == THREE)
	{
		for (Myint iy=0; iy<ny; iy++)
		{
			for (Myint ix=0; ix<nx; ix++)
			{
				for (Myint iz=0; iz<nz; iz++)
				{
					out_file << zcoord->pArray[iz] << " " << xcoord->pArray[ix] << " " << ycoord->pArray[iy] << "\n";
				}
			}
		}
	}

	out_file.close() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Model::initialize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Model::add_perturbation(Myint iz, Myint iy, Myint ix, Myfloat perturb)
{
	print_error(" Model::add_perturbation should not be called") ;
	return(RTN_CODE_KO) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Model::read()
{
	print_error(" Model::read should not be called") ;
	return(RTN_CODE_KO) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Model::write()
{
	print_error(" Model::write should not be called") ;
	return(RTN_CODE_KO) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Model::write(Myint iter)
{
	print_error(" Model::write should not be called") ;
	return(RTN_CODE_KO) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Model::info()
{
	print_debug(ALL, LIGHT_DEBUG, "IN Model::info");

	print_info(MASTER, "") ;
	print_line2() ;
	print_info(MASTER, " MODEL PARAMETERS") ;
	print_info(MASTER, "") ;

	// model dimension
	switch(dim)
	{
	case ONE:
		print_info(MASTER, " Dimension\t", "1D") ;
		break ;
	case TWO:
		print_info(MASTER, " Dimension\t", "2D") ;
		break ;
	case THREE:
		print_info(MASTER, " Dimension\t", "3D") ;
		break ;
	default:
		print_error(" Invalid model dimension", dim) ;
		return(RTN_CODE_KO) ;
	}

	// model type
	switch(type)
	{
	case GRID:
		print_info(MASTER, " Type\t\t", "GRID") ;
		break ;
	default:
		print_error(" Invalid model type", type) ;
		return(RTN_CODE_KO) ;
	}

	// model sub-type
	switch(sub_type)
	{
	case REGULAR:
		print_info(MASTER, " Sub-type\t", "REGULAR") ;
		break ;
	case UNREGULAR:
		print_info(MASTER, " Sub-type\t", "UNREGULAR") ;
		break ;
	default:
		print_error(" Invalid model sub-type", sub_type) ;
		return(RTN_CODE_KO) ;
	}

	// size
	if (dim >= ONE)
	{
		print_info(MASTER, " Nz\t\t", nz) ;
	}
	if (dim >= TWO)
	{
		print_info(MASTER, " Nx\t\t", nx) ;
	}
	if (dim >= THREE)
	{
		print_info(MASTER, " Ny\t\t", ny) ;
	}

	// sampling
	if (dim >= ONE)
	{
		if (sub_type == REGULAR)
		{
			print_info(MASTER, " Dz\t\t", dz) ;
		}
		else
		{
			// retrieve grid min and max sampling
			Myfloat min_hh = zcoord->pArray[1] - zcoord->pArray[0] ;
			Myfloat max_hh = min_hh ;
			for (Myint ii=0; ii<nz-1; ii++)
			{
				// current grid sampling
				Myfloat hh = zcoord->pArray[ii+1] - zcoord->pArray[ii] ;
				min_hh = min(min_hh, hh) ;
				max_hh = max(max_hh, hh) ;
			}
			print_info(MASTER, " Min. Dz\t", min_hh) ;
			print_info(MASTER, " Max. Dz\t", max_hh) ;
		}
	}
	if (dim >= TWO)
	{
		if (sub_type == REGULAR)
		{
			print_info(MASTER, " Dx\t\t", dx) ;
		}
		else
		{
			// retrieve grid min and max sampling
			Myfloat min_hh = xcoord->pArray[1] - xcoord->pArray[0] ;
			Myfloat max_hh = min_hh ;
			for (Myint ii=0; ii<nx-1; ii++)
			{
				// current grid sampling
				Myfloat hh = xcoord->pArray[ii+1] - xcoord->pArray[ii] ;
				min_hh = min(min_hh, hh) ;
				max_hh = max(max_hh, hh) ;
			}
			print_info(MASTER, " Min. Dx\t", min_hh) ;
			print_info(MASTER, " Max. Dx\t", max_hh) ;
		}
	}
	if (dim >= THREE)
	{
		if (sub_type == REGULAR)
		{
			print_info(MASTER, " Dy\t\t", dy) ;
		}
		else
		{
			// retrieve grid min and max sampling
			Myfloat min_hh = ycoord->pArray[1] - ycoord->pArray[0] ;
			Myfloat max_hh = min_hh ;
			for (Myint ii=0; ii<ny-1; ii++)
			{
				// current grid sampling
				Myfloat hh = ycoord->pArray[ii+1] - ycoord->pArray[ii] ;
				min_hh = min(min_hh, hh) ;
				max_hh = max(max_hh, hh) ;
			}
			print_info(MASTER, " Min. Dy\t", min_hh) ;
			print_info(MASTER, " Max. Dy\t", max_hh) ;
		}
	}

	// physical size
	if (dim >= ONE)
	{
		print_info(MASTER, " Size z (m)\t", zcoord->get_max()) ;
	}
	if (dim >= TWO)
	{
		print_info(MASTER, " Size x (m)\t", xcoord->get_max()) ;
	}
	if (dim >= THREE)
	{
		print_info(MASTER, " Size y (m)\t", ycoord->get_max()) ;
	}
	print_info(MASTER, "") ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Model::info");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Model::reset()
{
	print_error(" Model::reset should not be called") ;
	return(RTN_CODE_KO) ;
}

} // namespace django

