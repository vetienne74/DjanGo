//-------------------------------------------------------------------------------------------------------
//
// GRID 2D COMPLEX
//
//-------------------------------------------------------------------------------------------------------

#include "grid_2D_complex.h"

#include <cfloat>
#include <fstream>

#include "allocate_array.h"
#include "constant.h"
#include "output_report.h"
#include "type_def.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------

//#pragma ivdep not working with complex !

Rtn_code Grid_2D_complex::reset()
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_2D_complex::reset()");

	for (Myint ix=0; ix<nx; ix++)
	{
		for (Myint iz=0; iz<nz; iz++)
		{
			pArray[ix][iz] = ZERO_CMPLX ;
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_2D_complex::reset()");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

//#pragma ivdep not working with complex !

Rtn_code Grid_2D_complex::reset(Myfloat val)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_2D_complex::reset(Myfloat val)");

	for (Myint ix=0; ix<nx; ix++)
	{
		for (Myint iz=0; iz<nz; iz++)
		{
			pArray[ix][iz] = val ;
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_2D_complex::reset(Myfloat val)");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Grid_2D_complex::Grid_2D_complex(Myint nz2, Myint nx2, Myfloat dz2, Myfloat dx2)
{
	nz = nz2 ;
	nx = nx2 ;
	dz = dz2 ;
	dx = dx2 ;
	pArray = allocate_array<Mycomplex>(nx, nz) ;
}


//-------------------------------------------------------------------------------------------------------

Myfloat Grid_2D_complex::get_min()
{

	Myfloat min_val = FLT_MAX ;

	for (Myint ix=0; ix<nx; ix++)
	{
		for (Myint iz=0; iz<nz; iz++)
		{
			if (abs(pArray[ix][iz]) < min_val)
			{
				min_val = abs(pArray[ix][iz]) ;
			}
		}
	}
	return min_val ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid_2D_complex::get_max()
{

	Myfloat max_val = -FLT_MAX ;

	for (Myint ix=0; ix<nx; ix++)
	{
		for (Myint iz=0; iz<nz; iz++)
		{
			if (abs(pArray[ix][iz]) > max_val)
			{
				max_val = abs(pArray[ix][iz]) ;
			}
		}
	}
	return max_val ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_2D_complex::write_on_disk(string file_name)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_2D_complex::write_on_disk");
	print_debug(ALL, LIGHT_DEBUG, "write file", file_name.c_str());
	ofstream out_file ;

	out_file.open(file_name.c_str(), ios::binary | ios::app) ;
	assert(out_file.is_open());
	out_file.write((char*) *pArray, get_nb_grid_point() * sizeof(Mycomplex)) ;
	out_file.close() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_2D_complex::write_on_disk");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

void Grid_2D_complex::info(void)
{
	Grid_2D::info() ;
	print_info(ALL, "min", get_min() ) ;
	print_info(ALL, "max", get_max() ) ;
}

} // namespace django
