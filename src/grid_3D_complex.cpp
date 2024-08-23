//-------------------------------------------------------------------------------------------------------
//
// GRID 3D COMPLEX
//
//-------------------------------------------------------------------------------------------------------

#include "grid_3D_complex.h"

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

Rtn_code Grid_3D_complex::reset()
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_3D_complex::reset");

	for (Myint ix=0; ix<nx; ix++)
	{
		for (Myint iy=0; iy<ny; iy++)
		{
			for (Myint iz=0; iz<nz; iz++)
			{
				pArray[ix][iy][iz] = ZERO_CMPLX ;
			}
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_3D_complex::reset");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

//#pragma ivdep not working with complex !

Rtn_code Grid_3D_complex::reset(Myfloat val)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_3D_complex::reset(Myfloat val)");

	for (Myint ix=0; ix<nx; ix++)
	{
		for (Myint iy=0; iy<ny; iy++)
		{
			for (Myint iz=0; iz<nz; iz++)
			{
				pArray[ix][iy][iz] = val ;
			}
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_3D_complex::reset(Myfloat val)");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Grid_3D_complex::Grid_3D_complex(Myint nz2, Myint nx2, Myint ny2, Myfloat dz2, Myfloat dx2, Myfloat dy2)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_3D_complex::Grid_3D_complex");
	nz = nz2 ;
	nx = nx2 ;
	ny = ny2 ;
	dz = dz2 ;
	dx = dx2 ;
	dy = dy2 ;
	pArray = allocate_array<Mycomplex>(nx, ny, nz) ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_3D_complex::Grid_3D_complex");
}

//-------------------------------------------------------------------------------------------------------

Grid_3D_complex::~Grid_3D_complex(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_3D_complex::~Grid_3D_complex");
	deallocate_array<Mycomplex>(pArray, nx, ny, nz) ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_3D_complex::~Grid_3D_complex");
}


//-------------------------------------------------------------------------------------------------------

Myfloat Grid_3D_complex::get_min()
{

	Myfloat min_val = FLT_MAX ;

	for (Myint ix=0; ix<nx; ix++)
	{
		for (Myint iy=0; iy<ny; iy++)
		{
			for (Myint iz=0; iz<nz; iz++)
			{
				if (abs(pArray[ix][iy][iz]) < min_val)
				{
					min_val = abs(pArray[ix][iy][iz]) ;
				}
			}
		}
	}
	return min_val ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid_3D_complex::get_max()
{

	Myfloat max_val = -FLT_MAX ;

	for (Myint ix=0; ix<nx; ix++)
	{
		for (Myint iy=0; iy<ny; iy++)
		{
			for (Myint iz=0; iz<nz; iz++)
			{
				if (abs(pArray[ix][iy][iz]) > max_val)
				{
					max_val = abs(pArray[ix][iy][iz]) ;
				}
			}
		}
	}
	return max_val ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_3D_complex::write_on_disk(string file_name)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_3D_complex::write_on_disk");
	print_debug(ALL, LIGHT_DEBUG, "write file", file_name.c_str());
	ofstream out_file ;

	out_file.open(file_name.c_str(), ios::binary) ;
	assert(out_file.is_open());
	out_file.write((char*) **pArray, get_nb_grid_point() * sizeof(Mycomplex)) ;
	out_file.close() ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_3D_complex::write_on_disk");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

void Grid_3D_complex::info(void)
{
	Grid_3D::info() ;
	print_info(ALL, "min", get_min() ) ;
	print_info(ALL, "max", get_max() ) ;
}  

} // namespace django
