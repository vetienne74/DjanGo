//-------------------------------------------------------------------------------------------------------
//
// GRID 3D FLOAT
//
//-------------------------------------------------------------------------------------------------------

#include "grid_3D_float.h"

#include <cfloat>
#include <fstream>

#include "allocate_array.h"
#include "type_def.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------

Grid_3D_float::~Grid_3D_float(void)
{  
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_3D_float::~Grid_3D_float");
	deallocate_array<Myfloat>(pArray, nx, ny, nz) ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_3D_float::~Grid_3D_float");
}

//-------------------------------------------------------------------------------------------------------

Grid_3D_float::Grid_3D_float(Myint nz2, Myint ny2, Myint nx2, Myfloat dz2, Myfloat dy2, Myfloat dx2)
{
	nz = nz2 ;
	nx = nx2 ;
	ny = ny2 ;
	dz = dz2 ;
	dx = dx2 ;
	dy = dy2 ;
	pArray = allocate_array<Myfloat>(nx, ny, nz) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid_3D_float::get_min()
{

	Myfloat min_val = FLT_MAX ;

	for (Myint ix=0; ix<nx; ix++)
	{
		for (Myint iy=0; iy<ny; iy++)
		{
			for (Myint iz=0; iz<nz; iz++)
			{
				if (pArray[ix][iy][iz] < min_val)
				{
					min_val = pArray[ix][iy][iz] ;
				}
			}
		}
	}
	return min_val ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid_3D_float::get_max()
{

	Myfloat max_val = -FLT_MAX ;

	for (Myint ix=0; ix<nx; ix++)
	{
		for (Myint iy=0; iy<ny; iy++)
		{
			for (Myint iz=0; iz<nz; iz++)
			{
				if (pArray[ix][iy][iz] > max_val)
				{
					max_val = pArray[ix][iy][iz] ;
				}
			}
		}
	}
	return max_val ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_3D_float::read_from_disk(string file_name)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_3D_float::read_from_disk");

	ifstream in_file ;
	in_file.open(file_name.c_str(), ios::binary) ;
	assert(in_file.is_open());

	print_debug(ALL, LIGHT_DEBUG, " reading from file #values", get_nb_grid_point()) ;

	//Myfloat32 pArray_tmp[get_nb_grid_point()] ;
	Myfloat32* pArray_tmp = allocate_array<Myfloat32>(get_nb_grid_point()) ;

	in_file.read((char*) pArray_tmp, get_nb_grid_point() * sizeof(Myfloat32)) ;
	in_file.close() ;

	Myfloat* pArray_tmp2 = (Myfloat*) (**pArray) ;
#pragma ivdep
	for (Myint ii=0; ii < get_nb_grid_point(); ii++)
	{
		pArray_tmp2[ii] = pArray_tmp[ii] ;
	}

	deallocate_array<Myfloat32>(pArray_tmp, get_nb_grid_point() ) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_3D_float::read_from_disk");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_3D_float::write_on_disk(string file_name)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_3D_float::write_on_disk");
	print_debug(ALL, LIGHT_DEBUG, "write file", file_name.c_str());
	ofstream out_file ;

	out_file.open(file_name.c_str(), ios::binary) ;
	assert(out_file.is_open());
	out_file.write((char*) **pArray, get_nb_grid_point() * sizeof(Myfloat)) ;
	out_file.close() ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_3D_float::write_on_disk");
	return(RTN_CODE_OK) ;
}


//-------------------------------------------------------------------------------------------------------

void Grid_3D_float::info(void)
{
	Grid_3D::info() ;
	print_info(ALL, "min", get_min() ) ;
	print_info(ALL, "max", get_max() ) ;
}    

} // namespace django
