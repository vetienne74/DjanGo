//-------------------------------------------------------------------------------------------------------
//
// GRID 2D INT
//
//-------------------------------------------------------------------------------------------------------

#include "grid_2D_int.h"

#include <climits>
#include <fstream>

#include "allocate_array.h"
#include "type_def.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------

Grid_2D_int::~Grid_2D_int(void)
{  
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_2D_int::~Grid_2D_int");
	deallocate_array<Myint>(pArray, nx, nz) ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_2D_int::~Grid_2D_int");
}

//-------------------------------------------------------------------------------------------------------

Grid_2D_int::Grid_2D_int(Myint nz2, Myint nx2, Myfloat dz2, Myfloat dx2)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_2D_int::Grid_2D_int");
	nz = nz2 ;
	nx = nx2 ;
	dz = dz2 ;
	dx = dx2 ;
	pArray = allocate_array<Myint>(nx, nz) ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_2D_int::Grid_2D_int");
}

//-------------------------------------------------------------------------------------------------------

Myint Grid_2D_int::get_min()
{

	Myint min_val = INT_MAX ;

	for (Myint ix=0; ix<nx; ix++)
	{
		for (Myint iz=0; iz<nz; iz++)
		{
			if (pArray[ix][iz] < min_val)
			{
				min_val = pArray[ix][iz] ;
			}
		}
	}
	return min_val ;
}

//-------------------------------------------------------------------------------------------------------

Myint Grid_2D_int::get_max()
{

	Myint max_val = -INT_MAX ;

	for (Myint ix=0; ix<nx; ix++)
	{
		for (Myint iz=0; iz<nz; iz++)
		{
			if (pArray[ix][iz] > max_val)
			{
				max_val = pArray[ix][iz] ;
			}
		}
	}
	return max_val ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_2D_int::read_from_disk(string file_name)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_2D_int::read_from_disk");

	ifstream in_file ;
	in_file.open(file_name.c_str(), ios::binary) ;
	assert(in_file.is_open());

	Myint32 pArray_tmp[get_nb_grid_point()] ;
	in_file.read((char*) pArray_tmp, get_nb_grid_point() * sizeof(Myint32)) ;
	in_file.close() ;

	Myint* pArray_tmp2 = (Myint*) (*pArray) ;
#pragma ivdep
	for (Myint ii=0; ii < get_nb_grid_point(); ii++)
	{
		pArray_tmp2[ii] = pArray_tmp[ii] ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_2D_int::read_from_disk");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_2D_int::write_on_disk(string file_name)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_2D_int::write_on_disk");
	print_debug(ALL, LIGHT_DEBUG, "write file", file_name.c_str());
	ofstream out_file ;

	out_file.open(file_name.c_str(), ios::binary | ios::app) ;
	assert(out_file.is_open());
	out_file.write((char*) *pArray, get_nb_grid_point() * sizeof(Myint)) ;
	out_file.close() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_2D_int::write_on_disk");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_2D_int::reset()
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_2D_int::reset()");
	for (Myint ix=0; ix<nx; ix++)
	{
#pragma ivdep
		for (Myint iz=0; iz<nz; iz++)
		{
			pArray[ix][iz] = 0. ;
		}
	}
	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_2D_int::reset()");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_2D_int::reset(Myfloat val)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_2D_int::reset()");
	for (Myint ix=0; ix<nx; ix++)
	{
#pragma ivdep
		for (Myint iz=0; iz<nz; iz++)
		{
			pArray[ix][iz] = val ;
		}
	}
	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_2D_int::reset()");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

void Grid_2D_int::info(void)
{
	Grid_2D::info() ;
	print_info(ALL, "min", get_min() ) ;
	print_info(ALL, "max", get_max() ) ;
}

} // namespace django

