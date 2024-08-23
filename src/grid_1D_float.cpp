//-------------------------------------------------------------------------------------------------------
//
// GRID 1D FLOAT
//
//-------------------------------------------------------------------------------------------------------

#include "grid_1D_float.h"

#include <cfloat>
#include <fstream>

#include "allocate_array.h"
#include "type_def.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------

Grid_1D_float::~Grid_1D_float(void)
{  
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_1D_float::~Grid_1D_float");
	deallocate_array<Myfloat>(pArray, nz) ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_1D_float::~Grid_1D_float");
}

//-------------------------------------------------------------------------------------------------------

Grid_1D_float::Grid_1D_float(Myint nz2, Myfloat dz2)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_1D_float::Grid_1D_float");
	nz = nz2 ;
	dz = dz2 ;
	pArray = allocate_array<Myfloat>(nz) ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_1D_float::Grid_1D_float");
}

//-------------------------------------------------------------------------------------------------------

Grid_1D_float::Grid_1D_float(Myint nz2, Myfloat start, Myfloat end)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_1D_float::Grid_1D_float");
	nz = nz2 ;
	dz = (end - start) / (nz2-1) ;
	pArray = allocate_array<Myfloat>(nz) ;

	for (Myint iz=0; iz<nz; iz++)
	{
		pArray[iz] = start +iz * dz ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_1D_float::Grid_1D_float");
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid_1D_float::get_min()
{

	Myfloat min_val = FLT_MAX ;

	for (Myint iz=0; iz<nz; iz++)
	{
		if (pArray[iz] < min_val)
		{
			min_val = pArray[iz] ;
		}
	}
	return min_val ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid_1D_float::get_max()
{

	Myfloat max_val = -FLT_MAX ;

	for (Myint iz=0; iz<nz; iz++)
	{
		if (pArray[iz] > max_val)
		{
			max_val = pArray[iz] ;
		}
	}
	return max_val ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_1D_float::read_from_disk(string file_name)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_1D_float::read_from_disk");

	ifstream in_file ;
	in_file.open(file_name.c_str(), ios::binary) ;
	assert(in_file.is_open());

	Myfloat32 pArray_tmp[get_nb_grid_point()] ;
	in_file.read((char*) pArray_tmp, get_nb_grid_point() * sizeof(Myfloat32)) ;
	in_file.close() ;

#pragma ivdep
	for (Myint ii=0; ii < get_nb_grid_point(); ii++)
	{
		pArray[ii] = pArray_tmp[ii] ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_1D_float::read_from_disk");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_1D_float::write_on_disk(string file_name)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_1D_float::write_on_disk");
	ofstream out_file ;

	Myint64 size = get_nb_grid_point() * sizeof(Myfloat) ;
	out_file.open(file_name.c_str(), ios::binary | ios::app) ;
	assert(out_file.is_open());
	out_file.write((char*) pArray, size) ;
	out_file.close() ;

	print_debug(ALL, MID_DEBUG, " snapshot file", file_name) ;
	print_debug(ALL, MID_DEBUG, " snapshot size", size) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_1D_float::write_on_disk");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_1D_float::reset()
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_1D_float::reset()");
#pragma ivdep
	for (Myint iz=0; iz<nz; iz++)
	{
		pArray[iz] = 0. ;
	}
	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_1D_float::reset()");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_1D_float::reset(Myfloat val)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_1D_float::reset(Myfloat val)");
#pragma ivdep
	for (Myint iz=0; iz<nz; iz++)
	{
		pArray[iz] = val ;
	}
	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_1D_float::reset(Myfloat val)");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

void Grid_1D_float::info(void)
{
	Grid_1D::info() ;
	print_info(ALL, "min", get_min() ) ;
	print_info(ALL, "max", get_max() ) ;
}

} // namespace django
