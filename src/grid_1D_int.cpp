//-------------------------------------------------------------------------------------------------------
//
// GRID 1D INT
//
//-------------------------------------------------------------------------------------------------------

#include "grid_1D_int.h"

#include <algorithm>
#include <climits>
#include <fstream>

#include "allocate_array.h"
#include "type_def.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------

Grid_1D_int::~Grid_1D_int(void)
{  
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_1D_int::~Grid_1D_int");
	deallocate_array<Myint>(pArray, nz) ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_1D_int::~Grid_1D_int");
}

//-------------------------------------------------------------------------------------------------------

Grid_1D_int::Grid_1D_int(Myint nz2, Myfloat dz2)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_1D_int::Grid_1D_int");
	nz = nz2 ;
	dz = dz2 ;
	pArray = allocate_array<Myint>(nz) ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_1D_int::Grid_1D_int");
}

//-------------------------------------------------------------------------------------------------------

Myint Grid_1D_int::get_min()
{

	Myint min_val = INT_MAX ;

	for (Myint iz=0; iz<nz; iz++)
	{
		if (pArray[iz] < min_val)
		{
			min_val = min(min_val, pArray[iz]) ;
		}
	}
	return min_val ;
}

//-------------------------------------------------------------------------------------------------------

Myint Grid_1D_int::get_max()
{

	Myint max_val = INT_MIN ;

	for (Myint iz=0; iz<nz; iz++)
	{
		if (pArray[iz] > max_val)
		{
			max_val = max(max_val, pArray[iz]) ;
		}
	}
	return max_val ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_1D_int::read_from_disk(string file_name)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_1D_int::read_from_disk");

	ifstream in_file ;
	in_file.open(file_name.c_str(), ios::binary) ;
	assert(in_file.is_open());

	Myint32 pArray_tmp[get_nb_grid_point()] ;
	in_file.read((char*) pArray_tmp, get_nb_grid_point() * sizeof(Myint32)) ;
	in_file.close() ;

#pragma ivdep
	for (Myint ii=0; ii < get_nb_grid_point(); ii++)
	{
		pArray[ii] = pArray_tmp[ii] ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_1D_int::read_from_disk");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_1D_int::write_on_disk(string file_name)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_1D_int::write_on_disk");
	ofstream out_file ;

	Myint64 size = get_nb_grid_point() * sizeof(Myint32) ;
	out_file.open(file_name.c_str(), ios::binary | ios::app) ;
	assert(out_file.is_open());
	out_file.write((char*) pArray, size) ;
	out_file.close() ;

	print_debug(ALL, MID_DEBUG, " snapshot file", file_name) ;
	print_debug(ALL, MID_DEBUG, " snapshot size", size) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_1D_int::write_on_disk");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_1D_int::reset()
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_1D_int::reset()");
#pragma ivdep
	for (Myint iz=0; iz<nz; iz++)
	{
		pArray[iz] = 0. ;
	}
	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_1D_int::reset()");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_1D_int::reset(Myfloat val)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_1D_int::reset(Myint val)");
#pragma ivdep
	for (Myint iz=0; iz<nz; iz++)
	{
		pArray[iz] = val ;
	}
	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_1D_int::reset(Myint val)");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

void Grid_1D_int::info(void)
{
	Grid_1D::info() ;
	print_info(ALL, "min", get_min() ) ;
	print_info(ALL, "max", get_max() ) ;
}

} // namespace django
