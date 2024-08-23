//-------------------------------------------------------------------------------------------------------
//
// GRID 1D COMPLEX
//
//-------------------------------------------------------------------------------------------------------

#include "grid_1D_complex.h"

#include <cfloat>
#include <fstream>

#include "allocate_array.h"
#include "type_def.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------

Grid_1D_complex::Grid_1D_complex(Myint nz2, Myfloat dz2)
{
	nz = nz2 ;
	dz = dz2 ;
	pArray = allocate_array<Mycomplex>(nz) ;
}


//-------------------------------------------------------------------------------------------------------

Myfloat Grid_1D_complex::get_min()
{

	Myfloat min_val = FLT_MAX ;

	for (Myint iz=0; iz<nz; iz++)
	{
		if (abs(pArray[iz]) < min_val)
		{
			min_val = abs(pArray[iz]) ;
		}
	}
	return min_val ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid_1D_complex::get_max()
{

	Myfloat max_val = -FLT_MAX ;

	for (Myint iz=0; iz<nz; iz++)
	{
		if (abs(pArray[iz]) > max_val)
		{
			max_val = abs(pArray[iz]) ;
		}
	}
	return max_val ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_1D_complex::write_on_disk(string file_name)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Grid_1D_complex::write_on_disk");
	ofstream out_file ;

	out_file.open(file_name.c_str(), ios::binary) ;
	assert(out_file.is_open());
	out_file.write((char*) pArray, get_nb_grid_point() * sizeof(Mycomplex)) ;
	out_file.close() ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Grid_1D_complex::write_on_disk");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

void Grid_1D_complex::info(void)
{
	Grid_1D::info() ;
	print_info(ALL, "min", get_min() ) ;
	print_info(ALL, "max", get_max() ) ;
}

} // namespace django
