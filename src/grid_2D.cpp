//-------------------------------------------------------------------------------------------------------
//
// GRID 2D
//
//-------------------------------------------------------------------------------------------------------

#include "grid_2D.h"

#include "output_report.h"
#include "type_def.h"

using namespace std;

namespace django {

Myint64 Grid_2D::get_nb_grid_point()
{
	return(nz*nx) ;
}


//-------------------------------------------------------------------------------------------------------

void Grid_2D::info()
{
	Grid_1D::info() ;
	print_info(ALL, "nx", nx);
	print_info(ALL, "dx", dx);
}

} // namespace django
