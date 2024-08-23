//-------------------------------------------------------------------------------------------------------
//
// GRID 3D
//
//-------------------------------------------------------------------------------------------------------

#include "grid_3D.h"

#include "output_report.h"
#include "type_def.h"

using namespace std;

namespace django {

Myint64 Grid_3D::get_nb_grid_point()
{
	return(nz*nx*ny) ;
}

//-------------------------------------------------------------------------------------------------------

void Grid_3D::info()
{
	Grid_2D::info() ;
	print_info(ALL, "ny", ny);
	print_info(ALL, "dy", dy);
}

} // namespace django

