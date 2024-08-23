//-------------------------------------------------------------------------------------------------------
//
// PARENT CLASS FOR ALL SNAPSHOTS
//
//-------------------------------------------------------------------------------------------------------

#include "snapshot.h"

#include "constant.h"
#include "output_report.h"
#include "type_def.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------

Snapshot::Snapshot(void)
{  
	print_debug(ALL, LIGHT_DEBUG, "IN Snapshot::Snapshot");
	x_coord = NULL ;
	z_coord = NULL ;
	tmax    = 0.0 ;
	tmin    = 0.0 ;
	dt      = 0.0 ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Snapshot::Snapshot");
}

//-------------------------------------------------------------------------------------------------------

Snapshot::~Snapshot(void)
{  
	print_debug(ALL, LIGHT_DEBUG, "IN Snapshot::~Snapshot");
	delete(x_coord) ;
	delete(z_coord) ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Snapshot::~Snapshot");
}

//-------------------------------------------------------------------------------------------------------
Rtn_code Snapshot::set_time_param(Myfloat tmin_in, Myfloat tmax_in, Myfloat dt_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Snapshot::set_time_param");
	tmin = tmin_in ;
	tmax = tmax_in ;
	dt   = dt_in ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Snapshot::set_time_param");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code Snapshot::info(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Snapshot::info");

	print_info(MASTER, " Snapshot tmin\t", tmin) ;
	print_info(MASTER, " Snapshot tmax\t", tmax) ;
	print_info(MASTER, " Snapshot dt\t", dt) ;
	Myint nx = 1 ; // TO DO there should be a test on dim = 1, else nx = 0
	if (x_coord != NULL)
	{
		print_info(MASTER, " Snapshot xmin\t", x_coord->get_min()) ;
		print_info(MASTER, " Snapshot xmax\t", x_coord->get_max()) ;
		nx = x_coord->nz ;
		print_info(MASTER, " Snapshot nx\t", nx) ;
	}

	Myint nz = 0 ;
	if (z_coord != NULL)
	{
		print_info(MASTER, " Snapshot zmin\t", z_coord->get_min()) ;
		print_info(MASTER, " Snapshot zmax\t", z_coord->get_max()) ;
		nz = z_coord->nz ;
		print_info(MASTER, " Snapshot nz\t", nz) ;
	}

	Myint npixel = nx * nz ;
	if (npixel == 0)
	{
		print_warning(" Snapshot has no pixel") ;
	}


	print_info(MASTER, " Snapshot size (KB)", (Myfloat)(npixel * sizeof(Myfloat)/1000.)) ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Snapshot::info");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code Snapshot::set_pixel_param(Myfloat xmin, Myfloat xmax, Myint nx,
		Myfloat zmin, Myfloat zmax, Myint nz)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Snapshot::set_pixel_param");

	if (nx != NOT_SPECIFIED) x_coord = new Grid_1D_float(nx, xmin, xmax) ;
	if (nz != NOT_SPECIFIED) z_coord = new Grid_1D_float(nz, zmin, zmax) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Snapshot::set_pixel_param");
	return(RTN_CODE_OK) ;
}

} // namespace django

