//-------------------------------------------------------------------------------------------------------
//
// CLASS SCHEME_FACTORY
//
//-------------------------------------------------------------------------------------------------------

#include "scheme_factory.h"

#include <cstdlib>

#include "allocate_array.h"
#include "FDM_1D_ac_iso_1st_schemes.h"
#include "FDM_1D_ac_iso_2nd_schemes.h"
#include "FDM_1D_ac_lossy_1st_schemes.h"
#include "FDM_1D_ac_lossy_2nd_schemes.h"

#include "FDM_2D_ac_iso_1st_schemes.h"
#include "FDM_2D_ac_iso_2nd_schemes.h"
#include "FDM_2D_ac_lossy_1st_schemes.h"
#include "FDM_2D_ac_lossy_2nd_schemes.h"
#include "FDM_2D_el_iso_1st_schemes.h"

#include "FEM_1D_1st_ac_iso.h"
#include "FEM_1D_1st_ac_lossy.h"
#include "FEM_1D_2nd_ac_iso.h"

#include "FEM_2D_1st_ac_iso.h"
#include "FEM_2D_1st_ac_lossy.h"
#include "FEM_2D_2nd_ac_iso.h"
#include "FEM_2D_1st_el_iso.h"

#include "output_report.h"
#include "type_def.h"
#include "singleton.h"

namespace django {

Scheme_Factory::Scheme_Factory(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme_Factory::Scheme_Factory");

	pModel    = NULL ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme_Factory::Scheme_Factory");
}

//-------------------------------------------------------------------------------------------------------

Scheme* Scheme_Factory::create()
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme_Factory::create");

	Scheme* new_scheme = NULL ;
	if (pModel == NULL)
	{
		print_error(" pModel is NULL, Can not create scheme in Scheme_Factory::create") ;
	}
	else
	{
		// set dim
		Space_dim dim = pModel->get_dim() ;

		Eq_type eq_type      = Singleton::Instance()->eq_type ;
		Eq_order eq_order    = Singleton::Instance()->eq_order ;
		Myint space_order    = Singleton::Instance()->space_order ;
		Myint time_order     = Singleton::Instance()->time_order ;
		Scheme_method method = Singleton::Instance()->method ;
		Scheme_type type     = Singleton::Instance()->type ;

		if (method == SCHEME_FDM)
		{
			if (type == SCHEME_STAGGERED)
			{
				if (dim == ONE)
				{
					if (eq_type == ACOUSTIC)
					{
						if (eq_order == ORDER_1ST)
						{
							if (time_order == 2)
							{
								if (space_order == 2)  new_scheme = new FDM_1D_ac_iso_1st_2_2() ;
								if (space_order == 4)  new_scheme = new FDM_1D_ac_iso_1st_4_2() ;
								if (space_order == 8)  new_scheme = new FDM_1D_ac_iso_1st_8_2() ;
								if (space_order == 12) new_scheme = new FDM_1D_ac_iso_1st_12_2() ;
								if (space_order == 16) new_scheme = new FDM_1D_ac_iso_1st_16_2() ;
							}
						}
						else if (eq_order == ORDER_2ND)
						{
							if (time_order == 2)
							{
								if (space_order == 2)  new_scheme = new FDM_1D_ac_iso_2nd_2_2() ;
								if (space_order == 4)  new_scheme = new FDM_1D_ac_iso_2nd_4_2() ;
								if (space_order == 8)  new_scheme = new FDM_1D_ac_iso_2nd_8_2() ;
								if (space_order == 12) new_scheme = new FDM_1D_ac_iso_2nd_12_2() ;
								if (space_order == 16) new_scheme = new FDM_1D_ac_iso_2nd_16_2() ;
							}
						}
					}
					else if (eq_type == AC_LOSSY)
					{
						if (eq_order == ORDER_1ST)
						{
							if (time_order == 2)
							{
								if (space_order == 2)  new_scheme = new FDM_1D_ac_lossy_1st_2_2() ;
								if (space_order == 4)  new_scheme = new FDM_1D_ac_lossy_1st_4_2() ;
								if (space_order == 8)  new_scheme = new FDM_1D_ac_lossy_1st_8_2() ;
								if (space_order == 12) new_scheme = new FDM_1D_ac_lossy_1st_12_2() ;
								if (space_order == 16) new_scheme = new FDM_1D_ac_lossy_1st_16_2() ;
							}
						}
						else if (eq_order == ORDER_2ND)
						{
							if (time_order == 2)
							{
								if (space_order == 2)  new_scheme = new FDM_1D_ac_lossy_2nd_2_2() ;
								if (space_order == 4)  new_scheme = new FDM_1D_ac_lossy_2nd_4_2() ;
								if (space_order == 8)  new_scheme = new FDM_1D_ac_lossy_2nd_8_2() ;
								if (space_order == 12) new_scheme = new FDM_1D_ac_lossy_2nd_12_2() ;
								if (space_order == 16) new_scheme = new FDM_1D_ac_lossy_2nd_16_2() ;
							}
						}
					}
				}
				else if (dim == TWO)
				{
					if (eq_type == ACOUSTIC)
					{
						if (eq_order == ORDER_1ST)
						{
							if (time_order == 2)
							{
								if (space_order == 2)  new_scheme = new FDM_2D_ac_iso_1st_2_2() ;
								if (space_order == 4)  new_scheme = new FDM_2D_ac_iso_1st_4_2() ;
								if (space_order == 8)  new_scheme = new FDM_2D_ac_iso_1st_8_2() ;
								if (space_order == 12) new_scheme = new FDM_2D_ac_iso_1st_12_2() ;
								if (space_order == 16) new_scheme = new FDM_2D_ac_iso_1st_16_2() ;
							}
						}
						else if (eq_order == ORDER_2ND)
						{
							if (time_order == 2)
							{
								if (space_order == 2)  new_scheme = new FDM_2D_ac_iso_2nd_2_2() ;
								if (space_order == 4)  new_scheme = new FDM_2D_ac_iso_2nd_4_2() ;
								if (space_order == 8)  new_scheme = new FDM_2D_ac_iso_2nd_8_2() ;
								if (space_order == 12) new_scheme = new FDM_2D_ac_iso_2nd_12_2() ;
								if (space_order == 16) new_scheme = new FDM_2D_ac_iso_2nd_16_2() ;
							}
						}
					}
					else if (eq_type == AC_LOSSY)
					{
						if (eq_order == ORDER_1ST)
						{
							if (time_order == 2)
							{
								if (space_order == 2)  new_scheme = new FDM_2D_ac_lossy_1st_2_2() ;
								if (space_order == 4)  new_scheme = new FDM_2D_ac_lossy_1st_4_2() ;
								if (space_order == 8)  new_scheme = new FDM_2D_ac_lossy_1st_8_2() ;
								if (space_order == 12) new_scheme = new FDM_2D_ac_lossy_1st_12_2() ;
								if (space_order == 16) new_scheme = new FDM_2D_ac_lossy_1st_16_2() ;
							}
						}
						else if (eq_order == ORDER_2ND)
						{
							if (time_order == 2)
							{
								if (space_order == 2)  new_scheme = new FDM_2D_ac_lossy_2nd_2_2() ;
								if (space_order == 4)  new_scheme = new FDM_2D_ac_lossy_2nd_4_2() ;
								if (space_order == 8)  new_scheme = new FDM_2D_ac_lossy_2nd_8_2() ;
								if (space_order == 12) new_scheme = new FDM_2D_ac_lossy_2nd_12_2() ;
								if (space_order == 16) new_scheme = new FDM_2D_ac_lossy_2nd_16_2() ;
							}
						}
					}
					else if (eq_type == ELASTIC)
					{
						if (eq_order == ORDER_1ST)
						{
							if (time_order == 2)
							{
								if (space_order == 2)  new_scheme = new FDM_2D_el_iso_1st_2_2() ;
								if (space_order == 4)  new_scheme = new FDM_2D_el_iso_1st_4_2() ;
								if (space_order == 8)  new_scheme = new FDM_2D_el_iso_1st_8_2() ;
								if (space_order == 12) new_scheme = new FDM_2D_el_iso_1st_12_2() ;
								if (space_order == 16) new_scheme = new FDM_2D_el_iso_1st_16_2() ;
							}
						}
						else if (eq_order == ORDER_2ND)
						{
							//new_scheme = new FDM_2D_el_2nd() ;
						}
					}
				}
				else if (dim == THREE)
				{
					//new_scheme = new FDM_3D() ;
				}
			}
		}

		else if (method == SCHEME_FEM)
		{
			if (dim == ONE)
			{
				if (eq_type == ACOUSTIC)
				{
					if (eq_order == ORDER_1ST)
					{
						new_scheme = new FEM_1D_1st_ac_iso() ;
					}
					else if (eq_order == ORDER_2ND)
					{
						new_scheme = new FEM_1D_2nd_ac_iso() ;
					}
				}
				else if (eq_type == AC_LOSSY)
				{
					if (eq_order == ORDER_1ST)
					{
						new_scheme = new FEM_1D_1st_ac_lossy() ;
					}
				}
			}
			else if (dim == TWO)
			{
				if (eq_type == ACOUSTIC)
				{
					if (eq_order == ORDER_1ST)
					{
						new_scheme = new FEM_2D_1st_ac_iso() ;
					}
					else if (eq_order == ORDER_2ND)
					{
						new_scheme = new FEM_2D_2nd_ac_iso() ;
					}
				}
				else if (eq_type == AC_LOSSY)
				{
					if (eq_order == ORDER_1ST)
					{
						new_scheme = new FEM_2D_1st_ac_lossy() ;
					}
				}
				else if (eq_type == ELASTIC)
				{
					if (eq_order == ORDER_1ST)
					{
						new_scheme = new FEM_2D_1st_el_iso() ;
					}
				}
			}
		}
	}

	// check scheme has been instanciated
	if (new_scheme == NULL)
	{
		print_error(" Can not create scheme in Scheme_Factory::create") ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme_Factory::create");
	return(new_scheme) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Scheme_Factory::set_model(Model* pModel_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme_Factory::set_model");
	if (pModel_in == NULL)
	{
		print_error(" Invalid Model in Scheme_Factory::set_model") ;
		return(RTN_CODE_KO) ;
	}
	else
	{
		pModel = pModel_in ;
	}
	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme_Factory::set_model");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Scheme_Factory::info(void)
{
	print_info(MASTER, "") ;
	print_line2() ;
	print_info(MASTER, " SCHEME_FACTORY PARAMETERS") ;
	print_info(MASTER, "") ;

	// method
	switch(Singleton::Instance()->method)
	{
	case SCHEME_FDM:
		print_info(MASTER, " Method\t\t", "FDM") ;
		break ;
	case SCHEME_FEM:
		print_info(MASTER, " Method\t\t", "FEM") ;
		break ;
	default:
		print_error(" Invalid numerical method", Singleton::Instance()->method) ;
		return(RTN_CODE_KO) ;
	}

	// type
	switch(Singleton::Instance()->type)
	{
	case SCHEME_STAGGERED:
		print_info(MASTER, " Type\t\t", "STAGGERED") ;
		break ;
	case SCHEME_CGM:
		print_info(MASTER, " Type\t\t", "CGM") ;
		break ;
	case SCHEME_DGM:
		print_info(MASTER, " Type\t\t", "DGM") ;
		break ;
	case SCHEME_MGM:
		print_info(MASTER, " Type\t\t", "MGM") ;
		break ;
	default:
		print_error(" Invalid numerical type", Singleton::Instance()->type) ;
		return(RTN_CODE_KO) ;
	}

	// dim
	switch(Singleton::Instance()->dim)
	{
	case ONE:
		print_info(MASTER, " Dimension\t", "1D") ;
		break ;
	case TWO:
		print_info(MASTER, " Dimension\t", "2D") ;
		break ;
	case THREE:
		print_info(MASTER, " Dimension\t", "3D") ;
		break ;
	default:
		print_error(" Invalid dimension", Singleton::Instance()->dim) ;
		return(RTN_CODE_KO) ;
	}

	// equation type
	switch(Singleton::Instance()->eq_type)
	{
	case ACOUSTIC:
		print_info(MASTER, " Equation\t", "ACOUSTIC") ;
		break ;
	case ELASTIC:
		print_info(MASTER, " Equation\t", "ELASTIC") ;
		break ;
	case AC_LOSSY:
		print_info(MASTER, " Equation\t", "ACOUSTIC LOSSY") ;
		break ;
	default:
		print_error(" Invalid equation type", Singleton::Instance()->eq_type) ;
		return(RTN_CODE_KO) ;
	}

	// equation order
	switch(Singleton::Instance()->eq_order)
	{
	case ORDER_1ST:
		print_info(MASTER, " Equation\t", "1ST ORDER") ;
		break ;
	case ORDER_2ND:
		print_info(MASTER, " Equation\t", "2ND ORDER") ;
		break ;
	default:
		print_error(" Invalid equation order", Singleton::Instance()->eq_order) ;
		return(RTN_CODE_KO) ;
	}

	// space order
	Myint space_order = Singleton::Instance()->space_order ;
	if (space_order < 0)
	{
		print_error(" Invalid space order", space_order) ;
		return(RTN_CODE_KO) ;
	}
	else
	{
		print_info(MASTER, " Spatial order\t", space_order) ;
	}

	// time order
	Myint time_order = Singleton::Instance()->time_order ;
	if (time_order < 0)
	{
		print_error(" Invalid time order", time_order) ;
		return(RTN_CODE_KO) ;
	}
	else
	{
		print_info(MASTER, " Time order\t", time_order) ;
	}

	return(RTN_CODE_OK) ;
}

} // namespace django
