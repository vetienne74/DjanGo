//-------------------------------------------------------------------------------------------------------
//
// MANAGEMENT OF ACQUISITION
//
// Input acquisition file
// Order of coordinates:
// 1D config: Z
// 2D config: Z, X
// 3D config: Z, X, Y
//
// Output VTK files
// Order of coordinates: X, Y, Z (for 1D, 2D, and 3D)
//
//-------------------------------------------------------------------------------------------------------

#include "acquisition.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <string>

#include "allocate_array.h"
#include "constant.h"
#include "output_report.h"
#include "type_def.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------

Acquisition::~Acquisition(void) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Acquisition::~Acquisition");

	switch(acquisition_dim)
	{
	case ONE:
	{
		deallocate_array<Myfloat>(zrec, nrec) ;
		break ;
	}
	case TWO:
	{
		deallocate_array<Myfloat>(zrec, nrec) ;
		deallocate_array<Myfloat>(xrec, nrec) ;
		break ;
	}
	case THREE:
	{
		deallocate_array<Myfloat>(zrec, nrec) ;
		deallocate_array<Myfloat>(xrec, nrec) ;
		deallocate_array<Myfloat>(yrec, nrec) ;
		break ;
	}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Acquisition::~Acquisition");
}

//-------------------------------------------------------------------------------------------------------

Acquisition::Acquisition(string file_name, Space_dim dim) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Acquisition::Acquisition");
	acquisition_file_name = file_name ;
	acquisition_dim       = dim ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Acquisition::Acquisition");
}

//-------------------------------------------------------------------------------------------------------


Rtn_code Acquisition::initialize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Acquisition::initialize");

	Myint   int_tmp ;
	Myfloat float_tmp ;

	// open acquisition file

	print_info(MASTER, " acquisition file:", acquisition_file_name.c_str()) ;
	ifstream acquisition_file (acquisition_file_name.c_str()) ;
	assert(acquisition_file.is_open()) ;

	//-------------------------------------------------------------------------------------------------------
	// read header of acquisition file
	//-------------------------------------------------------------------------------------------------------

	// get file type

	acquisition_file >> int_tmp ;

	if (int_tmp == SEISMIC_FIXED)
	{
		acquisition_type = SEISMIC_FIXED ;
		print_info(MASTER, " acquisition type: \tSEISMIC FIXED") ;
	}
	else if (int_tmp == SEISMIC_STREAMER)
	{
		acquisition_type = SEISMIC_STREAMER ;
		print_info(MASTER, " acquisition type: \tSEISMIC STREAMER") ;
		print_error(" acquisition file type is not supported") ;
		return(RTN_CODE_KO) ;
	}
	else
	{
		print_error(" acquisition file type is not supported") ;
		return(RTN_CODE_KO) ;
	}

	// space dimension

	acquisition_file >> int_tmp ;

	if (int_tmp != acquisition_dim)
	{
		print_error(" dim in acqui file is incompatible with dim in config", int_tmp, acquisition_dim) ;
		return(RTN_CODE_KO) ;
	}

	// get nb of sources

	acquisition_file >> nsrc ;

	if (nsrc >= 1)
	{
		print_info(MASTER, " nb sources: \t", nsrc) ;
	}
	else
	{
		print_error(" nb of shots should be >= 1") ;
		return(RTN_CODE_KO) ;
	}

	//-------------------------------------------------------------------------------------------------------
	// specific initialisations for FIXED geometry
	//-------------------------------------------------------------------------------------------------------

	if (acquisition_type == SEISMIC_FIXED)

	{

		// loop on the sources
		// temporary: only last src coordinates is stored
		// the acquisition file is re-read at each source computation
		switch(acquisition_dim)
		{
		case ONE:
		{
			for (Myint isrc=1; isrc<=nsrc; isrc++)
			{
				acquisition_file >> int_tmp >> zsrc ;
			}
			break ;
		}
		case TWO:
		{
			for (Myint isrc=1; isrc<=nsrc; isrc++)
			{
				acquisition_file >> int_tmp >> zsrc >> xsrc ;
			}
			break ;
		}
		case THREE:
		{
			for (Myint isrc=1; isrc<=nsrc; isrc++)
			{
				acquisition_file >> int_tmp >> zsrc >> xsrc >> ysrc ;
			}
			break ;
		}
		}

		// get nb of receivers
		acquisition_file >> nrec ;
		print_info(MASTER, " nb receivers: \t", nrec) ;

		// read coordinates of the current source
		switch(acquisition_dim)
		{
		case ONE:
		{
			zrec = allocate_array<Myfloat>(nrec) ;
			for (Myint irec=0; irec<nrec; irec++)
			{
				acquisition_file >> int_tmp >> zrec[irec] ;
			}
			break ;
		}
		case TWO:
		{
			zrec = allocate_array<Myfloat>(nrec) ;
			xrec = allocate_array<Myfloat>(nrec) ;
			for (Myint irec=0; irec<nrec; irec++)
			{
				acquisition_file >> int_tmp >> zrec[irec] >> xrec[irec] ;
			}
			break ;
		}
		case THREE:
		{
			zrec = allocate_array<Myfloat>(nrec) ;
			xrec = allocate_array<Myfloat>(nrec) ;
			yrec = allocate_array<Myfloat>(nrec) ;
			for (Myint irec=0; irec<nrec; irec++)
			{
				acquisition_file >> int_tmp >> zrec[irec] >> xrec[irec] >> yrec[irec] ;
			}
			break ;
		}
		}
	}

	// intialise current shot id
	current_src = 1 ;

	// write acqui in VTK format
	write_src_VTK() ;
	write_rec_VTK() ;

	print_info(MASTER, "") ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Acquisition::initialize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Acquisition::get_position_for_current_shot(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Acquisition::get_position_for_current_shot");
	Myint int_tmp ;

	if ((current_src <= 0) || (current_src > nsrc))
	{
		print_error(" get_position_for_current_shot with invalid current_src") ;
		return(RTN_CODE_KO) ;
	}

	// open acquisition file
	ifstream acquisition_file (acquisition_file_name.c_str()) ;
	assert(acquisition_file.is_open()) ;

	// get file type
	acquisition_file >> int_tmp ;

	// space dimension
	acquisition_file >> int_tmp ;

	// get nb of sources
	acquisition_file >> int_tmp ;

	// read coordinates of the current source
	switch(acquisition_dim)
	{
	case ONE:
	{
		for (Myint isrc=1; isrc<=nsrc; isrc++)
		{
			acquisition_file >> int_tmp >> zsrc ;
			if (isrc == current_src)
			{
				break ;
			}
		}
		break ;
	}
	case TWO:
	{
		for (Myint isrc=1; isrc<=nsrc; isrc++)
		{
			acquisition_file >> int_tmp >> zsrc >> xsrc ;
			if (isrc == current_src)
			{
				break ;
			}
		}
		break ;
	}
	case THREE:
	{
		for (Myint isrc=1; isrc<=nsrc; isrc++)
		{
			acquisition_file >> int_tmp >> zsrc >> xsrc >> ysrc ;
			if (isrc == current_src)
			{
				break ;
			}
		}
		break ;
	}
	}

	//cout << "current_src " << current_src << " xsrc " << xsrc << "\n" ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Acquisition::get_position_for_current_shot");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

bool Acquisition::remain_shot_to_compute(void)
{
	if (current_src <= nsrc)
	{
		return true ;
	}
	else
	{
		return false ;
	}
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Acquisition::move_to_next_shot(void)
{
	current_src++ ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Acquisition::move_to_first_shot(void)
{
	current_src = 1 ;
	return(RTN_CODE_OK) ;
}

//-----------------------------------------------------------------------------------------
Rtn_code Acquisition::write_rec_VTK(void)
{   
	print_debug(ALL, LIGHT_DEBUG, "IN Acquisition::write_rec_VTK");

	ofstream out_file ;
	out_file.open(REC_VTK_OUT_FILE, ios::trunc | ios::out) ;
	out_file << "# vtk DataFile Version 2.0\n" ;
	out_file << "MESH FILE DJANGO\n" ;
	out_file << "ASCII\n" ;
	out_file << "DATASET POLYDATA\n" ;
	out_file << "POINTS " << nrec << " float\n" ;

	switch(acquisition_dim)
	{
	case ONE:
	{
		for (Myint irec = 0; irec < nrec; irec++)
		{
			out_file << "0.0 0.0 " << zrec[irec] << "\n" ;
		}
		break ;
	}
	case TWO:
	{
		for (Myint irec = 0; irec < nrec; irec++)
		{
			out_file << xrec[irec] << " 0.0 " << zrec[irec] << "\n" ;
		}
		break ;
	}
	case THREE:
	{
		for (Myint irec = 0; irec < nrec; irec++)
		{
			out_file << xrec[irec] << " " << yrec[irec] << " " << zrec[irec] << "\n" ;
		}

		break ;
	}
	}

	out_file.close() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Acquisition::write_rec_VTK");
	return(RTN_CODE_OK) ;
}

//-----------------------------------------------------------------------------------------
Rtn_code Acquisition::write_src_VTK(void)
{   
	print_debug(ALL, LIGHT_DEBUG, "IN Acquisition::write_src_VTK");

	ofstream out_file ;
	out_file.open(SRC_VTK_OUT_FILE, ios::trunc | ios::out) ;
	out_file << "# vtk DataFile Version 2.0\n" ;
	out_file << "MESH FILE DJANGO\n" ;
	out_file << "ASCII\n" ;
	out_file << "DATASET POLYDATA\n" ;
	out_file << "POINTS " << nsrc << " float\n" ;

	switch(acquisition_dim)
	{
	case ONE:
	{
		out_file << "0.0 0.0 " << zsrc << "\n" ;
		break ;
	}
	case TWO:
	{
		out_file << xsrc << " 0.0 " << zsrc << "\n" ;
		break ;
	}
	case THREE:
	{
		out_file << xsrc << " " << ysrc << " " << zsrc << "\n" ;
		break ;
	}
	}

	out_file.close() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Acquisition::write_src_VTK");
	return(RTN_CODE_OK) ;
}

} // namespace django


