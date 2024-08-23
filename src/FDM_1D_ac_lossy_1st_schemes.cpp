//------------------------------------------------------------------------------------------------------
//
// MODELLING WITH FDM IN 1D
//
//  * ACOUSTIC ISOTROPIC MEDIA 
//  * 1ST ORDER WAVE EQUATION (PRESSURE / VELOCITY):
//
// PARENT CLASS: Scheme
//   DERIVED CLASS: FDM
//     DERIVED CLASS: FDM_1D
//       DERIVED CLASS: FDM_1D_ac_lossy
//         DERIVED CLASS: FDM_1D_ac_lossy_1st
//           DERIVED CLASS: FDM_1D_ac_lossy_1st_XX_2
//
//-------------------------------------------------------------------------------------------------------

#include "FDM_1D_ac_lossy_1st_schemes.h"

#include <iostream>

#include "output_report.h"
#include "singleton.h"

using namespace std;

namespace django {

//=======================================================================================================
//
// CONSTRUCTOR
//
//=======================================================================================================

// 2ND SPATIAL ORDER
FDM_1D_ac_lossy_1st_2_2::FDM_1D_ac_lossy_1st_2_2(void)
{
	space_order = 2 ;
	time_order  = 2 ;
	lstencil = FDM_O1_2_2_LSTENCIL ;
	nb_op_d  = NB_OP_O1_2 ;
	CFL      = FDM_O1_2_2_CFL ;
}

// 4TH SPATIAL ORDER
FDM_1D_ac_lossy_1st_4_2::FDM_1D_ac_lossy_1st_4_2(void)
{
	space_order = 4 ;
	time_order  = 2 ;
	lstencil = FDM_O1_4_2_LSTENCIL ;
	nb_op_d  = NB_OP_O1_4 ;
	CFL      = FDM_O1_4_2_CFL ;
}

// 8TH SPATIAL ORDER
FDM_1D_ac_lossy_1st_8_2::FDM_1D_ac_lossy_1st_8_2(void)
{
	space_order = 8 ;
	time_order  = 2 ;
	lstencil = FDM_O1_8_2_LSTENCIL ;
	nb_op_d  = NB_OP_O1_8 ;
	CFL      = FDM_O1_8_2_CFL ;
}

// 12TH SPATIAL ORDER
FDM_1D_ac_lossy_1st_12_2::FDM_1D_ac_lossy_1st_12_2(void)
{
	space_order = 12 ;
	time_order  = 2 ;
	lstencil = FDM_O1_12_2_LSTENCIL ;
	nb_op_d  = NB_OP_O1_12 ;
	CFL      = FDM_O1_12_2_CFL ;
}

// 16TH SPATIAL ORDER
FDM_1D_ac_lossy_1st_16_2::FDM_1D_ac_lossy_1st_16_2(void)
{
	space_order = 16 ;
	time_order  = 2 ;
	lstencil = FDM_O1_16_2_LSTENCIL ;
	nb_op_d  = NB_OP_O1_16 ;
	CFL      = FDM_O1_16_2_CFL ;
}


//=======================================================================================================
//
// FDM OPERATOR
//
//=======================================================================================================

// 2ND SPATIAL ORDER
inline Myfloat FDM_1D_ac_lossy_1st_2_2::D_Z(Myfloat* U, Myint iz) 
{
	return ((U[iz] - U[iz-1]) * inv_dz) ;
}

// 4TH SPATIAL ORDER
inline Myfloat FDM_1D_ac_lossy_1st_4_2::D_Z(Myfloat* U, Myint iz) 
{
	return ((FDM_O1_4_2_A1 * (U[iz] - U[iz-1]) + FDM_O1_4_2_A2 * (U[iz+1] - U[iz-2])) * inv_dz) ;
}

// 8TH SPATIAL ORDER
inline Myfloat FDM_1D_ac_lossy_1st_8_2::D_Z(Myfloat* U, Myint iz) 
{
	return ((FDM_O1_8_2_A1 * (U[iz] - U[iz-1]) + FDM_O1_8_2_A2 * (U[iz+1] - U[iz-2])
			+ FDM_O1_8_2_A3 * (U[iz+2] - U[iz-3]) + FDM_O1_8_2_A4 * (U[iz+3] - U[iz-4])) * inv_dz) ;
}

// 12TH SPATIAL ORDER
inline Myfloat FDM_1D_ac_lossy_1st_12_2::D_Z(Myfloat* U, Myint iz) 
{
	return ((FDM_O1_12_2_A1 * (U[iz] - U[iz-1]) + FDM_O1_12_2_A2 * (U[iz+1] - U[iz-2])
			+ FDM_O1_12_2_A3 * (U[iz+2] - U[iz-3]) + FDM_O1_12_2_A4 * (U[iz+3] - U[iz-4])
			+ FDM_O1_12_2_A5 * (U[iz+4] - U[iz-5]) + FDM_O1_12_2_A6 * (U[iz+5] - U[iz-6])) * inv_dz) ;
}

// 16TH SPATIAL ORDER
inline Myfloat FDM_1D_ac_lossy_1st_16_2::D_Z(Myfloat* U, Myint iz) 
{
	return ((FDM_O1_16_2_A1 * (U[iz] - U[iz-1]) + FDM_O1_16_2_A2 * (U[iz+1] - U[iz-2])
			+ FDM_O1_16_2_A3 * (U[iz+2] - U[iz-3]) + FDM_O1_16_2_A4 * (U[iz+3] - U[iz-4])
			+ FDM_O1_16_2_A5 * (U[iz+4] - U[iz-5]) + FDM_O1_16_2_A6 * (U[iz+5] - U[iz-6])
			+ FDM_O1_16_2_A7 * (U[iz+6] - U[iz-7]) + FDM_O1_16_2_A8 * (U[iz+7] - U[iz-8])) * inv_dz) ;
}

//=======================================================================================================
//
// COMPUTE PRESSURE COMPONENT
//
//=======================================================================================================

// 2ND SPATIAL ORDER
void FDM_1D_ac_lossy_1st_2_2::compute_pressure()
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_1D_ac_lossy_1st_2_2::compute_pressure");

#include <FDM_1D_ac_lossy_1st_compute_pressure.inc>

	print_debug(ALL, FULL_DEBUG, "OUT FDM_1D_ac_lossy_1st_2_2::compute_pressure");
}

// 4TH SPATIAL ORDER
void FDM_1D_ac_lossy_1st_4_2::compute_pressure()
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_1D_ac_lossy_1st_4_2::compute_pressure");

#include <FDM_1D_ac_lossy_1st_compute_pressure.inc>

	print_debug(ALL, FULL_DEBUG, "OUT FDM_1D_ac_lossy_1st_4_2::compute_pressure");
}

// 8TH SPATIAL ORDER
void FDM_1D_ac_lossy_1st_8_2::compute_pressure()
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_1D_ac_lossy_1st_8_2::compute_pressure");

#include <FDM_1D_ac_lossy_1st_compute_pressure.inc>

	print_debug(ALL, FULL_DEBUG, "OUT FDM_1D_ac_lossy_1st_8_2::compute_pressure");
}

// 12TH SPATIAL ORDER
void FDM_1D_ac_lossy_1st_12_2::compute_pressure()
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_1D_ac_lossy_1st_12_2::compute_pressure");

#include <FDM_1D_ac_lossy_1st_compute_pressure.inc>

	print_debug(ALL, FULL_DEBUG, "OUT FDM_1D_ac_lossy_1st_12_2::compute_pressure");
}

// 16TH SPATIAL ORDER
void FDM_1D_ac_lossy_1st_16_2::compute_pressure()
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_1D_ac_lossy_1st_16_2::compute_pressure");

#include <FDM_1D_ac_lossy_1st_compute_pressure.inc>

	print_debug(ALL, FULL_DEBUG, "OUT FDM_1D_ac_lossy_1st_16_2::compute_pressure");
}


//=======================================================================================================
//
// COMPUTE VELOCITY COMPONENT
//
//=======================================================================================================

// 2ND SPATIAL ORDER
void FDM_1D_ac_lossy_1st_2_2::compute_velocity()
{      
	print_debug(ALL, FULL_DEBUG, "IN FDM_1D_ac_lossy_1st_2_2::compute_velocity");

#include <FDM_1D_ac_lossy_1st_compute_velocity.inc>  

	print_debug(ALL, FULL_DEBUG, "OUT FDM_1D_ac_lossy_1st_2_2::compute_velocity");
} 

// 4TH SPATIAL ORDER
void FDM_1D_ac_lossy_1st_4_2::compute_velocity()
{     

	print_debug(ALL, FULL_DEBUG, "IN FDM_1D_ac_lossy_1st_4_2::compute_velocity");

#include <FDM_1D_ac_lossy_1st_compute_velocity.inc>  

	print_debug(ALL, FULL_DEBUG, "OUT FDM_1D_ac_lossy_1st_4_2::compute_velocity");
} 

// 8TH SPATIAL ORDER
void FDM_1D_ac_lossy_1st_8_2::compute_velocity()
{      
	print_debug(ALL, FULL_DEBUG, "IN FDM_1D_ac_lossy_1st_8_2::compute_velocity");

#include <FDM_1D_ac_lossy_1st_compute_velocity.inc>  

	print_debug(ALL, FULL_DEBUG, "OUT FDM_1D_ac_lossy_1st_8_2::compute_velocity");
} 

// 12TH SPATIAL ORDER
void FDM_1D_ac_lossy_1st_12_2::compute_velocity()
{      
	print_debug(ALL, FULL_DEBUG, "IN FDM_1D_ac_lossy_1st_12_2::compute_velocity");

#include <FDM_1D_ac_lossy_1st_compute_velocity.inc>  

	print_debug(ALL, FULL_DEBUG, "OUT FDM_1D_ac_lossy_1st_12_2::compute_velocity");
} 

// 16TH SPATIAL ORDER
void FDM_1D_ac_lossy_1st_16_2::compute_velocity()
{      
	print_debug(ALL, FULL_DEBUG, "IN FDM_1D_ac_lossy_1st_16_2::compute_velocity");

#include <FDM_1D_ac_lossy_1st_compute_velocity.inc>  

	print_debug(ALL, FULL_DEBUG, "OUT FDM_1D_ac_lossy_1st_16_2::compute_velocity");
} 

} // namespace django
