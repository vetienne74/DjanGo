//------------------------------------------------------------------------------------------------------
//
// MODELLING WITH FDM IN 2D
//
//  * ACOUSTIC ISOTROPIC MEDIA 
//  * 1ST ORDER WAVE EQUATION (PRESSURE / VELOCITY):
//
// PARENT CLASS: Scheme
//   DERIVED CLASS: FDM
//     DERIVED CLASS: FDM_1D
//       DERIVED CLASS: FDM_2D
//         DERIVED CLASS: FDM_2D_ac_iso
//           DERIVED CLASS: FDM_2D_ac_iso_1st
//             DERIVED CLASS: FDM_2D_ac_iso_1st_XX_2
//
//-------------------------------------------------------------------------------------------------------

#include "FDM_2D_ac_iso_1st_schemes.h"

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
FDM_2D_ac_iso_1st_2_2::FDM_2D_ac_iso_1st_2_2(void)
{
	space_order = 2 ;
	time_order  = 2 ;
	lstencil = FDM_O1_2_2_LSTENCIL ;
	nb_op_d  = NB_OP_O1_2 ;
	CFL      = FDM_O1_2_2_CFL ;
}

// 4TH SPATIAL ORDER
FDM_2D_ac_iso_1st_4_2::FDM_2D_ac_iso_1st_4_2(void)
{
	space_order = 4 ;
	time_order  = 2 ;
	lstencil = FDM_O1_4_2_LSTENCIL ;
	nb_op_d  = NB_OP_O1_4 ;
	CFL      = FDM_O1_4_2_CFL ;
}

// 8TH SPATIAL ORDER
FDM_2D_ac_iso_1st_8_2::FDM_2D_ac_iso_1st_8_2(void)
{
	space_order = 8 ;
	time_order  = 2 ;
	lstencil = FDM_O1_8_2_LSTENCIL ;
	nb_op_d  = NB_OP_O1_8 ;
	CFL      = FDM_O1_8_2_CFL ;
}

// 12TH SPATIAL ORDER
FDM_2D_ac_iso_1st_12_2::FDM_2D_ac_iso_1st_12_2(void)
{
	space_order = 12 ;
	time_order  = 2 ;
	lstencil = FDM_O1_12_2_LSTENCIL ;
	nb_op_d  = NB_OP_O1_12 ;
	CFL      = FDM_O1_12_2_CFL ;
}

// 16TH SPATIAL ORDER
FDM_2D_ac_iso_1st_16_2::FDM_2D_ac_iso_1st_16_2(void)
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
inline Myfloat FDM_2D_ac_iso_1st_2_2::D_Z(Myfloat** U, Myint iz, Myint ix) 
{
	return ((U[ix][iz] - U[ix][iz-1]) * inv_dz) ;
}
inline Myfloat FDM_2D_ac_iso_1st_2_2::D_X(Myfloat** U, Myint iz, Myint ix) 
{
	return ((U[ix][iz] - U[ix-1][iz]) * inv_dx) ;
}

# define D_Z_2D_O4(U, iz, ix) ((FDM_O1_4_2_A1 * (U[ix][iz] - U[ix][iz-1]) + FDM_O1_4_2_A2 * (U[ix][iz+1] - U[ix][iz-2])) * inv_dz)
# define D_X_2D_O4(U, iz, ix) ((FDM_O1_4_2_A1 * (U[ix][iz] - U[ix-1][iz]) + FDM_O1_4_2_A2 * (U[ix+1][iz] - U[ix-2][iz])) * inv_dx) 

// 4TH SPATIAL ORDER
inline Myfloat FDM_2D_ac_iso_1st_4_2::D_Z(Myfloat** U, Myint iz, Myint ix) 
{
	//return ((FDM_O1_4_2_A1 * (U[ix][iz] - U[ix][iz-1]) + FDM_O1_4_2_A2 * (U[ix][iz+1] - U[ix][iz-2])) * inv_dz) ;
	return D_Z_2D_O4(U, iz, ix) ;
}
inline Myfloat FDM_2D_ac_iso_1st_4_2::D_X(Myfloat** U, Myint iz, Myint ix) 
{
	//return ((FDM_O1_4_2_A1 * (U[ix][iz] - U[ix-1][iz]) + FDM_O1_4_2_A2 * (U[ix+1][iz] - U[ix-2][iz])) * inv_dx) ;
	return D_X_2D_O4(U, iz, ix) ;
}

// 8TH SPATIAL ORDER
inline Myfloat FDM_2D_ac_iso_1st_8_2::D_Z(Myfloat** U, Myint iz, Myint ix) 
{
	return ((FDM_O1_8_2_A1 * (U[ix][iz] - U[ix][iz-1]) + FDM_O1_8_2_A2 * (U[ix][iz+1] - U[ix][iz-2])
			+ FDM_O1_8_2_A3 * (U[ix][iz+2] - U[ix][iz-3]) + FDM_O1_8_2_A4 * (U[ix][iz+3] - U[ix][iz-4])) * inv_dz) ;
}
inline Myfloat FDM_2D_ac_iso_1st_8_2::D_X(Myfloat** U, Myint iz, Myint ix) 
{
	return ((FDM_O1_8_2_A1 * (U[ix][iz] - U[ix-1][iz]) + FDM_O1_8_2_A2 * (U[ix+1][iz] - U[ix-2][iz])
			+ FDM_O1_8_2_A3 * (U[ix+2][iz] - U[ix-3][iz]) + FDM_O1_8_2_A4 * (U[ix+3][iz] - U[ix-4][iz])) * inv_dx) ;
}

// 12TH SPATIAL ORDER
inline Myfloat FDM_2D_ac_iso_1st_12_2::D_Z(Myfloat** U, Myint iz, Myint ix)
{
	return ((FDM_O1_12_2_A1 * (U[ix][iz] - U[ix][iz-1]) + FDM_O1_12_2_A2 * (U[ix][iz+1] - U[ix][iz-2])
			+ FDM_O1_12_2_A3 * (U[ix][iz+2] - U[ix][iz-3]) + FDM_O1_12_2_A4 * (U[ix][iz+3] - U[ix][iz-4])
			+ FDM_O1_12_2_A5 * (U[ix][iz+4] - U[ix][iz-5]) + FDM_O1_12_2_A6 * (U[ix][iz+5] - U[ix][iz-6])) * inv_dz) ;
}
inline Myfloat FDM_2D_ac_iso_1st_12_2::D_X(Myfloat** U, Myint iz, Myint ix)
{
	return ((FDM_O1_12_2_A1 * (U[ix][iz] - U[ix-1][iz]) + FDM_O1_12_2_A2 * (U[ix+1][iz] - U[ix-2][iz])
			+ FDM_O1_12_2_A3 * (U[ix+2][iz] - U[ix-3][iz]) + FDM_O1_12_2_A4 * (U[ix+3][iz] - U[ix-4][iz])
			+ FDM_O1_12_2_A5 * (U[ix+4][iz] - U[ix-5][iz]) + FDM_O1_12_2_A6 * (U[ix+5][iz] - U[ix-6][iz])) * inv_dx) ;
}

// 16TH SPATIAL ORDER
inline Myfloat FDM_2D_ac_iso_1st_16_2::D_Z(Myfloat** U, Myint iz, Myint ix)
{
	return ((FDM_O1_16_2_A1 * (U[ix][iz] - U[ix][iz-1]) + FDM_O1_16_2_A2 * (U[ix][iz+1] - U[ix][iz-2])
			+ FDM_O1_16_2_A3 * (U[ix][iz+2] - U[ix][iz-3]) + FDM_O1_16_2_A4 * (U[ix][iz+3] - U[ix][iz-4])
			+ FDM_O1_16_2_A5 * (U[ix][iz+4] - U[ix][iz-5]) + FDM_O1_16_2_A6 * (U[ix][iz+5] - U[ix][iz-6])
			+ FDM_O1_16_2_A7 * (U[ix][iz+6] - U[ix][iz-7]) + FDM_O1_16_2_A8 * (U[ix][iz+7] - U[ix][iz-8])) * inv_dz) ;
}
inline Myfloat FDM_2D_ac_iso_1st_16_2::D_X(Myfloat** U, Myint iz, Myint ix)
{
	return ((FDM_O1_16_2_A1 * (U[ix][iz] - U[ix-1][iz]) + FDM_O1_16_2_A2 * (U[ix+1][iz] - U[ix-2][iz])
			+ FDM_O1_16_2_A3 * (U[ix+2][iz] - U[ix-3][iz]) + FDM_O1_16_2_A4 * (U[ix+3][iz] - U[ix-4][iz])
			+ FDM_O1_16_2_A5 * (U[ix+4][iz] - U[ix-5][iz]) + FDM_O1_16_2_A6 * (U[ix+5][iz] - U[ix-6][iz])
			+ FDM_O1_16_2_A7 * (U[ix+6][iz] - U[ix-7][iz]) + FDM_O1_16_2_A8 * (U[ix+7][iz] - U[ix-8][iz])) * inv_dx) ;
}

//=======================================================================================================
//
// COMPUTE PRESSURE COMPONENT
//
//=======================================================================================================

// 2ND SPATIAL ORDER
void FDM_2D_ac_iso_1st_2_2::compute_pressure()
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_2D_ac_iso_1st_2_2::compute_pressure");

#include <FDM_2D_ac_iso_1st_compute_pressure.inc>

	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D_ac_iso_1st_2_2::compute_pressure");
}

// 4TH SPATIAL ORDER
void FDM_2D_ac_iso_1st_4_2::compute_pressure()
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_2D_ac_iso_1st_4_2::compute_pressure");

#include <FDM_2D_ac_iso_1st_compute_pressure.inc>

	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D_ac_iso_1st_4_2::compute_pressure");
}

// 8TH SPATIAL ORDER
void FDM_2D_ac_iso_1st_8_2::compute_pressure()
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_2D_ac_iso_1st_8_2::compute_pressure");

#include <FDM_2D_ac_iso_1st_compute_pressure.inc>

	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D_ac_iso_1st_8_2::compute_pressure");
}

// 12TH SPATIAL ORDER
void FDM_2D_ac_iso_1st_12_2::compute_pressure()
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_2D_ac_iso_1st_12_2::compute_pressure");

#include <FDM_2D_ac_iso_1st_compute_pressure.inc>

	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D_ac_iso_1st_12_2::compute_pressure");
}

// 16TH SPATIAL ORDER
void FDM_2D_ac_iso_1st_16_2::compute_pressure()
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_2D_ac_iso_1st_16_2::compute_pressure");

#include <FDM_2D_ac_iso_1st_compute_pressure.inc>

	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D_ac_iso_1st_16_2::compute_pressure");
}


//=======================================================================================================
//
// COMPUTE VELOCITY COMPONENT
//
//=======================================================================================================

// 2ND SPATIAL ORDER
void FDM_2D_ac_iso_1st_2_2::compute_velocity()
{      
	print_debug(ALL, FULL_DEBUG, "IN FDM_2D_ac_iso_1st_2_2::compute_velocity");

#include <FDM_2D_ac_iso_1st_compute_velocity.inc>  

	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D_ac_iso_1st_2_2::compute_velocity");
} 

// 4TH SPATIAL ORDER
void FDM_2D_ac_iso_1st_4_2::compute_velocity()
{     

	print_debug(ALL, FULL_DEBUG, "IN FDM_2D_ac_iso_1st_4_2::compute_velocity");

#include <FDM_2D_ac_iso_1st_compute_velocity.inc>  

	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D_ac_iso_1st_4_2::compute_velocity");
} 

// 8TH SPATIAL ORDER
void FDM_2D_ac_iso_1st_8_2::compute_velocity()
{      
	print_debug(ALL, FULL_DEBUG, "IN FDM_2D_ac_iso_1st_8_2::compute_velocity");

#include <FDM_2D_ac_iso_1st_compute_velocity.inc>  

	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D_ac_iso_1st_8_2::compute_velocity");
} 

// 12TH SPATIAL ORDER
void FDM_2D_ac_iso_1st_12_2::compute_velocity()
{      
	print_debug(ALL, FULL_DEBUG, "IN FDM_2D_ac_iso_1st_12_2::compute_velocity");

#include <FDM_2D_ac_iso_1st_compute_velocity.inc>  

	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D_ac_iso_1st_12_2::compute_velocity");
} 

// 16TH SPATIAL ORDER
void FDM_2D_ac_iso_1st_16_2::compute_velocity()
{      
	print_debug(ALL, FULL_DEBUG, "IN FDM_2D_ac_iso_1st_16_2::compute_velocity");

#include <FDM_2D_ac_iso_1st_compute_velocity.inc>  

	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D_ac_iso_1st_16_2::compute_velocity");
} 

} // namespace django
