//------------------------------------------------------------------------------------------------------
//
// MODELLING WITH FDM IN 2D
//
//  * ACOUSTIC ISOTROPIC MEDIA 
//  * 2ND ORDER WAVE EQUATION (PRESSURE):
//
// PARENT CLASS: Schemme
//   DERIVED CLASS: FDM
//     DERIVED CLASS: FDM_1D
//       DERIVED CLASS: FDM_2D
//         DERIVED CLASS: FDM_2D_ac_iso
//           DERIVED CLASS: FDM_2D_ac_iso_2nd
//             DERIVED CLASS: FDM_2D_ac_iso_2nd_XX_2
//
//-------------------------------------------------------------------------------------------------------

#include "FDM_2D_ac_iso_2nd_schemes.h"

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
FDM_2D_ac_iso_2nd_2_2::FDM_2D_ac_iso_2nd_2_2(void)
{
	space_order = 2 ;
	time_order  = 2 ;
	lstencil = FDM_O2_2_2_LSTENCIL ;
	nb_op_d  = NB_OP_O1_2 ;
	nb_op_d2 = NB_OP_O2_2 ;
	CFL      = FDM_O2_2_2_CFL ;
}

// 4TH SPATIAL ORDER
FDM_2D_ac_iso_2nd_4_2::FDM_2D_ac_iso_2nd_4_2(void)
{
	space_order = 4 ;
	time_order  = 2 ;
	lstencil = FDM_O2_4_2_LSTENCIL ;
	nb_op_d  = NB_OP_O1_4 ;
	nb_op_d2 = NB_OP_O2_4 ;
	CFL      = FDM_O2_4_2_CFL ;
}

// 8TH SPATIAL ORDER
FDM_2D_ac_iso_2nd_8_2::FDM_2D_ac_iso_2nd_8_2(void)
{
	space_order = 8 ;
	time_order  = 2 ;
	lstencil = FDM_O2_8_2_LSTENCIL ;
	nb_op_d  = NB_OP_O1_8 ;
	nb_op_d2 = NB_OP_O2_8 ;
	CFL      = FDM_O2_8_2_CFL ;
}

// 12TH SPATIAL ORDER
FDM_2D_ac_iso_2nd_12_2::FDM_2D_ac_iso_2nd_12_2(void)
{
	space_order = 12 ;
	time_order  = 2 ;
	lstencil = FDM_O2_12_2_LSTENCIL ;
	nb_op_d  = NB_OP_O1_12 ;
	nb_op_d2 = NB_OP_O2_12 ;
	CFL      = FDM_O2_12_2_CFL ;
}

// 16TH SPATIAL ORDER
FDM_2D_ac_iso_2nd_16_2::FDM_2D_ac_iso_2nd_16_2(void)
{
	space_order = 16 ;
	time_order  = 2 ;
	lstencil = FDM_O2_16_2_LSTENCIL ;
	nb_op_d  = NB_OP_O1_16 ;
	nb_op_d2 = NB_OP_O2_16 ;
	CFL      = FDM_O2_16_2_CFL ;
}


//=======================================================================================================
//
// FDM OPERATOR
//
//=======================================================================================================

//-------------------------------------------------------------------------------------------------------
// 1ST DERIVATIVE
//-------------------------------------------------------------------------------------------------------

// 2ND SPATIAL ORDER
inline Myfloat FDM_2D_ac_iso_2nd_2_2::D_Z(Myfloat** U, Myint iz, Myint ix) 
{
	return ((U[ix][iz] - U[ix][iz-1]) * inv_dz) ;
}
inline Myfloat FDM_2D_ac_iso_2nd_2_2::D_X(Myfloat** U, Myint iz, Myint ix) 
{
	return ((U[ix][iz] - U[ix-1][iz]) * inv_dx) ;
}

// 4TH SPATIAL ORDER
inline Myfloat FDM_2D_ac_iso_2nd_4_2::D_Z(Myfloat** U, Myint iz, Myint ix) 
{
	return ((FDM_O1_4_2_A1 * (U[ix][iz] - U[ix][iz-1]) + FDM_O1_4_2_A2 * (U[ix][iz+1] - U[ix][iz-2])) * inv_dz) ;
}
inline Myfloat FDM_2D_ac_iso_2nd_4_2::D_X(Myfloat** U, Myint iz, Myint ix) 
{
	return ((FDM_O1_4_2_A1 * (U[ix][iz] - U[ix-1][iz]) + FDM_O1_4_2_A2 * (U[ix+1][iz] - U[ix-2][iz])) * inv_dx) ;
}

// 8TH SPATIAL ORDER
inline Myfloat FDM_2D_ac_iso_2nd_8_2::D_Z(Myfloat** U, Myint iz, Myint ix) 
{
	return ((FDM_O1_8_2_A1 * (U[ix][iz] - U[ix][iz-1]) + FDM_O1_8_2_A2 * (U[ix][iz+1] - U[ix][iz-2])
			+ FDM_O1_8_2_A3 * (U[ix][iz+2] - U[ix][iz-3]) + FDM_O1_8_2_A4 * (U[ix][iz+3] - U[ix][iz-4])) * inv_dz) ;
}
inline Myfloat FDM_2D_ac_iso_2nd_8_2::D_X(Myfloat** U, Myint iz, Myint ix) 
{
	return ((FDM_O1_8_2_A1 * (U[ix][iz] - U[ix-1][iz]) + FDM_O1_8_2_A2 * (U[ix+1][iz] - U[ix-2][iz])
			+ FDM_O1_8_2_A3 * (U[ix+2][iz] - U[ix-3][iz]) + FDM_O1_8_2_A4 * (U[ix+3][iz] - U[ix-4][iz])) * inv_dx) ;
}

// 12TH SPATIAL ORDER
inline Myfloat FDM_2D_ac_iso_2nd_12_2::D_Z(Myfloat** U, Myint iz, Myint ix)
{
	return ((FDM_O1_12_2_A1 * (U[ix][iz] - U[ix][iz-1]) + FDM_O1_12_2_A2 * (U[ix][iz+1] - U[ix][iz-2])
			+ FDM_O1_12_2_A3 * (U[ix][iz+2] - U[ix][iz-3]) + FDM_O1_12_2_A4 * (U[ix][iz+3] - U[ix][iz-4])
			+ FDM_O1_12_2_A5 * (U[ix][iz+4] - U[ix][iz-5]) + FDM_O1_12_2_A6 * (U[ix][iz+5] - U[ix][iz-6])) * inv_dz) ;
}
inline Myfloat FDM_2D_ac_iso_2nd_12_2::D_X(Myfloat** U, Myint iz, Myint ix)
{
	return ((FDM_O1_12_2_A1 * (U[ix][iz] - U[ix-1][iz]) + FDM_O1_12_2_A2 * (U[ix+1][iz] - U[ix-2][iz])
			+ FDM_O1_12_2_A3 * (U[ix+2][iz] - U[ix-3][iz]) + FDM_O1_12_2_A4 * (U[ix+3][iz] - U[ix-4][iz])
			+ FDM_O1_12_2_A5 * (U[ix+4][iz] - U[ix-5][iz]) + FDM_O1_12_2_A6 * (U[ix+5][iz] - U[ix-6][iz])) * inv_dx) ;
}

// 16TH SPATIAL ORDER
inline Myfloat FDM_2D_ac_iso_2nd_16_2::D_Z(Myfloat** U, Myint iz, Myint ix)
{
	return ((FDM_O1_16_2_A1 * (U[ix][iz] - U[ix][iz-1]) + FDM_O1_16_2_A2 * (U[ix][iz+1] - U[ix][iz-2])
			+ FDM_O1_16_2_A3 * (U[ix][iz+2] - U[ix][iz-3]) + FDM_O1_16_2_A4 * (U[ix][iz+3] - U[ix][iz-4])
			+ FDM_O1_16_2_A5 * (U[ix][iz+4] - U[ix][iz-5]) + FDM_O1_16_2_A6 * (U[ix][iz+5] - U[ix][iz-6])
			+ FDM_O1_16_2_A7 * (U[ix][iz+6] - U[ix][iz-7]) + FDM_O1_16_2_A8 * (U[ix][iz+7] - U[ix][iz-8])) * inv_dz) ;
}
inline Myfloat FDM_2D_ac_iso_2nd_16_2::D_X(Myfloat** U, Myint iz, Myint ix)
{
	return ((FDM_O1_16_2_A1 * (U[ix][iz] - U[ix-1][iz]) + FDM_O1_16_2_A2 * (U[ix+1][iz] - U[ix-2][iz])
			+ FDM_O1_16_2_A3 * (U[ix+2][iz] - U[ix-3][iz]) + FDM_O1_16_2_A4 * (U[ix+3][iz] - U[ix-4][iz])
			+ FDM_O1_16_2_A5 * (U[ix+4][iz] - U[ix-5][iz]) + FDM_O1_16_2_A6 * (U[ix+5][iz] - U[ix-6][iz])
			+ FDM_O1_16_2_A7 * (U[ix+6][iz] - U[ix-7][iz]) + FDM_O1_16_2_A8 * (U[ix+7][iz] - U[ix-8][iz])) * inv_dx) ;
}

//-------------------------------------------------------------------------------------------------------
// 2ND DERIVATIVE
//-------------------------------------------------------------------------------------------------------

// 2ND SPATIAL ORDER
inline Myfloat FDM_2D_ac_iso_2nd_2_2::D2_Z(Myfloat** U, Myint iz, Myint ix) 
{
	return ((-2 * U[ix][iz] + U[ix][iz-1] + U[ix][iz+1]) * inv_dz2) ;
}
inline Myfloat FDM_2D_ac_iso_2nd_2_2::D2_X(Myfloat** U, Myint iz, Myint ix) 
{
	return ((-2 * U[ix][iz] + U[ix-1][iz] + U[ix+1][iz]) * inv_dx2) ;
}

// 4TH SPATIAL ORDER
inline Myfloat FDM_2D_ac_iso_2nd_4_2::D2_Z(Myfloat** U, Myint iz, Myint ix) 
{
	return ((FDM_O2_4_2_A0 * U[ix][iz] + FDM_O2_4_2_A1 * (U[ix][iz-1] + U[ix][iz+1])
			+ FDM_O2_4_2_A2 * (U[ix][iz-2] + U[ix][iz+2])) * inv_dz2) ;
}
inline Myfloat FDM_2D_ac_iso_2nd_4_2::D2_X(Myfloat** U, Myint iz, Myint ix) 
{
	return ((FDM_O2_4_2_A0 * U[ix][iz] + FDM_O2_4_2_A1 * (U[ix-1][iz] + U[ix+1][iz])
			+ FDM_O2_4_2_A2 * (U[ix-2][iz] + U[ix+2][iz])) * inv_dx2) ;
}

// 8TH SPATIAL ORDER
inline Myfloat FDM_2D_ac_iso_2nd_8_2::D2_Z(Myfloat** U, Myint iz, Myint ix) 
{
	return ((FDM_O2_8_2_A0 * U[ix][iz] + FDM_O2_8_2_A1 * (U[ix][iz-1] + U[ix][iz+1])
			+ FDM_O2_8_2_A2 * (U[ix][iz-2] + U[ix][iz+2]) + FDM_O2_8_2_A3 * (U[ix][iz-3] + U[ix][iz+3])
			+ FDM_O2_8_2_A4 * (U[ix][iz-4] + U[ix][iz+4])) * inv_dz2) ;
}
inline Myfloat FDM_2D_ac_iso_2nd_8_2::D2_X(Myfloat** U, Myint iz, Myint ix) 
{
	return ((FDM_O2_8_2_A0 * U[ix][iz] + FDM_O2_8_2_A1 * (U[ix-1][iz] + U[ix+1][iz])
			+ FDM_O2_8_2_A2 * (U[ix-2][iz] + U[ix+2][iz]) + FDM_O2_8_2_A3 * (U[ix-3][iz] + U[ix+3][iz])
			+ FDM_O2_8_2_A4 * (U[ix-4][iz] + U[ix+4][iz])) * inv_dx2) ;
}

// 12TH SPATIAL ORDER
inline Myfloat FDM_2D_ac_iso_2nd_12_2::D2_Z(Myfloat** U, Myint iz, Myint ix)
{
	return ((FDM_O2_12_2_A0 * U[ix][iz] + FDM_O2_12_2_A1 * (U[ix][iz-1] + U[ix][iz+1])
			+ FDM_O2_12_2_A2 * (U[ix][iz-2] + U[ix][iz+2]) + FDM_O2_12_2_A3 * (U[ix][iz-3] + U[ix][iz+3])
			+ FDM_O2_12_2_A4 * (U[ix][iz-4] + U[ix][iz+4]) + FDM_O2_12_2_A5 * (U[ix][iz-5] + U[ix][iz+5])
			+ FDM_O2_12_2_A6 * (U[ix][iz-6] + U[ix][iz+6])) * inv_dz2) ;
}
inline Myfloat FDM_2D_ac_iso_2nd_12_2::D2_X(Myfloat** U, Myint iz, Myint ix)
{
	return ((FDM_O2_12_2_A0 * U[ix][iz] + FDM_O2_12_2_A1 * (U[ix-1][iz] + U[ix+1][iz])
			+ FDM_O2_12_2_A2 * (U[ix-2][iz] + U[ix+2][iz]) + FDM_O2_12_2_A3 * (U[ix-3][iz] + U[ix+3][iz])
			+ FDM_O2_12_2_A4 * (U[ix-4][iz] + U[ix+4][iz]) + FDM_O2_12_2_A5 * (U[ix-5][iz] + U[ix+5][iz])
			+ FDM_O2_12_2_A6 * (U[ix-6][iz] + U[ix+6][iz])) * inv_dx2) ;
}

// 16TH SPATIAL ORDER
inline Myfloat FDM_2D_ac_iso_2nd_16_2::D2_Z(Myfloat** U, Myint iz, Myint ix)
{
	return ((FDM_O2_16_2_A0 * U[ix][iz] + FDM_O2_16_2_A1 * (U[ix][iz-1] + U[ix][iz+1])
			+ FDM_O2_16_2_A2 * (U[ix][iz-2] + U[ix][iz+2]) + FDM_O2_16_2_A3 * (U[ix][iz-3] + U[ix][iz+3])
			+ FDM_O2_16_2_A4 * (U[ix][iz-4] + U[ix][iz+4]) + FDM_O2_16_2_A5 * (U[ix][iz-5] + U[ix][iz+5])
			+ FDM_O2_16_2_A6 * (U[ix][iz-6] + U[ix][iz+6]) + FDM_O2_16_2_A7 * (U[ix][iz-7] + U[ix][iz+7])
			+ FDM_O2_16_2_A8 * (U[ix][iz-8] + U[ix][iz+8])) * inv_dz2) ;
}
inline Myfloat FDM_2D_ac_iso_2nd_16_2::D2_X(Myfloat** U, Myint iz, Myint ix)
{
	return ((FDM_O2_16_2_A0 * U[ix][iz] + FDM_O2_16_2_A1 * (U[ix-1][iz] + U[ix+1][iz])
			+ FDM_O2_16_2_A2 * (U[ix-2][iz] + U[ix+2][iz]) + FDM_O2_16_2_A3 * (U[ix-3][iz] + U[ix+3][iz])
			+ FDM_O2_16_2_A4 * (U[ix-4][iz] + U[ix+4][iz]) + FDM_O2_16_2_A5 * (U[ix-5][iz] + U[ix+5][iz])
			+ FDM_O2_16_2_A6 * (U[ix-6][iz] + U[ix+6][iz]) + FDM_O2_16_2_A7 * (U[ix-7][iz] + U[ix+7][iz])
			+ FDM_O2_16_2_A8 * (U[ix-8][iz] + U[ix+8][iz])) * inv_dx2) ;
}

//=======================================================================================================
//
// COMPUTE PRESSURE COMPONENT
//
//=======================================================================================================

// 2ND SPATIAL ORDER
void FDM_2D_ac_iso_2nd_2_2::compute_pressure()
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_2D_ac_iso_2nd_2_2::compute_pressure");

#include <FDM_2D_ac_iso_2nd_compute_pressure.inc>

	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D_ac_iso_2nd_2_2::compute_pressure");
}

// 4TH SPATIAL ORDER
void FDM_2D_ac_iso_2nd_4_2::compute_pressure()
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_2D_ac_iso_2nd_4_2::compute_pressure");

#include <FDM_2D_ac_iso_2nd_compute_pressure.inc>

	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D_ac_iso_2nd_4_2::compute_pressure");
}

// 8TH SPATIAL ORDER
void FDM_2D_ac_iso_2nd_8_2::compute_pressure()
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_2D_ac_iso_2nd_8_2::compute_pressure");

#include <FDM_2D_ac_iso_2nd_compute_pressure.inc>

	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D_ac_iso_2nd_8_2::compute_pressure");
}

// 12TH SPATIAL ORDER
void FDM_2D_ac_iso_2nd_12_2::compute_pressure()
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_2D_ac_iso_2nd_12_2::compute_pressure");

#include <FDM_2D_ac_iso_2nd_compute_pressure.inc>

	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D_ac_iso_2nd_12_2::compute_pressure");
}

// 16TH SPATIAL ORDER
void FDM_2D_ac_iso_2nd_16_2::compute_pressure()
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_2D_ac_iso_2nd_16_2::compute_pressure");

#include <FDM_2D_ac_iso_2nd_compute_pressure.inc>

	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D_ac_iso_2nd_16_2::compute_pressure");
}

} // namespace django
