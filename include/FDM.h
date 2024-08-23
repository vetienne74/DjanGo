#ifndef DJANGO_FDM_H_
#define DJANGO_FDM_H_

#include "scheme.h"

#include "type_def.h"
#include "snapshot.h"

namespace django {

//------------------------------------------------------------------------------------
// Constants definitions
//------------------------------------------------------------------------------------

//
// computation of FISRT order derivative on STAGGERED grid
//
// cf. B. Fornberg "Generation of finite-difference formulas on arbitrarily spaced grids"
//     in Mathematics of computation, vol 51, nb 184, october 1988, pages 699-706
//

// 2nd order in space / 2nd order in time
const Myint   FDM_O1_2_2_LSTENCIL = 1 ;
const Myint   FDM_O1_2_2_A1       = 1. ;
const Myfloat FDM_O1_2_2_CFL      = 1.000 ;
const Myint   NB_OP_O1_2          = 2 ;

// 4th order in space / 2nd order in time
const Myint   FDM_O1_4_2_LSTENCIL = 2 ;
const Myfloat FDM_O1_4_2_A1       = 1.125 ;
const Myfloat FDM_O1_4_2_A2       = -1./24. ;
const Myfloat FDM_O1_4_2_CFL      = 0.857 ;
const Myint   NB_OP_O1_4          = 6 ;

// 8th order in space / 2nd order in time
const Myint   FDM_O1_8_2_LSTENCIL = 4 ;
const Myfloat FDM_O1_8_2_A1       = 1225./1024. ;
const Myfloat FDM_O1_8_2_A2       = -245./3072. ;
const Myfloat FDM_O1_8_2_A3       = 49./5120. ;
const Myfloat FDM_O1_8_2_A4       = -5./7168. ;
const Myfloat FDM_O1_8_2_CFL      = 0.777 ;
const Myint   NB_OP_O1_8          = 12 ;

// 12th order in space / 2nd order in time
const Myint   FDM_O1_12_2_LSTENCIL = 6 ;
const Myfloat FDM_O1_12_2_A1       = 1225./1003. ;
const Myfloat FDM_O1_12_2_A2       = -338./3487. ;
const Myfloat FDM_O1_12_2_A3       = 35./2006. ;
const Myfloat FDM_O1_12_2_A4       = -127./42800. ;
const Myfloat FDM_O1_12_2_A5       = 19./52924. ;
const Myfloat FDM_O1_12_2_A6       = -6./274627. ;
const Myfloat FDM_O1_12_2_CFL      = 0.746 ;
const Myint   NB_OP_O1_12          = 18 ;

// 16th order in space / 2nd order in time
const Myint   FDM_O1_16_2_LSTENCIL = 8 ;
const Myfloat FDM_O1_16_2_A1       = 543./440. ;
const Myfloat FDM_O1_16_2_A2       = -85./797. ;
const Myfloat FDM_O1_16_2_A3       = 83./3603. ;
const Myfloat FDM_O1_16_2_A4       = -192./35939. ;
const Myfloat FDM_O1_16_2_A5       = 70./64979. ;
const Myfloat FDM_O1_16_2_A6       = -15./90134. ;
const Myfloat FDM_O1_16_2_A7       = 2./117497.;
const Myfloat FDM_O1_16_2_A8       = -0.000000852346420 ;
const Myfloat FDM_O1_16_2_CFL      = 0.729 ;
const Myint   NB_OP_O1_16          = 24 ;

//
// computation of SECOND order derivative on FULL grid
//
// cf. B. Fornberg "Generation of finite-difference formulas on arbitrarily spaced grids"
//     in Mathematics of computation, vol 51, nb 184, october 1988, pages 699-706
//

// 2nd order in space / 2nd order in time
const Myint   FDM_O2_2_2_LSTENCIL = 1 ;
const Myfloat FDM_O2_2_2_CFL      = FDM_O1_2_2_CFL ;
const Myint   NB_OP_O2_2          = 4 ;

// 4th order in space / 2nd order in time
const Myint   FDM_O2_4_2_LSTENCIL = 2 ;
const Myfloat FDM_O2_4_2_A0       = -5./2. ;
const Myfloat FDM_O2_4_2_A1       = 4./3. ;
const Myfloat FDM_O2_4_2_A2       = -1./12. ;
const Myfloat FDM_O2_4_2_CFL      = FDM_O1_4_2_CFL ;
const Myint   NB_OP_O2_4          = 8 ;

// 8th order in space / 2nd order in time
const Myint   FDM_O2_8_2_LSTENCIL = 4 ;
const Myfloat FDM_O2_8_2_A0       = -205./72. ;
const Myfloat FDM_O2_8_2_A1       = 8./5. ;
const Myfloat FDM_O2_8_2_A2       = -1/5. ;
const Myfloat FDM_O2_8_2_A3       = 8./315. ;
const Myfloat FDM_O2_8_2_A4       = -1/560. ;
const Myfloat FDM_O2_8_2_CFL      = FDM_O1_8_2_CFL ;
const Myint   NB_OP_O2_8          = 14 ;

// 12th order in space / 2nd order in time
const Myint   FDM_O2_12_2_LSTENCIL = 6 ;
const Myfloat FDM_O2_12_2_A0       = -2598./871. ;
const Myfloat FDM_O2_12_2_A1       = 12./7. ;
const Myfloat FDM_O2_12_2_A2       = -15./56. ;
const Myfloat FDM_O2_12_2_A3       = 10./189. ;
const Myfloat FDM_O2_12_2_A4       = -1./112. ;
const Myfloat FDM_O2_12_2_A5       = 2./1925. ;
const Myfloat FDM_O2_12_2_A6       = -1./16632. ;
const Myfloat FDM_O2_12_2_CFL      = FDM_O1_12_2_CFL ;
const Myint   NB_OP_O2_12          = 20 ;

// 16th order in space / 2nd order in time
const Myint   FDM_O2_16_2_LSTENCIL = 8 ;
const Myfloat FDM_O2_16_2_A0       = -1671./547. ;
const Myfloat FDM_O2_16_2_A1       = 16./9. ;
const Myfloat FDM_O2_16_2_A2       = -14./45. ;
const Myfloat FDM_O2_16_2_A3       = 112./1485. ;
const Myfloat FDM_O2_16_2_A4       = -7./396. ;
const Myfloat FDM_O2_16_2_A5       = 65./18673. ;
const Myfloat FDM_O2_16_2_A6       = -2./3861. ;
const Myfloat FDM_O2_16_2_A7       = 5./98536. ;
const Myfloat FDM_O2_16_2_A8       = -1./411840. ;
const Myfloat FDM_O2_16_2_CFL      = FDM_O1_16_2_CFL ;
const Myint   NB_OP_O2_16          = 26 ;


//------------------------------------------------------------------------------------

class FDM: public Scheme

{
public:

	FDM(void) ;

	// initialize
	virtual Rtn_code initialize(void) ;

	// finalize
	virtual Rtn_code finalize(void) ;

	// print FDM info
	virtual Rtn_code info(void) ;

protected:

	// space order
	Myint space_order ;

	// time order
	Myint time_order ;

	// width of the FDM stencil (half of the number of points used for derivative computation)
	Myint lstencil ;

	// total # grid points
	Myint npoint ;

	// # grid points in layers
	Myint npoint_lay ;

	// # grid points in medium
	Myint npoint_med ;

	// locate src and rec in the grid
	virtual Rtn_code locate_src_and_rec_in_grid(Acquisition*) = 0 ;

	// locate snapshot pixel in the grid
	virtual Rtn_code locate_pixel_in_grid(Snapshot*) = 0 ;

} ;

} // namespace django

#endif
