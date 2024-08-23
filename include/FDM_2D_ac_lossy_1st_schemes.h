#ifndef DJANGO_FDM_2D_AC_LOSSY_1ST_SCHEMES_H_
#define DJANGO_FDM_2D_AC_LOSSY_1ST_SCHEMES_H_

#include "FDM_2D_ac_lossy_1st.h"

#include "type_def.h"

namespace django {

// 2ND SPATIAL ORDER
class FDM_2D_ac_lossy_1st_2_2: public FDM_2D_ac_lossy_1st
{

public:

	// constructor
	FDM_2D_ac_lossy_1st_2_2(void) ;

protected:

	// compute pressure component
	void compute_pressure(void) ;

	// compute velocity component
	void compute_velocity(void) ;

	// FDM operator
	inline Myfloat D_Z(Myfloat** U, Myint iz, Myint ix) ;
	inline Myfloat D_X(Myfloat** U, Myint iz, Myint ix) ;

} ;

// 4TH SPATIAL ORDER
class FDM_2D_ac_lossy_1st_4_2: public FDM_2D_ac_lossy_1st
{

public:

	// constructor
	FDM_2D_ac_lossy_1st_4_2(void) ;

protected:

	// compute pressure component
	void compute_pressure(void) ;

	// compute velocity component
	void compute_velocity(void) ;

	// FDM operator
	inline Myfloat D_Z(Myfloat** U, Myint iz, Myint ix) ;
	inline Myfloat D_X(Myfloat** U, Myint iz, Myint ix) ;

} ;

// 8TH SPATIAL ORDER
class FDM_2D_ac_lossy_1st_8_2: public FDM_2D_ac_lossy_1st
{

public:

	// constructor
	FDM_2D_ac_lossy_1st_8_2(void) ;

protected:

	// compute pressure component
	void compute_pressure(void) ;

	// compute velocity component
	void compute_velocity(void) ;

	// FDM operator
	inline Myfloat D_Z(Myfloat** U, Myint iz, Myint ix) ;
	inline Myfloat D_X(Myfloat** U, Myint iz, Myint ix) ;

} ;

// 12TH SPATIAL ORDER
class FDM_2D_ac_lossy_1st_12_2: public FDM_2D_ac_lossy_1st
{

public:

	// constructor
	FDM_2D_ac_lossy_1st_12_2(void) ;

protected:

	// compute pressure component
	void compute_pressure(void) ;

	// compute velocity component
	void compute_velocity(void) ;

	// FDM operator
	inline Myfloat D_Z(Myfloat** U, Myint iz, Myint ix) ;
	inline Myfloat D_X(Myfloat** U, Myint iz, Myint ix) ;

} ;

// 16TH SPATIAL ORDER
class FDM_2D_ac_lossy_1st_16_2: public FDM_2D_ac_lossy_1st
{

public:

	// constructor
	FDM_2D_ac_lossy_1st_16_2(void) ;

protected:

	// compute pressure component
	void compute_pressure(void) ;

	// compute velocity component
	void compute_velocity(void) ;

	// FDM operator
	inline Myfloat D_Z(Myfloat** U, Myint iz, Myint ix) ;
	inline Myfloat D_X(Myfloat** U, Myint iz, Myint ix) ;

} ;

} // namespace django

#endif
