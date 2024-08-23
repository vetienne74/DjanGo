#ifndef DJANGO_FDM_2D_AC_ISO_2ND_SCHEMES_H_
#define DJANGO_FDM_2D_AC_ISO_2ND_SCHEMES_H_

#include "FDM_2D_ac_iso_2nd.h"

#include "type_def.h"

namespace django {

// 2ND SPATIAL ORDER
class FDM_2D_ac_iso_2nd_2_2: public FDM_2D_ac_iso_2nd
{

public:

	// constructor
	FDM_2D_ac_iso_2nd_2_2(void) ;

protected:

	// compute pressure component
	void compute_pressure(void) ;

	// FDM operator
	inline Myfloat D_Z(Myfloat** U, Myint iz, Myint ix) ;
	inline Myfloat D2_Z(Myfloat** U, Myint iz, Myint ix) ;
	inline Myfloat D_X(Myfloat** U, Myint iz, Myint ix) ;
	inline Myfloat D2_X(Myfloat** U, Myint iz, Myint ix) ;

} ;

// 4TH SPATIAL ORDER
class FDM_2D_ac_iso_2nd_4_2: public FDM_2D_ac_iso_2nd
{

public:

	// constructor
	FDM_2D_ac_iso_2nd_4_2(void) ;

protected:

	// compute pressure component
	void compute_pressure(void) ;

	// FDM operator
	inline Myfloat D_Z(Myfloat** U, Myint iz, Myint ix) ;
	inline Myfloat D2_Z(Myfloat** U, Myint iz, Myint ix) ;
	inline Myfloat D_X(Myfloat** U, Myint iz, Myint ix) ;
	inline Myfloat D2_X(Myfloat** U, Myint iz, Myint ix) ;

} ;

// 8TH SPATIAL ORDER
class FDM_2D_ac_iso_2nd_8_2: public FDM_2D_ac_iso_2nd
{

public:

	// constructor
	FDM_2D_ac_iso_2nd_8_2(void) ;

protected:

	// compute pressure component
	void compute_pressure(void) ;

	// FDM operator
	inline Myfloat D_Z(Myfloat** U, Myint iz, Myint ix) ;
	inline Myfloat D2_Z(Myfloat** U, Myint iz, Myint ix) ;
	inline Myfloat D_X(Myfloat** U, Myint iz, Myint ix) ;
	inline Myfloat D2_X(Myfloat** U, Myint iz, Myint ix) ;

} ;

// 12TH SPATIAL ORDER
class FDM_2D_ac_iso_2nd_12_2: public FDM_2D_ac_iso_2nd
{

public:

	// constructor
	FDM_2D_ac_iso_2nd_12_2(void) ;

protected:

	// compute pressure component
	void compute_pressure(void) ;

	// FDM operator
	inline Myfloat D_Z(Myfloat** U, Myint iz, Myint ix) ;
	inline Myfloat D2_Z(Myfloat** U, Myint iz, Myint ix) ;
	inline Myfloat D_X(Myfloat** U, Myint iz, Myint ix) ;
	inline Myfloat D2_X(Myfloat** U, Myint iz, Myint ix) ;

} ;

// 16TH SPATIAL ORDER
class FDM_2D_ac_iso_2nd_16_2: public FDM_2D_ac_iso_2nd
{

public:

	// constructor
	FDM_2D_ac_iso_2nd_16_2(void) ;

protected:

	// compute pressure component
	void compute_pressure(void) ;

	// FDM operator
	inline Myfloat D_Z(Myfloat** U, Myint iz, Myint ix) ;
	inline Myfloat D2_Z(Myfloat** U, Myint iz, Myint ix) ;
	inline Myfloat D_X(Myfloat** U, Myint iz, Myint ix) ;
	inline Myfloat D2_X(Myfloat** U, Myint iz, Myint ix) ;

} ;

} // namespace django

#endif
