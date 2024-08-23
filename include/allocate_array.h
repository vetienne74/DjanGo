#ifndef DJANGO_ALLOCATE_ARRAY_H_
#define DJANGO_ALLOCATE_ARRAY_H_

#include <cassert>
#include <iostream>

#include "global.h"
#include "output_report.h"
#include "type_def.h"

namespace django {

template<class T> T* allocate_array(int mm) ;
template<class T> T** allocate_array(int mm, int nn) ;
template<class T> T*** allocate_array(int mm, int nn, int oo) ;

//-------------------------------------------------------------------------------------------------------

template<class T> T* allocate_array(int mm)
{
	print_debug(ALL, MID_DEBUG, "IN allocate_array 1D");
	print_debug(ALL, MID_DEBUG, "n1", mm);

	if (mm < 0)
	{
		print_error(" Try to allocate array with negative dim") ;
		return(NULL) ;
	}

	T* A = new T[mm] ;
	assert(A != NULL) ;

	// memory management
	Myint64 val = mm * sizeof(T) ;
	print_debug(ALL, MID_DEBUG, "alloc.   1D array", val);
	current_mem += val ;
	print_debug(ALL, MID_DEBUG, "current mem", current_mem);
	if (current_mem > max_mem) max_mem = current_mem ;

	print_debug(ALL, MID_DEBUG, "OUT allocate_array 1D");
	return(A) ;
} 

//-------------------------------------------------------------------------------------------------------

template<class T> void deallocate_array(T* A, int mm)
{
	print_debug(ALL, MID_DEBUG, "IN deallocate_array 1D");
	print_debug(ALL, MID_DEBUG, "n1", mm);

	if (mm < 0)
	{
		print_error(" Try to deallocate array with negative dim") ;
	}

	if (A != NULL)
	{
		delete[] A ;

		// memory management
		Myint64 val = mm * sizeof(T) ;
		print_debug(ALL, MID_DEBUG, "dealloc. 1D array", val);
		current_mem -= val ;
		print_debug(ALL, MID_DEBUG, "current mem", current_mem);
		if (current_mem > max_mem) max_mem = current_mem ;
	}

	print_debug(ALL, MID_DEBUG, "OUT deallocate_array 1D");
} 

//-------------------------------------------------------------------------------------------------------

template<class T> T** allocate_array(int mm, int nn)
{
	print_debug(ALL, MID_DEBUG, "IN allocate_array 2D");
	print_debug(ALL, MID_DEBUG, "n1", nn);
	print_debug(ALL, MID_DEBUG, "n2", mm);

	if ((mm < 0) || (nn < 0))
	{
		print_error(" Try to allocate array with negative dim") ;
		return(NULL) ;
	}

	T** A = new T*[mm] ;
	assert(A != NULL) ;
	T*  B = new T[mm*nn] ;
	assert(B != NULL) ;
	for (Myint ii=0; ii<mm; ii++)
	{
		A[ii] = B + ii * nn ;
	}

	// memory management
	Myint64 val = mm*nn * sizeof(T) + mm * sizeof(T*) ;
	print_debug(ALL, MID_DEBUG, "alloc.   2D array", val);
	current_mem += val ;
	print_debug(ALL, MID_DEBUG, "current mem", current_mem);
	if (current_mem > max_mem) max_mem = current_mem ;

	print_debug(ALL, MID_DEBUG, "OUT allocate_array 2D");
	return(A) ;
} 

//-------------------------------------------------------------------------------------------------------

template<class T> void deallocate_array(T** A, int mm, int nn)
{
	print_debug(ALL, MID_DEBUG, "IN deallocate_array 2D");
	print_debug(ALL, MID_DEBUG, "n1", nn);
	print_debug(ALL, MID_DEBUG, "n2", mm);

	if ((mm < 0) || (nn < 0))
	{
		print_error(" Try to deallocate array with negative dim") ;
	}

	delete A[0] ;
	delete A ;

	// memory management
	Myint64 val = mm*nn * sizeof(T) + mm * sizeof(T*) ;
	print_debug(ALL, MID_DEBUG, "dealloc. 2D array", val);
	current_mem -= val ;
	print_debug(ALL, MID_DEBUG, "current mem", current_mem);
	if (current_mem > max_mem) max_mem = current_mem ;

	print_debug(ALL, MID_DEBUG, "OUT deallocate_array 2D");
} 

//-------------------------------------------------------------------------------------------------------

template<class T> T*** allocate_array(int mm, int nn, int oo)
{
	print_debug(ALL, MID_DEBUG, "IN allocate_array 3D");
	print_debug(ALL, MID_DEBUG, "n1", oo);
	print_debug(ALL, MID_DEBUG, "n2", nn);
	print_debug(ALL, MID_DEBUG, "n3", mm);

	if ((mm < 0) || (nn < 0) || (oo < 0))
	{
		print_error(" Try to allocate array with negative dim") ;
		return(NULL) ;
	}

	T*** A = new T**[mm] ;
	assert(A != NULL) ;
	T**  B = new T*[mm*nn] ;
	assert(B != NULL) ;
	T* C = new T[mm*nn*oo] ;
	assert(C != NULL) ;
	for (Myint ii=0; ii<mm; ii++)
	{
		for (Myint jj=0; jj<nn; jj++)
		{
			B[nn*ii+jj] = C + (nn*ii+jj)*oo ;
		}
		A[ii] = B + nn*ii ;
	}

	// memory management
	Myint64 val = mm*nn*oo * sizeof(T) + mm*nn * sizeof(T*) + mm * sizeof(T**) ;
	print_debug(ALL, MID_DEBUG, "alloc.   3D array", val);
	current_mem += val ;
	print_debug(ALL, MID_DEBUG, "current mem", current_mem);
	if (current_mem > max_mem) max_mem = current_mem ;

	print_debug(ALL, MID_DEBUG, "OUT allocate_array 3D");
	return(A) ;
} 

//-------------------------------------------------------------------------------------------------------

template<class T> void deallocate_array(T*** A, int mm, int nn, int oo)
{
	print_debug(ALL, MID_DEBUG, "IN deallocate_array 3D");
	print_debug(ALL, MID_DEBUG, "n1", oo);
	print_debug(ALL, MID_DEBUG, "n2", nn);
	print_debug(ALL, MID_DEBUG, "n3", mm);

	if ((mm < 0) || (nn < 0) || (oo < 0))
	{
		print_error(" Try to deallocate array with negative dim") ;
	}

	delete A[0][0] ;
	delete A[0] ;
	delete A ;

	// memory management
	Myint64 val = mm*nn*oo * sizeof(T) + mm*nn * sizeof(T*) + mm * sizeof(T**) ;
	print_debug(ALL, MID_DEBUG, "dealloc. 3D array", val);
	current_mem -= val ;
	print_debug(ALL, MID_DEBUG, "current mem", current_mem);
	if (current_mem > max_mem) max_mem = current_mem ;

	print_debug(ALL, MID_DEBUG, "OUT deallocate_array 3D");
} 

} // namespace django

#endif
