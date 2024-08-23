#ifndef DJANGO_VARIABLE_H_
#define DJANGO_VARIABLE_H_

#include <string>

#include "grid.h"
#include "grid_2D_float.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Variable

{
public:

	// constructor
	Variable(const char*) ;
	Variable(Var_type, string) ;

	// destructor
	~Variable(void) ;

	// allocate 1d grid
	//template <class Type> Rtn_code allocate_grid(Myint, Myfloat) ;
	Rtn_code allocate_grid(Myint, Myfloat) ;
	Rtn_code allocate_grid(Myint, Myint, Myfloat, Myfloat) ;


	// initialize grid
	Rtn_code initialize_grid(Myint, Myfloat) ;
	Rtn_code initialize_grid(Myint, Myint, Myfloat, Myfloat) ;

	// getters
	Grid* get_grid(void) ;
	Var_type get_type(void) ;
	string get_name(void) ;
	Myint get_id(void) ;

	// reset
	Rtn_code reset_grid() ;

	// setters
	Rtn_code set_type(Var_type) ;
	Rtn_code set_grid(Grid* pGrid) ;
	Rtn_code set_self_init_mode(Self_init_type) ;
	Rtn_code set_self_init_const(Myfloat) ;
	Rtn_code set_self_init_file(string) ;
	Rtn_code set_id(Myint) ;

	// display info
	void info(void) ;

	// read from disk
	Rtn_code read(void) ;

	// write on disk
	Rtn_code write(void) ;
	Rtn_code write(string) ;
	Rtn_code write_rec(Myint nrec, Myint* izrec) ;
	Rtn_code write_rec(Myint nrec, Myint* izrec, Myint* ixrec) ;
	Rtn_code flush_rec(void) ;

private:

	// variable id
	Myint id ;

	// variable type
	Var_type type ;

	// variable name
	string name ;

	// self-init mode
	Self_init_type self_init ;
	string self_init_file ;
	Myfloat self_init_const ;

	// grid
	Grid* pGrid ;

	// IO buffer
	Grid_2D_float* pIO_buffer ;
	Myint idx_buffer ;
	Myint max_buffer ;
	Myint nrec ;
	Myfloat nb_byte_write ;
	double time_in_write ;

} ;

} // namespace django

#endif
