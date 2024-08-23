#ifndef DJANGO_MODEL_H_
#define DJANGO_MODEL_H_

#include "grid_1D_float.h"
#include "variable.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Model

{
public:

	// constructor
	Model() ;
	Model(Space_dim, Model_type) ;
	Model(const Model&, const Model&, Myfloat) ;

	// destructor
	~Model() ;

	// get access to array
	Variable* get_parameter(Var_type) ;

	// register parameter
	Variable* register_parameter(Var_type, string) ;

	// add perturbation on one grid point
	virtual Rtn_code add_perturbation(Myint, Myint, Myint, Myfloat) ;

	// read model
	virtual Rtn_code read(void) ;

	// initialize model from xml
	Rtn_code initialize(void) ;

	// write model
	virtual Rtn_code write(void) ;
	virtual Rtn_code write(Myint) ;

	// print info
	Rtn_code info(void) ;

	// reset model
	virtual Rtn_code reset(void) ;

	// setters
	Rtn_code set_size(Myint, Myint, Myint) ;
	Rtn_code set_sampling(Myfloat, Myfloat, Myfloat) ;
	Rtn_code set_sampling_random(Myfloat, Myfloat, Myfloat) ;

	// getters
	Space_dim get_dim(void) ;
	Myint get_nx(void) ;
	Myint get_ny(void) ;
	Myint get_nz(void) ;
	Myfloat get_dx(void) ;
	Myfloat get_dy(void) ;
	Myfloat get_dz(void) ;
	Model_type get_type(void) ;
	Model_sub_type get_sub_type(void) ;
	Grid_1D_float* get_xcoord(void) ;
	Grid_1D_float* get_zcoord(void) ;
	Grid_1D_float* get_ycoord(void) ;

protected:

	// model dimension
	Space_dim dim ;

	// model type
	Model_type type ;

	// model sub-type
	Model_sub_type sub_type ;

	// size
	Myint nx, ny, nz ;

	// sampling
	Myfloat dx, dy, dz ;

	// sampling random
	Myfloat dxrandom, dyrandom, dzrandom ;

	// coordinates of grid points
	// along x,y and z directions
	Grid_1D_float *xcoord ;
	Grid_1D_float *ycoord ;
	Grid_1D_float *zcoord ;

} ;

} // namespace django

#endif
