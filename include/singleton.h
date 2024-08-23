
// this class is a singleton pattern

#ifndef DJANGO_SINGLETON_H_
#define DJANGO_SINGLETON_H_

#include <string>

#include "program.h"
#include "variable.h"
#include "type_def.h"

namespace django {

const Myint MAX_VARIABLE = 100 ;

class Singleton
{
public:

	static Singleton* Instance() ;

	// xml version
	Myfloat xml_version ;

	// cache blocking size
	Myint cbx, cby, cbz ;

	// io buffer size
	Myint iobufsize ;

	// subdomain decomposition
	Myint nsubx, nsuby, nsubz ;

	// xml config file
	string xml_config_file ;

	// dry run
	bool dryrun ;

	// objects
	Program* pProgram ;

	// equation
	Eq_type eq_type ;

	// equation order
	Eq_order eq_order ;

	// space order
	Myint space_order ;

	// time order
	Myint time_order ;

	// dt for computation
	Myfloat dt ;

	// ratio cfl
	Myfloat ratio_cfl ;

	// node type
	Node_type node_type ;

	// node distribution
	Node_distrib node_distrib ;

	// node integration
	Node_integ node_integ ;

	// flux type
	Flux_type flux_type ;

	// scheme type
	Scheme_type scheme_type ;

	// min. polynomial
	Myint pmin ;

	// max polynomial
	Myint pmax ;

	// no. element along axis in medium
	Myint nelem_x_med, nelem_z_med ;

	// physical properties within elements
	Prop_type prop ;

	// numerical method
	Scheme_method method ;

	// type of numerical method
	Scheme_type type ;

	// dimension
	Space_dim dim ;

	// front type
	Front_type front_type ;

	// model
	Model* pModel ;

	// adaptivity parameter
	Myfloat fmax0 ;
	Adapt_type adaptType ;
	vector<Myfloat> adaptTmin ;
	vector<Myfloat> adaptFmax ;

	// read django configuration parameters
	Rtn_code initialize(void) ;

	// print config info
	Rtn_code info(void) ;

	// get variable
	Variable* get_variable(Var_type) ;
	Variable* get_variable(Myint var_id) ;

	// register one variable
	Variable* register_variable(Var_type, string) ;

	// swap two variables
	Rtn_code swap_variable(Var_type, Var_type) ;

	// delete all variables
	Rtn_code delete_all_variable(void) ;

	// delete variable
	Rtn_code delete_variable(Var_type) ;
	Rtn_code delete_variable(Myint var_id) ;

	// print all variables
	Rtn_code print_all_variable(void) ;

	// boundary conditions
	Myint nBoundary ;
	Boundary* pBoundary[MAX_BOUNDARY] ;

	// create boundary
	Rtn_code create_boundary(Edge_type, Boundary_type, Myint, Myfloat) ;

	// get boundary
	Boundary* get_boundary(Edge_type edge_in) ;

	// get boundary size
	Myint get_boundary_width(Edge_type edge_in) ;

	// get boundary coef
	Myfloat get_boundary_coef(Edge_type edge_in) ;

	// get boundary type
	Boundary_type get_boundary_type(Edge_type edge_in) ;

	// delete all boundary
	Rtn_code delete_all_boundary(void) ;

private:

	// private so that it can  not be called
	Singleton(void) ;

	// copy constructor is private
	Singleton(Singleton const&){} ;

	// assignment operator is private
	Singleton& operator=(Singleton const&){} ;

	// singleton
	static Singleton* m_pInstance ;

	// variable
	Myint nVariable ;
	Variable* pVariable[MAX_VARIABLE] ;

} ;

} // namespace django

#endif
