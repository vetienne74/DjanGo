#ifndef DJANGO_SCHEME_H_
#define DJANGO_SCHEME_H_

#include <vector>

#include "acquisition.h"
#include "boundary.h"
#include "constant.h"
#include "data.h"
#include "freq_group.h"
#include "grid.h"
#include "model.h"
#include "snapshot.h"
#include "type_def.h"
#include "variable.h"

namespace django {

//------------------------------------------------------------------------------------
// Constant definitions
//------------------------------------------------------------------------------------

const Myfloat MAX_EIGEN_VAL    = 10.0 ;
const Myfloat CPML_FREQ        = 5. ;
const Myfloat CPML_POW         = 2. ;
const Myfloat RAND_BOUND_COEF1 = 0.75 ;
const Myfloat SPONGE_COEF      = 0.99 ;

const Myfloat TEMP_VP_CONST    = 4000. ;
const Myfloat TEMP_RHO_CONST   = 1. ;

//------------------------------------------------------------------------------------

class Scheme

{
public:

	// constructor
	Scheme(void) ;

	// initialize
	Rtn_code initialize(void) ;
	virtual Rtn_code initialize(Model*, Myfloat fmax) = 0 ;

	// finalize
	virtual Rtn_code finalize(void) ;

	// print scheme info
	virtual Rtn_code info(void) ;

	// setters
	Rtn_code set_tmax(Myfloat) ;
	Rtn_code set_dt_out(Myfloat) ;
	Rtn_code set_dt(Myfloat dt) ;
	Rtn_code set_energy_param(bool) ;
	Rtn_code set_ratio_cfl(Myfloat ratio_cfl) ;
	Rtn_code set_nt_and_dt(Myfloat optimal_dt) ;
	Rtn_code set_src_param(Src_type src_type_in, Src_stype src_stype_in, Myfloat src_sigma) ;
	Rtn_code set_snapshot_param(Myfloat, Myfloat, Myfloat) ;

	// getters
	Rtn_code get_src_time_function(void) ;
	Myfloat get_dt(void) ;
	Myint get_nt(void) ;

	// solve forward problem for one shot
	virtual Rtn_code solve_current_shot(Acquisition*, Data*, Wavefield_type, Snapshot*) = 0 ;

	// display statistics, computation time and others
	Rtn_code display_stat(void) ;

	// write statistics on disk, computation time and others
	Rtn_code log_stat(void) ;

	// create boundary
	Rtn_code create_boundary(Edge_type, Boundary_type, Myint, Myfloat) ;

	// reset scheme
	virtual Rtn_code reset(void) ;

	// write energy vs time
	Rtn_code write_energy(void) ;

	// write stat for cuurent time step
	Rtn_code write_stat_current_step(void) ;

protected:

	// numerical method
	Scheme_method method ;

	// scheme type
	Scheme_type type ;

	// dimension
	Space_dim dim ;

	// equation
	Eq_type eq_type ;

	// equation order
	Eq_order eq_order ;

	// CFL associated to the FEM stencil
	Myfloat CFL ;

	// src time function
	Myfloat* src_time_function ;

	// modelling max time
	Myfloat tmax ;

	// number of time steps for modelling
	Myint nt ;

	// number of time steps for output
	Myint nt_out ;

	// $$$ obsolete to be removed
	Freq_group* pFreq_group ;

	// modelling time step
	Myfloat dt ;

	// current time step index
	Myint it ;

	// ratio cfl
	Myfloat ratio_cfl ;

	// dynamic front
	Front_type front_type ;

	// output time step
	Myfloat dt_out ;

	// snapshot
	bool output_snapshot ;
	Myfloat tmin_snapshot ;
	Myfloat tmax_snapshot ;
	Myfloat dt_snapshot ;

	// output frequency flag
	bool output_frequency_flag ;

	// compute energy flag
	bool compute_energy_flag ;

	// kinetic, potential, total and max energy
	Myfloat energy_kin, energy_pot, energy_tot, energy_tot_max ;

	// perc_completion
	Myint perc_completion ;

	// Nb lattice update
	Myint64 ncell ;

	// Nb math operation in kernel
	Myint64 nb_op_kernel ;

	// Nb math operation in boundary condition
	Myint64 nb_op_bound ;

	// nb math operation (1st derivative)
	Myfloat nb_op_d ;

	// nb math operation (2nd derivative)
	Myfloat nb_op_d2 ;

	// computing time (kernel)
	double time_in_kernel ;

	// time spent for freq. extraction
	double time_freq_extract ;

	// BOUNDARY
	//---------

	// delete all boundary
	Rtn_code delete_all_boundary(void) ;

	// boundary conditions
	Myint nBoundary ;
	Boundary* pBoundary[MAX_BOUNDARY] ;

	// get boundary
	Boundary* get_boundary(Edge_type edge_in) ;

	// get boundary size
	Myint get_boundary_width(Edge_type edge_in) ;

	// get boundary coef
	Myfloat get_boundary_coef(Edge_type edge_in) ;

	// get boundary type
	Boundary_type get_boundary_type(Edge_type edge_in) ;

	// check if boundary type is active
	bool is_there_boundary_type(Boundary_type type_in) ;

	// display remaining computation time
	Rtn_code display_remaining_computation_time(Myint it, Myfloat tt) ;

	// free source time function
	Rtn_code free_src_time_function(void) ;

	// Source type
	Src_type src_type ;

	// Source subtype
	Src_stype src_stype ;

	// Source gaussian smoothing
	Myfloat src_sigma ;

	// Source file name
	string src_file ;

	// Write snapshot
	virtual Rtn_code write_snapshot(Variable*, Myint) ;

	// number of receivers
	Myint nrec ;

	// number of pixels for the snapshot
	Myint npixel ;

	// error in eigen mode at node positions
	Myfloat eigen_error_node ;

	// error in eigen mode at rec positions
	Myfloat eigen_error_rec_l2 ;
	Myfloat eigen_error_rec_l1 ;
	Myfloat eigen_error_rec_sum2 ;

	// number of modes
	Myint eigen_nmode ;

	// adaptivity
	Myfloat fmax ;
	Myfloat fmax0 ;
	Adapt_type adaptType ;
	vector<Myfloat> adaptTmin ;
	vector<Myfloat> adaptFmax ;
	Myint adaptTimeIdx ;

} ;

} // namespace django

#endif
