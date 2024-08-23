//-------------------------------------------------------------------------------------------------------
//
// CLASS SCHEME
//
//-------------------------------------------------------------------------------------------------------

#include "scheme.h"

#include <cstdlib>
#include <iostream>

#include "allocate_array.h"
#include "output_report.h"
#include "singleton.h"
#include "type_def.h"

using namespace std;

namespace django {

static const Myint NT_MAX = 100000000 ;

//-------------------------------------------------------------------------------------------------------

Scheme::Scheme(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::Scheme");

	method            = Singleton::Instance()->method ;
	type              = Singleton::Instance()->type ;
	dim               = Singleton::Instance()->dim ;
	eq_type           = Singleton::Instance()->eq_type ;
	eq_order          = Singleton::Instance()->eq_order ;

	tmax              = 0.0 ;
	dt                = Singleton::Instance()->dt ;
	dt_out            = 0.0 ;
	nt                = 0 ;
	nt_out            = 0 ;
	ratio_cfl         = Singleton::Instance()->ratio_cfl ;
	front_type        = Singleton::Instance()->front_type ;
	CFL               = 0.0 ;

	ncell             = 0.0 ;
	nb_op_kernel      = 0.0 ;
	nb_op_bound       = 0.0 ;
	time_in_kernel    = 0.0 ;

	nb_op_d           = 0.0 ;
	nb_op_d2          = 0.0 ;

	output_snapshot   = false ;
	tmin_snapshot     = 0.0 ;
	tmax_snapshot     = 0.0 ;
	dt_snapshot       = 0.0 ;

	output_frequency_flag = false ;
	time_freq_extract = 0.0 ;

	src_type          = NO_SRC_TYPE ;
	src_stype         = NO_SRC_STYPE ;

	src_time_function = NULL ;

	eigen_error_node     = 0.0 ;
	eigen_error_rec_l1   = 0.0 ;
	eigen_error_rec_l2   = 0.0 ;
	eigen_error_rec_sum2 = 0.0 ;

	eigen_nmode = 0 ;

	perc_completion = 10 ;

	compute_energy_flag = false ;
	energy_kin = 0.0 ;
	energy_pot = 0.0 ;
	energy_tot = 0.0 ;
	energy_tot_max = 0.0 ;

	npixel = 0 ;

	it     = 0 ;
	pFreq_group = NULL ;

	nrec = 0 ;
	src_sigma = 0.0 ;

	fmax         = 0.0 ;
	fmax0        = Singleton::Instance()->fmax0 ;
	adaptType    = Singleton::Instance()->adaptType ;
	adaptTimeIdx = 0 ;
	adaptTmin    = Singleton::Instance()->adaptTmin ;
	adaptFmax    = Singleton::Instance()->adaptFmax ;

	// initialize boundaries
	nBoundary         = 0 ;

	// add boundary to the scheme
	for (Myint iBoundary = 1; iBoundary<=Singleton::Instance()->nBoundary; iBoundary++)
	{
		Boundary_type bound_type = Singleton::Instance()->pBoundary[iBoundary-1]->get_type() ;
		Myint width = Singleton::Instance()->pBoundary[iBoundary-1]->get_width() ;
		Myfloat coef = Singleton::Instance()->pBoundary[iBoundary-1]->get_coef() ;
		Edge_type edge = Singleton::Instance()->pBoundary[iBoundary-1]->get_edge() ;
		this->create_boundary(edge, bound_type, width, coef) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme::Scheme");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Scheme::log_stat(void)
{
	print_debug(ALL, FULL_DEBUG, "IN Scheme::log_stat");

	// tmp
	Myfloat c_algo = 0.0 ;

	// log perf in ouput file
	// ofstream perf_file ;
	// perf_file.open(PERF_OUT_FILE, ios::app) ;
	// assert(perf_file.is_open());
	// perf_file << c_algo << " " << fdm_cbx << " " << fdm_cby << " " << fdm_cbz << " " << nproc_world << " " <<
	//   omp_get_max_threads() << " " << time_in_kernel << "\n" ;
	// perf_file.close() ;

	print_debug(ALL, FULL_DEBUG, "OUT Scheme::log_stat");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Scheme::display_stat(void)
{
	print_debug(ALL, FULL_DEBUG, "IN Scheme::display_stat");

	print_info(ALL, " Computation time:", (Myfloat) time_in_kernel) ;
	print_info(ALL, " Freq. extract. time:", (Myfloat) time_freq_extract) ;

	Myfloat nb_tot_op = nb_op_kernel + nb_op_bound ;
	print_info(ALL, "");
	print_info(ALL, " # Flop:\t ", nb_tot_op) ;

	Myfloat c_algo ;
	if (time_in_kernel != 0.)
	{
		c_algo = (Myfloat) (nb_tot_op/(1.e+9*time_in_kernel)) ;
	}
	else
	{
		c_algo = (Myfloat) -999 ;
	}
	print_info(ALL, " Speed (GFlop/s): ", c_algo) ;

	Myfloat lup_algo ;
	Myfloat inv_lup_algo ;
	Myfloat flop_per_cell ;
	if (time_in_kernel != 0.)
	{
		lup_algo = (Myfloat) (ncell/(1.e+9*time_in_kernel)) ;
		inv_lup_algo = 1/lup_algo ;
		flop_per_cell = nb_tot_op/ncell ;
	}
	else
	{
		lup_algo = (Myfloat) -999 ;
		inv_lup_algo = -999 ;
		flop_per_cell = -999 ;
	}
	print_info(ALL, " cell:\t\t", ncell) ;
	print_info(ALL, " Speed (Gcell/s):", lup_algo) ;
	print_info(ALL, " s/cell (ns):\t", inv_lup_algo) ;
	print_info(ALL, " avg. flop/cell:", flop_per_cell) ;

	print_debug(ALL, FULL_DEBUG, "OUT Scheme::display_stat");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Scheme::display_remaining_computation_time(Myint it, Myfloat tt)
{
	print_debug(ALL, FULL_DEBUG, "IN Scheme::display_remaing_computation_time");

	if ((perc_completion != 0) && (nt/perc_completion != 0))
	{
		if ((it == 0) || ( it%(nt/perc_completion) == 0))
		{
			Myint perc_comp   = Myfloat(it+1) / Myfloat(nt+1) * 100. ;
			print_info(MASTER, " % of completion:", perc_comp) ;

			Myfloat remain_time =  (Myfloat(nt+1) / Myfloat(it+1) * tt) - tt ;
			if (remain_time < 60.0)
			{
				print_info(MASTER, " > Remaining (s):", remain_time) ;
			}
			else
			{
				print_info(MASTER, " > Remaining (m):", Myfloat(remain_time/60.0)) ;
			}

			system("grep -i mhz /proc/cpuinfo >> system.info.django.out.ascii") ;
		}
	}

	print_debug(ALL, FULL_DEBUG, "OUT Scheme::display_remaing_computation_time");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Boundary* Scheme::get_boundary(Edge_type edge_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::get_boundary");

	Boundary* boundary_tmp = NULL ;
	for (Myint ib=1; ib<=nBoundary; ib++)
	{
		if (pBoundary[ib-1]->get_edge() == edge_in) boundary_tmp = pBoundary[ib-1] ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme::get_boundary");
	return boundary_tmp ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Scheme::set_src_param(Src_type src_type_in, Src_stype src_stype_in, Myfloat src_sigma_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::set_src_param");

	// Source type ;
	src_type = src_type_in ;

	// Source subtype ;
	src_stype = src_stype_in ;

	// Source gaussian smoothing ;
	src_sigma = src_sigma_in ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme::set_src_param");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Scheme::set_snapshot_param(Myfloat tmin, Myfloat tmax, Myfloat dt)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::set_snapshot_param");

	output_snapshot = true ;
	tmin_snapshot = tmin ;
	tmax_snapshot = tmax ;
	dt_snapshot   = dt ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme::set_snapshot_param");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code Scheme::set_tmax(Myfloat tmax_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::set_tmax");

	if (tmax_in < 0)
	{
		print_error(" Scheme::set_tmax, tmax < 0") ;
		return(RTN_CODE_KO) ;
	}
	tmax = tmax_in ;

	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::set_tmax");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code Scheme::set_dt_out(Myfloat dt_out_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::set_dt_out");

	if (dt_out_in < 0)
	{
		print_error(" Scheme::set_dt_out, dt_out < 0") ;
		return(RTN_CODE_KO) ;
	}
	dt_out = dt_out_in ;

	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::set_dt_out");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code Scheme::set_nt_and_dt(Myfloat optim_dt)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::set_nt_and_dt");

	//----------------------------------
	// check user dt against optimal dt
	//----------------------------------

	if (optim_dt <= 0.)
	{
		print_error(" Optimal dt is invalid") ;
		return(RTN_CODE_KO) ;
	}

	// user defined dt > optimal dt
	// warning and set dt to optimal dt
	if (dt > optim_dt)
	{
		print_warning(" User defined dt is too large ! \n --> Reset to optimal dt") ;
		dt = optim_dt ;
	}
	// set to optimal dt
	else if  (dt == 0.)
	{
		dt = optim_dt ;

	}
	// user defined dt < optimal dt
	else if  (dt < optim_dt)
	{
		print_warning(" User defined dt is smaller than optimal dt") ;
	}

	//--------------------------------------------------------
	// find appropriate dt such as dt_out is a multiple of dt
	//--------------------------------------------------------
	Myint decim = 1 ;
	if (dt_out > 0)
	{
		decim = ceil(dt_out / (dt*(1.0 + MY_EPSILON)) ) ;
		//cout << "dt_out " << dt_out << "\n" ;
		//cout << "dt " << dt << "\n" ;
		//cout << "ceil " << ceil((Myfloat64)dt_out / (Myfloat64)dt) << "\n" ;
		//cout << "round " << round((Myfloat64)dt_out / (Myfloat64)dt) << "\n" ;
		dt = (Myfloat64)dt_out / (Myfloat64)decim ;
	}
	print_info(MASTER, " Effective dt (s):", dt) ;

	// compute nt for modelling
	if (dt_out > 0)
	{
		// take into account nt output should reach max time
		nt_out = ceil(tmax / dt_out + 1) ;
		nt = (nt_out - 1) * decim + 1 ;
	}
	else
	{
		nt = ceil(tmax / dt) + 1 ;
	}
	print_info(MASTER, " No. time step:\t", nt) ;

	// check nt
	if (nt <= 0)
	{
		print_error(" Nt cannot be <= 0") ;
		return(RTN_CODE_KO) ;
	}
	else if (nt > NT_MAX)
	{
		print_error(" Nt larger than allowed") ;
		return(RTN_CODE_KO) ;
	}

	// compute nt for output
	if (dt_out > 0)
	{
		nt_out = (nt - 1) / decim + 1 ;
		print_info(MASTER, " Nt for output:\t", nt_out) ;
	}

	// set nt and dt in modelling class
	Singleton::Instance()->pProgram->pModelling->set_nt(nt) ;
	Singleton::Instance()->pProgram->pModelling->set_dt(dt) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme::set_nt_and_dt");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code Scheme::set_dt(Myfloat dt_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::set_dt");

	if (dt < 0.0)
	{
		print_error(" Scheme::set_dt, dt < 0.0") ;
		return(RTN_CODE_KO) ;
	}
	else
	{
		dt = dt_in ;
	}
	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme::set_dt");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code Scheme::set_energy_param(bool compute_energy_flag_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::set_energy_param");
	compute_energy_flag = compute_energy_flag_in ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme::set_energy_param");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code Scheme::set_ratio_cfl(Myfloat ratio_cfl_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::set_ratio_cfl");

	if (ratio_cfl_in <= 0.0)
	{
		print_error(" Scheme::set_ratio_cfl, ratio_cfl < 0.0") ;
		return(RTN_CODE_KO) ;
	}
	else
	{
		ratio_cfl = ratio_cfl_in ;
	}
	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme::set_ratio_cfl");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Scheme::get_src_time_function()
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::get_src_time_function");

	// src_time_function is NULL if no source function has been defined
	// this can happen for special test cases like eigen mode
	src_time_function = Singleton::Instance()->pProgram->pModelling->get_src_time_function() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme::get_src_time_function");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Scheme::get_dt()
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::get_dt");
	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme::get_dt");
	return(dt) ;
}

//-------------------------------------------------------------------------------------------------------

Myint Scheme::get_nt()
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::get_nt");
	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme::get_nt");
	return(nt) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Scheme::free_src_time_function()
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::free_src_time_function");
	deallocate_array<Myfloat>(src_time_function, nt+1) ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme::free_src_time_function");
	return(RTN_CODE_OK) ;
}


//-------------------------------------------------------------------------------------------------------
Myint Scheme::get_boundary_width(Edge_type edge_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::get_boundary_width");

	Myint width_tmp = 0 ;
	for (Myint ib=1; ib<=nBoundary; ib++)
	{
		if (pBoundary[ib-1]->get_edge() == edge_in) width_tmp = pBoundary[ib-1]->get_width() ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme::get_boundary_width");
	return width_tmp ;
}

//-------------------------------------------------------------------------------------------------------
Myfloat Scheme::get_boundary_coef(Edge_type edge_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::get_boundary_coef");

	Myfloat coef_tmp = 0 ;
	for (Myint ib=1; ib<=nBoundary; ib++)
	{
		if (pBoundary[ib-1]->get_edge() == edge_in) coef_tmp = pBoundary[ib-1]->get_coef() ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme::get_boundary_coef");
	return coef_tmp ;
}

//-------------------------------------------------------------------------------------------------------
Boundary_type Scheme::get_boundary_type(Edge_type edge_in)
{
	print_debug(ALL, FULL_DEBUG, "IN Scheme::get_boundary_type");

	Boundary_type type_tmp = NO_BOUNDARY ;
	for (Myint ib=1; ib<=nBoundary; ib++)
	{
		if (pBoundary[ib-1]->get_edge() == edge_in) type_tmp = pBoundary[ib-1]->get_type() ;
	}

	print_debug(ALL, FULL_DEBUG, "OUT Scheme::get_boundary_type");
	return type_tmp ;
}

//-------------------------------------------------------------------------------------------------------
bool Scheme::is_there_boundary_type(Boundary_type type_in)
{
	print_debug(ALL, FULL_DEBUG, "IN Scheme::is_there_boundary_type");

	for (Myint ib=1; ib<=nBoundary; ib++)
	{
		if (pBoundary[ib-1]->get_type() == type_in)
		{
			print_debug(ALL, FULL_DEBUG, "OUT Scheme::is_there_boundary_type");
			return true ;
		}
	}

	print_debug(ALL, FULL_DEBUG, "OUT Scheme::is_there_boundary_type");
	return false ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Scheme::initialize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::initialize(void)");

	// get source time function
	Rtn_code rtn_code = this->get_src_time_function() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme::initialize(void)");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Scheme::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::finalize");

	// free source function
	(*this).free_src_time_function() ;

	// test if scheme is stable for eigen mode
	if(Singleton::Instance()->pProgram->pModelling->get_case() == MODELLING_EIGEN)
	{
		cout << " EIGEN MODE L2 ERROR AT NODES " << eigen_error_node << "\n" ;
		if (eigen_error_rec_sum2 != 0)
		{
			cout << " EIGEN MODE NRMS ERROR AT RECEIVERS " << sqrt(eigen_error_rec_l2 / eigen_error_rec_sum2) << "\n" ;
		}
		else
		{
			cout << " EIGEN MODE NRMS ERROR AT RECEIVERS " << 9999 << "\n" ;
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme::finalize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Scheme::info(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::info");

	print_line2() ;
	print_info(MASTER, " SCHEME PARAMETERS") ;
	print_info(MASTER, "") ;

	// method
	switch(method)
	{
	case SCHEME_FDM:
		print_info(MASTER, " Method\t\t", "FDM") ;
		break ;
	case SCHEME_FEM:
		print_info(MASTER, " Method\t\t", "FEM") ;
		break ;
	default:
		print_error(" Invalid numerical method", method) ;
		return(RTN_CODE_KO) ;
	}

	// type
	switch(type)
	{
	case SCHEME_STAGGERED:
		print_info(MASTER, " Type\t\t", "STAGGERED") ;
		break ;
	case SCHEME_CGM:
		print_info(MASTER, " Type\t\t", "CONTINUOUS") ;
		break ;
	case SCHEME_DGM:
		print_info(MASTER, " Type\t\t", "DISCONTINUOUS") ;
		break ;
	case SCHEME_MGM:
		print_info(MASTER, " Type\t\t", "MIXED D/C GALERKIN") ;
		break ;
	default:
		print_error(" Invalid numerical type", type) ;
		return(RTN_CODE_KO) ;
	}

	// dim
	switch(dim)
	{
	case ONE:
		print_info(MASTER, " Dimension\t", "1D") ;
		break ;
	case TWO:
		print_info(MASTER, " Dimension\t", "2D") ;
		break ;
	case THREE:
		print_info(MASTER, " Dimension\t", "3D") ;
		break ;
	default:
		print_error(" Invalid dimension", dim) ;
		return(RTN_CODE_KO) ;
	}

	// equation type
	switch(eq_type)
	{
	case ACOUSTIC:
		print_info(MASTER, " Equation\t", "ACOUSTIC") ;
		break ;
	case ELASTIC:
		print_info(MASTER, " Equation\t", "ELASTIC") ;
		break ;
	case AC_LOSSY:
		print_info(MASTER, " Equation\t", "ACOUSTIC LOSSY") ;
		break ;
	default:
		print_error(" Invalid equation", eq_type) ;
		return(RTN_CODE_KO) ;
	}

	// equation order
	switch(eq_order)
	{
	case ORDER_1ST:
		print_info(MASTER, " Equation\t", "1ST ORDER") ;
		break ;
	case ORDER_2ND:
		print_info(MASTER, " Equation\t", "2ND ORDER") ;
		break ;
	default:
		print_error(" Invalid equation order", eq_order) ;
		return(RTN_CODE_KO) ;
	}

	// dynamic
	switch(front_type)
	{
	case FRONT_STATIC:
		print_info(MASTER, " Front\t\t", "STATIC") ;
		break ;
	case FRONT_DYN_VMAX:
		print_info(MASTER, " Front\t\t", "DYNAMIC VMAX") ;
		break ;
	default:
		print_error(" Invalid front type", front_type) ;
		return(RTN_CODE_KO) ;
	}

	// print boundary info
	print_info(MASTER, "") ;
	print_info(MASTER, " # Boundaries\t", nBoundary) ;
	for (Myint ib = 0; ib < nBoundary; ib++)
	{
		print_info(MASTER, "") ;
		pBoundary[ib]->info() ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme::info");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Scheme::create_boundary(Edge_type edge, Boundary_type type, Myint width, Myfloat coef)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::create_boundary");

	if ((nBoundary < 0) || (nBoundary >= MAX_BOUNDARY))
	{
		print_error(" Invalid nBoundary", nBoundary) ;
		return(RTN_CODE_KO) ;
	}

	pBoundary[nBoundary] = new Boundary(edge, type, width, coef) ;

	nBoundary++ ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme::create_boundary");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Scheme::delete_all_boundary(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::delete_all_boundary");
	for (Myint ibound=1; ibound <= nBoundary; ibound++)
	{
		delete(pBoundary[ibound-1]) ;
	}
	nBoundary = 0 ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme::delete_all_boundary");
	return(RTN_CODE_OK);
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Scheme::write_snapshot(Variable* pVar, Myint it)
{
	print_debug(ALL, FULL_DEBUG, "IN Scheme::write_snapshot");

	// snapshot disabled
	if (!output_snapshot)
	{
		return(RTN_CODE_OK) ;
	}

	// write snapshot
	Myfloat tcur = it * dt ;

	if ( (tcur >= (tmin_snapshot-dt)) && (tcur <= (tmax_snapshot+dt)) )
	{
		Myfloat dtt = tcur - tmin_snapshot ;
		Myint   itt = round(dtt / dt_snapshot) ;

		if (abs(itt*dt_snapshot - dtt) < dt/2.)
		{
			print_debug(ALL, LIGHT_DEBUG, " * snapshot taken (s)", tcur) ;
			string filename = pVar->get_name() + TIME_SNAPSHOT_OUT_FILE ;
			Rtn_code rtn_code = pVar->write(filename) ;
			if (rtn_code != RTN_CODE_OK) return (rtn_code) ;
		}
	}
	print_debug(ALL, FULL_DEBUG, "OUT Scheme::write_snapshot");
	return(RTN_CODE_OK);
}

//-------------------------------------------------------------------------------------------------------
Rtn_code Scheme::reset(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::reset");

	// reset counters
	ncell         = 0.0 ;
	nb_op_kernel   = 0.0 ;
	nb_op_bound    = 0.0 ;
	time_in_kernel = 0.0 ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme::reset");
	return(RTN_CODE_OK) ;

} ;

//-------------------------------------------------------------------------------------------------------
Rtn_code Scheme::write_energy(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::write_energy");

	if (compute_energy_flag)
	{
		// check if ouput is required at the current time step
		Myint decim = round(dt_out / dt) ;

		if (it%decim != 0)
		{
			print_debug(ALL, FULL_DEBUG, "OUT Scheme::write_energy");
			return(RTN_CODE_OK) ;
		}

		// open, write and close energy
		ofstream pFile ;
		pFile.open(ENERGY_TIME_REC_OUT_FILE, ios::app | ios::out) ;
		assert(pFile.is_open());
		pFile << it * dt << " " << energy_kin << " " << energy_pot << " " << energy_tot << "\n" ;
		pFile.close() ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme::write_energy");
	return(RTN_CODE_OK) ;

}

//-------------------------------------------------------------------------------------------------------
Rtn_code Scheme::write_stat_current_step(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Scheme::write_stat_current_step");

	// check if ouput is required at the current time step
	Myint decim = round(dt_out / dt) ;

	if (it%decim != 0)
	{
		print_debug(ALL, FULL_DEBUG, "OUT Scheme::write_stat_current_step");
		return(RTN_CODE_OK) ;
	}

	// log stat in ouput file
	ofstream stat_file ;
	stat_file.open(STAT_TIME_STEP_OUT_FILE, ios::app) ;
	assert(stat_file.is_open());
	stat_file << it*dt << " " << nb_op_kernel << " " << time_in_kernel << "\n" ;
	stat_file.close() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Scheme::write_stat_current_step");
	return(RTN_CODE_OK) ;
}

} // namespace django
