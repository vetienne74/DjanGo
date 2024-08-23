//-------------------------------------------------------------------------------------------------------
//
// MODELLING WITH FDM IN 2D
//
// PARENT CLASS: Scheme
//   DERIVED CLASS: FDM
//     DERIVED CLASS: FDM_1D
//       DERIVED CLASS: FDM_2D
//
//-------------------------------------------------------------------------------------------------------

#include "FDM_2D.h"

#include <cassert>
#include <fstream>
#include <iostream>

#include "acquisition.h"
#include "allocate_array.h"
#include "data_std.h"
#include "grid_3D_complex.h"
#include "output_report.h"
#include "type_def.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------

FDM_2D::FDM_2D(void) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D::FDM_2D");
	dim = TWO ;
	ix_rec   = NULL ;
	iz_rec   = NULL ;
	ix_src   = NULL ;
	iz_src   = NULL ;
	weight_src = NULL ;
	ix_pixel = NULL ;
	iz_pixel = NULL ;
	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D::FDM_2D");
}

//-------------------------------------------------------------------------------------------------------

// initialize scheme geometry based on input Grid

Rtn_code FDM_2D::initialize(Grid* pGrid_in, Myfloat fmax)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D::initialize(Grid* pGrid_in)");

	// check input grid is Grid_2D_float
	Grid_2D_float *pGrid = dynamic_cast<Grid_2D_float*>(pGrid_in) ;
	if (pGrid == NULL)
	{
		print_error("IN FDM_2D::initialize --> input grid is not Grid_2D_float");
		return(RTN_CODE_KO) ;
	}

	// set grid size and sampling
	nz = pGrid->nz ;
	dz = pGrid->dz ;
	inv_dz = 1. / dz ;
	inv_dz2 = inv_dz * inv_dz ;

	nx = pGrid->nx ;
	dx = pGrid->dx ;
	inv_dx = 1. / dx ;
	inv_dx2 = inv_dx * inv_dx ;

	// compute optimal time step
	Myfloat optimal_dt = compute_optimal_time_step(pGrid) ;

	// set appropriate nt and dt
	if (set_nt_and_dt(optimal_dt) == RTN_CODE_KO) return(RTN_CODE_KO) ;

	// set boundary width
	nlayer_zBeg = get_boundary_width(ZBEG) ;
	nlayer_zEnd = get_boundary_width(ZEND) ;
	nlayer_xBeg = get_boundary_width(XBEG) ;
	nlayer_xEnd = get_boundary_width(XEND) ;

	// usefull indices
	izBeg  = 0 ;
	izBeg1 = izBeg + lstencil ;
	izBeg2 = izBeg1 + nlayer_zBeg ;
	izEnd2 = izBeg2 + nz  ;
	izEnd1 = izEnd2 + nlayer_zEnd ;
	izEnd  = izEnd1 + lstencil -1 ;

	ixBeg  = 0 ;
	ixBeg1 = ixBeg + lstencil ;
	ixBeg2 = ixBeg1 + nlayer_xBeg ;
	ixEnd2 = ixBeg2 + nx  ;
	ixEnd1 = ixEnd2 + nlayer_xEnd ;
	ixEnd  = ixEnd1 + lstencil -1 ;

	// update total and layer grid points numbers
	npoint_med = nz * nx ;
	npoint     = (nx + nlayer_xBeg + nlayer_xEnd) * (nz + nlayer_zBeg + nlayer_zEnd) ;
	npoint_lay = npoint - npoint_med ;

	// call parent initialization
	Rtn_code rtn_code = FDM::initialize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D::initialize(Grid* pGrid_in)");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_2D::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D::finalize");

	// call parent finalize
	Rtn_code rtn_code = FDM::finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D::finalize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat FDM_2D::compute_optimal_time_step(Grid_2D_float *vp)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D::compute_optimal_time_step");

	// evaluate optimal time step
	Myfloat max_vel = vp->get_max() ;
	Myfloat optim_dt = CFL / (max_vel * sqrt(1./(dx*dx) + 1./(dz*dz))) ;
	print_info(MASTER, " Optimal dt (s):", optim_dt) ;

	if (ratio_cfl != 1.0)
	{
		optim_dt *= ratio_cfl ;
		print_info(MASTER, " Ratio CFL:\t", ratio_cfl) ;
		print_info(MASTER, " Modified dt (s):", optim_dt) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D::compute_optimal_time_step");
	return(optim_dt) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_2D::source_excitation(Myfloat** time_u, Myint it, Wavefield_type wtype, Data* pData)
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_2D:source_excitation(Myfloat** time_u, Myint it, Wavefield_type wtype, Data* pData)");

	if (wtype == INCIDENT)
	{
		// apply source function in time domain
		// excitation at the source position

		for (Myint ip=0; ip<npoint_src; ip++)
		{
			if ((iz_src[ip] != NOT_FOUND) && (ix_src[ip] != NOT_FOUND))
			{
				if (it >= 0) time_u[ix_src[ip]][iz_src[ip]] += weight_src[ip] * src_time_function[it] ;
			}
		}
	}
	else if (wtype == ADJOINT)
	{

		// build excitation in time domain via an inverse DFT of the frequency adjoint terms
		// excitation at the receivers positions

		// loop over frequencies
		Mycomplex **data_adj = ((Data_std*) pData)->pr_freq_adj_rec ;
		Myfloat coef2 = TEMP_RHO_CONST * pow(TEMP_VP_CONST, 2) * (2. / (nt * dz*dx)) ;

		// taper
		//Myfloat coef3 = min((Myfloat)it / 10., 1.) ;
		//if (it < 20) cout << coef3 << "\n" ;
		for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)
		{
			Myfloat coef_tmp = 2. * PI * pFreq_group->pFreq_list[ifreq] * (nt - it -1) * dt;
			//Myfloat coef_tmp = 2. * PI * pFreq_group->pFreq_list[ifreq] * (nt - it) * dt;
			//Myfloat coef_tmp = 2. * PI * pFreq_group->pFreq_list[ifreq] * (nt - it) * float(nt-1) / float(nt) * dt;
			Mycomplex coef (cos(coef_tmp), sin(coef_tmp)) ;

			for (Myint irec=0; irec<nrec; irec++)
			{
				time_u[ix_rec[irec]][iz_rec[irec]] += -1. * (data_adj[irec][ifreq] * coef).real() * coef2 ;
				//time_u[ix_rec[irec]][iz_rec[irec]] += (data_adj[irec][ifreq] * coef).real() * coef2 * coef3 ;
			}
		}
	}

	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D::source_excitation(Myfloat** time_u, Myint it, Wavefield_type wtype, Data* pData)");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_2D::source_excitation(Myfloat** time_u, Myint it, Wavefield_type wtype, Data* pData, Myfloat src_factor)
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_2D:source_excitation(Myfloat** time_u, Myint it, Wavefield_type wtype, Data* pData, Myfloat src_factor)");

	if (wtype == INCIDENT)
	{
		// apply source function in time domain
		// excitation at the source position

		print_debug(ALL, FULL_DEBUG, "src it", it);
		print_debug(ALL, FULL_DEBUG, "src amp",  src_time_function[it] * src_factor);

		for (Myint ip=0; ip<npoint_src; ip++)
		{
			if ((iz_src[ip] != NOT_FOUND) && (ix_src[ip] != NOT_FOUND))
			{
				if (it >= 0) time_u[ix_src[ip]][iz_src[ip]] += weight_src[ip] * src_time_function[it] * src_factor ;
			}
		}
	}
	else if (wtype == ADJOINT)
	{

		// build excitation in time domain via an inverse DFT of the frequency adjoint terms
		// excitation at the receivers positions

		// loop over frequencies
		Mycomplex **data_adj = ((Data_std*) pData)->pr_freq_adj_rec ;
		Myfloat coef2 = TEMP_RHO_CONST * pow(TEMP_VP_CONST, 2) * (2. / (nt * dz*dx)) ;

		// taper
		//Myfloat coef3 = min((Myfloat)it / 10., 1.) ;
		//if (it < 20) cout << coef3 << "\n" ;
		for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)
		{
			Myfloat coef_tmp = 2. * PI * pFreq_group->pFreq_list[ifreq] * (nt - it -1) * dt;
			//Myfloat coef_tmp = 2. * PI * pFreq_group->pFreq_list[ifreq] * (nt - it) * dt;
			//Myfloat coef_tmp = 2. * PI * pFreq_group->pFreq_list[ifreq] * (nt - it) * float(nt-1) / float(nt) * dt;
			Mycomplex coef (cos(coef_tmp), sin(coef_tmp)) ;

			for (Myint irec=0; irec<nrec; irec++)
			{
				time_u[ix_rec[irec]][iz_rec[irec]] += -1. * (data_adj[irec][ifreq] * coef).real() * coef2 ;
				//time_u[ix_rec[irec]][iz_rec[irec]] += (data_adj[irec][ifreq] * coef).real() * coef2 * coef3 ;
			}
		}
	}

	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D::source_excitation(Myfloat** time_u, Myint it, Wavefield_type wtype, Data* pData, Myfloat src_factor)");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_2D::update_freq_sol_grid(Myfloat** time_u, Mycomplex*** freq_u, Myint it)
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_2D::update_freq_sol_grid");

	// freq. extraction not required
	if ((!output_frequency_flag) || (freq_u == NULL))
	{
		print_debug(ALL, FULL_DEBUG, "OUT FDM_2D::update_freq_sol_grid");
		return(RTN_CODE_OK) ;
	}

	Myint decim = 13 ;
	if ((it % decim) != 0)
	{
		return(RTN_CODE_OK) ;
	}

	Myint nz1 = 0 ;
	Myint nz2 = nz ;
	Myint nx1 = 0 ;
	Myint nx2 = nx ;

	// coef for DFT
	for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)
	{

		// // determine decim
		// cout << " pFreq_group->pFreq_list[ifreq] " << pFreq_group->pFreq_list[ifreq] << "\n" ;
		// //Myfloat dt_nyquist = 0.5 / pFreq_group->pFreq_list[ifreq] ;
		// Myfloat dt_nyquist = 0.5 / 37.5 ;
		// cout << " dt_nyquist " << dt_nyquist << "\n" ;
		// Myint   decim  = floor(dt_nyquist / dt) ;
		// cout << " decim " << decim << "\n" ;
		// if ((it % decim) != 0) break ;

		Myfloat coef_tmp = -2. * PI * pFreq_group->pFreq_list[ifreq] * it * dt ;
		//Myfloat coef_tmp = -2. * PI * pFreq_group->pFreq_list[ifreq] * it * float(nt-1)/float(nt) * dt ;

		Mycomplex coef (cos(coef_tmp)*dt*decim, sin(coef_tmp)*dt*decim) ;

		//Mycomplex **freq_u2 = freq_u[ifreq] ;

#pragma omp parallel for
		for (Myint ix = nx1; ix < nx2; ix++)
		{
#pragma ivdep
			for (Myint iz = nz1; iz < nz2; iz++)
			{
				freq_u[ifreq][ix][iz] += time_u[ix+ixBeg2][iz+izBeg2] * coef ;
				//freq_u2[ix][iz] += time_u[ix+ixBeg2][iz+izBeg2] * coef ;
			}
		}
	}

	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D::update_freq_sol_grid");
	return(RTN_CODE_OK) ;

}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_2D::update_freq_sol_rec(Myfloat** time_u, Data* pData, Myint it)
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_2D::update_freq_sol_rec");

	Mycomplex **pSol = ((Data_std*) pData)->pr_freq_sol_rec ;

	// coef for DFT
	for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)
	{
		Myfloat coef_tmp = -2. * PI * pFreq_group->pFreq_list[ifreq] * it * dt ;
		//Myfloat coef_tmp = -2. * PI * pFreq_group->pFreq_list[ifreq] * it * float(nt-1)/float(nt) * dt ;

		Mycomplex coef (cos(coef_tmp)*dt, sin(coef_tmp)*dt) ;

		for (Myint irec=0; irec<nrec; irec++)
		{
			pSol[irec][ifreq] += time_u[ix_rec[irec]][iz_rec[irec]] * coef ;
		}
	}

	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D::update_freq_sol_rec");
	return(RTN_CODE_OK) ;

}

//-------------------------------------------------------------------------------------------------------

void FDM_2D::compute_energy(Myfloat** time_u)
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_2D::compute_energy");

	if (compute_energy_flag)
	{
		// compute energy
		Myint nz1 = izBeg2 ;
		Myint nz2 = izEnd2 ;
		Myint nx1 = ixBeg2 ;
		Myint nx2 = ixEnd2 ;

		Myfloat energy = 0. ;

		for (Myint ix = nx1; ix < nx2; ix++)
		{
			for (Myint iz = nz1; iz < nz2; iz++)
			{
				energy += time_u[ix][iz] * time_u[ix][iz] ;
			}
		}

		// write energy
		write_energy() ;
	}

	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D::compute_energy");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_2D::compute_eigen_error(Myfloat** pr, Myint it, Acquisition* pAcquisition) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D::compute_eigen_error");

	// compute error only at output time steps
	Myint decim = round(dt_out / dt) ;
	if (it%decim != 0) return(RTN_CODE_OK) ;

	Myfloat ttime = dt * it ;
	Myfloat sq2 = sqrt(2.0) ;

	// compute error at the rec positions
	//-----------------------------------

	//cout << "*** ttime " << ttime << "\n" ;

	// loop on receivers
	for (Myint irec = 0; irec < nrec; irec++)
	{
		if ((ix_rec[irec] == NOT_FOUND) && (iz_rec[irec] == NOT_FOUND)) continue ;

		Myfloat znode = pAcquisition->zrec[irec] ;
		Myfloat xnode = pAcquisition->xrec[irec] ;
		Myfloat eigen_sol = -sq2 * sin(M_PI*xnode * eigen_nmode)
		* sin(M_PI*znode * eigen_nmode) * sin(sq2*M_PI*ttime * eigen_nmode) ;

		Myfloat pr_rec = pr[ ix_rec[irec] ][ iz_rec[irec] ] ;

		// update L2 norm
		eigen_error_rec_l2 += pow(eigen_sol - pr_rec, 2) ;

		// update L1 norm
		eigen_error_rec_l1 += abs(eigen_sol - pr_rec) ;

		// update sum of square of ref_sol
		eigen_error_rec_sum2 += pow(eigen_sol, 2) ;
	}

	// open, write and close
	Myfloat nb_tot_op = nb_op_kernel + nb_op_bound ;
	ofstream pFile ;
	pFile.open(EIGEN_ERROR_REC_OUT_FILE, ios::app | ios::out) ;
	assert(pFile.is_open());
	pFile << ttime << " "
			<< nb_tot_op << " "
			<< nx * nz << " "
			<< eigen_error_rec_l1 << " "
			<< sqrt(eigen_error_rec_l2) << " "
			<< sqrt(eigen_error_rec_l2 / eigen_error_rec_sum2) << "\n" ;
	pFile.close() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D::compute_eigen_error");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_2D::write_trace(Variable* pVar, Myint it)
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_2D::write_trace");

	if (dt_out != 0)
	{
		// check if ouput is required at the current time step
		Myint decim = round(dt_out / dt) ;
		string file_name = pVar->get_name() + TIME_REC_OUT_FILE ;

		if (it%decim != 0)
		{
			print_debug(ALL, FULL_DEBUG, "OUT FDM_2D::write_trace");
			return(RTN_CODE_OK) ;
		}
		else
		{
			print_debug(ALL, FULL_DEBUG, " write_trace in ", file_name) ;
			print_debug(ALL, FULL_DEBUG, " output trace (s) ", it*dt) ;
		}

		// write traces
		Rtn_code rtn_code = pVar->write_rec(nrec, iz_rec, ix_rec) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	}
	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D::write_trace");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_2D::write_snapshot(Variable* pVar, Myint it)
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_2D::write_snapshot");

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
			string file_name = pVar->get_name() + TIME_SNAPSHOT_OUT_FILE ;
			print_debug(ALL, FULL_DEBUG, " write_snapshot in ", file_name) ;
			print_info(ALL, " snapshot captured (s)", it * dt) ;
			print_info(ALL, " snapshot componenent", pVar->get_name()) ;

			// retrieve grid
			Grid_2D_float *val_grid = (Grid_2D_float*) (pVar->get_grid()) ;
			Myfloat** const val = val_grid->pArray ;

			// allocate snapshot
			Myfloat* snap = allocate_array<Myfloat>(npixel) ;

			// extract snapshot
			for (Myint ipixel=0; ipixel<npixel; ipixel++)
			{
				if ((ix_pixel[ipixel] == NOT_FOUND) || (iz_pixel[ipixel] == NOT_FOUND))
				{
					snap[ipixel] = 0.0 ;
				}
				else
				{
					snap[ipixel] = val[ix_pixel[ipixel]][iz_pixel[ipixel]] ;
				}
			}

			// write snapshot
			ofstream pFile ;
			pFile.open(file_name.c_str(), ios::binary | ios::app | ios::out) ;
			pFile.write((char*) &(snap[0]), npixel * sizeof(Myfloat)) ;
			pFile.close() ;

			// deallocate snapshot
			deallocate_array<Myfloat>(snap, npixel) ;
		}
	}
	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D::write_snapshot");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_2D::write_time_sol_snapshot(Myfloat** time_u, const char* file_name, Myint it)
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_2D::write_time_sol_snapshot");

	// snapshot disabled
	if (dt_snapshot == 0.)
	{
		return(RTN_CODE_OK) ;
	}

	// write snapshot
	else if (it % Myint(dt_snapshot / dt) == 0)
	{
		print_debug(ALL, FULL_DEBUG, " write_time_sol_snapshot in ", file_name) ;
		print_info(ALL, " snapshot captured (s)", it * dt) ;
		ofstream pFile ;
		pFile.open(file_name, ios::binary | ios::app | ios::out) ;
		write_time_sol_grid(time_u, &pFile) ;
		pFile.close() ;
	}

	print_debug(ALL, FULL_DEBUG, "OUT FDM_2D::write_time_sol_snapshot");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_2D::write_freq_sol_rec(Data* pData)
{
	print_debug(ALL, MID_DEBUG, "IN FDM_2D::write_freq_sol_rec");

	// freq. extraction not required
	if (!output_frequency_flag)
	{
		print_debug(ALL, MID_DEBUG, "OUT FDM_2D::write_freq_sol_rec");
		return(RTN_CODE_OK) ;
	}

	Mycomplex **pSol = ((Data_std*) pData)->pr_freq_sol_rec ;
	ofstream pFile ;
	pFile.open(PR_FREQ_REC_OUT_FILE, ios::binary | ios::app | ios::out) ;
	assert(pFile.is_open());
	pFile.write((char*)&pSol[0][0], nrec * pFreq_group->nb_freq * sizeof(complex<float>)) ;
	pFile.close() ;

	print_debug(ALL, MID_DEBUG, "OUT FDM_2D::write_freq_sol_rec");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_2D::write_freq_sol_grid(Mycomplex*** freq_u)
{
	print_debug(ALL, MID_DEBUG, "IN FDM_2D::write_freq_sol_grid");

	// freq. extraction not required
	if (freq_u == NULL)
	{
		print_debug(ALL, FULL_DEBUG, "OUT FDM_2D::write_freq_sol_grid");
		return(RTN_CODE_OK) ;
	}

	// freq. extraction not required
	if (freq_u == NULL)
	{
		print_debug(ALL, FULL_DEBUG, "OUT FDM_2D::write_freq_sol_grid");
		return(RTN_CODE_OK) ;
	}

	// loop over frequencies
	ofstream pFile ;
	pFile.open(PR_FREQ_GRID_OUT_FILE, ios::binary | ios::app | ios::out) ;
	assert(pFile.is_open());
	for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)
	{
		pFile.write((char*) &(freq_u[ifreq][0][0]), nz * nx * sizeof(complex<float>)) ;
	}
	pFile.close() ;

	print_debug(ALL, MID_DEBUG, "OUT FDM_2D::write_freq_sol_grid");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_2D::write_time_sol_grid(Myfloat** time_u, ofstream* pFile)
{
	print_debug(ALL, MID_DEBUG, "IN FDM_2D::write_time_sol_grid");

	Myint nz2 = izEnd+1 ;
	Myint nx2 = ixEnd+1 ;

	pFile->write((char*) &(time_u[0][0]), nz2*nx2 * sizeof(float)) ;

	print_debug(ALL, MID_DEBUG, "OUT FDM_2D::write_time_sol_grid");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_2D::free_position_arrays(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D::free_position_arrays");

	if (ix_rec != NULL)     deallocate_array<Myint>(ix_rec, nrec) ;
	if (iz_rec != NULL)     deallocate_array<Myint>(iz_rec, nrec) ;

	if (ix_src != NULL)     deallocate_array<Myint>(ix_src, npoint_src) ;
	if (iz_src != NULL)     deallocate_array<Myint>(iz_src, npoint_src) ;
	if (weight_src != NULL) deallocate_array<Myfloat>(weight_src, npoint_src) ;

	if (ix_pixel != NULL)   deallocate_array<Myint>(ix_pixel, npixel) ;
	if (iz_pixel != NULL)   deallocate_array<Myint>(iz_pixel, npixel) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D::free_position_arrays");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_2D::locate_pixel_in_grid(Snapshot *pSnapshot)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D::locate_pixel_in_grid");

	if (output_snapshot)
	{

		// retrieve number of pixel
		Myint npixel_z = 0 ;
		if (pSnapshot->z_coord != NULL) npixel_z = pSnapshot->z_coord->nz ;
		if (npixel_z <= 0)
		{
			print_error(" Error in FDM_2D::locate_pixel_in_grid, npixel_z <= 0") ;
			return(RTN_CODE_KO) ;
		}
		Myint npixel_x = 0 ;
		if (pSnapshot->x_coord != NULL) npixel_x = pSnapshot->x_coord->nz ;
		if (npixel_x <= 0)
		{
			print_error(" Error in FDM_2D::locate_pixel_in_grid, npixel_x <= 0") ;
			return(RTN_CODE_KO) ;
		}
		npixel = npixel_z * npixel_x ;

		// retrieve position of receivers
		ix_pixel = allocate_array<Myint>(npixel) ;
		iz_pixel = allocate_array<Myint>(npixel) ;

		Myint ipixel = 0 ;
		Myint npixel_not_found = 0 ;
		for (Myint izp=0; izp<npixel_z; izp++)
		{
			for (Myint ixp=0; ixp<npixel_x; ixp++)
			{
				ix_pixel[ipixel] = (Myint)(ixBeg2 + pSnapshot->x_coord->pArray[ixp] / dx) ;
				if (ix_pixel[ipixel] < ixBeg2)
				{
					ix_pixel[ipixel] = NOT_FOUND ;
				}
				else if (ix_pixel[ipixel] > ixEnd2-1)
				{
					ix_pixel[ipixel] = NOT_FOUND ;
				}

				iz_pixel[ipixel] = (Myint)(izBeg2 + pSnapshot->z_coord->pArray[izp] / dz) ;
				if (iz_pixel[ipixel] < izBeg2)
				{
					iz_pixel[ipixel] = NOT_FOUND ;
				}
				else if (iz_pixel[ipixel] > izEnd2-1)
				{
					iz_pixel[ipixel] = NOT_FOUND ;
				}

				// update number of not found
				if ((ix_pixel[ipixel] == NOT_FOUND) || (iz_pixel[ipixel] == NOT_FOUND)) npixel_not_found++ ;

				ipixel++ ;
			}
		}

		//print_info(ALL, " Number of pixel", npixel) ;
		if (npixel_not_found > 0)
		{
			print_warning(" number of pixel not found", npixel_not_found) ;
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D::locate_pixel_in_grid");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_2D::locate_src_and_rec_in_grid(Acquisition *acquisition)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D::locate_src_and_rec_in_grid");

	// retrieve points for the source excitation
	//#########################################

	npoint_src = 0 ;
	Myfloat sum_src = 0.0 ;

	// NO source smoothing
	// excitation on single node
	//==============================

	if (src_type == SRC_POINT)
	{
		iz_src = allocate_array<Myint>(1) ;
		ix_src = allocate_array<Myint>(1) ;
		weight_src = allocate_array<Myfloat>(1) ;
		npoint_src    = 1 ;
		weight_src[0] = 1.0 ;
		sum_src       = 1.0 ;

		// retrieve position of source
		// --> nearest grid point
		ix_src[0] = (Myint)(ixBeg2 + (*acquisition).xsrc / dx) ;
		if (ix_src[0] < ixBeg2)
		{
			ix_src[0] = NOT_FOUND ;
		}
		else if (ix_src[0] > ixEnd2-1)
		{
			ix_src[0] = NOT_FOUND ;
		}

		iz_src[0] = (Myint)(izBeg2 + (*acquisition).zsrc / dz) ;
		if (iz_src[0] < izBeg2)
		{
			iz_src[0] = NOT_FOUND ;
		}
		else if (iz_src[0] > izEnd2-1)
		{
			iz_src[0] = NOT_FOUND ;
		}

		if ((iz_src[0] == NOT_FOUND) || (ix_src[0] == NOT_FOUND))
		{
			print_warning(" source not in domain") ;
		}
	}

	// WITH source gaussian
	// excitation on several node
	//================================
	else if (src_type == SRC_GAUSSIAN)

	{
		const Myfloat radius = src_sigma*3 ;

		// 1st loop to get number of points
		Myint izMin = izBeg2 ;
		Myint izMax = izEnd2 ;
		Myint ixMin = ixBeg2 ;
		Myint ixMax = ixEnd2 ;

		for (Myint ix = ixMin; ix < ixMax; ix++)
		{
			for (Myint iz = izMin; iz < izMax; iz++)
			{
				Myfloat xp = (ix-ixMin) * dx ;
				Myfloat zp = (iz-izMin) * dz ;
				Myfloat dist =  sqrt( pow(acquisition->zsrc - zp, 2) +
						pow(acquisition->xsrc - xp, 2) ) ;
				if (dist <= radius)
				{
					npoint_src++ ;
				}
			}
		}

		if (npoint_src == 0)
		{
			print_warning(" Source is not in domain") ;
		}
		else
		{
			ix_src = allocate_array<Myint>(npoint_src) ;
			iz_src = allocate_array<Myint>(npoint_src) ;
			weight_src = allocate_array<Myfloat>(npoint_src) ;
			npoint_src = 0 ;

			// 2nd loop to store the points and weights

			for (Myint ix = ixMin; ix < ixMax; ix++)
			{
				for (Myint iz = izMin; iz < izMax; iz++)
				{
					Myfloat xp = (ix-ixMin) * dx ;
					Myfloat zp = (iz-izMin) * dz ;
					Myfloat dist =  sqrt( pow(acquisition->zsrc - zp, 2) +
							pow(acquisition->xsrc - xp, 2) ) ;
					if (dist <= radius)
					{
						ix_src[npoint_src] = ix ;
						iz_src[npoint_src] = iz ;
						weight_src[npoint_src] = exp(-dist*dist/(2.0*src_sigma*src_sigma)) ;
						sum_src += weight_src[npoint_src] ;
						npoint_src++ ;
					}
				}
			}
		}
	}

	// normalisation
	print_info(MASTER, " # npoints for src ", npoint_src) ;
	for (Myint ii=0; ii<npoint_src; ii++)
	{
		weight_src[ii] /= sum_src ;
	}

	// retrieve position of receivers
	nrec = (*acquisition).nrec ;
	ix_rec = allocate_array<Myint>(nrec) ;
	iz_rec = allocate_array<Myint>(nrec) ;

	for (Myint irec=0; irec<nrec; irec ++)
	{
		ix_rec[irec] = (Myint)(ixBeg2 + (*acquisition).xrec[irec] / dx) ;
		if (ix_rec[irec] < ixBeg2)
		{
			ix_rec[irec] = NOT_FOUND ;
		}
		else if (ix_rec[irec] > ixEnd2-1)
		{
			ix_rec[irec] = NOT_FOUND ;
		}

		iz_rec[irec] = (Myint)(izBeg2 + (*acquisition).zrec[irec] / dz) ;
		if (iz_rec[irec] < izBeg2)
		{
			iz_rec[irec] = NOT_FOUND ;
		}
		else if (iz_rec[irec] > izEnd2-1)
		{
			iz_rec[irec] = NOT_FOUND ;
		}

		if ((iz_rec[irec] == NOT_FOUND) || (ix_rec[irec] == NOT_FOUND))
		{
			print_warning(" receiver not in domain with idx=", irec) ;
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D::locate_src_and_rec_in_grid");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FDM_2D::allocate_grid_for_freq_wavefield(Grid** pGrid)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D::allocate_grid_for_freq_wavefield");

	*pGrid = new Grid_3D_complex(nz, pFreq_group->nb_freq, nx, dz, 1., dx) ;
	Rtn_code rtn_code = (*pGrid)->reset() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D::allocate_grid_for_freq_wavefield");
	return(RTN_CODE_OK) ;
} 

//-------------------------------------------------------------------------------------------------------
Rtn_code FDM_2D::deallocate_grid_for_freq_wavefield(Grid* pGrid)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_2D::deallocate_grid_for_freq_wavefield");

	Grid_3D_complex *pGrid2 = dynamic_cast<Grid_3D_complex*>(pGrid) ;
	if (pGrid2 == NULL)
	{
		print_error(" Error in FDM_2D::deallocate_grid_for_freq_wavefield") ;
		return(RTN_CODE_KO) ;
	}


	delete(pGrid2) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_2D::deallocate_grid_for_freq_wavefield");
	return(RTN_CODE_OK) ;
} 

} // namespace django
