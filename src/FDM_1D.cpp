//-------------------------------------------------------------------------------------------------------
//
// MODELLING WITH FDM IN 1D
//
// PARENT CLASS: Scheme
//   DERIVED CLASS: FDM
//     DERIVED CLASS: FDM_1D
//
//-------------------------------------------------------------------------------------------------------

#include "FDM_1D.h"

#include <cassert>
#include <fstream>
#include <iostream>

#include "acquisition.h"
#include "allocate_array.h"
#include "data_std.h"
#include "grid_1D_float.h"
#include "grid_2D_complex.h"
#include "output_report.h"
#include "type_def.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------

FDM_1D::FDM_1D(void) : FDM()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D::FDM_1D");
	dim = ONE ;
	iz_rec   = NULL ;
	iz_src   = NULL ;
	weight_src = NULL ;
	iz_pixel = NULL ;
	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D::FDM_1D");
}

//-------------------------------------------------------------------------------------------------------

// initialize scheme geometry based on input Grid

Rtn_code FDM_1D::initialize(Grid* pGrid_in, Myfloat fmax)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D::initialize(Grid* pGrid_in)");

	// check input grid is Grid_1D_float
	Grid_1D_float *pGrid = dynamic_cast<Grid_1D_float*>(pGrid_in) ;
	if (pGrid == NULL)
	{
		print_error("IN FDM_1D::initialize --> input grid is not Grid_1D_float");
		return(RTN_CODE_KO) ;
	}

	// set grid size and sampling
	nz = pGrid->nz ;
	dz = pGrid->dz ;
	inv_dz = 1. / dz ;
	inv_dz2 = inv_dz * inv_dz ;

	// compute optimal time step
	Myfloat optimal_dt = compute_optimal_time_step(pGrid) ;

	// set appropriate nt and dt
	if (set_nt_and_dt(optimal_dt) == RTN_CODE_KO) return(RTN_CODE_KO) ;

	// set boundary width
	nlayer_zBeg = get_boundary_width(ZBEG) ;
	nlayer_zEnd = get_boundary_width(ZEND) ;

	// usefull indices
	izBeg  = 0 ;
	izBeg1 = izBeg + lstencil ;
	izBeg2 = izBeg1 + nlayer_zBeg ;
	izEnd2 = izBeg2 + nz  ;
	izEnd1 = izEnd2 + nlayer_zEnd ;
	izEnd  = izEnd1 + lstencil -1 ;

	// update total and layer grid points numbers
	npoint_med = nz ;
	npoint_lay = nlayer_zBeg + nlayer_zEnd ;
	npoint     = npoint_med + npoint_lay ;

	// call parent initialization
	Rtn_code rtn_code = FDM::initialize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D::initialize(Grid* pGrid_in)");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_1D::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D::finalize");

	// call parent finalize
	Rtn_code rtn_code = FDM::finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D::finalize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat FDM_1D::compute_optimal_time_step(Grid_1D_float *vp)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D::compute_optimal_time_step");

	// evaluate optimal time step
	Myfloat max_vel = vp->get_max() ;
	Myfloat optim_dt = CFL * dz / max_vel ;
	print_info(MASTER, " Optimal dt (s):", optim_dt) ;

	if (ratio_cfl != 1.0)
	{
		optim_dt *= ratio_cfl ;
		print_info(MASTER, " Ratio CFL:\t", ratio_cfl) ;
		print_info(MASTER, " Modified dt (s):", optim_dt) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D::compute_optimal_time_step");
	return(optim_dt) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_1D::source_excitation(Myfloat* time_u, Myint it, Wavefield_type wtype, Data* pData)
{  
	print_debug(ALL, FULL_DEBUG, "IN FDM_1D::source_excitation(Myfloat* time_u, Myint it, Wavefield_type wtype, Data* pData)");

	if (wtype == INCIDENT)
	{

		// apply source function in time domain
		// excitation at the source position

		for (Myint ip=0; ip<npoint_src; ip++)
		{
			if (iz_src[ip] != NOT_FOUND)
			{
				if (it >= 0) time_u[iz_src[ip]] += weight_src[ip] * src_time_function[it] ;
			}
		}

	}
	else if (wtype == ADJOINT)
	{

		// build excitation in time domain via an inverse DFT of the frequency adjoint terms
		// excitation at the receivers positions

		// loop over frequencies
		Mycomplex **data_adj = ((Data_std*) pData)->pr_freq_adj_rec ;
		Myfloat coef2 = TEMP_RHO_CONST * pow(TEMP_VP_CONST, 2) * (2. / (nt * dz)) ;

		for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)
		{
			Myfloat coef_tmp = 2. * PI * pFreq_group->pFreq_list[ifreq] * (nt - it + 1) * dt ;
			Mycomplex coef (cos(coef_tmp), sin(coef_tmp)) ;

			for (Myint irec=0; irec<nrec; irec++)
			{
				time_u[iz_rec[irec]] += -1. * (data_adj[irec][ifreq] * coef).real() * coef2 ;
			}
		}
	}

	print_debug(ALL, FULL_DEBUG, "OUT FDM_1D::source_excitation(Myfloat* time_u, Myint it, Wavefield_type wtype, Data* pData)");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_1D::source_excitation(Myfloat* time_u, Myint it, Wavefield_type wtype, Data* pData, Myfloat src_factor)
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_1D::source_excitation(Myfloat* time_u, Myint it, Wavefield_type wtype, Data* pData, Myfloat src_factor)");

	if (wtype == INCIDENT)
	{

		// apply source function in time domain
		// excitation at the source postion

		print_debug(ALL, FULL_DEBUG, "src it", it);
		print_debug(ALL, FULL_DEBUG, "src amp",  src_time_function[it] * src_factor);

		for (Myint ip=0; ip<npoint_src; ip++)
		{
			if (iz_src[ip] != NOT_FOUND)
			{
				if (it >= 0) time_u[iz_src[ip]] += weight_src[ip] * src_time_function[it] * src_factor ;
			}
		}
	}
	else if (wtype == ADJOINT)
	{

		// build excitation in time domain via an inverse DFT of the frequency adjoint terms
		// excitation at the receivers positions

		// loop over frequencies
		Mycomplex **data_adj = ((Data_std*) pData)->pr_freq_adj_rec ;
		Myfloat coef2 = TEMP_RHO_CONST * pow(TEMP_VP_CONST, 2) * (2. / (nt * dz)) ;

		for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)
		{
			Myfloat coef_tmp = 2. * PI * pFreq_group->pFreq_list[ifreq] * (nt - it + 1) * dt ;
			Mycomplex coef (cos(coef_tmp), sin(coef_tmp)) ;

			for (Myint irec=0; irec<nrec; irec++)
			{
				time_u[iz_rec[irec]] += -1. * (data_adj[irec][ifreq] * coef).real() * coef2 ;
			}
		}
	}

	print_debug(ALL, FULL_DEBUG, "OUT FDM_1D::source_excitation");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_1D::update_freq_sol_grid(Myfloat* time_u, Mycomplex** freq_u, Myint it)
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_1D::update_freq_sol_grid");

	// freq. extraction not required
	if ((!output_frequency_flag) || (freq_u == NULL))
	{
		print_debug(ALL, FULL_DEBUG, "OUT FDM_1D::update_freq_sol_grid");
		return(RTN_CODE_OK) ;
	}

	Myint nz1 = 0 ;
	Myint nz2 = nz ;

	// coef for DFT
	for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)
	{
		Myfloat coef_tmp = -2. * PI * pFreq_group->pFreq_list[ifreq] * it * dt ;
		//Myfloat coef_tmp = -2. * PI * pFreq_group->pFreq_list[ifreq] * it * float(nt-1)/float(nt) * dt ;

		Mycomplex coef (cos(coef_tmp)*dt, sin(coef_tmp)*dt) ;

#pragma omp parallel for      
		for (Myint iz = nz1; iz < nz2; iz++)
		{
			freq_u[ifreq][iz] += time_u[iz+izBeg2] * coef ;
		}

	}

	print_debug(ALL, FULL_DEBUG, "OUT FDM_1D::update_freq_sol_grid");
	return(RTN_CODE_OK) ;

}

//-------------------------------------------------------------------------------------------------------

void FDM_1D::compute_energy(Myfloat* time_u)
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_1D::compute_energy");

	if (compute_energy_flag)
	{
		// compute energy
		Myint nz1 = izBeg2 ;
		Myint nz2 = izEnd2 ;

		Myfloat energy = 0. ;

		for (Myint iz = nz1; iz < nz2; iz++)
		{
			energy += time_u[iz] * time_u[iz] ;
		}

		// write energy
		write_energy() ;
	}

	print_debug(ALL, FULL_DEBUG, "OUT FDM_1D::compute_energy");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_1D::write_trace(Variable* pVar, Myint it)
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_1D::write_trace");

	if (dt_out != 0)
	{
		// check if ouput is required at the current time step
		Myint decim = round(dt_out / dt) ;
		string file_name = pVar->get_name() + TIME_REC_OUT_FILE ;

		if (it%decim != 0)
		{
			print_debug(ALL, FULL_DEBUG, "OUT FDM_1D::write_trace");
			return(RTN_CODE_OK) ;
		}
		else
		{
			print_debug(ALL, FULL_DEBUG, " write_trace in ", file_name) ;
			print_debug(ALL, FULL_DEBUG, " output trace (s) ", it*dt) ;
		}

		// write traces
		Rtn_code rtn_code = pVar->write_rec(nrec, iz_rec) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	}
	print_debug(ALL, FULL_DEBUG, "OUT FDM_1D::write_trace");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_1D::write_time_sol_snapshot(Myfloat* time_u, const char* file_name, Myint it)
{
	print_debug(ALL, FULL_DEBUG, "IN FDM_1D::write_time_sol_snapshot");

	// snapshot disabled
	if (!output_snapshot)
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

	print_debug(ALL, FULL_DEBUG, "OUT FDM_1D::write_time_sol_snapshot");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_1D::write_time_sol_grid(Myfloat* time_u, ofstream* pFile)
{
	print_debug(ALL, MID_DEBUG, "IN FDM_1D::write_time_sol_grid");

	Myint nz2 = izEnd+1 ;

	pFile->write((char*) &(time_u[izBeg2]), nz2 * sizeof(float)) ;

	print_debug(ALL, MID_DEBUG, "OUT FDM_1D::write_time_sol_grid");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_1D::write_freq_sol_rec(Mycomplex** freq_u)
{
	print_debug(ALL, MID_DEBUG, "IN FDM_1D::write_freq_sol_rec");

	// freq. extraction not required
	if (freq_u == NULL)
	{
		print_debug(ALL, MID_DEBUG, "OUT FDM_1D::write_freq_sol_rec");
		return(RTN_CODE_OK) ;
	}

	Mycomplex u_tmp[nrec*pFreq_group->nb_freq] ;

	// loop over receivers
	Myint ipos = 0 ;
	for (Myint irec=0; irec<nrec; irec++)
	{
		// loop over frequencies
		for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)
		{
			u_tmp[ipos] = freq_u[ifreq][iz_rec[irec]-izBeg2] ;
			ipos++ ;
		}
	}

	ofstream pFile ;
	pFile.open(PR_FREQ_REC_OUT_FILE, ios::binary | ios::app | ios::out) ;
	assert(pFile.is_open());
	pFile.write((char*)u_tmp, nrec * pFreq_group->nb_freq * sizeof(complex<float>)) ;
	pFile.close() ;

	print_debug(ALL, MID_DEBUG, "OUT FDM_1D::write_freq_sol_rec");
	return(RTN_CODE_OK) ;

}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_1D::store_freq_sol_rec(Mycomplex** freq_u, Data* pData)
{
	print_debug(ALL, MID_DEBUG, "IN FDM_1D::store_freq_sol_rec");

	// freq. extraction not required
	if (freq_u == NULL)
	{
		print_debug(ALL, MID_DEBUG, "OUT FDM_1D::store_freq_sol_rec");
		return(RTN_CODE_OK) ;
	}

	Mycomplex **pSol = ((Data_std*) pData)->pr_freq_sol_rec ;

	// loop over frequencies
	for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)
	{
		for (Myint irec=0; irec<nrec; irec++)
		{
			pSol[irec][ifreq] = freq_u[ifreq][iz_rec[irec]-izBeg2] ;
		}
	}

	print_debug(ALL, MID_DEBUG, "OUT FDM_1D::store_freq_sol_rec");
	return(RTN_CODE_OK) ;

}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_1D::write_freq_sol_grid(Mycomplex** freq_u)
{
	print_debug(ALL, MID_DEBUG, "IN FDM_1D::write_freq_sol_grid");

	// freq. extraction not required
	if (freq_u == NULL)
	{
		print_debug(ALL, FULL_DEBUG, "OUT FDM_1D::write_freq_sol_grid");
		return(RTN_CODE_OK) ;
	}

	// loop over frequencies
	ofstream pFile ;
	pFile.open(PR_FREQ_GRID_OUT_FILE, ios::binary | ios::app | ios::out) ;
	assert(pFile.is_open());
	for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)
	{
		pFile.write((char*) &(freq_u[ifreq][0]), nz * sizeof(complex<float>)) ;
	}
	pFile.close() ;

	print_debug(ALL, MID_DEBUG, "OUT FDM_1D::write_freq_sol_grid");
	return(RTN_CODE_OK) ;

}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_1D::free_position_arrays(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D::free_position_arrays");

	if (iz_rec != NULL)     deallocate_array<Myint>(iz_rec, nrec) ;

	if (iz_src != NULL)     deallocate_array<Myint>(iz_src, npoint_src) ;

	if (weight_src != NULL) deallocate_array<Myfloat>(weight_src, npoint_src) ;

	if (iz_pixel != NULL)   deallocate_array<Myint>(iz_pixel, npixel) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D::free_position_arrays");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_1D::locate_pixel_in_grid(Snapshot *pSnapshot)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D::locate_pixel_in_grid");

	print_warning(" locate_pixel_in_grid not yet implemented") ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D::locate_pixel_in_grid");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FDM_1D::locate_src_and_rec_in_grid(Acquisition *acquisition)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D::locate_src_and_rec_in_grid");

	// retrieve points for the source excitation
	//#########################################

	npoint_src = 0 ;
	Myfloat sum_src = 0.0 ;

	// case one, NO source smoothing
	// excitation on single node
	//==============================

	if (src_type == SRC_POINT)
	{
		iz_src = allocate_array<Myint>(1) ;
		weight_src = allocate_array<Myfloat>(1) ;
		npoint_src    = 1 ;
		weight_src[0] = 1.0 ;
		sum_src       = 1.0 ;

		// retrieve position of source
		// --> nearest grid point
		iz_src[0] = (Myint)(izBeg2 + acquisition->zsrc / dz) ;
		if (iz_src[0] < izBeg2)
		{
			iz_src[0] = NOT_FOUND ;
		}
		else if (iz_src[0] > izEnd2-1)
		{
			iz_src[0] = NOT_FOUND ;
		}

		if (iz_src[0] == NOT_FOUND)
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

		for (Myint iz = izMin; iz < izMax; iz++)
		{
			Myfloat zp = (iz-izMin) * dz ;
			Myfloat dist = abs(acquisition->zsrc - zp) ;
			if (dist <= radius)
			{
				npoint_src++ ;
			}
		}

		if (npoint_src == 0)
		{
			print_warning(" Source is not in domain") ;
		}
		else
		{
			iz_src = allocate_array<Myint>(npoint_src) ;
			weight_src = allocate_array<Myfloat>(npoint_src) ;
			npoint_src = 0 ;

			// 2nd loop to store the points and weights

			for (Myint iz = izMin; iz < izMax; iz++)
			{
				Myfloat zp = (iz-izMin) * dz ;
				Myfloat dist =  abs(acquisition->zsrc - zp) ;
				if (dist <= radius)
				{
					iz_src[npoint_src] = iz ;
					weight_src[npoint_src] = exp(-dist*dist/(2.0*src_sigma*src_sigma)) ;
					sum_src += weight_src[npoint_src] ;
					npoint_src++ ;
				}
			}
		}
	}

	print_info(MASTER, " # points for src ", npoint_src) ;
	for (Myint ii=0; ii<npoint_src; ii++)
	{
		weight_src[ii] /= sum_src ;
	}

	// retrieve positions of receivers
	nrec = acquisition->nrec ;
	iz_rec = allocate_array<Myint>(nrec) ;

	for (Myint irec=0; irec<nrec; irec ++)
	{
		iz_rec[irec] = (Myint)(izBeg2 + acquisition->zrec[irec] / dz) ;
		if (iz_rec[irec] < izBeg2)
		{
			iz_rec[irec] = NOT_FOUND ;
		}
		else if (iz_rec[irec] > izEnd2-1)
		{
			iz_rec[irec] = NOT_FOUND ;
		}

		if (iz_rec[irec] == NOT_FOUND)
		{
			print_warning(" receiver not in domain with idx=", irec) ;
		}

	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D::locate_src_and_rec_in_grid");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FDM_1D::allocate_grid_for_freq_wavefield(Grid** pGrid)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D::allocate_grid_for_freq_wavefield");

	*pGrid = new Grid_2D_complex(nz, pFreq_group->nb_freq, dz, 1.) ;
	Rtn_code rtn_code = (*pGrid)->reset() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D::allocate_grid_for_freq_wavefield");
	return(RTN_CODE_OK) ;

}

//-------------------------------------------------------------------------------------------------------
Rtn_code FDM_1D::deallocate_grid_for_freq_wavefield(Grid* pGrid)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FDM_1D::deallocate_grid_for_freq_wavefield");

	Grid_2D_complex *pGrid2 = dynamic_cast<Grid_2D_complex*>(pGrid) ;
	if (pGrid2 == NULL)
	{
		print_error(" Error in FDM_1D::deallocate_grid_for_freq_wavefield") ;
		return(RTN_CODE_KO) ;
	}

	delete(pGrid2) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FDM_1D::deallocate_grid_for_freq_wavefield");
	return(RTN_CODE_OK) ;
}

} // namespace django
