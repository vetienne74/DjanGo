//-------------------------------------------------------------------------------------------------------
//
// GRADIENT OF THE OBJECTIVE FUNCTION
//
//-------------------------------------------------------------------------------------------------------

#include "gradient.h"

#include <iostream>

#include "mpi.h"

#include "acquisition.h"
#include "allocate_array.h"
#include "data.h"
#include "data_std.h"
#include "freq_group.h"
#include "grid.h"
#include "grid_2D_complex.h"
#include "grid_3D_complex.h"
#include "model.h"
#include "modelling.h"
#include "output_report.h"
#include "snapshot.h"

using namespace std;

namespace django {

static const Myfloat TEMP_VP_PERTURB = 10. ;

//-------------------------------------------------------------------------------------------------------
//
// CONSTRUCTORS
//
//-------------------------------------------------------------------------------------------------------

Gradient::Gradient(Acquisition* pAcquisition, Model* pModel, Freq_group* pFreq_group_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Gradient::Gradient");

	pFreq_group = pFreq_group_in ;
	src_coef = allocate_array<Mycomplex>(pAcquisition->nsrc, pFreq_group->nb_freq) ;

	this->pAcquisition = pAcquisition ;

	//-----------------------
	// Input Model_2D_ac_iso
	//-----------------------
	// Model_2D_ac_iso *pModel2D = dynamic_cast<Model_2D_ac_iso*>(pModel) ;
	// if (pModel2D != NULL)
	//   {
	//     Myint nz2   = pModel2D->vp->nz ;
	//     Myint nx2   = pModel2D->vp->nx ;
	//     Myfloat dz2 = pModel2D->vp->dz ;
	//     Myfloat dx2 = pModel2D->vp->dx ;

	//     this->pGrad_model = new Model_2D_ac_iso ;
	//     (((Model_2D_ac_iso*) pGrad_model))->vp = new Grid_2D_float(nz2, nx2, dz2, dx2) ;
	//     (((Model_2D_ac_iso*) pGrad_model))->vp_file_name = GRADIENT_VP_OUT_FILE ;

	//     this->pPrecond = new Model_2D_ac_iso ;
	//     (((Model_2D_ac_iso*) pPrecond))->vp = new Grid_2D_float(nz2, nx2, dz2, dx2) ;
	//     (((Model_2D_ac_iso*) pPrecond))->vp_file_name = PRECOND_VP_OUT_FILE ;

	//     this->pGrad_precond = new Model_2D_ac_iso ;
	//     (((Model_2D_ac_iso*) pGrad_precond))->vp = new Grid_2D_float(nz2, nx2, dz2, dx2) ;
	//     (((Model_2D_ac_iso*) pGrad_precond))->vp_file_name = GRAD_PRECOND_VP_OUT_FILE ;
	//   }

	//-----------------------
	// Input Model_1D_ac_iso
	//-----------------------
	// Model_1D_ac_iso *pModel1D = dynamic_cast<Model_1D_ac_iso*>(pModel) ;
	// if (pModel1D != NULL)
	//   {
	//     Myint nz2   = pModel1D->vp->nz ;
	//     Myfloat dz2 = pModel1D->vp->dz ;

	//     this->pGrad_model = new Model_1D_ac_iso ;
	//     (((Model_1D_ac_iso*) pGrad_model))->vp = new Grid_1D_float(nz2, dz2) ;
	//     (((Model_1D_ac_iso*) pGrad_model))->vp_file_name = GRADIENT_VP_OUT_FILE ;

	//     this->pPrecond = new Model_1D_ac_iso ;
	//     (((Model_1D_ac_iso*) pPrecond))->vp = new Grid_1D_float(nz2, dz2) ;
	//     (((Model_1D_ac_iso*) pPrecond))->vp_file_name = PRECOND_VP_OUT_FILE ;

	//     this->pGrad_precond = new Model_1D_ac_iso ;
	//     (((Model_1D_ac_iso*) pGrad_precond))->vp = new Grid_1D_float(nz2, dz2) ;
	//     (((Model_1D_ac_iso*) pGrad_precond))->vp_file_name = GRAD_PRECOND_VP_OUT_FILE ;
	//   }

	//if ((pModel2D == NULL) && (pModel1D == NULL))
	if (true)
	{
		print_error(" Cannot instanciate gradient") ;
	}
	else
	{
		//-------------------------------------------------------------------------------------------------------
		// read gradient configuration file
		//-------------------------------------------------------------------------------------------------------
		Rtn_code rtn_code = (*this).read_config() ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Gradient::Gradient");
}

//-------------------------------------------------------------------------------------------------------
//
// DESTRUCTOR
//
//-------------------------------------------------------------------------------------------------------

Gradient::~Gradient(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Gradient::~Gradient");

	deallocate_array<Mycomplex>(src_coef, pAcquisition->nsrc, pFreq_group->nb_freq) ;

	// Model_2D_ac_iso* pGrad_2D = dynamic_cast<Model_2D_ac_iso*>(pGrad_model) ;
	// Model_1D_ac_iso* pGrad_1D = dynamic_cast<Model_1D_ac_iso*>(pGrad_model) ;

	// if (pGrad_2D != NULL)
	//   {
	//     delete( (Model_2D_ac_iso*) pGrad_model) ;
	//     delete( (Model_2D_ac_iso*) pPrecond) ;
	//     delete( (Model_2D_ac_iso*) pGrad_precond) ;
	//   }

	// if (pGrad_1D != NULL)
	//   {
	//     delete( (Model_1D_ac_iso*) pGrad_model) ;
	//     delete( (Model_1D_ac_iso*) pPrecond) ;
	//     delete( (Model_1D_ac_iso*) pGrad_precond) ;
	//   }

	print_debug(ALL, LIGHT_DEBUG, "OUT Gradient::~Gradient");
}

//-------------------------------------------------------------------------------------------------------
void Gradient::open_config_file(const char* file_name) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Gradient::open_config_file");
	config_file.open(file_name) ;
	assert(config_file.is_open());
	print_debug(ALL, LIGHT_DEBUG, "OUT Gradient::open_config_file");

}

//-------------------------------------------------------------------------------------------------------
void Gradient::close_config_file(void) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Gradient::close_config_file");
	config_file.close() ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Gradient::close_config_file");

}

//-------------------------------------------------------------------------------------------------------

Rtn_code Gradient::read_config()
{
	print_debug(ALL, LIGHT_DEBUG, "IN Gradient::read_config");

	//-------------------------------------------------------------------------------------------------------
	// open configuration file
	//-------------------------------------------------------------------------------------------------------

	open_config_file(GRAD_CONFIG_IN_FILE) ;

	print_info(MASTER, "") ;
	print_line2() ;
	print_info(MASTER, " GRADIENT PARAMETERS") ;
	print_info(MASTER, "");

	Myint int_tmp ;

	// get data type
	config_file >> int_tmp ;

	switch(int_tmp)
	{
	case TIME:
		domain = TIME ;
		print_info(MASTER, " Gradient domain: \tTIME DOMAIN") ;
		break ;
	case FREQ:
		domain = FREQ ;
		print_info(MASTER, " Gradient domain: \tFREQUENCY DOMAIN") ;
		break ;
	default:
		print_error(" Invalid gradient domain") ;
		return(RTN_CODE_KO) ;
	}

	// gradient computation type
	config_file >> int_tmp ;

	if (int_tmp == ADJOINT_METHOD)
	{
		print_info(MASTER, " Grad. computation:", "ADJOINT METHOD") ;
		type = ADJOINT_METHOD ;
	}
	else if (int_tmp == FD_FCOST)
	{
		print_info(MASTER, " Grad. computation:", "FD OF OBJECTIVE FUNCTION") ;
		type = FD_FCOST ;
	}
	else
	{
		print_error(" Invalid parameter for grad. computation type") ;
		return(RTN_CODE_KO) ;
	}

	// fcost function
	config_file >> int_tmp ;

	if (int_tmp == DIFF_L2)
	{
		print_info(MASTER, " Objective function:", "DIFF. L2") ;
		fcost_type = DIFF_L2 ;
	}
	else if (int_tmp == DIFF_L1)
	{
		print_info(MASTER, " Objective function:", "DIFF. L1") ;
		fcost_type = DIFF_L1 ;
	}
	else if (int_tmp == LOG_L2)
	{
		print_info(MASTER, " Objective function:", "LOG. L2") ;
		fcost_type = LOG_L2 ;
	}
	else if (int_tmp == LOG_L2_PHASE)
	{
		print_info(MASTER, " Objective function:", "LOG. L2 PHASE ONLY") ;
		fcost_type = LOG_L2_PHASE ;
	}
	else if (int_tmp == LOG_L2_AMP)
	{
		print_info(MASTER, " Objective function:", "LOG. L2 AMPLITUDE ONLY") ;
		fcost_type = LOG_L2_AMP ;
	}
	else
	{
		print_error(" Invalid parameter for objective function type") ;
		return(RTN_CODE_KO) ;
	}

	// source estimation
	config_file >> int_tmp ;

	if (int_tmp == 0)
	{
		print_info(MASTER, " Source estimation:", "NO") ;
		src_estim_flag = false ;
	}
	else if (int_tmp == 1)
	{
		print_info(MASTER, " Source estimation:", "YES") ;
		src_estim_flag = true ;
	}
	else
	{
		print_error(" Invalid parameter for source estimation") ;
		return(RTN_CODE_KO) ;
	}

	// min and max offset for source estimation
	config_file >> src_estim_offset_min ;
	if (src_estim_offset_min >= 0.)
	{
		print_info(MASTER, " Min offset/src estim:", src_estim_offset_min) ;
	}
	else
	{
		print_error(" Min offset/src estim should be >= 0.") ;
		return(RTN_CODE_KO) ;
	}
	config_file >> src_estim_offset_max ;
	if (src_estim_offset_max >= 0.)
	{
		print_info(MASTER, " Max offset/src estim:", src_estim_offset_max) ;
		if (src_estim_offset_min > src_estim_offset_max)
		{
			print_error(" Max offset/src estim should be >= Min offset/src estim") ;
			return(RTN_CODE_KO) ;
		}
	}
	else
	{
		print_error("  Max offset/src estim should be >= 0.") ;
		return(RTN_CODE_KO) ;
	}

	// min and max offset for gradient computation
	config_file >> grad_offset_min ;
	if (grad_offset_min >= 0.)
	{
		print_info(MASTER, " Min offset/gradient:", grad_offset_min) ;
	}
	else
	{
		print_error(" Min offset/grad. comp. should be >= 0.") ;
		return(RTN_CODE_KO) ;
	}
	config_file >> grad_offset_max ;
	if (grad_offset_max >= 0.)
	{
		print_info(MASTER, " Max offset/gradient:", grad_offset_max) ;
		if (grad_offset_min > grad_offset_max)
		{
			print_error(" Max offset/grad. comp. should be >= Min offset/grad. comp.") ;
			return(RTN_CODE_KO) ;
		}
	}
	else
	{
		print_error(" Max offset/grad. comp. should be >= 0.") ;
		return(RTN_CODE_KO) ;
	}

	// coef for gradient preconditionning
	config_file >> coef_grad_precond ;

	if (coef_grad_precond >= 0.)
	{
		print_info(MASTER, " Coef. grad. precond.:", coef_grad_precond) ;
	}
	else
	{
		print_error(" Coef. grad. precond. should be >= 0.") ;
		return(RTN_CODE_KO) ;
	}

	print_line2() ;

	//-------------------------------------------------------------------------------------------------------
	// close configuration file
	//-------------------------------------------------------------------------------------------------------
	close_config_file() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Gradient::read_config");
	return(RTN_CODE_OK) ;
}


//-------------------------------------------------------------------------------------------------------

Rtn_code Gradient::initialize() 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Gradient::initialize");

	if (pGrad_model->reset()   != RTN_CODE_OK) return(RTN_CODE_KO) ;
	if (pPrecond->reset()      != RTN_CODE_OK) return(RTN_CODE_KO) ;
	if (pGrad_precond->reset() != RTN_CODE_OK) return(RTN_CODE_KO) ;

	for (Myint isrc = 0; isrc < pAcquisition->nsrc; isrc++)
	{
		for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)
		{
			src_coef[isrc][ifreq]  = 1 ;
		}
	}
	print_debug(ALL, LIGHT_DEBUG, "OUT Gradient::initialize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Gradient::store_value(Myint iz, Myint ix, Myint iy, Myfloat grad_val)
{

	print_debug(ALL, LIGHT_DEBUG, "IN Gradient::store_value");

	//-----------------------
	// Input Model_2D_ac_iso
	//-----------------------
	// Model_2D_ac_iso *pModel2D = dynamic_cast<Model_2D_ac_iso*>(pGrad_model) ;
	// if (pModel2D != NULL)
	//   {
	//     ((Model_2D_ac_iso*) pGrad_model)->vp->pArray[ix][iz] = grad_val ;
	//   }

	//-----------------------
	// Input Model_1D_ac_iso
	//-----------------------
	// Model_1D_ac_iso *pModel1D = dynamic_cast<Model_1D_ac_iso*>(pGrad_model) ;
	// if (pModel1D != NULL)
	//   {
	//     ((Model_1D_ac_iso*) pGrad_model)->vp->pArray[iz] = grad_val ;
	//   }

	//if ((pModel2D == NULL) && (pModel1D == NULL))
	if (true)
	{
		print_error(" Cannot store value in gradient") ;
		return(RTN_CODE_KO) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Gradient::store_value");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Gradient::update(Grid* pIncident, Grid* pAdjoint, Model* pModel, Modelling* pModelling)
{

	print_debug(ALL, LIGHT_DEBUG, "IN Gradient::update");

	Freq_group* pFreq_group = pModelling->get_freq_group() ;

	//-----------------------
	// Input Model_2D_ac_iso
	//-----------------------
	// Model_2D_ac_iso *pModel2D = dynamic_cast<Model_2D_ac_iso*>(pModel) ;
	// if (pModel2D != NULL)

	//   {

	//     // gradient formula:
	//     // coef_i = 2 / (rho_i * vp_i**3) * -i * freq
	//     // grad = real(coef_i * incident_i * conj(adjoint_i) )

	//     Grid_3D_complex* pIncident2 = (Grid_3D_complex*) pIncident ;
	//     Grid_3D_complex* pAdjoint2  = (Grid_3D_complex*) pAdjoint ;
	//     Model_2D_ac_iso* pModel2    = (Model_2D_ac_iso*) pModel ;

	//     Myint nz2 = pIncident2->nz ;
	//     Myint nx2 = pIncident2->ny ;
	//     Myfloat **vp2 = pModel2->vp->pArray ;

	//     // loop over frequencies
	//     for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)
	// 	{

	// 	  Mycomplex*** pIncident3 = pIncident2->pArray ;
	// 	  Mycomplex*** pAdjoint3  = pAdjoint2->pArray ;

	// 	  // compute coef = 2iw (1st order wave eq.)
	// 	  Mycomplex coef = (Mycomplex)I_CMPLX * (Myfloat) (2. * 2. * PI * pFreq_group->pFreq_list[ifreq]) ;

	// 	  // compute coef = 2 w^2 (2nd order wave eq.)
	// 	  //Mycomplex coef = 2. * pow(2. * PI * pFreq_group->pFreq_list[ifreq], 2.) ;

	// 	  // loop over the grid points
	// 	  for (Myint ix = 0; ix < nx2; ix++)
	// 	    {
	// 	      for (Myint iz = 0; iz < nz2; iz++)
	// 		{
	// 		  // dA/dm
	// 		  Mycomplex coef1 = coef / (TEMP_RHO_CONST * pow(vp2[ix][iz],3)) ;

	// 		  // preconditionner (pseudo-hessian)
	// 		  ((Model_2D_ac_iso*) pPrecond)->vp->pArray[ix][iz] += real(coef1 * pIncident3[ifreq][ix][iz]
	// 									    * conj(coef1 * pIncident3[ifreq][ix][iz])) ;

	// 		  // gradient
	// 		  //((Model_2D_ac_iso*) pGrad_model)->vp->pArray[ix][iz]
	// 		  //+= real(coef1* pIncident3[ifreq][ix][iz] * conj(pAdjoint3[ifreq][ix][iz])) ;
	// 		  ((Model_2D_ac_iso*) pGrad_model)->vp->pArray[ix][iz] += real(coef1 * pIncident3[ifreq][ix][iz]
	// 									       * pAdjoint3[ifreq][ix][iz] ) ;
	// 		}
	// 	    }
	// 	}
	//   }

	//-----------------------
	// Input Model_1D_ac_iso
	//-----------------------
	// Model_1D_ac_iso *pModel1D = dynamic_cast<Model_1D_ac_iso*>(pModel) ;
	// if (pModel1D != NULL)

	//   {

	//     // gradient formula:
	//     // coef_i = 2 / (rho_i * vp_i**3) * -i * freq
	//     // grad = real(coef_i * incident_i * conj(adjoint_i) )

	//     Grid_2D_complex* pIncident2 = (Grid_2D_complex*) pIncident ;
	//     Grid_2D_complex* pAdjoint2  = (Grid_2D_complex*) pAdjoint ;
	//     Model_1D_ac_iso* pModel2    = (Model_1D_ac_iso*) pModel ;

	//     Myint nz2 = pIncident2->nz ;
	//     Myfloat *vp2 = pModel2->vp->pArray ;

	//     // loop over frequencies
	//     for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)
	// 	{

	// 	  Mycomplex** pIncident3 = pIncident2->pArray ;
	// 	  Mycomplex** pAdjoint3  = pAdjoint2->pArray ;

	// 	  // compute coef = 2iw (1st order wave eq.)
	// 	  Mycomplex coef = (Mycomplex)I_CMPLX * (Myfloat)(2. * 2. * PI * pFreq_group->pFreq_list[ifreq]) ;

	// 	  // compute coef = 2 w^2 (2nd order wave eq.)
	// 	  //Mycomplex coef = 2. * pow(2. * PI * pFreq_group->pFreq_list[ifreq], 2.) ;

	// 	  // loop over the grid points
	// 	  for (Myint iz = 0; iz < nz2; iz++)
	// 	    {
	// 	      // dA/dm
	// 	      Mycomplex coef1 = coef / (TEMP_RHO_CONST * pow(vp2[iz],3)) ;

	// 	      // preconditionner (pseudo-hessian)
	// 	      ((Model_1D_ac_iso*) pPrecond)->vp->pArray[iz] += real(coef1 * pIncident3[ifreq][iz] * conj(coef1 * pIncident3[ifreq][iz])) ;

	// 	      // gradient
	// 	      //((Model_1D_ac_iso*) pGrad_model)->vp->pArray[iz] += real(coef1 * pIncident3[ifreq][iz] * conj(pAdjoint3[ifreq][iz])) ;
	// 	      ((Model_1D_ac_iso*) pGrad_model)->vp->pArray[iz] += real(coef1 * pIncident3[ifreq][iz] * pAdjoint3[ifreq][iz] ) ;
	// 	    }
	// 	}
	//   }

	//if ((pModel2D == NULL) && (pModel1D == NULL))
	if (true)
	{
		print_error(" Cannot update gradient") ;
		return(RTN_CODE_KO) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Gradient::update");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Gradient::gather(void)
{

	print_debug(ALL, LIGHT_DEBUG, "IN Gradient::gather");

	//-----------------------
	// Input Model_2D_ac_iso
	//-----------------------
	// Model_2D_ac_iso *pModel2D = dynamic_cast<Model_2D_ac_iso*>(pGrad_model) ;
	// if (pModel2D != NULL)

	//   {
	//     Myint nz2 = ((Model_2D_ac_iso*) pGrad_model)->vp->nz ;
	//     Myint nx2 = ((Model_2D_ac_iso*) pGrad_model)->vp->nx ;

	//     Myfloat prewhite = coef_grad_precond * max( abs(((Model_2D_ac_iso*) pPrecond)->vp->get_max()), abs(((Model_2D_ac_iso*) pPrecond)->vp->get_min()) ) ;

	//     // loop over the grid points
	//     for (Myint ix = 0; ix < nx2; ix++)
	// 	{
	// 	  for (Myint iz = 0; iz < nz2; iz++)
	// 	    {
	// 	      // divide gradient by preconditionner
	// 	      ((Model_2D_ac_iso*) pGrad_precond)->vp->pArray[ix][iz] = ((Model_2D_ac_iso*) pGrad_model)->vp->pArray[ix][iz]
	// 		/ ( ((Model_2D_ac_iso*) pPrecond)->vp->pArray[ix][iz] + prewhite) ;

	// 	      //((Model_2D_ac_iso*) pGrad_precond)->vp->pArray[ix][iz] = ((Model_2D_ac_iso*) pGrad_model)->vp->pArray[ix][iz] ;
	// 	    }
	// 	}
	//   }

	//-----------------------
	// Input Model_1D_ac_iso
	//-----------------------
	// Model_1D_ac_iso *pModel1D = dynamic_cast<Model_1D_ac_iso*>(pGrad_model) ;
	// if (pModel1D != NULL)

	//   {
	//     Myint nz2 = ((Model_1D_ac_iso*) pGrad_model)->vp->nz ;

	//     Myfloat prewhite = coef_grad_precond * max( abs(((Model_1D_ac_iso*) pPrecond)->vp->get_max()),
	// 						  abs(((Model_1D_ac_iso*) pPrecond)->vp->get_min()) ) ;

	//     // loop over the grid points
	//     for (Myint iz = 0; iz < nz2; iz++)
	// 	{
	// 	  // divide gradient by preconditionner
	// 	  ((Model_1D_ac_iso*) pGrad_precond)->vp->pArray[iz] = ((Model_1D_ac_iso*) pGrad_model)->vp->pArray[iz]
	// 	    / ( ((Model_1D_ac_iso*) pPrecond)->vp->pArray[iz] + prewhite) ;

	// 	  //((Model_2D_ac_iso*) pGrad_precond)->vp->pArray[ix][iz] = ((Model_2D_ac_iso*) pGrad_model)->vp->pArray[ix][iz] ;
	// 	}
	//   }

	//if ((pModel2D == NULL) && (pModel1D == NULL))
	if (true)
	{
		print_error(" Cannot gather gradient") ;
		return(RTN_CODE_KO) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT Gradient::gather");
	return(RTN_CODE_OK) ;
}


//-------------------------------------------------------------------------------------------------------

Rtn_code Gradient::evaluate_fcost_wavelet_and_adjoint_source(Acquisition* pAcquisition, Data* pData2, Modelling* pModelling, Myfloat* fcost) 

{

	print_debug(ALL, LIGHT_DEBUG, "IN Gradient::evaluate_fcost_wavelet_and_adjoint_source");

	Data_std *pData = dynamic_cast<Data_std*>(pData2) ;
	if (pData == NULL)
	{
		print_error(" Only seismic data are supported ") ;
		return(RTN_CODE_KO) ;
	}

	//-----------
	// read data
	//-----------
	Rtn_code rtn_code = pData->read(pAcquisition) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	Myint isrc_current = pAcquisition->current_src - 1;

	//-------------------
	// source estimation
	//-------------------
	print_debug(ALL, MID_DEBUG, "perform source estimation");

	Freq_group* pFreq_group = pModelling->get_freq_group() ;
	Mycomplex src_coef1, src_coef2 ;

	for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)
	{

		if (src_estim_flag)
		{
			src_coef1 = ZERO_CMPLX  ;
			src_coef2 = ZERO_CMPLX ;
			src_coef[isrc_current][ifreq] = ZERO_CMPLX ;

			for (Myint irec=0; irec<pAcquisition->nrec; irec++)
			{
				Myfloat dist = sqrt( pow(pAcquisition->xsrc-pAcquisition->xrec[irec], 2.) +  pow(pAcquisition->zsrc-pAcquisition->zrec[irec], 2.)) ;

				if ((dist < src_estim_offset_min) || (dist > src_estim_offset_max)) continue ;

				// estimate correction to be applied on computed data
				//src_coef1 +=  pData->pr_freq_obs_rec[irec][ifreq] * conj( pData->pr_freq_sol_rec[irec][ifreq] ) ;
				//src_coef2 +=  pData->pr_freq_sol_rec[irec][ifreq] * conj( pData->pr_freq_sol_rec[irec][ifreq] ) ;

				// estimate correction to be applied on observed data
				src_coef1 +=  pData->pr_freq_sol_rec[irec][ifreq] * conj( pData->pr_freq_obs_rec[irec][ifreq] ) ;
				src_coef2 +=  pData->pr_freq_obs_rec[irec][ifreq] * conj( pData->pr_freq_obs_rec[irec][ifreq] ) ;

			} // for (Myint irec=0; irec<pAcquisition->nrec; irec++)

			src_coef[isrc_current][ifreq] = src_coef1 / src_coef2 ;
		}
		else
		{
			src_coef[isrc_current][ifreq] = 1. ;
		}

		//cout << "src_coef1 " << src_coef1 << "\n" ;
		//cout << "src_coef2 " << src_coef2 << "\n" ;

		if ((isnan(real(src_coef[isrc_current][ifreq]))) || (isnan(imag(src_coef[isrc_current][ifreq]))))
		{
			cout << "src_coef1 " << src_coef1 << "\n" ;
			cout << "src_coef2 " << src_coef2 << "\n" ;
			cout << " src_coef[isrc_current][ifreq] "  << src_coef[isrc_current][ifreq] << "\n" ;
			print_error(" NAN detected in source estimation") ;
			return(RTN_CODE_KO) ;
		}

	} // for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)

	//---------------------------------
	// output estimated source wavelet
	//---------------------------------
	print_debug(ALL, MID_DEBUG, "output estimated source wavelet");
	//pModelling->output_estimated_src_wavelet(src_coef) ;

	//--------------------------
	// fcost and adjoint source
	//--------------------------
	print_debug(ALL, MID_DEBUG, "compute fcost and adjoint soure");
	*fcost = 0. ;
	for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)
	{
		print_debug(ALL, FULL_DEBUG, "ifreq", ifreq);
		for (Myint irec=0; irec<pAcquisition->nrec; irec++)
		{

			// correction of observed data
			pData->pr_freq_obs_rec[irec][ifreq] = pData->pr_freq_obs_rec[irec][ifreq] * src_coef[isrc_current][ifreq] ;

			//cout << "pData->pr_freq_obs_rec[irec][ifreq] " << pData->pr_freq_obs_rec[irec][ifreq] << "\n" ;
			//cout << "pData->pr_freq_sol_rec[irec][ifreq] " << pData->pr_freq_sol_rec[irec][ifreq] << "\n" ;

			// 1d
			//Myfloat dist = abs(pAcquisition->zsrc-pAcquisition->zrec[irec]) ;

			// 2d
			Myfloat dist = sqrt( pow(pAcquisition->xsrc-pAcquisition->xrec[irec], 2.) +
					pow(pAcquisition->zsrc-pAcquisition->zrec[irec], 2.)) ;

			if ((dist < grad_offset_min) || (dist > grad_offset_max)) {
				Mycomplex res = ZERO_CMPLX ;
				pData->pr_freq_adj_rec[irec][ifreq] = res ;
				*fcost += 0.5 * real( res * conj(res) ) ;
				continue ;
			}

			if (fcost_type == DIFF_L2)
			{
				//Mycomplex res = (pData->pr_freq_obs_rec[irec][ifreq] - pData->pr_freq_sol_rec[irec][ifreq] * src_coef[ifreq]) ;
				//pData->pr_freq_adj_rec[irec][ifreq] = res * conj(src_coef[ifreq]) ;

				Mycomplex res = (pData->pr_freq_obs_rec[irec][ifreq] - pData->pr_freq_sol_rec[irec][ifreq]) ;
				pData->pr_freq_adj_rec[irec][ifreq] = res ;

				*fcost += 0.5 * real( res * conj(res) ) ;
			}
			else if (fcost_type == DIFF_L1)
			{
				//Mycomplex res = (pData->pr_freq_obs_rec[irec][ifreq] - pData->pr_freq_sol_rec[irec][ifreq] * src_coef[ifreq]) ;
				//pData->pr_freq_adj_rec[irec][ifreq] = conj(src_coef[ifreq]) * ( res / abs(res) ) ;
				Mycomplex res = (pData->pr_freq_obs_rec[irec][ifreq] - pData->pr_freq_sol_rec[irec][ifreq]) ;
				pData->pr_freq_adj_rec[irec][ifreq] = ( res / abs(res) ) ;
				*fcost += abs(res) ;
			}
			else if (fcost_type == LOG_L2)
			{
				pData->pr_freq_adj_rec[irec][ifreq] = (log(pData->pr_freq_sol_rec[irec][ifreq] / pData->pr_freq_obs_rec[irec][ifreq]) / pData->pr_freq_sol_rec[irec][ifreq]) ;
				*fcost += 0.5 * real(log(pData->pr_freq_sol_rec[irec][ifreq] / pData->pr_freq_obs_rec[irec][ifreq]) * conj(log(pData->pr_freq_sol_rec[irec][ifreq] / pData->pr_freq_obs_rec[irec][ifreq]))) ;
			}
			else if (fcost_type == LOG_L2_AMP)
			{
				pData->pr_freq_adj_rec[irec][ifreq] = real(log(pData->pr_freq_sol_rec[irec][ifreq] / pData->pr_freq_obs_rec[irec][ifreq])) / pData->pr_freq_sol_rec[irec][ifreq] ;
				*fcost += 0.5 * pow(real(log(pData->pr_freq_sol_rec[irec][ifreq] / pData->pr_freq_obs_rec[irec][ifreq])), 2.) ;
			}
			else if (fcost_type == LOG_L2_PHASE)
			{
				pData->pr_freq_adj_rec[irec][ifreq] = - I_CMPLX * imag(log(pData->pr_freq_sol_rec[irec][ifreq] / pData->pr_freq_obs_rec[irec][ifreq])) / pData->pr_freq_sol_rec[irec][ifreq] ;
				*fcost += 0.5 * pow(imag(log(pData->pr_freq_sol_rec[irec][ifreq] / pData->pr_freq_obs_rec[irec][ifreq])), 2.) ;
			}

			// RTM
			//pr_freq_adj_rec[irec][ifreq] = pr_freq_obs_rec[irec][ifreq] ;

			//cout << " pData->pr_freq_obs_rec[irec][ifreq] "  << pData->pr_freq_obs_rec[irec][ifreq] << "\n" ;
			//cout << " pData->pr_freq_sol_rec[irec][ifreq] "  << pData->pr_freq_sol_rec[irec][ifreq] << "\n" ;
			//cout << " pData-> pr_freq_adj_rec[irec][ifreq] " << pData->pr_freq_adj_rec[irec][ifreq] << "\n" ;

		} // for (Myint irec=0; irec<pAcquisition->nrec; irec++)

	} // for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)

	// write adjoint source term
	print_debug(ALL, MID_DEBUG, "write adjoint source term");
	ofstream res_file ;
	res_file.open(PR_FREQ_ADJ_SRC_OUT_FILE, ios::binary| ios::app ) ;
	assert(res_file.is_open());
	res_file.write((char*) &(pData->pr_freq_adj_rec[0][0]), pAcquisition->nrec * pFreq_group->nb_freq * sizeof(Mycomplex)) ;
	res_file.close() ;

	// write observed data
	print_debug(ALL, MID_DEBUG, "write observed data");
	ofstream obs_file ;
	obs_file.open("OBSERVED_DATA.django.out", ios::binary | ios::app) ;
	assert(obs_file.is_open());
	obs_file.write((char*) &(pData->pr_freq_obs_rec[0][0]), pAcquisition->nrec * pFreq_group->nb_freq * sizeof(Mycomplex)) ;
	obs_file.close() ;

	// write computed data
	print_debug(ALL, MID_DEBUG, "write computed data");
	ofstream cal_file ;
	cal_file.open("COMPUTED_DATA.django.out", ios::binary | ios::app) ;
	assert(cal_file.is_open());
	cal_file.write((char*) &(pData->pr_freq_sol_rec[0][0]), pAcquisition->nrec * pFreq_group->nb_freq * sizeof(Mycomplex)) ;
	cal_file.close() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Gradient::evaluate_fcost_wavelet_and_adjoint_source");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Gradient::evaluate_fcost(Acquisition* pAcquisition, Data* pData2, Modelling* pModelling, Myfloat* fcost) 

{

	print_debug(ALL, LIGHT_DEBUG, "IN Gradient::evaluate_fcost");

	Data_std *pData = dynamic_cast<Data_std*>(pData2) ;
	if (pData == NULL)
	{
		print_error(" Only seismic data are supported ") ;
		return(RTN_CODE_KO) ;
	}

	//-----------
	// read data
	//-----------
	Rtn_code rtn_code = pData->read(pAcquisition) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	Myint isrc_current = pAcquisition->current_src - 1 ;

	//-----------
	// fcost
	//-----------
	*fcost = 0. ;
	Freq_group* pFreq_group = pModelling->get_freq_group() ;
	for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)
	{

		for (Myint irec=0; irec<pAcquisition->nrec; irec++)
		{

			//cout << "src_coef[isrc_current][ifreq] " << src_coef[isrc_current][ifreq] << "\n" ;

			// correction of observed data
			pData->pr_freq_obs_rec[irec][ifreq] = pData->pr_freq_obs_rec[irec][ifreq] * src_coef[isrc_current][ifreq] ;

			Myfloat dist = sqrt( pow(pAcquisition->xsrc-pAcquisition->xrec[irec], 2.) +  pow(pAcquisition->zsrc-pAcquisition->zrec[irec], 2.)) ;
			if ((dist < grad_offset_min) || (dist > grad_offset_max)) {
				Mycomplex res = ZERO_CMPLX ;
				//pData->pr_freq_adj_rec[irec][ifreq] = res ;
				*fcost += 0.5 * real( res * conj(res) ) ;
				continue ;
			}

			if (fcost_type == DIFF_L2)
			{
				//Mycomplex res = (pData->pr_freq_obs_rec[irec][ifreq] - pData->pr_freq_sol_rec[irec][ifreq] * src_coef[ifreq]) ;
				//pData->pr_freq_adj_rec[irec][ifreq] = res * conj(src_coef[ifreq]) ;
				Mycomplex res = (pData->pr_freq_obs_rec[irec][ifreq] - pData->pr_freq_sol_rec[irec][ifreq]) ;
				*fcost += 0.5 * real( res * conj(res) ) ;
			}
			else if (fcost_type == DIFF_L1)
			{
				//Mycomplex res = (pData->pr_freq_obs_rec[irec][ifreq] - pData->pr_freq_sol_rec[irec][ifreq] * src_coef[ifreq]) ;
				//pData->pr_freq_adj_rec[irec][ifreq] = conj(src_coef[ifreq]) * ( res / abs(res) ) ;
				Mycomplex res = (pData->pr_freq_obs_rec[irec][ifreq] - pData->pr_freq_sol_rec[irec][ifreq]) ;
				*fcost += abs(res) ;
			}
			else if (fcost_type == LOG_L2)
			{
				*fcost += 0.5 * real(log(pData->pr_freq_sol_rec[irec][ifreq] / pData->pr_freq_obs_rec[irec][ifreq]) * conj(log(pData->pr_freq_sol_rec[irec][ifreq] / pData->pr_freq_obs_rec[irec][ifreq]))) ;
			}
			else if (fcost_type == LOG_L2_AMP)
			{
				*fcost += 0.5 * pow(real(log(pData->pr_freq_sol_rec[irec][ifreq] / pData->pr_freq_obs_rec[irec][ifreq])), 2.) ;
			}
			else if (fcost_type == LOG_L2_PHASE)
			{
				*fcost += 0.5 * pow(imag(log(pData->pr_freq_sol_rec[irec][ifreq] / pData->pr_freq_obs_rec[irec][ifreq])), 2.) ;
			}

		} // for (Myint irec=0; irec<pAcquisition->nrec; irec++)

	} // for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)

	print_debug(ALL, LIGHT_DEBUG, "OUT Gradient::evaluate_fcost");
	return(RTN_CODE_OK) ;
}


//-------------------------------------------------------------------------------------------------------

Rtn_code Gradient::compute(Acquisition* pAcquisition, Data* pData, Modelling* pModelling, Model* pModel, Myfloat* fcost) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Gradient::compute");

	// return code
	Rtn_code rtn_code ;

	print_info(MASTER, "") ;
	print_line2() ;

	double t0 = MPI_Wtime() ;
	if (this->domain == TIME)
	{
		//=======================================================================================================
		//                                         TIME DOMAIN
		//=======================================================================================================

		print_info(MASTER, " START GRADIENT COMPUTATION (TIME DOMAIN)") ;
		print_info(MASTER, "") ;

		rtn_code = compute_gradient_time_domain(pAcquisition, pData, pModel, pModelling) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
	}
	else
	{
		//=======================================================================================================
		//                                         FREQ. DOMAIN
		//=======================================================================================================

		print_info(MASTER, " START GRADIENT COMPUTATION (FREQ. DOMAIN)") ;
		print_info(MASTER, "") ;

		if (type == ADJOINT_METHOD)
		{
			rtn_code = compute_gradient_freq_domain_with_adjoint_state_method(pAcquisition, pData, pModel, pModelling, fcost) ;
		}
		else if (type == FD_FCOST)
		{
			rtn_code = compute_gradient_freq_domain_with_fd_fcost(pAcquisition, pData, pModel, pModelling) ;
		}

		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
	}
	double t1 = MPI_Wtime() ;

	print_info(MASTER, " COMP. TIME GRADIENT: ", (float) (t1 - t0)) ;
	print_line2() ;

	print_info(MASTER, "") ;
	print_info(MASTER, " END OF GRADIENT COMPUTATION") ;
	print_line2() ;
	print_info(MASTER, "") ;

	// write gradient
	rtn_code = (this->pGrad_model)->write() ;
	rtn_code = (this->pPrecond)->write() ;
	rtn_code = (this->pGrad_precond)->write() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Gradient::compute");
	return(RTN_CODE_OK) ;
}

//=======================================================================================================
//
// COMPUTE FWI GRADIENT IN TIME DOMAIN
// WITH ADJOINT STATE METHOD
//
//=======================================================================================================

Rtn_code Gradient::compute_gradient_time_domain(Acquisition* pAcquisition, Data* pData, Model* pModel, Modelling* pModelling)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Gradient::compute_gradient_time_domain");

	// return code
	Rtn_code rtn_code ;

	print_info(MASTER, "") ;
	print_line2() ;
	print_info(MASTER, " START MODELLING") ;
	print_info(MASTER, "") ;

	double t0 = MPI_Wtime() ;

	// allocate grid for frequency wavefield extraction of incident wavefield
	Grid* pIncident ;
	//rtn_code = pModelling->allocate_grid_for_freq_wavefield(&pIncident) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// allocate grid for frequency wavefield extraction of adjoint wavefield
	Grid* pAdjoint ;
	//rtn_code = pModelling->allocate_grid_for_freq_wavefield(&pAdjoint) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// instanciate snapshot
	Snapshot* pSnapshot ;

	// initialize gradient
	rtn_code = initialize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// obj function
	Myfloat fcost ;

	//-------------------------------------------------------------------------------------------------------
	// start loop on shots until all have been computed
	//-------------------------------------------------------------------------------------------------------

	while ( pAcquisition->remain_shot_to_compute() )

	{

		//-------------------------------------------------------------------------------------------------------
		// retrieve src and rec positions the current shot
		//-------------------------------------------------------------------------------------------------------

		rtn_code = pAcquisition->get_position_for_current_shot() ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		//-------------------------------------------------------------------------------------------------------
		// reset data
		//-------------------------------------------------------------------------------------------------------

		rtn_code = pData->reset(pAcquisition, FWI_PROG) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		//-------------------------------------------------------------------------------------------------------
		// compute INCIDENT wavefield
		//-------------------------------------------------------------------------------------------------------

		print_info(ALL, " COMPUTE INCIDENT WAVEFIELD FOR SOURCE:", pAcquisition->current_src) ;
		double t1 = MPI_Wtime() ;

		rtn_code = pModelling->pScheme->solve_current_shot(pAcquisition, pData, INCIDENT, pSnapshot) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		//-------------------------------------------------------------------------------------------------------
		// evaluate cost function, src wavelet and build adjoint source
		//-------------------------------------------------------------------------------------------------------

		rtn_code = evaluate_fcost_wavelet_and_adjoint_source(pAcquisition, pData, pModelling, &fcost) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		//-------------------------------------------------------------------------------------------------------
		// compute ADJOINT wavefield
		//-------------------------------------------------------------------------------------------------------

		print_info(ALL, " COMPUTE ADJOINT WAVEFIELD FOR SOURCE:", pAcquisition->current_src) ;
		double t1bis = MPI_Wtime() ;

		rtn_code = pModelling->pScheme->solve_current_shot(pAcquisition, pData, ADJOINT, pSnapshot) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		//-------------------------------------------------------------------------------------------------------
		// update gradient
		//-------------------------------------------------------------------------------------------------------

		//rtn_code = pGradient->update(pIncident, pAdjoint, pModel) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		// increment shot counter

		rtn_code = pAcquisition->move_to_next_shot() ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	}

	//-------------------------------------------------------------------------------------------------------
	// end loop on shots until all have been computed
	//-------------------------------------------------------------------------------------------------------

	//-------------------------------------------------------------------------------------------------------
	// finalize modelling
	//-------------------------------------------------------------------------------------------------------

	rtn_code = pModelling->finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_info(MASTER, "") ;
	print_info(MASTER, " END MODELLING") ;
	print_info(MASTER, "") ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Gradient::compute_gradient_time_domain");
	return(RTN_CODE_OK) ;
} ;

//=======================================================================================================
//
// COMPUTE FWI GRADIENT IN FREQUENCY DOMAIN
// WITH FINITE DIFFERENCE OF THE OBJ FUNCTION
//
//=======================================================================================================


Rtn_code Gradient::compute_gradient_freq_domain_with_fd_fcost(Acquisition* pAcquisition, 
		Data* pData,
		Model* pModel,
		Modelling* pModelling)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Gradient::compute_gradient_freq_domain_with_fd_fcost");

	// return code
	Rtn_code rtn_code ;

	print_info(MASTER, "") ;
	print_line2() ;
	print_info(MASTER, " START MODELLING") ;
	print_info(MASTER, "") ;

	//-------------------------------------------------------------------------------------------------------
	// initialize modelling
	//-------------------------------------------------------------------------------------------------------

	rtn_code = pModelling->initialize(pModel) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	double t0 = MPI_Wtime() ;

	// allocate grid for frequency wavefield extraction of incident wavefield
	Grid* pIncident ;
	pIncident = NULL ;

	// initialize gradient
	rtn_code = initialize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// obj function
	Myfloat fcost0, fcost ;

	//-------------------------------------------------------------------------------------------------------
	// first step compute obj function in m0
	//-------------------------------------------------------------------------------------------------------
	fcost0 = 0. ;
	print_info(MASTER, " MODELLING IN M0") ;
	//rtn_code = compute_fcost_freq_domain(pAcquisition, pData, pModel, &fcost0) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
	print_info(MASTER, " OBJ IN M0", fcost0) ;
	print_info(MASTER, "") ;

	//-------------------------------------------------------------------------------------------------------
	// second step compute obj function at each grid point
	//-------------------------------------------------------------------------------------------------------

	// 2D
	Myint offset = 21 ;
	Myint iz = 200 + offset ;
	Myint iy = 0 ;
	Myint ixmin = 421 - offset ;
	Myint ixmax = 421 + offset ;

	// 1D
	// Myint offset = 21 ;
	// Myint ix = 0 ;
	// Myint iy = 0 ;
	// Myint izmin = 300 + offset ;
	// Myint izmax = 600 + offset ;

	// loop over the grid points
	for (Myint ix=ixmin; ix <=ixmax; ix++) // 2D
		// for (Myint iz=izmin; iz <=izmax; iz++) // 1D

	{

		// update model with perturbation on one grid point
		print_info(MASTER, " ADD PERTURBATION") ;
		print_info(MASTER, " IX ", ix) ;
		print_info(MASTER, " IY ", iy) ;
		print_info(MASTER, " IZ ", iz) ;

		rtn_code = pModel->add_perturbation(iz, ix, iy, TEMP_VP_PERTURB) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
		pModel->info() ;
		rtn_code = pModelling->initialize(pModel) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		// evaluate cost function
		fcost = 0. ;
		//rtn_code = compute_fcost_freq_domain(pAcquisition, pData, pModel, &fcost) ;

		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		print_info(MASTER, " OBJ DUE TO PERTURB.", fcost) ;

		// compute the gradient at the grid point
		Myfloat grad_val = (fcost - fcost0) / TEMP_VP_PERTURB ;
		print_info(MASTER, " GRADIENT AT THAT POINT", grad_val) ;
		print_info(MASTER, "") ;

		// store the value of the gradient
		rtn_code = store_value(iz-offset, ix-offset, iy-offset, grad_val) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		// undo the perturbation
		rtn_code = pModel->add_perturbation(iz, ix, iy, -TEMP_VP_PERTURB) ;
		//this->update_physical_coef(iz, ix, iy, -TEMP_VP_PERTURB) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	}

	// rtn_code = this->finalize() ;
	// if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_info(MASTER, "") ;
	print_info(MASTER, " END MODELLING") ;
	print_info(MASTER, "") ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Gradient::compute_gradient_freq_domain_with_fd_fcost");
	return(RTN_CODE_OK) ;
} ;


//=======================================================================================================
//
// COMPUTE FWI GRADIENT IN FREQUENCY DOMAIN
// WITH THE ADJOINT STATE METHOD
//
//=======================================================================================================


Rtn_code Gradient::compute_gradient_freq_domain_with_adjoint_state_method(Acquisition* pAcquisition, 
		Data* pData,
		Model* pModel,
		Modelling* pModelling,
		Myfloat* fcost)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Gradient::compute_gradient_freq_domain_with_adjoint_state_method");

	// return code
	Rtn_code rtn_code ;

	print_info(MASTER, "") ;
	print_line2() ;
	print_info(MASTER, " START MODELLING") ;
	print_info(MASTER, "") ;

	//-------------------------------------------------------------------------------------------------------
	// initialize modelling
	//-------------------------------------------------------------------------------------------------------

	//rtn_code = pModelling->initialize(pModel) ;
	//if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	double t0 = MPI_Wtime() ;

	// allocate grid for frequency wavefield extraction of incident wavefield
	Grid* pIncident ;
	//rtn_code = pModelling->allocate_grid_for_freq_wavefield(&pIncident) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// allocate grid for frequency wavefield extraction of adjoint wavefield
	Grid* pAdjoint ;
	//rtn_code = pModelling->allocate_grid_for_freq_wavefield(&pAdjoint) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// instanciate snapshot
	Snapshot *pSnapshot ;

	// initialize gradient
	rtn_code = initialize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// test print info
	(pGrad_model)->info() ;
	(pPrecond)->info() ;
	(pGrad_precond)->info() ;

	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	//-------------------------------------------------------------------------------------------------------
	// start loop on shots until all have been computed
	//-------------------------------------------------------------------------------------------------------

	pAcquisition->move_to_first_shot() ;
	*fcost = 0 ;

	while ( pAcquisition->remain_shot_to_compute() )

	{

		//-------------------------------------------------------------------------------------------------------
		// retrieve src and rec positions the current shot
		//-------------------------------------------------------------------------------------------------------

		rtn_code = pAcquisition->get_position_for_current_shot() ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		//-------------------------------------------------------------------------------------------------------
		// reset data
		//-------------------------------------------------------------------------------------------------------

		rtn_code = pData->reset(pAcquisition, FWI_PROG) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		//-------------------------------------------------------------------------------------------------------
		// compute INCIDENT wavefield
		//-------------------------------------------------------------------------------------------------------

		print_info(ALL, " COMPUTE INCIDENT WAVEFIELD FOR SOURCE:", pAcquisition->current_src) ;
		double t1 = MPI_Wtime() ;

		rtn_code = pModelling->pScheme->solve_current_shot(pAcquisition, pData, INCIDENT, pSnapshot) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		//-------------------------------------------------------------------------------------------------------
		// evaluate obj function and build adjoint source
		//-------------------------------------------------------------------------------------------------------
		Myfloat fcost_tmp = 0. ;
		rtn_code = evaluate_fcost_wavelet_and_adjoint_source(pAcquisition, pData, pModelling, &fcost_tmp) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
		*fcost += fcost_tmp ;

		//-------------------------------------------------------------------------------------------------------
		// compute ADJOINT wavefield
		//-------------------------------------------------------------------------------------------------------

		print_info(ALL, " COMPUTE ADJOINT WAVEFIELD FOR SOURCE:", pAcquisition->current_src) ;
		double t1bis = MPI_Wtime() ;

		rtn_code = pModelling->pScheme->solve_current_shot(pAcquisition, pData, ADJOINT, pSnapshot) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		//-------------------------------------------------------------------------------------------------------
		// update gradient
		//-------------------------------------------------------------------------------------------------------

		rtn_code = update(pIncident, pAdjoint, pModel, pModelling) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		// increment shot counter

		rtn_code = pAcquisition->move_to_next_shot() ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	}

	// deallocate grid for frequency wavefield extraction of incident wavefield
	//rtn_code = pModelling->deallocate_grid_for_freq_wavefield(pIncident) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// deallocate grid for frequency wavefield extraction of adjoint wavefield
	//rtn_code = pModelling->deallocate_grid_for_freq_wavefield(pAdjoint) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	//-------------------------------------------------------------------------------------------------------
	// end loop on shots until all have been computed
	//-------------------------------------------------------------------------------------------------------

	// sum the gradient contributions
	rtn_code = gather() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	//-------------------------------------------------------------------------------------------------------
	// finalize modelling
	//-------------------------------------------------------------------------------------------------------

	rtn_code = pModelling->finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_info(MASTER, "") ;
	print_info(MASTER, " END MODELLING") ;
	print_info(MASTER, "") ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Gradient::compute_gradient_freq_domain_with_adjoint_state_method");
	return(RTN_CODE_OK) ;
}

} // namespace django
