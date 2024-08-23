//-------------------------------------------------------------------------------------------------------
//
// MODELLING WITH FEM
//
// PARENT CLASS: Scheme
//   DERIVED CLASS: FEM
//
//-------------------------------------------------------------------------------------------------------

#include "FEM.h"

#include <cassert>
#include <fstream>
#include <iostream>

#include "allocate_array.h"
#include "constant.h"
#include "data_std.h"
#include "output_report.h"
#include "singleton.h"
#include "type_def.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------

FEM::FEM(void) : Scheme()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM::FEM");

	method = SCHEME_FEM ;

	nelem          = 0 ;
	nelem_med      = 0;
	nelem_lay      = 0 ;
	nelem_ghost    = 0 ;
	nelem_ac_iso   = 0 ;
	nelem_ac_lossy = 0 ;
	nelem_el_iso   = 0 ;
	nnode          = 0 ;

	for (Myint ii=0; ii<MAX_INV_MATRIX; ii++)
	{
		pMat_M_inv_glob_param[ii] = NULL ;
	}

	nnode_src   = 0 ;
	inode_src   = NULL ;
	wnode_src   = NULL ;

	ielem_rec   = NULL ;
	eta_rec     = NULL ;
	xi_rec      = NULL ;

	ielem_pixel = NULL ;
	eta_pixel   = NULL ;
	xi_pixel    = NULL ;

	min_z_mesh  = 0.0 ;
	max_z_mesh  = 0.0 ;
	min_x_mesh  = 0.0 ;
	max_x_mesh  = 0.0 ;

	min_dist_node = 0.0 ;
	max_dist_node = 0.0 ;

	for (Myint imat = 0; imat <= MAX_POLY_ORDER; imat++)
	{
		pMat_M_inv[imat] = NULL ;
	}

	pVec_k_glob     = NULL ;
	pVec_vp_glob    = NULL ;
	pVec_vs_glob    = NULL ;
	pVec_rho_glob   = NULL ;
	pMat_M_inv_glob = NULL ;

	node_type     = Singleton::Instance()->node_type ;
	node_distrib  = Singleton::Instance()->node_distrib ;
	node_integ    = Singleton::Instance()->node_integ ;
	flux_type     = Singleton::Instance()->flux_type ;

	pmin          = Singleton::Instance()->pmin ;
	pmax          = Singleton::Instance()->pmax ;

	nelem_x_med   = Singleton::Instance()->nelem_x_med ;
	nelem_z_med   = Singleton::Instance()->nelem_z_med ;
	prop          = Singleton::Instance()->prop ;

	nelem_x       = 0 ;
	nelem_z       = 0 ;

	pSponge_coef_node = NULL ;
	sigma_l           = 0.0 ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM::FEM");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FEM::initialize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM::initialize(void)");

	// call parent initialization
	Rtn_code rtn_code = Scheme::initialize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM::initialize(void)");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FEM::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM::finalize");

	// deallocate pMat_M_inv matrices
	for (Myint imat = 0; imat <= MAX_POLY_ORDER; imat++)
	{
		if (pMat_M_inv[imat] != NULL) delete(pMat_M_inv[imat]) ;
	}

	// delete pVec_k_glob
	if (pVec_k_glob != NULL) delete(pVec_k_glob) ;

	// delete global vectors of parameters
	if (pVec_vp_glob != NULL) delete(pVec_vp_glob) ;
	if (pVec_vs_glob != NULL) delete(pVec_vs_glob) ;
	if (pVec_rho_glob != NULL) delete(pVec_rho_glob) ;

	// delete global mass matrix
	if (pMat_M_inv_glob != NULL) delete(pMat_M_inv_glob) ;

	for (Myint ii=0; ii<MAX_INV_MATRIX; ii++)
	{
		if (pMat_M_inv_glob_param[ii] != NULL) delete(pMat_M_inv_glob_param[ii]) ;
	}

	// delete sponge coef
	if (pSponge_coef_node != NULL) delete(pSponge_coef_node) ;

	// call parent finialize
	Rtn_code rtn_code = Scheme::finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM::finalize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FEM::info(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM::info");

	Scheme::info() ;

	print_info(MASTER, "") ;
	switch(node_type)
	{
	case NODE_TYPE_EQD:
		print_info(MASTER, " Node type\t", "EQUIDISTANT") ;
		break ;
	case NODE_TYPE_GLL:
		print_info(MASTER, " Node type\t", "GAUSS-LEGENDRE-LOBATTO") ;
		break ;
	default:
		print_error(" FEM::info Invalid node type", node_type) ;
		return(RTN_CODE_KO) ;
	}

	switch(node_distrib)
	{
	case NODE_DISTRIB_UNIFORM:
		print_info(MASTER, " Node distribution", "UNIFORM") ;
		break ;
	case NODE_DISTRIB_RANDOM:
		print_info(MASTER, " Node distribution", "RANDOM") ;
		break ;
	case NODE_DISTRIB_MODEL:
		print_info(MASTER, " Node distribution", "FROM MODEL") ;
		break ;
	default:
		print_error(" FEM::info Invalid node distribution", node_distrib) ;
		return(RTN_CODE_KO) ;
	}

	switch(node_integ)
	{
	case NODE_INTEG_GL:
		print_info(MASTER, " Numerical integration", "GAUSS-LEGENDRE") ;
		break ;
	case NODE_INTEG_GLL:
		print_info(MASTER, " Numerical integration", "GAUSS-LEGENDRE-LOBATTO") ;
		break ;
	default:
		print_error(" FEM::info Invalid numerical integration", node_integ) ;
		return(RTN_CODE_KO) ;
	}

	switch(flux_type)
	{
	case FLUX_TYPE_CENTERED:
		print_info(MASTER, " Flux type\t", "CENTERED") ;
		break ;
	case FLUX_TYPE_UPWIND:
		print_info(MASTER, " Flux type\t", "UPWIND") ;
		break ;
	default:
		print_error(" FEM::info Invalid flux type", flux_type) ;
		return(RTN_CODE_KO) ;
	}

	// polynomial order
	Myint polynom_order = Singleton::Instance()->space_order ;
	if (polynom_order < 0)
	{
		print_error(" FEM::info Invalid polynomial order", polynom_order) ;
		return(RTN_CODE_KO) ;
	}
	else
	{
		print_info(MASTER, " Polynomial order", polynom_order) ;
	}

	// time order
	Myint time_order = Singleton::Instance()->time_order ;
	if (time_order < 0)
	{
		print_error(" FEM::info Invalid time order", time_order) ;
		return(RTN_CODE_KO) ;
	}
	else
	{
		print_info(MASTER, " Time order\t", time_order) ;
	}

	// adaptivity
	// only used for NODE_DISTRIB_MODEL
	if (node_distrib == NODE_DISTRIB_MODEL)
	{
		// adaptivity type
		if (adaptType == NO_ADAPT)
		{
			print_error(" FEM::info Invalid combination NODE_DISTRIB_MODEL & NO_ADAPT") ;
			return(RTN_CODE_KO) ;
		}
		else
		{
			if (adaptType == ADAPT_STATIC)
			{
				print_info(MASTER, " Adaptivity\t", "STATIC") ;
			}
			else if (adaptType == ADAPT_FREQTIME)
			{
				print_info(MASTER, " Adaptivity\t", "FREQ. VS TIME") ;

				// print list of freq.
				Myint nfreq = adaptTmin.size() ;
				print_info(MASTER, " # Frequencies\t", nfreq) ;

				if (nfreq <= 0.0)
				{
					print_error(" FEM::info nfreq should be at least >= 1") ;
					return(RTN_CODE_KO) ;
				}

				for (Myint ifreq = 0; ifreq < nfreq; ifreq++)
				{
					print_info(MASTER, " Tmin (s) / Freq. (Hz)",
							adaptTmin.at(ifreq), adaptFmax.at(ifreq)) ;

					// check frequency
					if (adaptFmax.at(ifreq) <= 0.0)
					{
						print_error(" FEM::info fmax <= 0") ;
						return(RTN_CODE_KO) ;
					}

					// check time in increasing order
					if (ifreq == 0)
					{
						if (adaptTmin.at(ifreq) <= 0.0)
						{
							print_error(" FEM::info tmin <= 0") ;
							return(RTN_CODE_KO) ;
						}
					}
					else
					{
						if (adaptTmin.at(ifreq) <= adaptTmin.at(ifreq-1))
						{
							print_error(" FEM::info tmin not in increasing order") ;
							return(RTN_CODE_KO) ;
						}
					}
				}
			}
		}

		// polynomial min and max
		if ((pmin < 0) || (pmin > MAX_POLY_ORDER))
		{
			print_error(" FEM::info Invalid pmin", pmin) ;
			return(RTN_CODE_KO) ;
		}
		else
		{
			print_info(MASTER, " Adaptivity pmin", pmin) ;
		}
		if ((pmax < 0) || (pmax > MAX_POLY_ORDER))
		{
			print_error(" FEM::info Invalid pmax", pmax) ;
			return(RTN_CODE_KO) ;
		}
		else
		{
			print_info(MASTER, " Adaptivity pmax", pmax) ;
		}
		if (pmin > pmax)
		{
			print_error(" FEM::info pmin cannot be > pmax") ;
			return(RTN_CODE_KO) ;
		}

		// fmax0
		if (fmax0 <= 0.0)
		{
			print_error(" FEM::info fmax0 <= 0") ;
			return(RTN_CODE_KO) ;
		}
		else
		{
			print_info(MASTER, " Adaptivity fmax0 (Hz)", fmax0) ;
		}
	}
	else
	{
		print_info(MASTER, " Adaptivity\t", "NONE") ;
	}

	// physical properties
	switch(prop)
	{
	case PROP_CONST:
		print_info(MASTER, " Properties/element", "CONSTANT") ;
		break ;
	case PROP_AVER:
		print_info(MASTER, " Properties/element", "AVERAGE") ;
		break ;
	case PROP_LINEAR:
		print_info(MASTER, " Properties/element", "LINEAR") ;
		break ;
	default:
		print_error(" FEM::info Invalid properties type", prop) ;
		return(RTN_CODE_KO) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM::info");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FEM::source_excitation(Myfloat* time_u, Myint it, Wavefield_type wtype, Data* pData, Myfloat src_factor)
{
	print_debug(ALL, FULL_DEBUG, "IN FEM::source_excitation");

	if (wtype == INCIDENT)
	{

		// apply source function in time domain
		// excitation at the source postion

		print_debug(ALL, FULL_DEBUG, "src it", it);
		print_debug(ALL, FULL_DEBUG, "src amp",  src_time_function[it]);

		// apply source excitation on the nodes
		for (Myint ii=0; ii<nnode_src; ii++)
		{
			time_u[inode_src[ii]] += wnode_src[ii] * src_time_function[it] * src_factor ;
		}
	}
	else if (wtype == ADJOINT)
	{

		// build excitation in time domain via an inverse DFT of the frequency adjoint terms
		// excitation at the receivers positions

		// loop over frequencies
		Mycomplex **data_adj = ((Data_std*) pData)->pr_freq_adj_rec ;
		//Myfloat coef2 = TEMP_RHO_CONST * pow(TEMP_VP_CONST, 2) * (2. / (nt * dz)) ;
		Myfloat coef2 = TEMP_RHO_CONST * pow(TEMP_VP_CONST, 2) * (2. / (nt)) ;

		for (Myint ifreq = 0; ifreq < pFreq_group->nb_freq; ifreq++)
		{
			Myfloat coef_tmp = 2. * PI * pFreq_group->pFreq_list[ifreq] * (nt - it + 1) * dt ;
			Mycomplex coef (cos(coef_tmp), sin(coef_tmp)) ;

			for (Myint irec=0; irec<nrec; irec++)
			{
				//time_u[inode_rec[irec]] += -1. * (data_adj[irec][ifreq] * coef).real() * coef2 ;
			}
		}
	}

	print_debug(ALL, FULL_DEBUG, "OUT FEM::source_excitation");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FEM::free_position_arrays(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM::free_position_arrays");

	deallocate_array<Myint>(inode_src, nnode_src) ;
	deallocate_array<Myfloat>(wnode_src, nnode_src) ;

	deallocate_array<Myint>(ielem_rec, nrec) ;
	deallocate_array<Myfloat>(eta_rec, nrec) ;
	deallocate_array<Myfloat>(xi_rec, nrec) ;

	deallocate_array<Myint>(ielem_pixel, npixel) ;
	deallocate_array<Myfloat>(eta_pixel, npixel) ;
	deallocate_array<Myfloat>(xi_pixel, npixel) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM::free_position_arrays");
	return(RTN_CODE_OK) ;
}

//-----------------------------------------------------------------------------------------
Rtn_code FEM::mesh_info(void)
{   
	print_debug(ALL, LIGHT_DEBUG, "IN FEM::mesh_info");

	print_info(MASTER, "") ;
	print_info(MASTER, " MESH STATISTICS") ;

	print_info(MASTER, " # total elements ", nelem) ;
	print_info(MASTER, " # medium elements", nelem_med) ;
	print_info(MASTER, " # layer elements ", nelem_lay) ;
	print_info(MASTER, " # ghost elements ", nelem_ghost) ;
	print_info(MASTER, " # acoustic iso. elem.", nelem_ac_iso) ;
	print_info(MASTER, " # acoustic lossy elem.", nelem_ac_lossy) ;
	print_info(MASTER, " # elastic iso. elem.", nelem_el_iso) ;
	print_info(MASTER, " Layer volume (%)", (Myfloat)((Myfloat)nelem_lay/nelem*100.)) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM::mesh_info");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FEM::initialize_node_ref(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM::initialize_node_ref");

	// loop on all possible polynomial order
	for (Myint ipol=0; ipol <= MAX_POLY_ORDER; ipol++)
	{

		// number of node = npoly + 1
		Myint ni = ipol+1 ;

		//--------------------------------------
		// compute coordinates of element nodes
		//--------------------------------------
		xnode_ref_elem_1D[ipol] = VecDoub(ni) ;

		//--------------------------------------
		// --> case equidistant node
		//--------------------------------------
		if (node_type == NODE_TYPE_EQD)
		{
			eqd(xnode_ref_elem_1D[ipol]) ;
		}
		//--------------------------------------
		// --> case GLL node
		//--------------------------------------
		else if (node_type == NODE_TYPE_GLL)
		{
			if (ni == 1)
			{
				xnode_ref_elem_1D[ipol][0] = 0.0 ;
			}
			else
			{
				gll(xnode_ref_elem_1D[ipol]) ;
			}
		}
		else
		{
			print_error(" Error in FEM::initialize_node_ref, node type not supported", node_integ) ;
			return(RTN_CODE_KO) ;
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM::initialize_node_ref");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code FEM::reset_scheme_from(FEM* schemeIn)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM::reset_scheme_from");

	// check schemes have same number of elements
	if (nelem != schemeIn->nelem)
	{
		print_error("FEM::reset_scheme_from, nelem != schemeIn.nelem") ;
		print_info(MASTER, "nelem", nelem) ;
		print_info(MASTER, "schemeIn.nelem", schemeIn->nelem) ;
		return(RTN_CODE_KO) ;
	}

	// check schemes have same dt
	if (abs(dt - schemeIn->dt) > dt * MY_EPSILON)
	{
		print_error("FEM::reset_scheme_from, dt != schemeIn.dt") ;
		print_info(MASTER, "dt", dt) ;
		print_info(MASTER, "schemeIn.dt", schemeIn->dt) ;
		return(RTN_CODE_KO) ;
	}

	//--------------------------------------------------------------------------
	// free unused data
	//--------------------------------------------------------------------------

	// free source function in new scheme
	schemeIn->free_src_time_function() ;

	// deallocate pMat_M_inv matrices in new scheme
	//---------------------------------------------
	for (Myint imat = 0; imat <= MAX_POLY_ORDER; imat++)
	{
		if (schemeIn->pMat_M_inv[imat] != NULL) delete(schemeIn->pMat_M_inv[imat]) ;
	}

	// delete global mass matrix in this scheme
	//-----------------------------------------
	if (pMat_M_inv_glob != NULL) delete(pMat_M_inv_glob) ;

	// delete pVec_k_glob in this scheme
	//----------------------------------
	if (pVec_k_glob != NULL) delete(pVec_k_glob) ;

	// delete global vectors of parameters in this scheme
	//---------------------------------------------------
	if (pVec_vp_glob != NULL) delete(pVec_vp_glob) ;
	if (pVec_vs_glob != NULL) delete(pVec_vs_glob) ;
	if (pVec_rho_glob != NULL) delete(pVec_rho_glob) ;

	for (Myint ii=0; ii<MAX_INV_MATRIX; ii++)
	{
		if (pMat_M_inv_glob_param[ii] != NULL) delete(pMat_M_inv_glob_param[ii]) ;
	}

	// delete sponge coef in this scheme
	if (pSponge_coef_node != NULL) delete(pSponge_coef_node) ;

	//--------------------------------------------------------------------------
	// reset scheme
	//--------------------------------------------------------------------------
	pMat_M_inv_glob   = schemeIn->pMat_M_inv_glob ;
	pVec_k_glob       = schemeIn->pVec_k_glob ;
	pVec_vp_glob      = schemeIn->pVec_vp_glob ;
	pVec_vs_glob      = schemeIn->pVec_vs_glob ;
	pVec_rho_glob     = schemeIn->pVec_rho_glob ;
	pSponge_coef_node = schemeIn->pSponge_coef_node ;
	for (Myint ii=0; ii<MAX_INV_MATRIX; ii++)
	{
		pMat_M_inv_glob_param[ii] = schemeIn->pMat_M_inv_glob_param[ii] ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM::reset_scheme_from");
	return(RTN_CODE_OK) ;
}

} // namespace django
