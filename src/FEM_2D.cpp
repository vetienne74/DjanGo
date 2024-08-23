//-------------------------------------------------------------------------------------------------------
//
// MODELLING WITH FEM IN 2D
//
// PARENT CLASS: Scheme
//   DERIVED CLASS: FEM
//     DERIVED CLASS: FEM_2D
//
//-------------------------------------------------------------------------------------------------------

#include "FEM_2D.h"

#include <algorithm>
#include <cassert>
#include <cfloat>
#include <climits>
#include <fstream>
#include <iostream>

#include "nr3.h"
#include "mpi.h"

#include "allocate_array.h"
#include "grid_2D_int.h"
#include "grid_2D_float.h"
#include "output_report.h"
#include "singleton.h"
#include "type_def.h"

//-------------------------------------------------------------------------------------------------------

using namespace std;

namespace django {

FEM_2D::FEM_2D(void)  : FEM()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D::FEM_2D");

	dim = TWO ;

	for (Myint imat = 0; imat <= MAX_POLY_ORDER; imat++)
	{
		pMat_Dz[imat] = NULL ;
		pMat_Dx[imat] = NULL ;
	}

	pElement = NULL ;
	pNode    = NULL ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D::FEM_2D");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FEM_2D::project_variable(FEM_2D& schemeSrc,
		Myint varSrcId, Myint varDestId)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D::project_variable");

	cout << "### varSrcId " << varSrcId << " varDestId " << varDestId << "\n" ;

	Variable* varSrc  = Singleton::Instance()->get_variable(varSrcId) ;
	if (varSrc == NULL)
	{
		print_error("FEM_2D::project_variable, varSrc == NULL") ;
		return(RTN_CODE_KO) ;
	}

	Variable* varDest = Singleton::Instance()->get_variable(varDestId) ;
	if (varDest == NULL)
	{
		print_error("FEM_2D::project_variable, varDest == NULL") ;
		return(RTN_CODE_KO) ;
	}

	// check src and dest variable type are identical
	if (varSrc->get_type() != varDest->get_type())
	{
		print_error("FEM_2D::project_variable, varSrc->get_type() != varDest->get_type()") ;
		return(RTN_CODE_KO) ;
	}

	// check number of elements in mesh src and dest are identical
	// TODO this is a current limitation
	// Only meshes with same geometry are supported
	// however, polynomials in elements can be different between src and dest
	if (this->nelem != schemeSrc.nelem)
	{
		print_error("FEM_2D::project_variable, this->nelem != schemeSrc.nelem") ;
		print_info(MASTER, "this->nelem", this->nelem) ;
		print_info(MASTER, "schemeSrc.nelem", schemeSrc.nelem) ;
		return(RTN_CODE_KO) ;
	}

	// get src grid
	Grid_1D_float* srcGrid = dynamic_cast<Grid_1D_float*>(varSrc->get_grid()) ;
	if (srcGrid == NULL)
	{
		print_error("FEM_2D::project_variable, srcGrid == NULL");
		return(RTN_CODE_KO) ;
	}
	Myfloat* srcArray = srcGrid->pArray ;
	//cout << "srcGrid->get_min()" << srcGrid->get_min() << "\n" ;
	//cout << "srcGrid->get_max()" << srcGrid->get_max() << "\n" ;

	// get dest grid
	Grid_1D_float* destGrid = dynamic_cast<Grid_1D_float*>(varDest->get_grid()) ;
	if (destGrid == NULL)
	{
		print_error("FEM_2D::project_variable, destGrid == NULL");
		return(RTN_CODE_KO) ;
	}
	Myfloat* destArray = destGrid->pArray ;
	//cout << "destGrid->get_min()" << destGrid->get_min() << "\n" ;
	//cout << "destGrid->get_max()" << destGrid->get_max() << "\n" ;

	// perform interpolation from src to dest mesh
	//--------------------------------------------

	// loop on elements in dest mesh
	for (Myint iel = 0; iel < nelem; iel++)
	{
		Myfloat znode_min = pNode[pElement[iel].inode[0]].zcoord ;
		Myfloat xnode_min = pNode[pElement[iel].inode[0]].xcoord ;

		// loop on nodes of element
		for (Myint inode=0; inode< pElement[iel].nnode ; inode++)
		{
			// get global node number
			Myint inodeGlob = pElement[iel].inode[inode] ;

			// get node coordinates
			Myfloat znode = pNode[inodeGlob].zcoord ;
			Myfloat xnode = pNode[inodeGlob].xcoord ;

			// convert to eta and xi coordinate in ref element
			Myfloat etaNode = -1.0 + 2.0*(znode - znode_min)/pElement[iel].size_z  ;
			Myfloat xiNode  = -1.0 + 2.0*(xnode - xnode_min)/pElement[iel].size_x  ;

			// call interpolate_variable
			Myfloat valNode = schemeSrc.interpolate_variable(srcArray, iel, etaNode, xiNode) ;

			// store value in destination variable
			destArray[inodeGlob] = valNode ;
		}
	}

	//cout << "after destGrid->get_min()" << destGrid->get_min() << "\n" ;
	//cout << "after destGrid->get_max()" << destGrid->get_max() << "\n" ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D::project_variable");
	return(RTN_CODE_OK) ;
}


//-------------------------------------------------------------------------------------------------------

Rtn_code FEM_2D::initialize_mesh(Model* pModel, Myfloat fmax)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D::initialize_mesh");

	// check model type is GRID
	//=========================

	if (pModel == NULL)
	{
		print_error("IN FEM_1D::initialize_mesh --> pModel == NULL");
		return(RTN_CODE_KO) ;
	}

	if (pModel->get_type() != GRID)
	{
		print_error("IN FEM_2D::initialize_mesh --> model type is not GRID");
		return(RTN_CODE_KO) ;
	}

	// get grid size and sampling
	//===========================

	Myint nz_model = pModel->get_nz() ;
	Myint nz ;
	if (nelem_z_med == 0)
	{
		// init from model
		nz = nz_model - 1 ; // nz grid points -> nz-1 elements
		nelem_z_med = nz ;
	}
	else
	{
		// init from XML
		nz = nelem_z_med ;
	}

	if (nz <= 0)
	{
		print_error("IN FEM_1D::initialize_mesh --> nz <= 0", nz);
		return(RTN_CODE_KO) ;
	}

	Myint nx_model = pModel->get_nx() ;
	Myint nx ;
	if (nelem_x_med == 0)
	{
		// init from model
		nx = nx_model - 1 ; // nx grid points -> nz-1 elements
		nelem_x_med = nx ;
	}
	else
	{
		// init from XML
		nx = nelem_x_med ;
	}

	if (nx <= 0)
	{
		print_error("IN FEM_1D::initialize_mesh --> nx <= 0", nx);
		return(RTN_CODE_KO) ;
	}

	// get boundary width
	//===================
	Myint nlayer_zBeg = get_boundary_width(ZBEG) ;
	Myint nlayer_zEnd = get_boundary_width(ZEND) ;
	Myint nlayer_xBeg = get_boundary_width(XBEG) ;
	Myint nlayer_xEnd = get_boundary_width(XEND) ;

	// set ghost elements
	//===================
	Myint nghost_scheme ;
	if (type == SCHEME_CGM)
	{
		nghost_scheme = 0 ;
	}
	else if (type == SCHEME_DGM)
	{
		nghost_scheme = 0 ;
	}
	else if (type == SCHEME_MGM)
	{
		nghost_scheme = 0 ;
	}
	Myint nghost_zBeg = nghost_scheme ;
	Myint nghost_zEnd = nghost_scheme ;
	Myint nghost_xBeg = nghost_scheme ;
	Myint nghost_xEnd = nghost_scheme ;

	// set number of elements
	//=======================
	nelem_z = nelem_z_med + nlayer_zBeg + nlayer_zEnd + nghost_zBeg + nghost_zEnd ;
	nelem_x = nelem_x_med + nlayer_xBeg + nlayer_xEnd + nghost_xBeg + nghost_xEnd ;
	nelem = nelem_z * nelem_x ;

	// allocate element array
	//=======================

	print_debug(ALL, MID_DEBUG, "allocate element array") ;
	if (nelem <= 0)
	{
		print_error("IN FEM_2D::initialize_mesh --> nelem <= 0");
		return(RTN_CODE_KO) ;
	}
	pElement = allocate_array<Quad_struct_type>(nelem) ;

	// retrieve model parameters
	//==========================

	// get VP model
	// needed for ACOUSTIC and ELASTIC
	Variable* vp_var = NULL ;
	Grid* vp_grid = NULL ;
	Grid_2D_float *vp_2D_grid = NULL ;
	vp_var = pModel->get_parameter(VP) ;
	if (vp_var == NULL)
	{
		print_error("IN FEM_2D::initialize_mesh --> VP model not found");
		return(RTN_CODE_KO) ;
	}
	vp_grid = vp_var->get_grid() ;
	if (vp_grid == NULL)
	{
		print_error("IN FEM_2D::initialize_mesh --> VP grid not initialized");
		return(RTN_CODE_KO) ;
	}
	vp_2D_grid = dynamic_cast<Grid_2D_float*>(vp_grid) ;
	if (vp_2D_grid == NULL)
	{
		print_error("IN FEM_2D::initialize_mesh --> VP grid is not Grid_2D_float");
		return(RTN_CODE_KO) ;
	}

	// get VS model
	// needed for ACOUSTIC and ELASTIC
	Variable* vs_var = NULL ;
	Grid* vs_grid = NULL ;
	Grid_2D_float *vs_2D_grid = NULL ;
	if (eq_type == ELASTIC)
	{
		vs_var = pModel->get_parameter(VS) ;
		if (vs_var == NULL)
		{
			print_error("IN FEM_2D::initialize_mesh --> VS model not found");
			return(RTN_CODE_KO) ;
		}
		vs_grid = vs_var->get_grid() ;
		if (vs_grid == NULL)
		{
			print_error("IN FEM_2D::initialize_mesh --> VS grid not initialized");
			return(RTN_CODE_KO) ;
		}
		vs_2D_grid = dynamic_cast<Grid_2D_float*>(vs_grid) ;
		if (vs_2D_grid == NULL)
		{
			print_error("IN FEM_2D::initialize_mesh --> VS grid is not Grid_2D_float");
			return(RTN_CODE_KO) ;
		}
	}

	// get RHO model
	// needed for ACOUSTIC and ELASTIC
	Variable* rho_var = NULL ;
	Grid* rho_grid = NULL ;
	Grid_2D_float *rho_2D_grid = NULL ;
	if (eq_type == ELASTIC)
	{
		rho_var = pModel->get_parameter(RHO) ;
		if (rho_var == NULL)
		{
			print_error("IN FEM_2D::initialize_mesh --> RHO model not found");
			return(RTN_CODE_KO) ;
		}
		rho_grid = rho_var->get_grid() ;
		if (rho_grid == NULL)
		{
			print_error("IN FEM_2D::initialize_mesh --> RHO grid not initialized");
			return(RTN_CODE_KO) ;
		}
		rho_2D_grid = dynamic_cast<Grid_2D_float*>(rho_grid) ;
		if (rho_2D_grid == NULL)
		{
			print_error("IN FEM_2D::initialize_mesh --> RHO grid is not Grid_2D_float");
			return(RTN_CODE_KO) ;
		}
	}

	// determine element region
	//=========================

	// $$$ tmp, does not allow variable grid spacing including random spacing
	Myfloat el_sizez = (pModel->get_zcoord()->get_max() - pModel->get_zcoord()->get_min()) / nz ;
	Myfloat el_sizex = (pModel->get_xcoord()->get_max() - pModel->get_xcoord()->get_min()) / nx ;

	Myfloat dz_model = pModel->get_dz() ;
	Myfloat dx_model = pModel->get_dx() ;

	print_debug(ALL, MID_DEBUG, "determine element region") ;
	Myint nelem_no_region = 0 ;

	Myint izGhostBeg = 0 ;
	Myint izLayerBeg = izGhostBeg + nghost_zBeg ;
	Myint izMediumBeg = izLayerBeg + nlayer_zBeg ;
	Myint izMediumEnd = izMediumBeg + nz-1 ;
	Myint izLayerEnd = izMediumEnd + nlayer_zEnd ;
	Myint izGhostEnd = izLayerEnd + nghost_zEnd ;
	Myint ixGhostBeg = 0 ;
	Myint ixLayerBeg = ixGhostBeg + nghost_xBeg ;
	Myint ixMediumBeg = ixLayerBeg + nlayer_xBeg ;
	Myint ixMediumEnd = ixMediumBeg + nx-1 ;
	Myint ixLayerEnd = ixMediumEnd + nlayer_xEnd ;
	Myint ixGhostEnd = ixLayerEnd + nghost_xEnd ;

	Myint iel = -1 ;
	for (Myint iz = 0; iz < nelem_z; iz++)
	{
		for (Myint ix = 0; ix < nelem_x; ix++)
		{
			iel++ ;

			// determine region
			if ((iz >= izMediumBeg) && (iz<= izMediumEnd) && (ix >= ixMediumBeg) && (ix<= ixMediumEnd))
			{
				nelem_med++ ;
				pElement[iel].region = MEDIUM ;
			}
			else if ((iz >= izLayerBeg) && (iz<= izLayerEnd) && (ix >= ixLayerBeg) && (ix<= ixLayerEnd))
			{
				nelem_lay++ ;
				pElement[iel].region = LAYER ;
			}
			else if ((iz >= izGhostBeg) && (iz<= izGhostEnd) && (ix >= ixGhostBeg) && (ix<= ixGhostEnd))
			{
				nelem_ghost++ ;
				pElement[iel].region = GHOST ;
			}
			else
			{
				nelem_no_region++ ;
				pElement[iel].region = NO_REGION ;
			}
		}
	}

	// check region distribution
	if (nelem != (nelem_ghost + nelem_lay + nelem_med))
	{
		print_error(" Error in FEM_2D::initialize_mesh, error in region distribution", nelem_no_region) ;
		return(RTN_CODE_KO) ;
	}

	// retrieve min velocity
	// this value will be used for p-adaptivity
	//=========================================
	print_debug(ALL, MID_DEBUG, "retrieve min velocity") ;

	// determine mesh origin coordinate
	// $$$ tmp, does not allow variable grid spacing including random spacing
	const Myfloat z_origin_mesh =  0.0 - el_sizez * (nlayer_zBeg + nghost_zBeg) ;
	const Myfloat x_origin_mesh =  0.0 - el_sizex * (nlayer_xBeg + nghost_xBeg) ;

	iel = -1 ;
	for (Myint iz = 0; iz < nelem_z; iz++)
	{
		for (Myint ix = 0; ix < nelem_x; ix++)
		{
			iel++ ;

			// determine min velocity per element
			//if (pElement[iel].region == MEDIUM)
			{
				Myint ix_model = (x_origin_mesh + ((ix + 0) * el_sizex)) / dx_model ;
				Myint iz_model = (z_origin_mesh + ((iz + 0) * el_sizez)) / dz_model ;
				ix_model = max(0,  ix_model) ;
				ix_model = min(nx_model-1, ix_model) ;
				iz_model = max(0,  iz_model) ;
				iz_model = min(nz_model-1, iz_model) ;
				Myfloat v1 = vp_2D_grid->pArray[ix_model][iz_model] ;

				ix_model = (x_origin_mesh + ((ix + 1) * el_sizex)) / dx_model ;
				iz_model = (z_origin_mesh + ((iz + 0) * el_sizez)) / dz_model ;
				ix_model = max(0,  ix_model) ;
				ix_model = min(nx_model-1, ix_model) ;
				iz_model = max(0,  iz_model) ;
				iz_model = min(nz_model-1, iz_model) ;
				Myfloat v2 = vp_2D_grid->pArray[ix_model][iz_model] ;

				ix_model = (x_origin_mesh + ((ix + 0) * el_sizex)) / dx_model ;
				iz_model = (z_origin_mesh + ((iz + 1) * el_sizez)) / dz_model ;
				ix_model = max(0,  ix_model) ;
				ix_model = min(nx_model-1, ix_model) ;
				iz_model = max(0,  iz_model) ;
				iz_model = min(nz_model-1, iz_model) ;
				Myfloat v3 = vp_2D_grid->pArray[ix_model][iz_model] ;

				ix_model = (x_origin_mesh + ((ix + 1) * el_sizex)) / dx_model ;
				iz_model = (z_origin_mesh + ((iz + 1) * el_sizez)) / dz_model ;
				ix_model = max(0,  ix_model) ;
				ix_model = min(nx_model-1, ix_model) ;
				iz_model = max(0,  iz_model) ;
				iz_model = min(nz_model-1, iz_model) ;
				Myfloat v4 = vp_2D_grid->pArray[ix_model][iz_model] ;

				ix_model = (x_origin_mesh + ((ix + 0.5) * el_sizex)) / dx_model ;
				iz_model = (z_origin_mesh + ((iz + 0.5) * el_sizez)) / dz_model ;
				ix_model = max(0,  ix_model) ;
				ix_model = min(nx_model-1, ix_model) ;
				iz_model = max(0,  iz_model) ;
				iz_model = min(nz_model-1, iz_model) ;
				Myfloat v0 = vp_2D_grid->pArray[ix_model][iz_model] ;

				// retrieve vp_min
				// ACOUSTIC and ELASTIC
				//---------------------

				Myfloat vp_min = 0.0 ;
				Myfloat vp_max = 0.0 ;

				// PROP_CONST
				if (prop == PROP_CONST)
				{
					// take velocity at cell barycenter
					vp_min = v0 ;
					vp_max = vp_min ;
				}

				// PROP_AVER
				else if (prop == PROP_AVER)
				{
					// average of values at 4 corners
					vp_min = (v1 + v2 + v3 + v4) / 4.0 ;
					vp_max = vp_min ;
				}

				// PROP_LINEAR
				else if (prop == PROP_LINEAR)
				{
					// take minimum of values at 4 corners
					vp_min = min(v1,v2) ;
					vp_min = min(vp_min,v3) ;
					vp_min = min(vp_min,v3) ;

					// take maximum of values at 4 corners
					vp_max = max(v1,v2) ;
					vp_max = max(vp_max,v3) ;
					vp_max = max(vp_max,v3) ;
				}

				// retrieve vs_min
				// ONLY for ELASTIC
				//-----------------

				Myfloat vs_min = 0.0 ;
				if (eq_type == ELASTIC)
				{
					Myint ix_model = (x_origin_mesh + ((ix + 0) * el_sizex)) / dx_model ;
					Myint iz_model = (z_origin_mesh + ((iz + 0) * el_sizez)) / dz_model ;
					ix_model = max(0,  ix_model) ;
					ix_model = min(nx_model-1, ix_model) ;
					iz_model = max(0,  iz_model) ;
					iz_model = min(nz_model-1, iz_model) ;
					Myfloat v1 = vs_2D_grid->pArray[ix_model][iz_model] ;

					ix_model = (x_origin_mesh + ((ix + 1) * el_sizex)) / dx_model ;
					iz_model = (z_origin_mesh + ((iz + 0) * el_sizez)) / dz_model ;
					ix_model = max(0,  ix_model) ;
					ix_model = min(nx_model-1, ix_model) ;
					iz_model = max(0,  iz_model) ;
					iz_model = min(nz_model-1, iz_model) ;

					Myfloat v2 = vs_2D_grid->pArray[ix_model][iz_model] ;
					ix_model = (x_origin_mesh + ((ix + 0) * el_sizex)) / dx_model ;
					iz_model = (z_origin_mesh + ((iz + 1) * el_sizez)) / dz_model ;
					ix_model = max(0,  ix_model) ;
					ix_model = min(nx_model-1, ix_model) ;
					iz_model = max(0,  iz_model) ;
					iz_model = min(nz_model-1, iz_model) ;
					Myfloat v3 = vs_2D_grid->pArray[ix_model][iz_model] ;

					ix_model = (x_origin_mesh + ((ix + 1) * el_sizex)) / dx_model ;
					iz_model = (x_origin_mesh + ((iz + 1) * el_sizez)) / dz_model ;
					ix_model = max(0,  ix_model) ;
					ix_model = min(nx_model-1, ix_model) ;
					iz_model = max(0,  iz_model) ;
					iz_model = min(nz_model-1, iz_model) ;
					Myfloat v4 = vs_2D_grid->pArray[ix_model][iz_model] ;

					ix_model = (x_origin_mesh + ((ix + 0.5) * el_sizex)) / dx_model ;
					iz_model = (z_origin_mesh + ((iz + 0.5) * el_sizez)) / dz_model ;
					ix_model = max(0,  ix_model) ;
					ix_model = min(nx_model-1, ix_model) ;
					iz_model = max(0,  iz_model) ;
					iz_model = min(nz_model-1, iz_model) ;
					Myfloat v0 = vs_2D_grid->pArray[ix_model][iz_model] ;

					// PROP_CONST
					if (prop == PROP_CONST)
					{
						// take velocity at cell barycenter
						vs_min = v0 ;
					}

					// PROP_AVER
					if (prop == PROP_AVER)
					{
						// average of values at 4 corners
						vs_min = (v1 + v2 + v3 + v4) / 4.0 ;
					}

					// PROP_LINEAR
					else if (prop == PROP_LINEAR)
					{
						// take minimum of values at 4 corners
						vs_min = min(v1,v2) ;
						vs_min = min(vs_min,v3) ;
						vs_min = min(vs_min,v3) ;
					}
				}

				pElement[iel].vmax = vp_max ;
				if (eq_type == ACOUSTIC)
				{
					nelem_ac_iso++ ;
					pElement[iel].vmin = vp_min ;
				}
				if (eq_type == AC_LOSSY)
				{
					nelem_ac_lossy++ ;
					pElement[iel].vmin = vp_min ;
				}
				else if (eq_type == ELASTIC)
				{
					// fluid region in elastic modelling
					if (vs_min == 0.0)
					{
						nelem_ac_iso++ ;
						pElement[iel].vmin = vp_min ;
					}
					else
					{
						nelem_el_iso++ ;
						pElement[iel].vmin = vs_min ;
					}
				}
			}
		}
	}

	// check ac/el distribution
	if (nelem != (nelem_ac_iso + nelem_ac_lossy + nelem_el_iso))
	{
		print_error(" Error in FEM_2D::initialize_mesh, error in region ac/el distribution") ;
		return(RTN_CODE_KO) ;
	}

	// determine neighbour elements
	//=============================

	print_debug(ALL, MID_DEBUG, "determine neighbour elements") ;
	iel = -1 ;
	for (Myint iz = 0; iz < nelem_z; iz++)
	{
		for (Myint ix = 0; ix < nelem_x; ix++)
		{
			iel++ ;

			if (iz == izGhostBeg)
			{
				pElement[iel].neigh[I_ZPREV] = NO_NEIGH ;
			}
			else
			{
				pElement[iel].neigh[I_ZPREV] = iel - nelem_x ;
			}

			if (iz == izGhostEnd)
			{
				pElement[iel].neigh[I_ZNEXT] = NO_NEIGH ;
			}
			else
			{
				pElement[iel].neigh[I_ZNEXT] = iel + nelem_x ;
			}

			if (ix == ixGhostBeg)
			{
				pElement[iel].neigh[I_XPREV] = NO_NEIGH ;
			}
			else
			{
				pElement[iel].neigh[I_XPREV] = iel - 1 ;
			}

			if (ix == ixGhostEnd)
			{
				pElement[iel].neigh[I_XNEXT] = NO_NEIGH ;
			}
			else
			{
				pElement[iel].neigh[I_XNEXT] = iel + 1 ;
			}
		}
	}

	// initialize size of elements
	//============================

	print_debug(ALL, MID_DEBUG, "initialize size of elements") ;

	// get pointer to xcoord array
	if (pModel->get_xcoord() == NULL)
	{
		print_error(" Error in FEM_2D::initialize_mesh, pModel->get_xcoord() is NULL") ;
		return(RTN_CODE_KO) ;
	}
	Myfloat* xcoord = pModel->get_xcoord()->pArray ;
	if (xcoord == NULL)
	{
		print_error(" Error in FEM_2D::initialize_mesh, xcoord is NULL") ;
		return(RTN_CODE_KO) ;
	}

	// get pointer to zcoord array
	if (pModel->get_zcoord() == NULL)
	{
		print_error(" Error in FEM_2D::initialize_mesh, pModel->get_zcoord() is NULL") ;
		return(RTN_CODE_KO) ;
	}
	Myfloat* zcoord = pModel->get_zcoord()->pArray ;
	if (zcoord == NULL)
	{
		print_error(" Error in FEM_2D::initialize_mesh, zcoord is NULL") ;
		return(RTN_CODE_KO) ;
	}

	// loop on all elements
	Myint ielem = -1 ;
	for (Myint iz = 0; iz < nelem_z; iz++)
	{
		for (Myint ix = 0; ix < nelem_x; ix++)
		{
			ielem++ ;
			// case model type is GRID
			// elements are defined by 2 consecutive points in the model
			//pElement[ielem].size_x = xcoord[ix+1] - xcoord[ix] ;
			//pElement[ielem].size_z = zcoord[iz+1] - zcoord[iz] ;

			pElement[ielem].size_x = el_sizex ;
			pElement[ielem].size_z = el_sizez ;
		}
	}

	// initialize number of nodes per element
	//=======================================

	print_debug(ALL, MID_DEBUG, "initialize number of nodes per element") ;

	// check node distribution
	if ((node_distrib != NODE_DISTRIB_UNIFORM) &&
			(node_distrib != NODE_DISTRIB_MODEL) &&
			(node_distrib != NODE_DISTRIB_RANDOM))
	{
		print_error(" Error in FEM_2D::initialize_mesh, unsupported node_distrib", node_distrib) ;
		return(RTN_CODE_KO) ;
	}

	// check number of node per element
	Myint polynom_order = Singleton::Instance()->space_order ;
	if ( (polynom_order < 0) || (polynom_order > MAX_POLY_ORDER) )
	{
		print_error(" Error in FEM_2D::initialize_mesh, invalid polynomial order ", polynom_order) ;
		return(RTN_CODE_KO) ;
	}

	// special case
	// for CGM, one node element is not allowed
	if ((type == SCHEME_CGM) && (polynom_order == 0))
	{
		print_error(" Error in FEM_2D::initialize_mesh, one node element not allowed for CGM") ;
		return(RTN_CODE_KO) ;
	}

	// special case
	// for CGM, random number of node element is not yet supported
	if ((type == SCHEME_CGM) && (node_distrib == NODE_DISTRIB_RANDOM))
	{
		print_error(" Error in FEM_2D::initialize_mesh, random node distribution not yet supported for CGM") ;
		return(RTN_CODE_KO) ;
	}

	// loop on all elements
	for (Myint ielem = 0; ielem < nelem; ielem++)
	{
		// uniform distribution
		if (node_distrib == NODE_DISTRIB_UNIFORM)
		{
			// number of nodes = polynom order + 1
			pElement[ielem].nnode_1D = polynom_order + 1 ;
		}

		// random distribution
		else if (node_distrib == NODE_DISTRIB_RANDOM)
		{
			// for CGM, one node element node is not allowed
			if (type == SCHEME_CGM)
			{
				// RANDOM not yet supported for 2D continuous formulation
				print_error(" Error in FEM_2D::initialize_mesh, RANDOM not suported for 2D CGM") ;
				return(RTN_CODE_KO) ;
			}
			else
			{
				Myint nnode_1D_min = pmin + 1 ;
				Myint nnode_1D_max = pmax + 1 ;
				Myint nnode_1D_delta  = nnode_1D_max - nnode_1D_min ;
				pElement[ielem].nnode_1D = rand() % (nnode_1D_delta+1) + nnode_1D_min ;
				if ((pElement[ielem].nnode_1D < 0) || (pElement[ielem].nnode_1D > MAX_POLY_ORDER+1))
				{
					print_error(" Error in FEM_2D::initialize_mesh, invalid nnode_1D", pElement[ielem].nnode_1D) ;
					return(RTN_CODE_KO) ;
				}
			}
		}

		// local adaptivity
		else if (node_distrib == NODE_DISTRIB_MODEL)
		{
			// determine nelem / lambda min
			Myfloat min_size = min(pElement[ielem].size_x, pElement[ielem].size_z) ;
			Myfloat nelperlamb = (pElement[ielem].vmin / fmax) / min_size ;

			if (nelperlamb >= P0_RATIO)
			{
				// P0
				pElement[ielem].nnode_1D = 0 + 1 ;
			}
			else if (nelperlamb >= P1_RATIO)
			{
				// P1
				pElement[ielem].nnode_1D = 1 + 1 ;
			}
			else if (nelperlamb >= P2_RATIO)
			{
				// P2
				pElement[ielem].nnode_1D = 2 + 1 ;
			}
			else if (nelperlamb >= P3_RATIO)
			{
				// P3
				pElement[ielem].nnode_1D = 3 + 1 ;
			}
			else if (nelperlamb >= P4_RATIO)
			{
				// P4
				pElement[ielem].nnode_1D = 4 + 1 ;
			}
			else if (nelperlamb >= P5_RATIO)
			{
				// P5
				pElement[ielem].nnode_1D = 5 + 1 ;
			}
			else
			{
				// P6
				pElement[ielem].nnode_1D = 6 + 1 ;
			}
		}

		// update total number of nodes
		// in 2D, nnode = nnode_1D * nnode_1D
		pElement[ielem].nnode = pElement[ielem].nnode_1D * pElement[ielem].nnode_1D ;
	}

	// determine total number of nodes in the mesh
	//============================================

	print_debug(ALL, MID_DEBUG, "determine total number of nodes in the mesh") ;
	nnode = 0 ;
	ielem = -1 ;

	for (Myint iz = 0; iz < nelem_z; iz++)
	{
		for (Myint ix = 0; ix < nelem_x; ix++)
		{
			ielem++ ;
			Myint nnode_x, nnode_z ;

			// DGM case
			// nodes are not shared between elements
			if ((type == SCHEME_DGM) || (type == SCHEME_MGM))
			{
				nnode_x = pElement[ielem].nnode_1D ;
				nnode_z = pElement[ielem].nnode_1D ;
			}

			// CGM case
			// nodes are shared between elements
			else if (type == SCHEME_CGM)
			{
				if (pElement[ielem].neigh[I_ZPREV] == NO_NEIGH)
				{
					nnode_z = pElement[ielem].nnode_1D ;
				}
				else
				{
					// one node is shared with previous element (-1)
					nnode_z = pElement[ielem].nnode_1D - 1 ;
				}
				if (pElement[ielem].neigh[I_XPREV] == NO_NEIGH)
				{
					nnode_x = pElement[ielem].nnode_1D ;
				}
				else
				{
					// one node is shared with previous element (-1)
					nnode_x = pElement[ielem].nnode_1D - 1 ;
				}
			}

			// add to global node number
			nnode += nnode_x * nnode_z ;
		}
	}

	// allocate global node array
	//===========================

	print_debug(ALL, MID_DEBUG, "allocate global node array") ;
	if (nnode <= 0)
	{
		print_error("IN FEM_2D::initialize_mesh --> nnode <= 0");
		return(RTN_CODE_KO) ;
	}
	pNode = allocate_array<Node_2D_struct_type>(nnode) ;

	// check node type
	if ((node_type != NODE_TYPE_EQD) &&
			(node_type != NODE_TYPE_GLL))
	{
		print_error(" Error in FEM_2D::initialize_mesh, unsupported node_type", node_type) ;
		return(RTN_CODE_KO) ;
	}

	// initialize global node array
	// compute node coordinates
	//=============================

	print_debug(ALL, MID_DEBUG, "initialize global node array") ;

	Myint inode_global = 0 ;
	Myint inode_local  = 0 ;
	ielem = -1 ;

	Myfloat z_origin_el = z_origin_mesh ;

	for (Myint iz = 0; iz < nelem_z; iz++)
	{
		Myfloat x_origin_el =  x_origin_mesh ;
		for (Myint ix = 0; ix < nelem_x; ix++)
		{
			ielem++ ;

			// store number of node along one direction
			Myint nnode_z = pElement[ielem].nnode_1D ;
			Myint nnode_x = nnode_z ;

			// allocate local node array
			pElement[ielem].inode = allocate_array<Myint>(pElement[ielem].nnode) ;

			// get coordinates in reference element
			Myint ng = nnode_z ;
			VecDoub xg = xnode_ref_elem_1D[ng-1] ;
			VecDoub zg = xnode_ref_elem_1D[ng-1] ;

			// assign global nodes to element
			//===============================

			inode_local = 0 ;
			for (Myint izn = 0; izn < nnode_z ; izn++)
			{

				// detect new node
				bool new_node_z = false ;

				if ((type == SCHEME_DGM) || (type == SCHEME_MGM))
				{
					// always new node
					new_node_z = true ;
				}
				else if (type == SCHEME_CGM)
				{
					if (izn == 0)
					{
						if (pElement[ielem].neigh[I_ZPREV] == NO_NEIGH)
						{
							new_node_z = true ;
						}
					}
					else
					{
						new_node_z = true ;
					}
				}

				for (Myint ixn = 0; ixn < nnode_x ; ixn++)
				{

					// detect new node
					bool new_node_x = false ;

					if ((type == SCHEME_DGM) || (type == SCHEME_MGM))
					{
						// always new node
						new_node_x = true ;
					}
					else if (type == SCHEME_CGM)
					{
						if (ixn == 0)
						{
							if (pElement[ielem].neigh[I_XPREV] == NO_NEIGH)
							{
								new_node_x = true ;
							}
						}
						else
						{
							new_node_x = true ;
						}
					}

					// assign new node
					if (new_node_x && new_node_z)
					{
						// new node
						if ((inode_global > nnode) || (inode_global < 0))
						{
							print_error("IN FEM_2D::initialize_mesh --> invalid inode_global", inode_global);
							return(RTN_CODE_KO) ;
						}

						// assign global node index
						pElement[ielem].inode[inode_local] = inode_global ;

						// node coordinates
						pNode[inode_global].zcoord = z_origin_el + (zg[izn] + 1.0) * pElement[ielem].size_z / 2. ;
						pNode[inode_global].xcoord = x_origin_el + (xg[ixn] + 1.0) * pElement[ielem].size_x / 2. ;

						// increment global node index
						inode_global++ ;
					}

					// shared node with neighbour element
					else
					{
						if (izn == 0)
						{
							// neighbour element
							Myint ielem_neigh = pElement[ielem].neigh[I_ZPREV] ;
							if (ielem_neigh != NO_NEIGH)
							{
								// corresponding node in neighbour element
								Myint ixn_neigh = ixn ;
								Myint izn_neigh = nnode_z - 1 ;
								Myint inode_neigh = izn_neigh * nnode_x + ixn_neigh ;
								pElement[ielem].inode[inode_local] = pElement[ielem_neigh].inode[inode_neigh] ;
							}
						}

						if (ixn == 0)
						{
							// neighbour element
							Myint ielem_neigh = pElement[ielem].neigh[I_XPREV] ;
							if (ielem_neigh != NO_NEIGH)
							{
								// corresponding node in neighbour element
								Myint ixn_neigh = nnode_x - 1 ;
								Myint izn_neigh = izn ;
								Myint inode_neigh = izn_neigh * nnode_x + ixn_neigh ;
								pElement[ielem].inode[inode_local] = pElement[ielem_neigh].inode[inode_neigh] ;
							}
						}
					}

					// increment local node index
					inode_local++ ;

				} // for (Myint ixn = 0; ixn < nnode_x ; ixn++)
			} // for (Myint izn = 0; izn < nnode_z ; izn++)

			// update element origin
			x_origin_el += el_sizex ;

		} // for (Myint ix = 0; ix < nelem_x; ix++)

		// update element origin
		z_origin_el += el_sizez ;

	} // for (Myint iz = 0; iz < nelem_z; iz++)

	// check all nodes have been assigned
	if (inode_global != nnode)
	{
		print_error(" Error in FEM_2D::initialize_mesh, inode_global != nnode", inode_global) ;
		return(RTN_CODE_KO) ;
	}

	// assign boundary type for each node
	//-----------------------------------
	print_debug(ALL, MID_DEBUG, "assign boundary type for each node") ;

	for (Myint ielem = 0; ielem < nelem; ielem++)
	{
		Myint nnode_z = pElement[ielem].nnode_1D ;
		Myint nnode_x = nnode_z ;
		Myint inode = 0 ;
		for (Myint izn = 0; izn < nnode_z ; izn++)
		{
			for (Myint ixn = 0; ixn < nnode_x ; ixn++)
			{

				// $$$ tmp: does not allow mix of boundaries
				// boundary: NO_BOUNDARY, PML, RANDOM, SPG, FREESURF, RIGID, INTERNAL
				// region: NO_REGION, GHOST, LAYER, MEDIUM
				Myint inode_global = pElement[ielem].inode[inode] ;

				// free surface is mandatory at the edges of the mesh
				// if not, scheme is unstable
				if ((pElement[ielem].neigh[I_ZPREV] == NO_NEIGH) && (izn == 0))
				{
					pNode[inode_global].boundary = FREESURF ;
				}
				else if ((pElement[ielem].neigh[I_ZNEXT] == NO_NEIGH) && (izn == nnode_z - 1))
				{
					pNode[inode_global].boundary = FREESURF ;
				}
				else if ((pElement[ielem].neigh[I_XPREV] == NO_NEIGH) && (ixn == 0))
				{
					pNode[inode_global].boundary = FREESURF ;
				}
				else if ((pElement[ielem].neigh[I_XNEXT] == NO_NEIGH) && (ixn == nnode_x - 1))
				{
					pNode[inode_global].boundary = FREESURF ;
				}

				// node is not at free surface
				else
				{
					if (pElement[ielem].region == MEDIUM)
					{
						pNode[inode_global].boundary = INTERNAL ;

					}
					else if (pElement[ielem].region == LAYER)
					{
						pNode[inode_global].boundary = SPG ;
					}
					else
					{
						pNode[inode_global].boundary = NO_BOUNDARY ;
					}
				}

				inode++ ;

			} // for (Myint ixn = 0; ixn < nnode_x ; ixn++)
		} // for (Myint izn = 0; izn < nnode_z ; izn++)
	} // for (Myint ielem = 0; ielem < nelem; ielem++)

	// flux type
	//----------
	for (Myint ielem = 0; ielem < nelem; ielem++)
	{
		for (Myint ineigh = 0; ineigh < 4; ineigh++)
		{
			if (type == SCHEME_CGM)
			{
				pElement[ielem].flux[ineigh] = NO_FLUX_TYPE ;
			}
			else
			{
				pElement[ielem].flux[ineigh] = FLUX_TYPE_CENTERED ;
			}
		}
	}

	// assign VP per node
	// for ACOUSTIC and ELASTIC
	//----------------------------------------------------------------------------------

	print_debug(ALL, MID_DEBUG, "assign VP per node") ;

	// initialize pVec_vp_glob
	pVec_vp_glob = new Grid_1D_float(nnode, 0.) ;
	Myfloat *vp_glob = pVec_vp_glob->pArray ;

	iel = -1 ;
	for (Myint iz = 0; iz < nelem_z; iz++)
	{
		for (Myint ix = 0; ix < nelem_x; ix++)
		{
			iel++ ;

			// determine vp
			//if (pElement[iel].region == MEDIUM)
			{
				// retrieve velocity at 4 corners of element
				// and take minimum value

				Myint ix_model = (x_origin_mesh + ((ix + 0) * el_sizex)) / dx_model ;
				Myint iz_model = (z_origin_mesh + ((iz + 0) * el_sizez)) / dz_model ;
				ix_model = max(0,  ix_model) ;
				ix_model = min(nx_model-1, ix_model) ;
				iz_model = max(0,  iz_model) ;
				iz_model = min(nz_model-1, iz_model) ;
				Myfloat v1 = vp_2D_grid->pArray[ix_model][iz_model] ;

				ix_model = (x_origin_mesh + ((ix + 1) * el_sizex)) / dx_model ;
				iz_model = (z_origin_mesh + ((iz + 0) * el_sizez)) / dz_model ;
				ix_model = max(0,  ix_model) ;
				ix_model = min(nx_model-1, ix_model) ;
				iz_model = max(0,  iz_model) ;
				iz_model = min(nz_model-1, iz_model) ;
				Myfloat v2 = vp_2D_grid->pArray[ix_model][iz_model] ;

				ix_model = (x_origin_mesh + ((ix + 0) * el_sizex)) / dx_model ;
				iz_model = (z_origin_mesh + ((iz + 1) * el_sizez)) / dz_model ;
				ix_model = max(0,  ix_model) ;
				ix_model = min(nx_model-1, ix_model) ;
				iz_model = max(0,  iz_model) ;
				iz_model = min(nz_model-1, iz_model) ;
				Myfloat v3 = vp_2D_grid->pArray[ix_model][iz_model] ;

				ix_model = (x_origin_mesh + ((ix + 1) * el_sizex)) / dx_model ;
				iz_model = (z_origin_mesh + ((iz + 1) * el_sizez)) / dz_model ;
				ix_model = max(0,  ix_model) ;
				ix_model = min(nx_model-1, ix_model) ;
				iz_model = max(0,  iz_model) ;
				iz_model = min(nz_model-1, iz_model) ;
				Myfloat v4 = vp_2D_grid->pArray[ix_model][iz_model] ;

				// assign vp to each node
				for (Myint inode = 0; inode < pElement[iel].nnode ; inode++)
				{
					// PROP_CONST
					if (prop == PROP_CONST)
					{
						ix_model = (x_origin_mesh + ((ix + 0.5) * el_sizex)) / dx_model ;
						iz_model = (z_origin_mesh + ((iz + 0.5) * el_sizez)) / dz_model ;
						ix_model = max(0,  ix_model) ;
						ix_model = min(nx_model-1, ix_model) ;
						iz_model = max(0,  iz_model) ;
						iz_model = min(nz_model-1, iz_model) ;
						Myfloat v0 = vp_2D_grid->pArray[ix_model][iz_model] ;
						vp_glob[ pElement[iel].inode[inode] ] = v0 ;
					}
					// PROP_AVER
					else if (prop == PROP_AVER)
					{
						vp_glob[ pElement[iel].inode[inode] ] = (v1 + v2 + v3 + v4) / 4.0 ;
					}
					// linear properties per element
					// PROP_LINEAR
					else if (prop == PROP_LINEAR)
					{
						Myfloat ix0 = x_origin_mesh + ix * el_sizex ;
						Myfloat iz0 = z_origin_mesh + iz * el_sizez ;
						Myfloat alpha_x = (pNode[ pElement[iel].inode[inode] ].xcoord - ix0) / el_sizex ;
						Myfloat alpha_z = (pNode[ pElement[iel].inode[inode] ].zcoord - iz0) / el_sizez ;
						Myfloat vx_down = v1 + alpha_x * (v2 - v1) ;
						Myfloat vx_up   = v3 + alpha_x * (v4 - v3) ;
						vp_glob[ pElement[iel].inode[inode] ] = vx_down + alpha_z * (vx_up - vx_down) ;
					}
				}
			}
		}
	}

	Myfloat minval = FLT_MAX ;
	Myfloat maxval = -FLT_MAX ;
	for (Myint inode = 0; inode < nnode; inode++)
	{
		minval = min(vp_glob[inode], minval) ;
		maxval = max(vp_glob[inode], maxval) ;
	}
	print_info(MASTER, " Min. Vp in mesh", minval) ;
	print_info(MASTER, " Max. Vp in mesh", maxval) ;
	if (minval < 0.0)
	{
		print_error(" Error in FEM_2D::initialize_mesh, min. vp < 0") ;
		return(RTN_CODE_KO) ;
	}
	if (maxval < 0.0)
	{
		print_error(" Error in FEM_2D::initialize_mesh, max. vp < 0") ;
		return(RTN_CODE_KO) ;
	}

	// assign VS per node
	// only for  ELASTIC
	//----------------------------------------------------------------------------------

	if (eq_type == ELASTIC)
	{
		print_debug(ALL, MID_DEBUG, "assign VS per node") ;

		// initialize pVec_vs_glob
		pVec_vs_glob = new Grid_1D_float(nnode, 0.) ;
		Myfloat *vs_glob = pVec_vs_glob->pArray ;

		iel = -1 ;
		for (Myint iz = 0; iz < nelem_z; iz++)
		{
			for (Myint ix = 0; ix < nelem_x; ix++)
			{
				iel++ ;

				// determine vs
				//if (pElement[iel].region == MEDIUM)
				{
					// retrieve velocity at 4 corners of element
					// and take minimum value

					Myint ix_model = (x_origin_mesh + ((ix + 0) * el_sizex)) / dx_model ;
					Myint iz_model = (z_origin_mesh + ((iz + 0) * el_sizez)) / dz_model ;
					ix_model = max(0,  ix_model) ;
					ix_model = min(nx_model-1, ix_model) ;
					iz_model = max(0,  iz_model) ;
					iz_model = min(nz_model-1, iz_model) ;
					Myfloat v1 = vs_2D_grid->pArray[ix_model][iz_model] ;

					ix_model = (x_origin_mesh + ((ix + 1) * el_sizex)) / dx_model ;
					iz_model = (z_origin_mesh + ((iz + 0) * el_sizez)) / dz_model ;
					ix_model = max(0,  ix_model) ;
					ix_model = min(nx_model-1, ix_model) ;
					iz_model = max(0,  iz_model) ;
					iz_model = min(nz_model-1, iz_model) ;
					Myfloat v2 = vs_2D_grid->pArray[ix_model][iz_model] ;

					ix_model = (x_origin_mesh + ((ix + 0) * el_sizex)) / dx_model ;
					iz_model = (z_origin_mesh + ((iz + 1) * el_sizez)) / dz_model ;
					ix_model = max(0,  ix_model) ;
					ix_model = min(nx_model-1, ix_model) ;
					iz_model = max(0,  iz_model) ;
					iz_model = min(nz_model-1, iz_model) ;
					Myfloat v3 = vs_2D_grid->pArray[ix_model][iz_model] ;

					ix_model = (x_origin_mesh + ((ix + 1) * el_sizex)) / dx_model ;
					iz_model = (z_origin_mesh + ((iz + 1) * el_sizez)) / dz_model ;
					ix_model = max(0,  ix_model) ;
					ix_model = min(nx_model-1, ix_model) ;
					iz_model = max(0,  iz_model) ;
					iz_model = min(nz_model-1, iz_model) ;
					Myfloat v4 = vs_2D_grid->pArray[ix_model][iz_model] ;

					// assign vs to each node
					for (Myint inode = 0; inode < pElement[iel].nnode ; inode++)
					{
						// PROP_CONST
						if (prop == PROP_CONST)
						{
							ix_model = (x_origin_mesh + ((ix + 0.5) * el_sizex)) / dx_model ;
							iz_model = (z_origin_mesh + ((iz + 0.5) * el_sizez)) / dz_model ;
							ix_model = max(0,  ix_model) ;
							ix_model = min(nx_model-1, ix_model) ;
							iz_model = max(0,  iz_model) ;
							iz_model = min(nz_model-1, iz_model) ;
							Myfloat v0 = vs_2D_grid->pArray[ix_model][iz_model] ;
							vs_glob[ pElement[iel].inode[inode] ] = v0 ;
						}
						// PROP_AVER
						else if (prop == PROP_AVER)
						{
							vs_glob[ pElement[iel].inode[inode] ] = (v1 + v2 + v3 + v4) / 4.0 ;
						}
						// linear properties per element
						// PROP_LINEAR
						else if (prop == PROP_LINEAR)
						{
							Myfloat ix0 = x_origin_mesh + ix * el_sizex ;
							Myfloat iz0 = z_origin_mesh + iz * el_sizez ;
							Myfloat alpha_x = (pNode[ pElement[iel].inode[inode] ].xcoord - ix0) / el_sizex ;
							Myfloat alpha_z = (pNode[ pElement[iel].inode[inode] ].zcoord - iz0) / el_sizez ;
							Myfloat vx_down = v1 + alpha_x * (v2 - v1) ;
							Myfloat vx_up   = v3 + alpha_x * (v4 - v3) ;
							vs_glob[ pElement[iel].inode[inode] ] = vx_down + alpha_z * (vx_up - vx_down) ;
						}
					}
				}
			}
		}

		Myfloat minval = FLT_MAX ;
		Myfloat maxval = -FLT_MAX ;
		for (Myint inode = 0; inode < nnode; inode++)
		{
			minval = min(vs_glob[inode], minval) ;
			maxval = max(vs_glob[inode], maxval) ;
		}
		print_info(MASTER, " Min. Vs in mesh", minval) ;
		print_info(MASTER, " Max. Vs in mesh", maxval) ;
		if (minval < 0.0)
		{
			print_error(" Error in FEM_2D::initialize_mesh, min. vs < 0") ;
			return(RTN_CODE_KO) ;
		}
		if (maxval < 0.0)
		{
			print_error(" Error in FEM_2D::initialize_mesh, max. vs < 0") ;
			return(RTN_CODE_KO) ;
		}
	}

	// assign RHO per node
	// only for  ELASTIC
	//----------------------------------------------------------------------------------

	if (eq_type == ELASTIC)
	{
		print_debug(ALL, MID_DEBUG, "assign RHO per node") ;

		// initialize pVec_rho_glob
		pVec_rho_glob = new Grid_1D_float(nnode, 0.) ;
		Myfloat *rho_glob = pVec_rho_glob->pArray ;

		iel = -1 ;
		for (Myint iz = 0; iz < nelem_z; iz++)
		{
			for (Myint ix = 0; ix < nelem_x; ix++)
			{
				iel++ ;

				// determine rho
				//if (pElement[iel].region == MEDIUM)
				{
					// retrieve velocity at 4 corners of element
					// and take minimum value

					Myint ix_model = (x_origin_mesh + ((ix + 0) * el_sizex)) / dx_model ;
					Myint iz_model = (z_origin_mesh + ((iz + 0) * el_sizez)) / dz_model ;
					ix_model = max(0,  ix_model) ;
					ix_model = min(nx_model-1, ix_model) ;
					iz_model = max(0,  iz_model) ;
					iz_model = min(nz_model-1, iz_model) ;
					Myfloat v1 = rho_2D_grid->pArray[ix_model][iz_model] ;

					ix_model = (x_origin_mesh + ((ix + 1) * el_sizex)) / dx_model ;
					iz_model = (z_origin_mesh + ((iz + 0) * el_sizez)) / dz_model ;
					ix_model = max(0,  ix_model) ;
					ix_model = min(nx_model-1, ix_model) ;
					iz_model = max(0,  iz_model) ;
					iz_model = min(nz_model-1, iz_model) ;
					Myfloat v2 = rho_2D_grid->pArray[ix_model][iz_model] ;

					ix_model = (x_origin_mesh + ((ix + 0) * el_sizex)) / dx_model ;
					iz_model = (z_origin_mesh + ((iz + 1) * el_sizez)) / dz_model ;
					ix_model = max(0,  ix_model) ;
					ix_model = min(nx_model-1, ix_model) ;
					iz_model = max(0,  iz_model) ;
					iz_model = min(nz_model-1, iz_model) ;
					Myfloat v3 = rho_2D_grid->pArray[ix_model][iz_model] ;

					ix_model = (x_origin_mesh + ((ix + 1) * el_sizex)) / dx_model ;
					iz_model = (z_origin_mesh + ((iz + 1) * el_sizez)) / dz_model ;
					ix_model = max(0,  ix_model) ;
					ix_model = min(nx_model-1, ix_model) ;
					iz_model = max(0,  iz_model) ;
					iz_model = min(nz_model-1, iz_model) ;
					Myfloat v4 = rho_2D_grid->pArray[ix_model][iz_model] ;

					// assign rho to each node
					for (Myint inode = 0; inode < pElement[iel].nnode ; inode++)
					{
						// PROP_CONST
						if (prop == PROP_CONST)
						{
							ix_model = (x_origin_mesh + ((ix + 0.5) * el_sizex)) / dx_model ;
							iz_model = (z_origin_mesh + ((iz + 0.5) * el_sizez)) / dz_model ;
							ix_model = max(0,  ix_model) ;
							ix_model = min(nx_model-1, ix_model) ;
							iz_model = max(0,  iz_model) ;
							iz_model = min(nz_model-1, iz_model) ;
							Myfloat v0 = rho_2D_grid->pArray[ix_model][iz_model] ;
							rho_glob[ pElement[iel].inode[inode] ] = v0 ;
						}
						// PROP_AVER
						else if (prop == PROP_AVER)
						{
							rho_glob[ pElement[iel].inode[inode] ] = (v1 + v2 + v3 + v4) / 4.0 ;
						}
						// linear properties per element
						// PROP_LINEAR
						else if (prop == PROP_LINEAR)
						{
							Myfloat ix0 = x_origin_mesh + ix * el_sizex ;
							Myfloat iz0 = z_origin_mesh + iz * el_sizez ;
							Myfloat alpha_x = (pNode[ pElement[iel].inode[inode] ].xcoord - ix0) / el_sizex ;
							Myfloat alpha_z = (pNode[ pElement[iel].inode[inode] ].zcoord - iz0) / el_sizez ;
							Myfloat vx_down = v1 + alpha_x * (v2 - v1) ;
							Myfloat vx_up   = v3 + alpha_x * (v4 - v3) ;
							rho_glob[ pElement[iel].inode[inode] ] = vx_down + alpha_z * (vx_up - vx_down) ;
						}
					}
				}
			}
		}

		Myfloat minval = FLT_MAX ;
		Myfloat maxval = -FLT_MAX ;
		for (Myint inode = 0; inode < nnode; inode++)
		{
			minval = min(rho_glob[inode], minval) ;
			maxval = max(rho_glob[inode], maxval) ;
		}
		print_info(MASTER, " Min. Rho in mesh", minval) ;
		print_info(MASTER, " Max. Rho in mesh", maxval) ;
		if (minval < 0.0)
		{
			print_error(" Error in FEM_2D::initialize_mesh, min. rho < 0") ;
			return(RTN_CODE_KO) ;
		}
		if (maxval < 0.0)
		{
			print_error(" Error in FEM_2D::initialize_mesh, max. rho < 0") ;
			return(RTN_CODE_KO) ;
		}
	}

	// initialize absorbing sponge
	//============================
	if (is_there_boundary_type(SPG))
	{
		print_debug(ALL, MID_DEBUG, "initialize absorbing sponge") ;

		// determine region coordinates
		Myfloat zGhostBeg  = z_origin_mesh + izGhostBeg  * el_sizez ;
		Myfloat zLayerBeg  = z_origin_mesh + izLayerBeg  * el_sizez ;
		Myfloat zMediumBeg = z_origin_mesh + izMediumBeg * el_sizez ;
		Myfloat zMediumEnd = z_origin_mesh + izMediumEnd * el_sizez ;
		Myfloat zLayerEnd  = z_origin_mesh + izLayerEnd  * el_sizez ;
		Myfloat zGhostEnd  = z_origin_mesh + izGhostEnd  * el_sizez ;

		Myfloat xGhostBeg  = x_origin_mesh + ixGhostBeg  * el_sizex ;
		Myfloat xLayerBeg  = x_origin_mesh + ixLayerBeg  * el_sizex ;
		Myfloat xMediumBeg = x_origin_mesh + ixMediumBeg * el_sizex ;
		Myfloat xMediumEnd = x_origin_mesh + ixMediumEnd * el_sizex ;
		Myfloat xLayerEnd  = x_origin_mesh + ixLayerEnd  * el_sizex ;
		Myfloat xGhostEnd  = x_origin_mesh + ixGhostEnd  * el_sizex ;

		// allocate array
		pSponge_coef_node = new Grid_1D_float(nnode, 0.) ;
		Myfloat* sponge_coef_node = pSponge_coef_node->pArray ;

		// compute sponge coef alpha
		// if coef is set to zero in xml, then constant coef is used in sponges
		Myfloat zLayerBeg_alpha = -1 ;
		if (get_boundary_coef(ZBEG) > 0) zLayerBeg_alpha = sqrt(-log(get_boundary_coef(ZBEG))) ;
		Myfloat zLayerEnd_alpha = -1 ;
		if (get_boundary_coef(ZEND) > 0) zLayerEnd_alpha = sqrt(-log(get_boundary_coef(ZEND))) ;

		Myfloat xLayerBeg_alpha = -1 ;
		if (get_boundary_coef(XBEG) > 0) xLayerBeg_alpha = sqrt(-log(get_boundary_coef(XBEG))) ;
		Myfloat xLayerEnd_alpha = -1 ;
		if (get_boundary_coef(XEND) > 0) xLayerEnd_alpha = sqrt(-log(get_boundary_coef(XEND))) ;

		// compute coef
		for (Myint inode=0; inode < nnode; inode++)
		{
			// set coef to 1
			sponge_coef_node[inode] = 1.0 ;

			// node in sponge
			if (pNode[inode].boundary == SPG)
			{
				// retrieve node coordinates
				Myfloat zcoord = pNode[inode].zcoord ;
				Myfloat xcoord = pNode[inode].xcoord ;

				// ZBEG
				if ((zcoord <= zMediumBeg) && (zcoord >= zLayerBeg))
				{
					if (zLayerBeg_alpha < 0)
					{
						sponge_coef_node[inode] *= SPONGE_COEF ;
					}
					else
					{
						Myfloat ratio = (zMediumBeg - zcoord) / (zMediumBeg - zLayerBeg) ;
						Myfloat gcoef = exp(-pow(zLayerBeg_alpha * ratio, 2.0)) ;
						//cout << zcoord << " " << ratio << " " << gcoef << "\n" ;
						sponge_coef_node[inode] *= gcoef ;
					}
				}

				// ZEND
				else if ((zcoord >= zMediumEnd) && (zcoord <= zLayerEnd))
				{
					if (zLayerEnd_alpha < 0)
					{
						sponge_coef_node[inode] *= SPONGE_COEF ;
					}
					else
					{
						Myfloat ratio = (zcoord - zMediumEnd) / (zLayerEnd - zMediumEnd) ;
						Myfloat gcoef = exp(-pow(zLayerEnd_alpha * ratio, 2.0)) ;
						//cout << zcoord << " " << ratio << " " << gcoef << "\n" ;
						sponge_coef_node[inode] *= gcoef ;
					}
				}

				// XBEG
				if ((xcoord <= xMediumBeg) && (xcoord >= xLayerBeg))
				{
					if (xLayerBeg_alpha < 0)
					{
						sponge_coef_node[inode] *= SPONGE_COEF ;
					}
					else
					{
						Myfloat ratio = (xMediumBeg - xcoord) / (xMediumBeg - xLayerBeg) ;
						Myfloat gcoef = exp(-pow(xLayerBeg_alpha * ratio, 2.0)) ;
						//cout << xcoord << " " << ratio << " " << gcoef << "\n" ;
						sponge_coef_node[inode] *= gcoef ;
					}
				}

				// XEND
				else if ((xcoord >= xMediumEnd) && (xcoord <= xLayerEnd))
				{
					if (xLayerEnd_alpha < 0)
					{
						sponge_coef_node[inode] *= SPONGE_COEF ;
					}
					else
					{
						Myfloat ratio = (xcoord - xMediumEnd) / (xLayerEnd - xMediumEnd) ;
						Myfloat gcoef = exp(-pow(xLayerEnd_alpha * ratio, 2.0)) ;
						//cout << xcoord << " " << ratio << " " << gcoef << "\n" ;
						sponge_coef_node[inode] *= gcoef ;
					}
				}

			}
		}
	}

	// save nodes coordinates on disk
	//-------------------------------

	print_debug(ALL, MID_DEBUG, "save nodes coordinates on disk") ;

	ofstream out_file ;
	out_file.open(NODE_COORD_OUT_FILE, ios::trunc | ios::out) ;
	for (Myint inode = 0; inode < nnode; inode++)
	{
		out_file << pNode[inode].zcoord << " " << pNode[inode].xcoord << "\n" ;
	}
	out_file.close() ;

	// display mesh info
	Rtn_code rtn_code = mesh_info() ;
	if (rtn_code != RTN_CODE_OK) return (rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D::initialize_mesh");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FEM_2D::initialize(Model* pModel, Myfloat fmax)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D::initialize");

	//======================================================
	// initialize node coordinates in ref element
	//======================================================
	Rtn_code rtn_code = initialize_node_ref() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	//======================================================
	// initialize mesh
	//======================================================
	rtn_code = initialize_mesh(pModel, fmax) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	//======================================================
	// compute elements matrices
	//======================================================
	rtn_code = compute_element_matrices() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	//======================================================
	// compute normal vectors
	//======================================================

	for (Myint ielem = 0; ielem < nelem; ielem++)
	{
		// $ tmp
		Myfloat hx = pElement[ielem].size_x ;
		Myfloat hz = pElement[ielem].size_z ;

		pElement[ielem].vnz[I_ZPREV] = -hx ;
		pElement[ielem].vnx[I_ZPREV] = 0.0 ;

		pElement[ielem].vnz[I_ZNEXT] = +hx ;
		pElement[ielem].vnx[I_ZNEXT] = 0.0 ;

		pElement[ielem].vnz[I_XPREV] = 0.0 ;
		pElement[ielem].vnx[I_XPREV] = -hz ;

		pElement[ielem].vnz[I_XNEXT] = 0.0 ;
		pElement[ielem].vnx[I_XNEXT] = +hz ;
	}

	// compute optimal time step
	Myfloat optimal_dt = compute_optimal_time_step() ;

	// set appropriate nt and dt
	if (set_nt_and_dt(optimal_dt) == RTN_CODE_KO) return(RTN_CODE_KO) ;

	// call parent initialization
	rtn_code = FEM::initialize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// initialize global vector pVec_k_glob
	//===============================================================
	pVec_k_glob = new Grid_1D_float(nnode, 0.) ;

	// initialize global invert mass matrix
	//===============================================================
	rtn_code = compute_global_inv_mass_matrix() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// save mesh in VTK format
	//write_mesh_VTK() ;
	//write_node_VTK() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D::initialize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FEM_2D::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D::finalize");

	// deallocate global node array
	deallocate_array<Node_2D_struct_type>(pNode, nnode) ;

	// deallocate local node array
	for (Myint ielem = 0; ielem < nelem; ielem++)
	{
		deallocate_array<Myint>(pElement[ielem].inode, pElement[ielem].nnode) ;
	}

	// deallocate element array
	deallocate_array<Quad_struct_type>(pElement, nelem) ;

	// call parent finalize
	Rtn_code rtn_code = FEM::finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// deallocate stifness matrices
	for (Myint imat = 0; imat <= MAX_POLY_ORDER; imat++)
	{
		// full matrix
		if (pMat_Dz[imat] != NULL) delete(pMat_Dz[imat]) ;
		if (pMat_Dx[imat] != NULL) delete(pMat_Dx[imat]) ;

		// sparse matrix
		if (pMat_Dz_i[imat] != NULL) delete(pMat_Dz_i[imat]) ;
		if (pMat_Dz_j[imat] != NULL) delete(pMat_Dz_j[imat]) ;
		if (pMat_Dz_c[imat] != NULL) delete(pMat_Dz_c[imat]) ;
		if (pMat_Dx_i[imat] != NULL) delete(pMat_Dx_i[imat]) ;
		if (pMat_Dx_j[imat] != NULL) delete(pMat_Dx_j[imat]) ;
		if (pMat_Dx_c[imat] != NULL) delete(pMat_Dx_c[imat]) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D::finalize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FEM_2D::compute_eigen_error(Myfloat* pr, Myint it, Acquisition* pAcquisition) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D::compute_eigen_error");

	// compute error only at output time steps
	Myint decim = round(dt_out / dt) ;
	if (it%decim != 0) return(RTN_CODE_OK) ;

	Myfloat ttime = dt * it ;
	Myfloat sq2 = sqrt(2.0) ;

	// compute error at the nodes positions
	//-------------------------------------
	{
		// loop on elements
		for (Myint iel = 0; iel < nelem; iel++)
		{
			// loop on nodes
			for (Myint ii=0; ii< pElement[iel].nnode ; ii++)
			{
				Myint inode = pElement[iel].inode[ii] ;
				Myfloat znode = pNode[inode].zcoord ;
				Myfloat xnode = pNode[inode].xcoord ;

				if (abs(pr[inode]) > MAX_EIGEN_VAL)
				{
					print_error(" MAX_EIGEN_VAL HAS BEEN REACHED") ;
					return(RTN_CODE_KO) ;
				}

				Myfloat eigen_sol = -sq2 * sin(M_PI*xnode * eigen_nmode)
				* sin(M_PI*znode * eigen_nmode) * sin(sq2*M_PI*ttime * eigen_nmode) ;
				eigen_error_node += pow(eigen_sol - pr[inode], 2) ;
			}
		}

		// open, write and close
		ofstream pFile ;
		pFile.open(EIGEN_ERROR_NODE_OUT_FILE, ios::app | ios::out) ;
		assert(pFile.is_open());
		pFile << ttime << " " << eigen_error_node << "\n" ;
		pFile.close() ;
	}

	// compute error at the rec positions
	//-----------------------------------
	{
		// loop on receivers
		for (Myint irec = 0; irec < nrec; irec++)
		{
			if (ielem_rec[irec] == NOT_FOUND) continue ;

			Myfloat znode = pAcquisition->zrec[irec] ;
			Myfloat xnode = pAcquisition->xrec[irec] ;
			Myfloat eigen_sol = -sq2 * sin(M_PI*xnode * eigen_nmode)
			* sin(M_PI*znode * eigen_nmode) * sin(sq2*M_PI*ttime * eigen_nmode) ;

			Myfloat pr_rec = interpolate_variable(pr, ielem_rec[irec], eta_rec[irec], xi_rec[irec]) ;

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
				<< nnode << " "
				<< eigen_error_rec_l1 << " "
				<< sqrt(eigen_error_rec_l2) << " "
				<< sqrt(eigen_error_rec_l2 / eigen_error_rec_sum2) << "\n" ;
		pFile.close() ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D::compute_eigen_error");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat FEM_2D::compute_optimal_time_step(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D::compute_optimal_time_step");

	Myfloat optim_dt = +9999 ;

	// global variables
	min_dist_node   = FLT_MAX ;
	max_dist_node   = FLT_MIN ;

	// loop over all elements
	for (Myint ielem = 0; ielem < nelem; ielem++)
	{

		// retrieve polynomial order
		Myint iorder = pElement[ielem].nnode_1D - 1 ;

		// retrieve minimal distance between nodes
		Myfloat min_dist_node_i = FLT_MAX ;
		Myfloat max_dist_node_i = FLT_MIN ;

		// retrieve min and max distance between node
		min_dist_node = FLT_MAX ;
		max_dist_node = FLT_MIN ;

		// P0 special case
		if (pElement[ielem].nnode_1D == 1)
		{
			Myfloat dist = min(pElement[ielem].size_x, pElement[ielem].size_z) ;
			max_dist_node_i = max(max_dist_node_i, dist) ;
			min_dist_node_i = min(min_dist_node_i, dist) ;
		}
		else
		{
			for (Myint ii=0; ii< pElement[ielem].nnode-1 ; ii++)
			{
				Myint inode1 = pElement[ielem].inode[ii] ;
				for (Myint jj=ii+1; jj< pElement[ielem].nnode ; jj++)
				{
					Myint inode2 = pElement[ielem].inode[jj] ;

					// distance between node1 and node2
					Myfloat dist = sqrt( pow(pNode[inode1].zcoord - pNode[inode2].zcoord, 2) +
							pow(pNode[inode1].xcoord - pNode[inode2].xcoord, 2) ) ;

					max_dist_node_i = max(max_dist_node_i, dist) ;
					min_dist_node_i = min(min_dist_node_i, dist) ;
				}
			}
		}

		// update global variables
		max_dist_node = max(max_dist_node_i, max_dist_node) ;
		min_dist_node = min(min_dist_node_i, min_dist_node) ;

		// Hesthaven formula
		//===================
		Myfloat dt_i = min_dist_node_i / pElement[ielem].vmax / sqrt(2.0) ;
		if (type == SCHEME_CGM)
		{
			dt_i *= CFL_CGSE[iorder] ;
		}
		else
		{
			dt_i *= CFL_DGSE[iorder] ;
		}

		optim_dt = min(optim_dt, dt_i) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "CFL", CFL);
	print_info(MASTER, " Optimal dt (s):", optim_dt) ;
	if (ratio_cfl != 1.0)
	{
		optim_dt  *= ratio_cfl ;
		print_info(MASTER, " Ratio CFL:\t", ratio_cfl) ;
		print_info(MASTER, " Modified dt (s):", optim_dt) ;
	}

	print_info(MASTER, " min dist. bw. node (m)", min_dist_node) ;
	print_info(MASTER, " max dist. bw. node (m)", max_dist_node) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D::compute_optimal_time_step");
	return(optim_dt) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FEM_2D::locate_pixel_in_mesh(Snapshot *pSnapshot)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D::locate_pixel_in_mesh");

	if (output_snapshot)
	{
		double t0 = MPI_Wtime() ;

		// retrieve number of pixel
		Myint npixel_z = 0 ;
		if (pSnapshot->z_coord != NULL) npixel_z = pSnapshot->z_coord->nz ;
		if (npixel_z <= 0)
		{
			print_error(" Error in FEM_2D::locate_pixel_in_grid, npixel_z <= 0") ;
			return(RTN_CODE_KO) ;
		}
		Myint npixel_x = 0 ;
		if (pSnapshot->x_coord != NULL) npixel_x = pSnapshot->x_coord->nz ;
		if (npixel_x <= 0)
		{
			print_error(" Error in FEM_2D::locate_pixel_in_grid, npixel_x <= 0") ;
			return(RTN_CODE_KO) ;
		}
		npixel = npixel_z * npixel_x ;

		// --> exact location with interpolation
		//----------------------------------------

		print_debug(ALL, MID_DEBUG, "receiver located at exact position with interpolation") ;
		ielem_pixel = allocate_array<Myint>(npixel) ;
		eta_pixel   = allocate_array<Myfloat>(npixel) ;
		xi_pixel    = allocate_array<Myfloat>(npixel) ;

		// define epsilon
		//Myfloat eps_dist = min_dist_node * MY_EPSILON ;
		Myfloat eps_dist = min_dist_node ;

		// regular or irregular mesh
		// TODO irregular mesh not supported
		bool regular_mesh = true ;

		// loop on all pixel
		Myint ipixel = 0 ;
		Myint npixel_not_found = 0 ;

		if (regular_mesh)
		{
			// get min mesh coordinates
			Myint ielem = 0 ;
			const Myfloat zmin = min_z_mesh ;
			const Myfloat xmin = min_x_mesh ;
			const Myfloat zmax = max_z_mesh ;
			const Myfloat xmax = max_x_mesh ;
			const Myfloat zsize = pElement[ielem].size_z ;
			const Myfloat xsize = pElement[ielem].size_x ;

			for (Myint izp=0; izp<npixel_z; izp++)
			{
				for (Myint ixp=0; ixp<npixel_x; ixp++)
				{
					ielem_pixel[ipixel] = NOT_FOUND ;
					Myfloat x_pixel = pSnapshot->x_coord->pArray[ixp] ;
					Myfloat z_pixel = pSnapshot->z_coord->pArray[izp] ;

					// retrieve element that contain the pixel
					Myint ielz ;
					// to avoid issues when z_pixel = max_z_mesh
					if (abs(z_pixel - max_z_mesh) < eps_dist)
					{
						ielz = nelem_z - 1 ;
					}
					else
					{
						ielz = floor((z_pixel - zmin) / zsize) ;
					}

					Myint ielx ;
					if (abs(x_pixel - max_x_mesh) < eps_dist)
					{
						ielx = nelem_x - 1 ;
					}
					else
					{
						ielx = floor((x_pixel - xmin) / xsize) ;
					}

					// pixel is in mesh
					if ((ielz >= 0) && (ielz < nelem_z) && (ielx >= 0) && (ielx < nelem_x))
					{
						ielem = ielx + ielz * nelem_x ;
						ielem_pixel[ipixel] = ielem ;

						// get extremums of element
						if (pElement[ielem].nnode > 1)
						{
							Myfloat znode_min = pNode[pElement[ielem].inode[0]].zcoord ;
							Myfloat xnode_min = pNode[pElement[ielem].inode[0]].xcoord ;
							eta_pixel[ipixel]   = -1.0 + 2.0*(z_pixel - znode_min)/pElement[ielem].size_z  ;
							xi_pixel[ipixel]    = -1.0 + 2.0*(x_pixel - xnode_min)/pElement[ielem].size_x  ;
						}
						// P0 element is special case
						else
						{
							Myfloat znode_min = pNode[pElement[ielem].inode[0]].zcoord - pElement[ielem].size_z / 2.0 ;
							Myfloat xnode_min = pNode[pElement[ielem].inode[0]].xcoord - pElement[ielem].size_x / 2.0 ;
							eta_pixel[ipixel]   = -1.0 + 2.0*(z_pixel - znode_min)/pElement[ielem].size_z  ;
							xi_pixel[ipixel]    = -1.0 + 2.0*(x_pixel - xnode_min)/pElement[ielem].size_x  ;
						}
					}

					// check pixel is in mesh
					if (ielem_pixel[ipixel] == NOT_FOUND)
					{
						npixel_not_found++ ;
					}

					ipixel++ ;
				}
			}
		}
		else
		{

			for (Myint izp=0; izp<npixel_z; izp++)
			{
				for (Myint ixp=0; ixp<npixel_x; ixp++)
				{
					ielem_pixel[ipixel] = NOT_FOUND ;
					Myfloat x_pixel = pSnapshot->x_coord->pArray[ixp] ;
					Myfloat z_pixel = pSnapshot->z_coord->pArray[izp] ;

					// loop on all elements
					Myfloat znode_min, znode_max ;
					Myfloat xnode_min, xnode_max ;
					for (Myint ielem = 0; ielem < nelem; ielem++)
					{

						// get extremums of element
						if (pElement[ielem].nnode > 1)
						{
							znode_min = pNode[pElement[ielem].inode[0]].zcoord ;
							znode_max = pNode[pElement[ielem].inode[ pElement[ielem].nnode-1]].zcoord ;
							xnode_min = pNode[pElement[ielem].inode[0]].xcoord ;
							xnode_max = pNode[pElement[ielem].inode[ pElement[ielem].nnode-1]].xcoord ;
						}
						// P0 element is special case
						else
						{
							znode_min = pNode[pElement[ielem].inode[0]].zcoord - pElement[ielem].size_z / 2.0 ;
							znode_max = pNode[pElement[ielem].inode[0]].zcoord + pElement[ielem].size_z / 2.0 ;
							xnode_min = pNode[pElement[ielem].inode[0]].xcoord - pElement[ielem].size_x / 2.0 ;
							xnode_max = pNode[pElement[ielem].inode[0]].xcoord + pElement[ielem].size_x / 2.0 ;
						}

						// check rec is within the extremums
						if ( (z_pixel >= (znode_min-eps_dist)) && (z_pixel <= (znode_max+eps_dist))
								&& (x_pixel >= (xnode_min-eps_dist)) && (x_pixel <= (xnode_max+eps_dist)) )
						{
							ielem_pixel[ipixel] = ielem ;
							eta_pixel[ipixel]   = -1.0 + 2.0*(z_pixel - znode_min)/pElement[ielem].size_z  ;
							xi_pixel[ipixel]    = -1.0 + 2.0*(x_pixel - xnode_min)/pElement[ielem].size_x  ;
							break ;
						}
					}

					// check pixel is in mesh
					if (ielem_pixel[ipixel] == NOT_FOUND) npixel_not_found++ ;

					ipixel++ ;
				}
			}

		} // if (regular_mesh)

		if (npixel_not_found > 0)
		{
			print_warning(" number of pixel not found", npixel_not_found) ;
		}

		//--------------------------------------
		// write mesh parameters in binary files
		//--------------------------------------

		{
			Grid_1D_float val_grid = Grid_1D_float(npixel, 0.0) ;
			Myfloat* val = val_grid.pArray ;

			// write elem idx
			{
				for (Myint ipixel = 0; ipixel < npixel ; ipixel++)
				{
					val[ipixel] = ielem_pixel[ipixel] ;
				}
				ofstream val_file(MESH_ELEM_OUT_FILE, ios::binary | ios::app | ios::out) ;
				val_file.write((char*)&(val[0]), npixel * sizeof(Myfloat)) ;
				val_file.close();
			}

			// write vp
			{
				if (pVec_vp_glob == NULL)
				{
					print_error(" Error in FEM_2D::write_node_VTK, pVec_vp_glob == NULL") ;
					return(RTN_CODE_KO) ;
				}
				else if (pVec_vp_glob->pArray == NULL)
				{
					print_error(" Error in FEM_2D::write_node_VTK, pVec_vp_glob->pArray == NULL") ;
					return(RTN_CODE_KO) ;
				}

				for (Myint ipixel = 0; ipixel < npixel ; ipixel++)
				{
					if (ielem_pixel[ipixel] == NOT_FOUND)
					{
						val[ipixel] = NOT_FOUND ;
					}
					else
					{
						val[ipixel] = interpolate_variable(pVec_vp_glob->pArray,
								ielem_pixel[ipixel],
								eta_pixel[ipixel],
								xi_pixel[ipixel]) ;
					}
				}
				ofstream val_file(MESH_VP_OUT_FILE, ios::binary | ios::app | ios::out) ;
				val_file.write((char*)&(val[0]), npixel * sizeof(Myfloat)) ;
				val_file.close();
			}

			// write vs
			if (eq_type == ELASTIC)
			{
				if (pVec_vs_glob == NULL)
				{
					print_error(" Error in FEM_2D::write_node_VTK, pVec_vs_glob == NULL") ;
					return(RTN_CODE_KO) ;
				}
				else if (pVec_vs_glob->pArray == NULL)
				{
					print_error(" Error in FEM_2D::write_node_VTK, pVec_vs_glob->pArray == NULL") ;
					return(RTN_CODE_KO) ;
				}

				for (Myint ipixel = 0; ipixel < npixel ; ipixel++)
				{
					if (ielem_pixel[ipixel] == NOT_FOUND)
					{
						val[ipixel] = NOT_FOUND ;
					}
					else
					{
						val[ipixel] = interpolate_variable(pVec_vs_glob->pArray,
								ielem_pixel[ipixel],
								eta_pixel[ipixel],
								xi_pixel[ipixel]) ;
					}
				}
				ofstream val_file(MESH_VS_OUT_FILE, ios::binary | ios::app | ios::out) ;
				val_file.write((char*)&(val[0]), npixel * sizeof(Myfloat)) ;
				val_file.close();
			}

			// write rho
			if (eq_type == ELASTIC)
			{
				if (pVec_rho_glob == NULL)
				{
					print_error(" Error in FEM_2D::write_node_VTK, pVec_rho_glob == NULL") ;
					return(RTN_CODE_KO) ;
				}
				else if (pVec_rho_glob->pArray == NULL)
				{
					print_error(" Error in FEM_2D::write_node_VTK, pVec_rho_glob->pArray == NULL") ;
					return(RTN_CODE_KO) ;
				}

				for (Myint ipixel = 0; ipixel < npixel ; ipixel++)
				{
					if (ielem_pixel[ipixel] == NOT_FOUND)
					{
						val[ipixel] = NOT_FOUND ;
					}
					else
					{
						val[ipixel] = interpolate_variable(pVec_rho_glob->pArray,
								ielem_pixel[ipixel],
								eta_pixel[ipixel],
								xi_pixel[ipixel]) ;
					}
				}
				ofstream val_file(MESH_RHO_OUT_FILE, ios::binary | ios::app | ios::out) ;
				val_file.write((char*)&(val[0]), npixel * sizeof(Myfloat)) ;
				val_file.close();
			}

			// write order
			{
				for (Myint ipixel = 0; ipixel < npixel ; ipixel++)
				{
					if (ielem_pixel[ipixel] == NOT_FOUND)
					{
						val[ipixel] = NOT_FOUND ;
					}
					else
					{
						val[ipixel] = pElement[ ielem_pixel[ipixel] ].nnode_1D - 1 ;
					}
				}
				ofstream val_file(MESH_ORDER_OUT_FILE, ios::binary | ios::app | ios::out) ;
				val_file.write((char*)&(val[0]), npixel * sizeof(Myfloat)) ;
				val_file.close();
			}

			// write tmin
			{
				for (Myint ipixel = 0; ipixel < npixel ; ipixel++)
				{
					if (ielem_pixel[ipixel] == NOT_FOUND)
					{
						val[ipixel] = NOT_FOUND ;
					}
					else
					{
						val[ipixel] = pElement[ ielem_pixel[ipixel] ].tmin ;
					}
				}
				ofstream val_file(MESH_TMIN_OUT_FILE, ios::binary | ios::app | ios::out) ;
				val_file.write((char*)&(val[0]), npixel * sizeof(Myfloat)) ;
				val_file.close();
			}

			// write nelem per lambda
			if (node_distrib == NODE_DISTRIB_MODEL)
			{
				for (Myint ipixel = 0; ipixel < npixel ; ipixel++)
				{
					if (ielem_pixel[ipixel] == NOT_FOUND)
					{
						val[ipixel] = NOT_FOUND ;
					}
					else
					{
						if (fmax != 0.0)
						{
							Myfloat min_size = min(pElement[ ielem_pixel[ipixel] ].size_x,
									pElement[ ielem_pixel[ipixel] ].size_z) ;
							val[ipixel] = (pElement[ ielem_pixel[ipixel] ].vmin / fmax) / min_size ;
						}
						else
						{
							val[ipixel] = NOT_FOUND ;
						}
					}
				}
				ofstream val_file(MESH_NELPERLAMBDA_OUT_FILE, ios::binary | ios::app | ios::out) ;
				val_file.write((char*)&(val[0]), npixel * sizeof(Myfloat)) ;
				val_file.close();
			}
		}

		double time_in_function = MPI_Wtime() - t0 ;
		print_info(MASTER, " Time locate pixel (s)", (Myfloat) time_in_function) ;

	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D::locate_pixel_in_mesh");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FEM_2D::compute_tmin(Acquisition *acquisition)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D::compute_tmin");

	if (front_type == FRONT_STATIC)
	{
		for (Myint ielem=0; ielem<nelem; ielem++)
		{
			pElement[ielem].tmin = 0.0 ;
		}
	}

	else if (front_type == FRONT_DYN_VMAX)
	{
		double t0 = MPI_Wtime() ;

		// first loop on elements
		// compute min distance from source for each element
		// retrieve vmax in mesh
		Myfloat max_vel = -FLT_MAX ;
		for (Myint ielem=0; ielem<nelem; ielem++)
		{
			Myfloat min_dist = FLT_MAX ;
			// loop on nodes
			for (Myint ii=0; ii< pElement[ielem].nnode; ii++)
			{
				// compute distance src - node
				Myfloat dist =  sqrt( pow(acquisition->zsrc - pNode[pElement[ielem].inode[ii]].zcoord, 2) +
						pow(acquisition->xsrc - pNode[pElement[ielem].inode[ii]].xcoord, 2) ) ;
				min_dist = min (min_dist, dist) ;
			}
			pElement[ielem].tmin = min_dist ;

			// retrieve max velocity
			max_vel = max(max_vel, pElement[ielem].vmax) ;
		}

		if (max_vel <= 0.0)
		{
			print_error("IN FEM_2D::FEM_2D::compute_tmin --> invalid vmax", max_vel);
			return(RTN_CODE_KO) ;
		}
		else
		{
			print_info(MASTER, " Dynamic front vmax", max_vel) ;
		}

		// second loop on elements
		// compute tmin = min distance / vmax

		// get size of element for additional time delay
		// needed to allow correct source excitation
		Myfloat size_elem = pElement[0].size_x ;
		if (pElement[0].size_z > pElement[0].size_x)
		{
			size_elem = pElement[0].size_z ;
		}
		Myfloat delay_dist = size_elem * 2 ;

		// loop on elements
		for (Myint ielem=0; ielem<nelem; ielem++)
		{
			//pElement[ielem].tmin /= max_vel ;
			// tmin - delay
			pElement[ielem].tmin = (pElement[ielem].tmin - delay_dist) / max_vel ;
		}

		double time_in_function = MPI_Wtime() - t0 ;
		print_info(MASTER, " Time compute tmin (s)", (Myfloat) time_in_function) ;
	}

	else
	{
		print_error("IN FEM_2D::FEM_2D::compute_tmin --> invalid front type", front_type);
		return(RTN_CODE_KO) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D::compute_tmin");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FEM_2D::locate_src_and_rec_in_mesh(Acquisition *acquisition)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D::locate_src_and_rec_in_mesh");

	double t0 = MPI_Wtime() ;

	// retrieve nodes for the source excitation
	//#########################################

	nnode_src = 0 ;
	Myfloat sum_src = 0.0 ;
	Myfloat* const Mat_M_inv_glob = pMat_M_inv_glob->pArray ;

	// NO source smoothing
	// excitation on single node
	//==============================

	if (src_type == SRC_POINT)
	{
		inode_src = allocate_array<Myint>(1) ;
		wnode_src = allocate_array<Myfloat>(1) ;

		// retrieve position of source
		// --> nearest node
		Myfloat min_dist = FLT_MAX ;
		for (Myint ielem = 0; ielem < nelem; ielem++)
		{
			for (Myint ii=0; ii< pElement[ielem].nnode; ii++)
			{
				Myfloat dist =  sqrt( pow(acquisition->zsrc - pNode[pElement[ielem].inode[ii]].zcoord, 2) +
						pow(acquisition->xsrc - pNode[pElement[ielem].inode[ii]].xcoord, 2) ) ;
				if (dist <= min_dist)
				{
					min_dist = dist ;
					inode_src[0] = pElement[ielem].inode[ii] ;
					wnode_src[0] = Mat_M_inv_glob[inode_src[0]] ;
					sum_src = 1.0 ;
				}
			}
		}

		// $$$ tmp, need to check source is in domain
		// if (min_dist > max_elem_size)
		// 	{
		// 	  print_warning(" Source is not in domain, min_dist", min_dist) ;
		// 	  inode_src[0] = NOT_FOUND ;
		// 	}
		// else
		{
			nnode_src = 1 ;
		}
	}

	// WITH source gaussian
	// excitation on several node
	//================================
	else if (src_type == SRC_GAUSSIAN)

	{
		const Myfloat radius = src_sigma*3 ;

		// 1st loop to get number of node
		for (Myint inode = 0; inode < nnode; inode++)
		{
			Myfloat dist =  sqrt( pow(acquisition->zsrc - pNode[inode].zcoord, 2) +
					pow(acquisition->xsrc - pNode[inode].xcoord, 2) ) ;
			if (dist <= radius)
			{
				nnode_src++ ;
			}
		}

		if (nnode_src == 0)
		{
			print_warning(" Source is not in domain") ;
		}
		else
		{
			inode_src = allocate_array<Myint>(nnode_src) ;
			wnode_src = allocate_array<Myfloat>(nnode_src) ;
			nnode_src = 0 ;

			// 2nd loop to store the nodes and weights
			for (Myint inode = 0; inode < nnode; inode++)
			{
				Myfloat dist =  sqrt( pow(acquisition->zsrc - pNode[inode].zcoord, 2) +
						pow(acquisition->xsrc - pNode[inode].xcoord, 2) ) ;
				if (dist <= radius)
				{
					inode_src[nnode_src] = inode ;
					wnode_src[nnode_src] = exp(-dist*dist/(2.0*src_sigma*src_sigma)) ;
					sum_src += wnode_src[nnode_src] / Mat_M_inv_glob[inode_src[nnode_src]] ;
					nnode_src++ ;
				}
			}
		}
	}

	// WITH source SINC
	// excitation on several node
	//================================
	else if (src_type == SRC_SINC)

	{
		const Myint nw = 3 ;
		const Myfloat radius = src_sigma*(1+nw) ;

		// 1st loop to get number of node
		for (Myint inode = 0; inode < nnode; inode++)
		{
			Myfloat dist =  sqrt( pow(acquisition->zsrc - pNode[inode].zcoord, 2) +
					pow(acquisition->xsrc - pNode[inode].xcoord, 2) ) ;
			if (dist <= radius)
			{
				nnode_src++ ;
			}
		}

		if (nnode_src == 0)
		{
			print_warning(" Source is not in domain") ;
		}
		else
		{
			inode_src = allocate_array<Myint>(nnode_src) ;
			wnode_src = allocate_array<Myfloat>(nnode_src) ;
			nnode_src = 0 ;

			// 2nd loop to store the nodes
			for (Myint inode = 0; inode < nnode; inode++)
			{
				Myfloat dist =  sqrt( pow(acquisition->zsrc - pNode[inode].zcoord, 2) +
						pow(acquisition->xsrc - pNode[inode].xcoord, 2) ) ;
				if (dist <= radius)
				{
					inode_src[nnode_src] = inode ;
					Myfloat zeta = dist / src_sigma ;
					if (zeta != 0.0)
					{
						wnode_src[nnode_src] = 0.5 * (1 + cos((PI*zeta)/(nw+1))) * sin(PI*zeta)/(PI*zeta) ;
					}
					else
					{
						wnode_src[nnode_src] = 1.0 ;
					}
					sum_src += wnode_src[nnode_src] / Mat_M_inv_glob[inode_src[nnode_src]] ;
					nnode_src++ ;
				}
			}
		}
	}

	// normalisation
	print_info(MASTER, " # nodes for src ", nnode_src) ;
	for (Myint ii=0; ii<nnode_src; ii++)
	{
		wnode_src[ii] /= sum_src ;
	}

	// retrieve positions of receivers
	//################################

	nrec = acquisition->nrec ;

	// --> exact location with interpolation
	//----------------------------------------

	print_debug(ALL, MID_DEBUG, "receiver located at exact position with interpolation") ;
	ielem_rec = allocate_array<Myint>(nrec) ;
	eta_rec   = allocate_array<Myfloat>(nrec) ;
	xi_rec    = allocate_array<Myfloat>(nrec) ;

	// define epsilon
	Myfloat eps_dist = min_dist_node / 1000.0 ;

	// loop on all receivers
	Myfloat znode_min, znode_max ;
	Myfloat xnode_min, xnode_max ;
	for (Myint irec=0; irec<nrec; irec ++)
	{
		ielem_rec[irec] = NOT_FOUND ;

		// loop on all elements
		for (Myint ielem = 0; ielem < nelem; ielem++)
		{

			// get extremums of element
			if (pElement[ielem].nnode > 1)
			{
				znode_min = pNode[pElement[ielem].inode[0]].zcoord ;
				znode_max = pNode[pElement[ielem].inode[ pElement[ielem].nnode-1]].zcoord ;
				xnode_min = pNode[pElement[ielem].inode[0]].xcoord ;
				xnode_max = pNode[pElement[ielem].inode[ pElement[ielem].nnode-1]].xcoord ;
			}
			// P0 element is special case
			else
			{
				znode_min = pNode[pElement[ielem].inode[0]].zcoord - pElement[ielem].size_z / 2.0 ;
				znode_max = pNode[pElement[ielem].inode[0]].zcoord + pElement[ielem].size_z / 2.0 ;
				xnode_min = pNode[pElement[ielem].inode[0]].xcoord - pElement[ielem].size_x / 2.0 ;
				xnode_max = pNode[pElement[ielem].inode[0]].xcoord + pElement[ielem].size_x / 2.0 ;
			}

			// check rec is within the extremums
			if ( (acquisition->zrec[irec] >= (znode_min-eps_dist)) && (acquisition->zrec[irec] <= (znode_max+eps_dist))
					&& (acquisition->xrec[irec] >= (xnode_min-eps_dist)) && (acquisition->xrec[irec] <= (xnode_max+eps_dist)) )
			{
				ielem_rec[irec] = ielem ;
				eta_rec[irec] = -1.0 + 2.0*(acquisition->zrec[irec] - znode_min)/pElement[ielem].size_z  ;
				xi_rec[irec] = -1.0 + 2.0*(acquisition->xrec[irec] - xnode_min)/pElement[ielem].size_x  ;
				break ;
			}

		}

		// check receiver is in mesh
		if (ielem_rec[irec] == NOT_FOUND)
		{
			print_warning(" receiver not in domain with idx=", irec) ;
		}
		else
		{
			print_debug(ALL, MID_DEBUG, "receiver found at element idx=", ielem_rec[irec]) ;
			print_debug(ALL, MID_DEBUG, "receiver             with eta=", eta_rec[irec]) ;
			print_debug(ALL, MID_DEBUG, "receiver              with xi=", xi_rec[irec]) ;
		}
	}

	double time_in_function = MPI_Wtime() - t0 ;
	print_info(MASTER, " Time src & rec (s)", (Myfloat) time_in_function) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D::locate_src_and_rec_in_mesh");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FEM_2D::write_trace(Variable* pVar, Myint it)
{
	print_debug(ALL, FULL_DEBUG, "IN FEM_2D::write_trace");

	// check if ouput is required at the current time step
	Myint decim = round(dt_out / dt) ;
	string file_name = pVar->get_name() + TIME_REC_OUT_FILE ;

	if (it%decim != 0)
	{
		print_debug(ALL, FULL_DEBUG, "OUT FEM_2D::write_trace");
		return(RTN_CODE_OK) ;
	}
	else
	{
		print_debug(ALL, FULL_DEBUG, " write_trace in ", file_name) ;
		print_debug(ALL, FULL_DEBUG, " output trace (s) ", it*dt) ;
	}

	// get access to variable's grid
	Grid_1D_float *pGrid = (Grid_1D_float*) pVar->get_grid() ;
	Myfloat * const time_u = pGrid->pArray ;

	Myfloat u_tmp[nrec] ;

	for (Myint irec=0; irec<nrec; irec++)
	{
		if (ielem_rec[irec] == NOT_FOUND)
		{
			u_tmp[irec] = 0.0 ;
		}
		else
		{
			u_tmp[irec] = interpolate_variable(time_u, ielem_rec[irec], eta_rec[irec], xi_rec[irec]) ;
		}
	}

	// open, write and close
	ofstream pFile ;
	pFile.open(file_name.c_str(), ios::binary | ios::app | ios::out) ;
	assert(pFile.is_open());
	pFile.write((char*)u_tmp, nrec * sizeof(Myfloat)) ;
	pFile.close() ;

	print_debug(ALL, FULL_DEBUG, "OUT FEM_2D::write_trace");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat FEM_2D::interpolate_variable(Myfloat* pVal, Myint ielem, Myfloat eta_coord, Myfloat xi_coord)
{
	print_debug(ALL, FULL_DEBUG, "IN FEM_2D::interpolate_variable");

	Myfloat return_val = 0.0 ;

	// get number of nodes in element
	Myint ni = pElement[ielem].nnode_1D ;

	// P0 element
	if (ni == 1)
	{
		print_debug(ALL, FULL_DEBUG, "OUT FEM_2D::interpolate_variable");
		return(pVal[pElement[ielem].inode[0]]) ;
	}

	// build vector of nodes
	VecDoub *xnode = &(xnode_ref_elem_1D[ni-1]) ;

	// loop over the nodes of the element
	Myint in = 0 ;
	for (Myint ieta_i=0; ieta_i<ni; ieta_i++)
	{
		for (Myint ixi_i=0; ixi_i<ni; ixi_i++)
		{
			Myfloat poly_val = lagran(*xnode, ieta_i, eta_coord) * lagran(*xnode, ixi_i, xi_coord) ;
			return_val += poly_val * pVal[pElement[ielem].inode[in]] ;
			in++ ;
		}
	}

	print_debug(ALL, FULL_DEBUG, "OUT FEM_2D::interpolate_variable");
	return(return_val) ;
}

//-----------------------------------------------------------------------------------------
// compute inverse of global mass matrix with parameters
// from element mass matrix pMat_M_inv
//
// INPUT var    = type of physical pamateer
//       coef   = index of the global matrix
//       pModel = physical model
//
// OUTPUT pMat_M_inv_glob_param[coef]
//-----------------------------------------------------------------------------------------
Rtn_code FEM_2D::compute_global_inv_mass_matrix(Model* pModel, Var_type var, Coef_type coef)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D::compute_global_inv_mass_matrix");

	// check coef
	if ((coef < 0) || (coef >= MAX_INV_MATRIX))
	{
		print_error("IN FEM_2D::compute_global_inv_mass_matrix --> invalid coef", coef);
		return(RTN_CODE_KO) ;
	}

	// initialize pMat_M_inv_glob_param[coef]
	pMat_M_inv_glob_param[coef] = new Grid_1D_float(nnode, 0.) ;

	// retrieve needed physical paramaters
	//------------------------------------

	// get VP model
	Variable* vp_var = NULL ;
	Grid* vp_grid = NULL ;
	Grid_2D_float *vp_2D_grid = NULL ;
	if ((var == INV_RHOVP2) || (var == INV_RHOVP2MINUSVS2))
	{
		vp_var = pModel->get_parameter(VP) ;
		if (vp_var == NULL)
		{
			print_error("IN FEM_2D::compute_global_inv_mass_matrix --> VP model not found");
			return(RTN_CODE_KO) ;
		}
		vp_grid = vp_var->get_grid() ;
		if (vp_grid == NULL)
		{
			print_error("IN FEM_2D::compute_global_inv_mass_matrix --> VP grid not initialized");
			return(RTN_CODE_KO) ;
		}
		vp_2D_grid = dynamic_cast<Grid_2D_float*>(vp_grid) ;
		if (vp_2D_grid == NULL)
		{
			print_error("IN FEM_2D::compute_global_inv_mass_matrix --> VP grid is not Grid_2D_float");
			return(RTN_CODE_KO) ;
		}
	}

	// get VS model
	Variable* vs_var = NULL ;
	Grid* vs_grid = NULL ;
	Grid_2D_float *vs_2D_grid = NULL ;
	if ((var == INV_RHOVS2) || (var == INV_RHOVP2MINUSVS2))
	{
		vs_var = pModel->get_parameter(VS) ;
		if (vs_var == NULL)
		{
			print_error("IN FEM_2D::compute_global_inv_mass_matrix --> VS model not found");
			return(RTN_CODE_KO) ;
		}
		vs_grid = vs_var->get_grid() ;
		if (vs_grid == NULL)
		{
			print_error("IN FEM_2D::compute_global_inv_mass_matrix --> VS grid not initialized");
			return(RTN_CODE_KO) ;
		}
		vs_2D_grid = dynamic_cast<Grid_2D_float*>(vs_grid) ;
		if (vs_2D_grid == NULL)
		{
			print_error("IN FEM_2D::compute_global_inv_mass_matrix --> VS grid is not Grid_2D_float");
			return(RTN_CODE_KO) ;
		}
	}

	// get RHO model
	Variable* rho_var = NULL ;
	Grid* rho_grid = NULL ;
	Grid_2D_float *rho_2D_grid = NULL ;
	if ((var == RHO) || (var == INV_RHOVP2) || (var == INV_RHOVS2) || (var == INV_RHOVP2MINUSVS2))
	{
		rho_var = pModel->get_parameter(RHO) ;
		if (rho_var == NULL)
		{
			print_error("IN FEM_2D::compute_global_inv_mass_matrix --> RHO model not found");
			return(RTN_CODE_KO) ;
		}
		rho_grid = rho_var->get_grid() ;
		if (rho_grid == NULL)
		{
			print_error("IN FEM_2D::compute_global_inv_mass_matrix --> RHO grid not initialized");
			return(RTN_CODE_KO) ;
		}
		rho_2D_grid = dynamic_cast<Grid_2D_float*>(rho_grid) ;
		if (rho_2D_grid == NULL)
		{
			print_error("IN FEM_2D::compute_global_inv_mass_matrix --> RHO grid is not Grid_2D_float");
			return(RTN_CODE_KO) ;
		}
	}

	// reset Mat_M_inv_glob
	//===============================================================

	Myfloat* const Mat_M_inv_glob = pMat_M_inv_glob_param[coef]->pArray ;
	for (Myint inode = 0; inode < nnode; inode++)
	{
		Mat_M_inv_glob[inode] = 0.0 ;
	}

	// assemble mass matrix
	// only diagonal terms are considered -> valid only for GLL nodes
	//===============================================================

	Myfloat vp_node, vp1, vp2, vp3, vp4 ;
	Myfloat vs_node, vs1, vs2, vs3, vs4 ;
	Myfloat rho_node, rho1, rho2, rho3, rho4 ;

	Myint nz_model = pModel->get_nz() ;
	Myint nz ;
	if (nelem_z_med == 0)
	{
		// init from model
		nz = nz_model - 1 ; // nz grid points -> nz-1 elements
	}
	else
	{
		// init from XML
		nz = nelem_z_med ;
	}

	if (nz <= 0)
	{
		print_error("IN FEM_1D::compute_global_inv_mass_matrix --> nz <= 0", nz);
		return(RTN_CODE_KO) ;
	}

	Myint nx_model = pModel->get_nx() ;
	Myint nx ;
	if (nelem_x_med == 0)
	{
		// init from model
		nx = nx_model - 1 ; // nx grid points -> nz-1 elements
	}
	else
	{
		// init from XML
		nx = nelem_x_med ;
	}

	if (nx <= 0)
	{
		print_error("IN FEM_1D::compute_global_inv_mass_matrix --> nx <= 0", nx);
		return(RTN_CODE_KO) ;
	}

	Myfloat dz_model = pModel->get_dz() ;
	Myfloat dx_model = pModel->get_dx() ;

	Myint nlayer_zBeg = get_boundary_width(ZBEG) ;
	Myint nlayer_zEnd = get_boundary_width(ZEND) ;
	Myint nlayer_xBeg = get_boundary_width(XBEG) ;
	Myint nlayer_xEnd = get_boundary_width(XEND) ;

	Myint nghost_scheme ;
	if (type == SCHEME_CGM)
	{
		nghost_scheme = 0 ;
	}
	else if (type == SCHEME_DGM)
	{
		nghost_scheme = 0 ;
	}
	else if (type == SCHEME_MGM)
	{
		nghost_scheme = 0 ;
	}
	Myint nghost_zBeg = nghost_scheme ;
	Myint nghost_zEnd = nghost_scheme ;
	Myint nghost_xBeg = nghost_scheme ;
	Myint nghost_xEnd = nghost_scheme ;

	nelem_z = nz + nlayer_zBeg + nlayer_zEnd + nghost_zBeg + nghost_zEnd ;
	nelem_x = nx + nlayer_xBeg + nlayer_xEnd + nghost_xBeg + nghost_xEnd ;

	Myfloat el_sizez = (pModel->get_zcoord()->get_max() - pModel->get_zcoord()->get_min()) / nz ;
	Myfloat el_sizex = (pModel->get_xcoord()->get_max() - pModel->get_xcoord()->get_min()) / nx ;

	const Myfloat z_origin_mesh =  0.0 - el_sizez * (nlayer_zBeg + nghost_zBeg) ;
	const Myfloat x_origin_mesh =  0.0 - el_sizex * (nlayer_xBeg + nghost_xBeg) ;

	// loop on all elements
	Myint iel = -1 ;
	for (Myint iz = 0; iz < nelem_z; iz++)
	{
		for (Myint ix = 0; ix < nelem_x; ix++)
		{
			iel++ ;

			// get local invert matrix
			Myint const mat_idx = pElement[iel].nnode_1D - 1 ;
			Myfloat** const Mat_M_inv = pMat_M_inv[mat_idx]->pArray ;
			Myfloat cvol = Myfloat(2.0)*Myfloat(2.0) / (pElement[iel].size_x*pElement[iel].size_z) ;

			// retrieve physical  at 4 corners of element
			// and take minimum value

			Myint ix_model = (x_origin_mesh + ((ix + 0) * el_sizex)) / dx_model ;
			Myint iz_model = (z_origin_mesh + ((iz + 0) * el_sizez)) / dz_model ;
			ix_model = max(0,  ix_model) ;
			ix_model = min(nx_model-1, ix_model) ;
			iz_model = max(0,  iz_model) ;
			iz_model = min(nz_model-1, iz_model) ;
			if ((var == RHO) || (var == INV_RHOVP2) || (var == INV_RHOVS2) || (var == INV_RHOVP2MINUSVS2))
				rho1 = rho_2D_grid->pArray[ix_model][iz_model] ;
			if ((var == INV_RHOVP2) || (var == INV_RHOVP2MINUSVS2)) vp1 = vp_2D_grid->pArray[ix_model][iz_model] ;
			if ((var == INV_RHOVS2) || (var == INV_RHOVP2MINUSVS2)) vs1 = vs_2D_grid->pArray[ix_model][iz_model] ;

			ix_model = (x_origin_mesh + ((ix + 1) * el_sizex)) / dx_model ;
			iz_model = (z_origin_mesh + ((iz + 0) * el_sizez)) / dz_model ;
			ix_model = max(0,  ix_model) ;
			ix_model = min(nx_model-1, ix_model) ;
			iz_model = max(0,  iz_model) ;
			iz_model = min(nz_model-1, iz_model) ;
			if ((var == RHO) || (var == INV_RHOVP2) || (var == INV_RHOVS2) || (var == INV_RHOVP2MINUSVS2))
				rho2 = rho_2D_grid->pArray[ix_model][iz_model] ;
			if ((var == INV_RHOVP2) || (var == INV_RHOVP2MINUSVS2)) vp2 = vp_2D_grid->pArray[ix_model][iz_model] ;
			if ((var == INV_RHOVS2) || (var == INV_RHOVP2MINUSVS2)) vs2 = vs_2D_grid->pArray[ix_model][iz_model] ;

			ix_model = (x_origin_mesh + ((ix + 0) * el_sizex)) / dx_model ;
			iz_model = (z_origin_mesh + ((iz + 1) * el_sizez)) / dz_model ;
			ix_model = max(0,  ix_model) ;
			ix_model = min(nx_model-1, ix_model) ;
			iz_model = max(0,  iz_model) ;
			iz_model = min(nz_model-1, iz_model) ;
			if ((var == RHO) || (var == INV_RHOVP2) || (var == INV_RHOVS2) || (var == INV_RHOVP2MINUSVS2))
				rho3 = rho_2D_grid->pArray[ix_model][iz_model] ;
			if ((var == INV_RHOVP2) || (var == INV_RHOVP2MINUSVS2)) vp3 = vp_2D_grid->pArray[ix_model][iz_model] ;
			if ((var == INV_RHOVS2) || (var == INV_RHOVP2MINUSVS2)) vs3 = vs_2D_grid->pArray[ix_model][iz_model] ;

			ix_model = (x_origin_mesh + ((ix + 1) * el_sizex)) / dx_model ;
			iz_model = (z_origin_mesh + ((iz + 1) * el_sizez)) / dz_model ;
			ix_model = max(0,  ix_model) ;
			ix_model = min(nx_model-1, ix_model) ;
			iz_model = max(0,  iz_model) ;
			iz_model = min(nz_model-1, iz_model) ;
			if ((var == RHO) || (var == INV_RHOVP2) || (var == INV_RHOVS2) || (var == INV_RHOVP2MINUSVS2))
				rho4 = rho_2D_grid->pArray[ix_model][iz_model] ;
			if ((var == INV_RHOVP2) || (var == INV_RHOVP2MINUSVS2)) vp4 = vp_2D_grid->pArray[ix_model][iz_model] ;
			if ((var == INV_RHOVS2) || (var == INV_RHOVP2MINUSVS2)) vs4 = vs_2D_grid->pArray[ix_model][iz_model] ;

			// loop on element nodes
			for (Myint inode = 0; inode < pElement[iel].nnode ; inode++)
			{
				// PROP_CONST
				if (prop == PROP_CONST)
				{
					ix_model = (x_origin_mesh + ((ix + 0.5) * el_sizex)) / dx_model ;
					iz_model = (z_origin_mesh + ((iz + 0.5) * el_sizez)) / dz_model ;
					ix_model = max(0,  ix_model) ;
					ix_model = min(nx_model-1, ix_model) ;
					iz_model = max(0,  iz_model) ;
					iz_model = min(nz_model-1, iz_model) ;
					if ((var == RHO) || (var == INV_RHOVP2) || (var == INV_RHOVS2) || (var == INV_RHOVP2MINUSVS2))
					{
						rho_node = rho_2D_grid->pArray[ix_model][iz_model] ;
					}
					if ((var == INV_RHOVP2) || (var == INV_RHOVP2MINUSVS2))
					{
						vp_node = vp_2D_grid->pArray[ix_model][iz_model] ;
					}
					if ((var == INV_RHOVS2) || (var == INV_RHOVP2MINUSVS2))
					{
						vs_node = vs_2D_grid->pArray[ix_model][iz_model] ;
					}
				}

				// PROP_AVER
				else if (prop == PROP_AVER)
				{
					if ((var == RHO) || (var == INV_RHOVP2) || (var == INV_RHOVS2) || (var == INV_RHOVP2MINUSVS2))
					{
						rho_node = (rho1 + rho2 + rho3 + rho4) / 4.0 ;
					}
					if ((var == INV_RHOVP2) || (var == INV_RHOVP2MINUSVS2))
					{
						vp_node = (vp1 + vp2 + vp3 + vp4) / 4.0 ;
					}
					if ((var == INV_RHOVS2) || (var == INV_RHOVP2MINUSVS2))
					{
						vs_node = (vs1 + vs2 + vs3 + vs4) / 4.0 ;
					}
				}

				// linear properties per element
				// PROP_LINEAR
				else if (prop == PROP_LINEAR)
				{
					Myfloat ix0 = x_origin_mesh + ix * el_sizex ;
					Myfloat iz0 = z_origin_mesh + iz * el_sizez ;
					Myfloat alpha_x = (pNode[ pElement[iel].inode[inode] ].xcoord - ix0) / el_sizex ;
					Myfloat alpha_z = (pNode[ pElement[iel].inode[inode] ].zcoord - iz0) / el_sizez ;

					if ((var == RHO) || (var == INV_RHOVP2) || (var == INV_RHOVS2) || (var == INV_RHOVP2MINUSVS2))
					{
						Myfloat vx_down = rho1 + alpha_x * (rho2 - rho1) ;
						Myfloat vx_up   = rho3 + alpha_x * (rho4 - rho3) ;
						rho_node = vx_down + alpha_z * (vx_up - vx_down) ;
					}
					if ((var == INV_RHOVP2) || (var == INV_RHOVP2MINUSVS2))
					{
						Myfloat vx_down = vp1 + alpha_x * (vp2 - vp1) ;
						Myfloat vx_up   = vp3 + alpha_x * (vp4 - vp3) ;
						vp_node = vx_down + alpha_z * (vx_up - vx_down) ;
					}
					if ((var == INV_RHOVS2) || (var == INV_RHOVP2MINUSVS2))
					{
						Myfloat vx_down = vs1 + alpha_x * (vs2 - vs1) ;
						Myfloat vx_up   = vs3 + alpha_x * (vs4 - vs3) ;
						vs_node = vx_down + alpha_z * (vx_up - vx_down) ;
					}
				}

				// compute coef
				Myfloat coef_node ;

				// no coef (param = 1)
				if (var == NO_VAR_TYPE)
				{
					coef_node = 1.0 ;
				}

				// rho
				else if (var == RHO)
				{
					coef_node = rho_node ;
				}

				// 1/(rho * vp^2)
				else if (var == INV_RHOVP2)
				{
					coef_node = 1.0 / (rho_node * vp_node * vp_node) ;
				}

				// 1/(rho * vs^2)
				else if (var == INV_RHOVS2)
				{
					coef_node = 1.0 / (rho_node * vs_node * vs_node) ;
				}

				// 1/(rho * (vp^2 - vs^2))
				else if (var == INV_RHOVP2MINUSVS2)
				{
					coef_node = 1.0 / (rho_node * (vp_node*vp_node - vs_node*vs_node)) ;
				}
				// unknown parameter
				else
				{
					print_error("IN FEM_2D::compute_global_inv_mass_matrix --> invalid var", var);
					return(RTN_CODE_KO) ;
				}

				// assemble coef in matrix
				Myint kk = pElement[iel].inode[inode] ;
				Mat_M_inv_glob[kk] += coef_node / (Mat_M_inv[inode][inode] * cvol) ;

			} // for (Myint inode = 0; inode < pElement[iel].nnode ; inode++)
		}
	}

	// compute inverse of  Mat_M_inv_glob
	for (Myint inode = 0; inode < nnode; inode++)
	{
		Mat_M_inv_glob[inode] = 1.0 / Mat_M_inv_glob[inode] ;
		print_debug(ALL, MID_DEBUG, " Mat_M_inv_glob ", Mat_M_inv_glob[inode]) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D::compute_global_inv_mass_matrix");
	return(RTN_CODE_OK) ;
}

//-----------------------------------------------------------------------------------------
// compute inverse of global mass matrix
//
// INPUT pMat_M_inv
//
// OUTPUT pMat_M_inv_glob
//-----------------------------------------------------------------------------------------
Rtn_code FEM_2D::compute_global_inv_mass_matrix(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D::compute_global_inv_mass_matrix");

	// initialize pMat_M_inv_glob
	pMat_M_inv_glob = new Grid_1D_float(nnode, 0.) ;

	// assemble mass matrix
	// only diagonal terms are considered -> valid only for GLL nodes
	//===============================================================

	// reset Mat_M_inv_glob
	Myfloat* const Mat_M_inv_glob = pMat_M_inv_glob->pArray ;
	for (Myint inode = 0; inode < nnode; inode++)
	{
		Mat_M_inv_glob[inode] = 0.0 ;
	}

	// assemble diagonal terms
	for (Myint iel = 0; iel < nelem; iel++)
	{
		Myint const ni = pElement[iel].nnode ;
		Myint const mat_idx = pElement[iel].nnode_1D - 1 ;
		Myfloat** const Mat_M_inv = pMat_M_inv[mat_idx]->pArray ;
		Myfloat cvol = Myfloat(2.0)*Myfloat(2.0) / (pElement[iel].size_x*pElement[iel].size_z) ;
		for (Myint ii = 0; ii < ni; ii++) {
			Myint jj = ii ;
			Myint kk = pElement[iel].inode[ii] ;
			Mat_M_inv_glob[kk] += 1.0 / (Mat_M_inv[ii][jj] * cvol) ;
		}
	}

	// compute inverse of  Mat_M_inv_glob
	for (Myint inode = 0; inode < nnode; inode++)
	{
		Mat_M_inv_glob[inode] = 1.0 / Mat_M_inv_glob[inode] ;
		print_debug(ALL, MID_DEBUG, " Mat_M_inv_glob ", Mat_M_inv_glob[inode]) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D::compute_global_inv_mass_matrix");
	return(RTN_CODE_OK) ;
}

//-----------------------------------------------------------------------------------------
// compute mass matrix with quadrature rule such as Gauss-Legendre
//
// (Mi)jr = Int_elem_i_vol (phi_ir x phi_ij) dx
//
// INPUT
// xg (VecDoub) = array of quadrature coordinates
// wg (VecDoub) = array of quadrature weights
// xnode (VecDoub) = array of node coordinates in element
//
// OUTPUT
// Mat_M (MatDoub) mass matrix
//-----------------------------------------------------------------------------------------
Rtn_code FEM_2D::compute_mass_matrix(VecDoub &xg, VecDoub& wg, VecDoub& xnode, MatDoub &Mat_M)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D::compute_mass_matrix");

	// # nodes in 1D
	Myint nn = xnode.size() ;
	if (nn == 0)
	{
		print_error(" Error in FEM_2D::compute_mass_matrix, xnode.size() = 0") ;
		return(RTN_CODE_KO) ;
	}

	// check size of matrix
	Myint ni1= sqrt(Mat_M.ncols()) ;
	if (nn != ni1)
	{
		print_error(" Error in FEM_2D::compute_mass_matrix, nn != ni1") ;
		return(RTN_CODE_KO) ;
	}
	Myint ni2= sqrt(Mat_M.nrows()) ;
	if (nn != ni2)
	{
		print_error(" Error in FEM_2D::compute_mass_matrix, nn != ni2") ;
		return(RTN_CODE_KO) ;
	}

	// # GL points
	Myint ng = xg.size() ;
	if (ng == 0)
	{
		print_error(" Error in FEM_2D::compute_mass_matrix, xg.size() = 0") ;
		return(RTN_CODE_KO) ;
	}

	// initialize mass matrix
	for (Myint ii = 0; ii < Mat_M.nrows(); ii++)
	{
		for (Myint jj = 0; jj < Mat_M.ncols(); jj++)
		{
			Mat_M[ii][jj] = 0.0 ;
		}
	}

	Myint nxi   = nn ;
	Myint neta  = nn ;
	Myint netag = ng ;
	Myint nxig  = ng ;

	// loop on element nodes i
	Myint ii = -1 ;
	for (Myint ieta_i=0; ieta_i<neta; ieta_i++)
	{
		for (Myint ixi_i=0; ixi_i<nxi; ixi_i++)
		{
			ii++ ;

			// loop on element nodes j
			Myint jj = -1 ;
			for (Myint ieta_j=0; ieta_j<neta; ieta_j++)
			{
				for (Myint ixi_j=0; ixi_j<nxi; ixi_j++)
				{
					jj++ ;

					// loop on quadrature points
					for (Myint ieta_g=0; ieta_g<netag; ieta_g++)
					{
						for (Myint ixi_g=0; ixi_g<nxig; ixi_g++)
						{

							// compute phi_i at [xg(ixi_g), xg(ieta_g)]
							//=========================================

							// phi_ieta_i at xg[ieta_g]
							Doub phi_ieta_i = lagran(xnode, ieta_i, xg[ieta_g]) ;

							// phi_ixi_i at xg[ixi_g]
							Doub phi_ixi_i = lagran(xnode, ixi_i, xg[ixi_g]) ;

							// compute phi_i
							Doub phi_i = phi_ieta_i * phi_ixi_i ;

							// compute phi_j at [xg(ixi_g), xg(ieta_g)]
							//=========================================

							// phi_ieta_j at xg[ieta_g]
							Doub phi_ieta_j = lagran(xnode, ieta_j, xg[ieta_g]) ;

							// phi_ixi_j at xg[ixi_g]
							Doub phi_ixi_j = lagran(xnode, ixi_j, xg[ixi_g]) ;

							// compute phi_j
							Doub phi_j = phi_ieta_j * phi_ixi_j ;

							// update mass matrix
							//===================
							Mat_M[ii][jj] = Mat_M[ii][jj] + wg[ixi_g] * wg[ieta_g] * phi_i * phi_j ;

						} // for (Myint ixi_g=0; ixi_g<nxig; ixi_g++)
					} // for (Myint ieta_g=0; ieta_g<netag; ieta_g++)

				} // for (Myint ixi_j=0; ixi_j<nxi; ixi_j++)
			} // for (Myint ieta_j=0; ieta_j<neta; ieta_j++)

		} // for (Myint ixi_i=0; ixi_i<nxi; ixi_i++)
	} // for (Myint ieta_i=0; ieta_i<neta; ieta_i++)

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D::compute_mass_matrix");
	return(RTN_CODE_OK) ;
}

//-----------------------------------------------------------------------------------------
// compute stiffness matrix with quadrature rule such as Gauss-Legendre
// in the reference element with xi in [-1, 1] and eta [-1, 1]
//
// For 1ST ORDER:
// (Dz_i)jr = Int_elem_i_vol (d_phi_ir/dz x phi_ij) dv
// (Dx_i)jr = Int_elem_i_vol (d_phi_ir/dx x phi_ij) dv
//
// For 2ND ORDER:
// (Dz_i)jr = Int_elem_i_vol (d_phi_ir/dz x d_phi_ij/dz) dv
// (Dx_i)jr = Int_elem_i_vol (d_phi_ir/dx x d_phi_ij/dx) dv
//
// INPUT
// xg (float*) = array of quadrature coordinates
// wg (float*) = array of quadrature weights
// xnode (float*) = array of node coordinates in element
// axis = AXIS_Z or AXIS_X
//
// OUTPUT
// Mat_K (MatDoub) stiffness matrix
//-----------------------------------------------------------------------------------------
Rtn_code FEM_2D::compute_stiffness_matrix(VecDoub &xg, VecDoub& wg, VecDoub& xnode, MatDoub &Mat_K, Axis_type axis)
{

	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D::compute_stiffness_matrix");

	// only Z_AXIS and X_AXIS are allowed in 2D
	if ((axis != Z_AXIS) && (axis != X_AXIS))
	{
		print_error(" Error in FEM_2D::compute_stiffness_matrix, (axis != Z_AXIS) && (axis != X_AXIS)", axis) ;
		return(RTN_CODE_KO) ;
	}

	// # nodes in 1D
	Myint nn = xnode.size() ;
	if (nn == 0)
	{
		print_error(" Error in FEM_2D::compute_stiffness_matrix, xnode.size() = 0") ;
		return(RTN_CODE_KO) ;
	}

	// check size of matrix
	Myint ni1= sqrt(Mat_K.ncols()) ;
	if (nn != ni1)
	{
		print_error(" Error in FEM_2D::compute_stiffness_matrix, nn != ni1") ;
		return(RTN_CODE_KO) ;
	}
	Myint ni2= sqrt(Mat_K.nrows()) ;
	if (nn != ni2)
	{
		print_error(" Error in FEM_2D::compute_stiffness_matrix, nn != ni2") ;
		return(RTN_CODE_KO) ;
	}

	// # GL points
	Myint ng = xg.size() ;
	if (ng == 0)
	{
		print_error(" Error in FEM_2D::compute_stiffness_matrix, xg.size() = 0") ;
		return(RTN_CODE_KO) ;
	}

	// initialize mass matrix
	for (Myint ii = 0; ii < Mat_K.nrows(); ii++)
	{
		for (Myint jj = 0; jj < Mat_K.ncols(); jj++)
		{
			Mat_K[ii][jj] = 0.0 ;
		}
	}

	Myint nxi   = nn ;
	Myint neta  = nn ;
	Myint netag = ng ;
	Myint nxig  = ng ;

	// loop on element nodes i
	Myint ii = -1 ;
	for (Myint ieta_i=0; ieta_i<neta; ieta_i++)
	{
		for (Myint ixi_i=0; ixi_i<nxi; ixi_i++)
		{
			ii++ ;

			// loop on element nodes j
			Myint jj = -1 ;
			for (Myint ieta_j=0; ieta_j<neta; ieta_j++)
			{
				for (Myint ixi_j=0; ixi_j<nxi; ixi_j++)
				{
					jj++ ;

					// loop on quadrature points
					for (Myint ieta_g=0; ieta_g<netag; ieta_g++)
					{
						for (Myint ixi_g=0; ixi_g<nxig; ixi_g++)
						{

							// compute d(phi_i)/dxi at [xg(ixi_g), xg(ieta_g)]
							//================================================
							Doub dphi_i = 0.0 ;
							if (axis == X_AXIS)
							{
								// phi_ieta_i at xg[ieta_g]
								Doub phi_ieta_i = lagran(xnode, ieta_i, xg[ieta_g]) ;

								// d(phi_ixi_i)/dxi at xg[ixi_g]
								Doub dphi_ixi_i = dlagran(xnode, ixi_i, xg[ixi_g]) ;

								// compute phi_i
								dphi_i = phi_ieta_i * dphi_ixi_i ;
							}

							// compute d(phi_i)/deta at [xg(ixi_g), xg(ieta_g)]
							//================================================
							else if (axis == Z_AXIS)
							{
								// d(phi_ieta_i)deta at xg[ieta_g]
								Doub dphi_ieta_i = dlagran(xnode, ieta_i, xg[ieta_g]) ;

								// phi_ixi_i at xg[ixi_g]
								Doub phi_ixi_i = lagran(xnode, ixi_i, xg[ixi_g]) ;

								// compute phi_i
								dphi_i = dphi_ieta_i * phi_ixi_i ;
							}

							// compute phi_j at [xg(ixi_g), xg(ieta_g)]
							//=========================================

							// phi_ieta_j at xg[ieta_g]
							Doub phi_ieta_j = lagran(xnode, ieta_j, xg[ieta_g]) ;

							// phi_ixi_j at xg[ixi_g]
							Doub phi_ixi_j = lagran(xnode, ixi_j, xg[ixi_g]) ;

							// compute phi_j
							Doub phi_j = phi_ieta_j * phi_ixi_j ;

							// compute d(phi_j)/dxi at [xg(ixi_g), xg(ieta_g)]
							//================================================
							Doub dphi_j = 0.0 ;
							if (axis == X_AXIS)
							{
								// phi_ieta_i at xg[ieta_g]
								Doub phi_ieta_j = lagran(xnode, ieta_j, xg[ieta_g]) ;

								// d(phi_ixi_j)/dxi at xg[ixi_g]
								Doub dphi_ixi_j = dlagran(xnode, ixi_j, xg[ixi_g]) ;

								// compute phi_i
								dphi_j = phi_ieta_j * dphi_ixi_j ;
							}

							// compute d(phi_j)/deta at [xg(ixi_g), xg(ieta_g)]
							//================================================
							else if (axis == Z_AXIS)
							{
								// d(phi_ieta_j)deta at xg[ieta_g]
								Doub dphi_ieta_j = dlagran(xnode, ieta_j, xg[ieta_g]) ;

								// phi_ixi_j at xg[ixi_g]
								Doub phi_ixi_j = lagran(xnode, ixi_j, xg[ixi_g]) ;

								// compute phi_i
								dphi_j = dphi_ieta_j * phi_ixi_j ;
							}

							// update mass matrix
							//===================
							if (eq_order == ORDER_1ST) Mat_K[ii][jj] = Mat_K[ii][jj] + wg[ixi_g] * wg[ieta_g] * dphi_i * phi_j ;
							if (eq_order == ORDER_2ND) Mat_K[ii][jj] = Mat_K[ii][jj] + wg[ixi_g] * wg[ieta_g] * dphi_i * dphi_j ;

						} // for (Myint ixi_g=0; ixi_g<nxig; ixi_g++)
					} // for (Myint ieta_g=0; ieta_g<netag; ieta_g++)

				} // for (Myint ixi_j=0; ixi_j<nxi; ixi_j++)
			} // for (Myint ieta_j=0; ieta_j<neta; ieta_j++)

		} // for (Myint ixi_i=0; ixi_i<nxi; ixi_i++)
	} // for (Myint ieta_i=0; ieta_i<neta; ieta_i++)

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D::compute_stiffness_matrix");
	return(RTN_CODE_OK) ;
}

//-----------------------------------------------------------------------------------------
Rtn_code FEM_2D::write_mesh_VTK(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D::write_mesh_VTK");

	double t0 = MPI_Wtime() ;

	Myfloat znode_min, znode_max, xnode_min, xnode_max ;

	//================================================================================================
	// VTK file
	//================================================================================================

	// VTK header
	ofstream out_file ;
	out_file.open(MESH_VTK_OUT_FILE, ios::trunc | ios::out) ;
	out_file << "# vtk DataFile Version 2.0\n" ;
	out_file << "MESH FILE DJANGO\n" ;
	out_file << "ASCII\n" ;
	out_file << "DATASET UNSTRUCTURED_GRID\n" ;

	// points
	out_file << "POINTS " << nnode << " float\n" ;
	for (Myint inode = 0; inode < nnode; inode++)
	{
		out_file << pNode[inode].xcoord << " 0.0 " << pNode[inode].zcoord << "\n" ;
	}

	// cells (quadrangle)
	out_file << "CELLS " << nelem << " " << nelem*5 << "\n" ;
	for (Myint ielem = 0; ielem < nelem; ielem++)
	{
		// number of nodes in quad
		out_file << 4 << " " ;
		// node down/left
		Myint node_tmp = 0 ;
		out_file << pElement[ielem].inode[node_tmp] << " " ;
		// node down/right
		node_tmp = pElement[ielem].nnode_1D-1 ;
		out_file << pElement[ielem].inode[node_tmp] << " " ;
		// node up/right
		node_tmp = pElement[ielem].nnode-1 ;
		out_file << pElement[ielem].inode[node_tmp] << " " ;
		// node up/left
		node_tmp = pElement[ielem].nnode-pElement[ielem].nnode_1D ;
		out_file << pElement[ielem].inode[node_tmp] << "\n" ;
	}

	// cell types (VTK_QUAD = 9)
	out_file << "CELL_TYPES " << nelem << "\n" ;
	for (Myint ielem = 0; ielem < nelem; ielem++)
	{
		out_file << "9 \n" ;
	}

	// polynomial order
	out_file << "CELL_DATA " << nelem << "\n" ;
	out_file << "SCALARS polynomial float 1 \n" ;
	out_file << "LOOKUP_TABLE default \n" ;
	for (Myint ielem = 0; ielem < nelem; ielem++)
	{
		out_file << pElement[ielem].nnode_1D - 1 << "\n" ;
	}

	// tmin (only for dynamic front)
	if (front_type != FRONT_STATIC)
	{
		out_file << "SCALARS tmin float 1 \n" ;
		out_file << "LOOKUP_TABLE default \n" ;
		for (Myint ielem = 0; ielem < nelem; ielem++)
		{
			out_file << pElement[ielem].tmin << "\n" ;
		}
	}

	out_file.close() ;

	double time_in_function = MPI_Wtime() - t0 ;
	print_info(MASTER, " Time write mesh (s)", (Myfloat) time_in_function) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D::write_mesh_VTK");
	return(RTN_CODE_OK) ;
}

//-----------------------------------------------------------------------------------------
Rtn_code FEM_2D::write_node_VTK(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D::write_node_VTK");

	ofstream out_file ;
	out_file.open(NODE_VTK_OUT_FILE, ios::trunc | ios::out) ;
	out_file << "# vtk DataFile Version 2.0\n" ;
	out_file << "MESH FILE DJANGO\n" ;
	out_file << "ASCII\n" ;
	out_file << "DATASET POLYDATA\n" ;
	out_file << "POINTS " << nnode << " float\n" ;
	for (Myint inode = 0; inode < nnode; inode++)
	{
		out_file << pNode[inode].xcoord << " 0.0 " << pNode[inode].zcoord << "\n" ;
	}
	out_file.close() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D::write_node_VTK");
	return(RTN_CODE_OK) ;
}

//-----------------------------------------------------------------------------------------
Rtn_code FEM_2D::mesh_info(void)
{

	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D::mesh_info");

	// call parent
	FEM::mesh_info() ;

	if (nelem <= 0)
	{
		print_warning("Mesh has no element!") ;
		return(RTN_CODE_OK) ;
	}

	// retrieve min and max z coordinate
	min_z_mesh = FLT_MAX ;
	max_z_mesh = FLT_MIN ;
	Myfloat min_z_med = FLT_MAX ;
	Myfloat max_z_med = FLT_MIN ;
	min_x_mesh = FLT_MAX ;
	max_x_mesh = FLT_MIN ;
	Myfloat min_x_med = FLT_MAX ;
	Myfloat max_x_med = FLT_MIN ;

	Myint min_nnode = pElement[0].nnode ;
	Myint max_nnode = min_nnode ;

	// loop on all element
	for (Myint ielem = 0; ielem < nelem; ielem++)
	{
		// loop on all nodes of current element
		for (Myint ii=0; ii<pElement[ielem].nnode; ii++)
		{
			min_z_mesh = min(min_z_mesh, pNode[pElement[ielem].inode[ii]].zcoord) ;
			max_z_mesh = max(max_z_mesh, pNode[pElement[ielem].inode[ii]].zcoord) ;
			min_x_mesh = min(min_x_mesh, pNode[pElement[ielem].inode[ii]].xcoord) ;
			max_x_mesh = max(max_x_mesh, pNode[pElement[ielem].inode[ii]].xcoord) ;
			if (pElement[ielem].region == MEDIUM)
			{
				min_z_med = min(min_z_med, pNode[pElement[ielem].inode[ii]].zcoord) ;
				max_z_med = max(max_z_med, pNode[pElement[ielem].inode[ii]].zcoord) ;
				min_x_med = min(min_x_med, pNode[pElement[ielem].inode[ii]].xcoord) ;
				max_x_med = max(max_x_med, pNode[pElement[ielem].inode[ii]].xcoord) ;
			}
		}

		min_nnode = min(min_nnode, pElement[ielem].nnode) ;
		max_nnode = max(max_nnode, pElement[ielem].nnode) ;

	}

	// retrieve min and max element size
	Myfloat min_elem_size_x = FLT_MAX ;
	Myfloat max_elem_size_x = FLT_MIN ;
	Myfloat min_elem_size_z = FLT_MAX ;
	Myfloat max_elem_size_z = FLT_MIN ;
	for (Myint ielem = 0; ielem < nelem; ielem++)
	{
		min_elem_size_x = min(min_elem_size_x, pElement[ielem].size_x) ;
		max_elem_size_x = max(max_elem_size_x, pElement[ielem].size_x) ;
		min_elem_size_z = min(min_elem_size_z, pElement[ielem].size_z) ;
		max_elem_size_z = max(max_elem_size_z, pElement[ielem].size_z) ;
	}

	print_info(MASTER, "") ;
	print_info(MASTER, " min z in mesh (m)", min_z_mesh) ;
	print_info(MASTER, " max z in mesh (m)", max_z_mesh) ;
	print_info(MASTER, " min x in mesh (m)", min_x_mesh) ;
	print_info(MASTER, " max x in mesh (m)", max_x_mesh) ;
	print_info(MASTER, "") ;
	print_info(MASTER, " min z medium reg. (m)", min_z_med) ;
	print_info(MASTER, " max z medium reg. (m)", max_z_med) ;
	print_info(MASTER, " min x medium reg. (m)", min_x_med) ;
	print_info(MASTER, " max x medium reg. (m)", max_x_med) ;
	print_info(MASTER, "") ;
	print_info(MASTER, " min elem size x (m)", min_elem_size_x) ;
	print_info(MASTER, " max elem size x (m)", max_elem_size_x) ;
	print_info(MASTER, " min elem size z (m)", min_elem_size_z) ;
	print_info(MASTER, " max elem size z (m)", max_elem_size_z) ;

	if (min_elem_size_x <= 0.0)
	{
		print_error(" Error in FEM_2D::mesh_info, min_elem_size_x <= 0.0", min_elem_size_x) ;
		return(RTN_CODE_KO) ;
	}
	if (max_elem_size_x <= 0.0)
	{
		print_error(" Error in FEM_2D::mesh_info, max_elem_size_x <= 0.0", max_elem_size_x) ;
		return(RTN_CODE_KO) ;
	}
	if (min_elem_size_z <= 0.0)
	{
		print_error(" Error in FEM_2D::mesh_info, min_elem_size_z <= 0.0", min_elem_size_z) ;
		return(RTN_CODE_KO) ;
	}
	if (max_elem_size_z <= 0.0)
	{
		print_error(" Error in FEM_2D::mesh_info, max_elem_size_z <= 0.0", max_elem_size_z) ;
		return(RTN_CODE_KO) ;
	}

	// neighbours
	Myint nneigh_min = INT_MAX ;
	Myint nneigh_max = INT_MIN ;

	// loop on all element
	for (Myint ielem = 0; ielem < nelem; ielem++)
	{
		// loop on all nighbours
		Myint nneigh = 0 ;
		for (Myint ii=0; ii<4; ii++)
		{
			if (pElement[ielem].neigh[ii] != NO_NEIGH) nneigh++ ;
		}
		nneigh_min = min(nneigh_min, nneigh) ;
		nneigh_max = max(nneigh_max, nneigh) ;
	}
	print_info(MASTER, "") ;
	print_info(MASTER, " min neigh. per element", nneigh_min) ;
	print_info(MASTER, " max neigh. per element", nneigh_max) ;

	print_info(MASTER, "") ;
	print_info(MASTER, " # total nodes\t", nnode) ;

	// type of node
	Myint nnode_no_boundary = 0 ;
	Myint nnode_pml         = 0 ;
	Myint nnode_random      = 0 ;
	Myint nnode_spg         = 0 ;
	Myint nnode_freesurf    = 0 ;
	Myint nnode_rigid       = 0 ;
	Myint nnode_internal    = 0 ;

	for (Myint inode = 0; inode < nnode; inode++)
	{
		if (pNode[inode].boundary == NO_BOUNDARY)
		{
			nnode_no_boundary++ ;
		}
		else if (pNode[inode].boundary == PML)
		{
			nnode_pml++ ;
		}
		else if (pNode[inode].boundary == RANDOM)
		{
			nnode_random++ ;
		}
		else if (pNode[inode].boundary == SPG)
		{
			nnode_spg++ ;
		}
		else if (pNode[inode].boundary == FREESURF)
		{
			nnode_freesurf++ ;
		}
		else if (pNode[inode].boundary == RIGID)
		{
			nnode_rigid++ ;
		}
		else if (pNode[inode].boundary == INTERNAL)
		{
			nnode_internal++ ;
		}
	}

	if (nnode_no_boundary > 0) print_info(MASTER, " # no bound. nodes", nnode_no_boundary) ;
	if (nnode_pml > 0)         print_info(MASTER, " # PML bound. nodes", nnode_pml) ;
	if (nnode_random > 0)      print_info(MASTER, " # rand. bound. nodes", nnode_random) ;
	if (nnode_spg > 0)         print_info(MASTER, " # spg. bound. nodes", nnode_spg) ;
	if (nnode_freesurf > 0)    print_info(MASTER, " # free surf. nodes", nnode_freesurf) ;
	if (nnode_rigid > 0)       print_info(MASTER, " # rigid surf. nodes", nnode_rigid) ;
	if (nnode_internal > 0)    print_info(MASTER, " # internal nodes", nnode_internal) ;

	print_info(MASTER, "") ;
	print_info(MASTER, " min node in element", min_nnode) ;
	print_info(MASTER, " max node in element", max_nnode) ;

	// flux type
	Myint nNO_FLUX_TYPE       = 0 ;
	Myint nFLUX_TYPE_CENTERED = 0 ;
	Myint nFLUX_TYPE_UPWIND   = 0 ;
	for (Myint ielem = 0; ielem < nelem; ielem++)
	{
		for (Myint ineigh = 0; ineigh < 4; ineigh++)
		{
			if (pElement[ielem].flux[ineigh] == NO_FLUX_TYPE) nNO_FLUX_TYPE++ ;
			if (pElement[ielem].flux[ineigh] == FLUX_TYPE_CENTERED) nFLUX_TYPE_CENTERED++ ;
			if (pElement[ielem].flux[ineigh] == FLUX_TYPE_UPWIND) nFLUX_TYPE_UPWIND++ ;
		}
	}
	print_info(MASTER, "") ;
	if (nNO_FLUX_TYPE > 0) print_info(MASTER, " # Face with NO flux", nNO_FLUX_TYPE) ;
	if (nFLUX_TYPE_CENTERED > 0) print_info(MASTER, " # Face CENTERED flux", nFLUX_TYPE_CENTERED) ;
	if (nFLUX_TYPE_UPWIND > 0) print_info(MASTER, " # Face UPWIND flux", nFLUX_TYPE_UPWIND) ;
	print_info(MASTER, "") ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D::mesh_info");
	return(RTN_CODE_OK) ;
}

//-----------------------------------------------------------------------------------------
Rtn_code FEM_2D::compute_element_matrices(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_2D::compute_element_matrices");

	// open output file with matrix layout
	ofstream out_file ;
	out_file.open(FEM_MATRIX_OUT_FILE, ios::trunc | ios::out) ;

	//-----------------------------------------------------------------------------------------
	// loop on mass matrices
	//-----------------------------------------------------------------------------------------

	for (Myint imat=0; imat <= MAX_POLY_ORDER; imat++)
	{
		Myint ni = imat+1 ;

		out_file << "\n====================================\n" ;
		out_file <<   "   MATRIX FOR POLYNOMIAL ORDER " << imat << endl ;
		out_file <<   "====================================\n" ;

		//------------------------------------
		// compute quadrature points
		// used for numerical integration
		// in FEM matrix computation
		//------------------------------------

		Myint ng = ni ;
		Doub x1 = -1.0 ;
		Doub x2 = +1.0 ;
		VecDoub xg(ng) ;
		VecDoub wg(ng) ;

		//------------------------------------
		// compute quadrature points with nr3
		// --> case Gauss-Legendre
		//------------------------------------
		if (node_integ == NODE_INTEG_GL)
		{
			gauleg(x1, x2, xg, wg) ;
			out_file << "\n*** GL quadrature ***\n" ;
		}
		//------------------------------------
		// compute quadrature points with gll.h
		// --> case Gauss-Lobatto-Legendre
		//------------------------------------
		else if (node_integ == NODE_INTEG_GLL)
		{
			if (ni == 1)
			{
				gauleg(x1, x2, xg, wg) ;
			}
			else
			{
				gll(xg, wg) ;
			}
			out_file << "\n*** GLL quadrature ***\n" ;
		}
		else
		{
			print_error(" Error in FEM_2D::compute_element_matrices, numerical integration not supported", node_integ) ;
			return(RTN_CODE_KO) ;
		}

		for (Myint ii=0; ii<xg.size(); ii++)
		{
			out_file << "xg[" << ii << "] " << xg[ii] << "\n" ;
			out_file << "wg[" << ii << "] " << wg[ii] << "\n" ;
		}

		//--------------------------------------
		// retrieve coordinates of element nodes
		//--------------------------------------
		VecDoub xnode = xnode_ref_elem_1D[imat] ;

		out_file << "\n*** Element nodes ***\n" ;
		for (Myint ii=0; ii<ni; ii++)
		{
			out_file << "xnode[" << ii << "] " << xnode[ii] << "\n" ;
		}

		//---------------------
		// compute mass matrix
		//---------------------
		MatDoub Mat_M2(ni*ni, ni*ni) ;
		Rtn_code rtn_code = compute_mass_matrix(xg, wg, xnode, Mat_M2) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		out_file << "\n*** Matrix Mat_M *** \n" ;
		for (int ii=0; ii< Mat_M2.nrows(); ii++)
		{
			for (int jj=0; jj< Mat_M2.ncols(); jj++)
			{
				out_file << Mat_M2[ii][jj] << " " ;
			}
			out_file << endl ;
		}

		//--------------------------------
		// compute inverse of mass matrix
		//--------------------------------
		pMat_M_inv[imat] = new Grid_2D_float(ni*ni, ni*ni, 1, 1) ;
		Myfloat** pMatM = pMat_M_inv[imat]->pArray ;
		gaussj(Mat_M2) ;
		out_file << "\n*** Invert Matrix Mat_M ***\n" ;

		for (int ii=0; ii< Mat_M2.nrows(); ii++)
		{
			for (int jj=0; jj< Mat_M2.ncols(); jj++)
			{
				pMatM[ii][jj] = Mat_M2[ii][jj] ;
				out_file << pMatM[ii][jj] << " " ;
			}
			out_file << endl ;
		}

		//-------------------------------------
		// compute stiffness matrices Dz and Dx
		//-------------------------------------

		// Dz
		{
			MatDoub Mat_K2(ni*ni, ni*ni) ;
			Axis_type axis = Z_AXIS ;
			rtn_code = compute_stiffness_matrix(xg, wg, xnode, Mat_K2, axis) ;
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

			out_file << "\n*** Matrix Mat_Dz *** \n" ;
			pMat_Dz[imat] = new Grid_2D_float(ni*ni, ni*ni, 1, 1) ;
			Myfloat** Mat_K = pMat_Dz[imat]->pArray ;

			for (int ii=0; ii< Mat_K2.nrows(); ii++)
			{
				for (int jj=0; jj< Mat_K2.ncols(); jj++)
				{
					Mat_K[ii][jj] = Mat_K2[ii][jj] ;
					out_file << Mat_K[ii][jj] << " " ;
				}
				out_file << endl ;
			}

			// convert full matrix to sparse matrix
			full2sparse(&pMat_Dz_i[imat], &pMat_Dz_j[imat], &pMat_Dz_c[imat], Mat_K2) ;

			// print sparse matrix
			Myint*   Mat_Dz_i = pMat_Dz_i[imat]->pArray ;
			Myint*   Mat_Dz_j = pMat_Dz_j[imat]->pArray ;
			Myfloat* Mat_Dz_c = pMat_Dz_c[imat]->pArray ;
			out_file << "\ni j c / sparse format\n" ;
			for (Myint ii=0; ii< pMat_Dz_i[imat]->nz; ii++)
			{
				out_file << Mat_Dz_i[ii] << " " << Mat_Dz_j[ii] << " " << Mat_Dz_c[ii] << endl ;
			}
		}

		// Dx
		{
			MatDoub Mat_K2(ni*ni, ni*ni) ;
			Axis_type axis = X_AXIS ;
			rtn_code = compute_stiffness_matrix(xg, wg, xnode, Mat_K2, axis) ;
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

			out_file << "\n*** Matrix Mat_Dx *** \n" ;
			pMat_Dx[imat] = new Grid_2D_float(ni*ni, ni*ni, 1, 1) ;
			Myfloat** Mat_K = pMat_Dx[imat]->pArray ;

			for (int ii=0; ii< Mat_K2.nrows(); ii++)
			{
				for (int jj=0; jj< Mat_K2.ncols(); jj++)
				{
					Mat_K[ii][jj] = Mat_K2[ii][jj] ;
					out_file << Mat_K[ii][jj] << " " ;
				}
				out_file << endl ;
			}

			// convert full matrix to sparse matrix
			full2sparse(&pMat_Dx_i[imat], &pMat_Dx_j[imat], &pMat_Dx_c[imat], Mat_K2) ;

			// print sparse matrix
			Myint*   Mat_Dx_i = pMat_Dx_i[imat]->pArray ;
			Myint*   Mat_Dx_j = pMat_Dx_j[imat]->pArray ;
			Myfloat* Mat_Dx_c = pMat_Dx_c[imat]->pArray ;
			out_file << "\ni j c / sparse format\n" ;
			for (Myint ii=0; ii< pMat_Dx_i[imat]->nz; ii++)
			{
				out_file << Mat_Dx_i[ii] << " " << Mat_Dx_j[ii] << " " << Mat_Dx_c[ii] << endl ;
			}
		}

		//----------------------
		// compute flux matrices
		//----------------------
		rtn_code = compute_flux_matrices(imat, &out_file, ni, xg, wg, xnode) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	} // for (Myint imat=0; imat <= MAX_POLY_ORDER; imat++)

	print_info(MASTER, " Compute FEM matrices", "OK") ;
	print_info(MASTER, " Max polynomial order", MAX_POLY_ORDER) ;
	out_file.close() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_2D::compute_element_matrices");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FEM_2D::write_snapshot(Variable* pVar, Myint it)
{
	print_debug(ALL, FULL_DEBUG, "IN FEM_2D::write_snapshot");

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
			Grid_1D_float *val_grid = (Grid_1D_float*) (pVar->get_grid()) ;
			Myfloat* const val = val_grid->pArray ;

			// allocate snapshot
			Myfloat* snap = allocate_array<Myfloat>(npixel) ;

			// extract snapshot
			for (Myint ipixel=0; ipixel<npixel; ipixel++)
			{
				if (ielem_pixel[ipixel] == NOT_FOUND)
				{
					snap[ipixel] = 0.0 ;
				}
				else
				{
					snap[ipixel] = interpolate_variable(val, ielem_pixel[ipixel], eta_pixel[ipixel], xi_pixel[ipixel]) ;
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

	print_debug(ALL, FULL_DEBUG, "OUT FEM_2D::write_snapshot");
	return(RTN_CODE_OK) ;
}

} // namespace django
