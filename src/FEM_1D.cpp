//-------------------------------------------------------------------------------------------------------
//
// MODELLING WITH FEM IN 1D
//
// PARENT CLASS: Scheme
//   DERIVED CLASS: FEM
//     DERIVED CLASS: FEM_1D
//
//-------------------------------------------------------------------------------------------------------

#include "FEM_1D.h"

#include <algorithm>
#include <cassert>
#include <cfloat>
#include <climits>
#include <fstream>
#include <iostream>

#include "nr3.h"
#include "bessel.h"

#include "allocate_array.h"
#include "grid_1D_float.h"
#include "output_report.h"
#include "singleton.h"
#include "type_def.h"

namespace django {

//-------------------------------------------------------------------------------------------------------

using namespace std;

FEM_1D::FEM_1D(void) : FEM()
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D::FEM_1D");

	dim = ONE ;

	for (Myint imat = 0; imat <= MAX_POLY_ORDER; imat++)
	{
		pMat_Dz[imat] = NULL ;
	}

	pElement = NULL ;
	pNode    = NULL ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D::FEM_1D");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FEM_1D::project_variable(FEM_1D& schemeSrc,
		Myint varSrcId, Myint varDestId)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D::project_variable");

	cout << "### varSrcId " << varSrcId << " varDestId " << varDestId << "\n" ;

	Variable* varSrc  = Singleton::Instance()->get_variable(varSrcId) ;
	if (varSrc == NULL)
	{
		print_error("FEM_1D::project_variable, varSrc == NULL") ;
		return(RTN_CODE_KO) ;
	}

	Variable* varDest = Singleton::Instance()->get_variable(varDestId) ;
	if (varDest == NULL)
	{
		print_error("FEM_1D::project_variable, varDest == NULL") ;
		return(RTN_CODE_KO) ;
	}

	// check src and dest variable type are identical
	if (varSrc->get_type() != varDest->get_type())
	{
		print_error("FEM_1D::project_variable, varSrc->get_type() != varDest->get_type()") ;
		return(RTN_CODE_KO) ;
	}

	// check number of elements in mesh src and dest are identical
	// TODO this is a current limitation
	// Only meshes with same geometry are supported
	// however, polynomials in elements can be different between src and dest
	if (this->nelem != schemeSrc.nelem)
	{
		print_error("FEM_1D::project_variable, this->nelem != schemeSrc,nelem") ;
		return(RTN_CODE_KO) ;
	}

	// get src grid
	Grid_1D_float* srcGrid = dynamic_cast<Grid_1D_float*>(varSrc->get_grid()) ;
	if (srcGrid == NULL)
	{
		print_error("FEM_1D::project_variable, srcGrid == NULL");
		return(RTN_CODE_KO) ;
	}
	Myfloat* srcArray = srcGrid->pArray ;
	//cout << "srcGrid->get_min()" << srcGrid->get_min() << "\n" ;
	//cout << "srcGrid->get_max()" << srcGrid->get_max() << "\n" ;

	// get dest grid
	Grid_1D_float* destGrid = dynamic_cast<Grid_1D_float*>(varDest->get_grid()) ;
	if (destGrid == NULL)
	{
		print_error("FEM_1D::project_variable, destGrid == NULL");
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

		// loop on nodes of element
		for (Myint inode=0; inode< pElement[iel].nnode ; inode++)
		{
			// get global node number
			Myint inodeGlob = pElement[iel].inode[inode] ;

			// get node coordinates
			Myfloat znode = pNode[inodeGlob].zcoord ;

			// convert to eta coordinate in ref element
			Myfloat etaNode = -1.0 + 2.0*(znode - znode_min)/pElement[iel].size  ;

			// call interpolate_variable
			Myfloat valNode = schemeSrc.interpolate_variable(srcArray, iel, etaNode) ;

			// store value in destination variable
			destArray[inodeGlob] = valNode ;
		}
	}

	//cout << "after destGrid->get_min()" << destGrid->get_min() << "\n" ;
	//cout << "after destGrid->get_max()" << destGrid->get_max() << "\n" ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D::project_variable");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FEM_1D::initialize_mesh(Model* pModel, Myfloat fmax)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D::initialize_mesh");

	// check model type is GRID
	//=========================

	if (pModel == NULL)
	{
		print_error("IN FEM_1D::initialize_mesh --> pModel == NULL");
		return(RTN_CODE_KO) ;
	}

	if (pModel->get_type() != GRID)
	{
		print_error("IN FEM_1D::initialize_mesh --> model type is not GRID");
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

	// get boundary width
	//===================
	Myint nlayer_zBeg = get_boundary_width(ZBEG) ;
	Myint nlayer_zEnd = get_boundary_width(ZEND) ;

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

	// set number of elements
	//=======================
	nelem_z = nelem_z_med + nlayer_zBeg + nlayer_zEnd + nghost_zBeg + nghost_zEnd ;
	nelem = nelem_z ;

	// allocate element array
	//=======================

	print_debug(ALL, MID_DEBUG, "allocate element array") ;
	if (nelem <= 0)
	{
		print_error("IN FEM_1D::initialize_mesh --> nelem <= 0");
		return(RTN_CODE_KO) ;
	}
	pElement = allocate_array<Segment_struct_type>(nelem) ;

	// retrieve model parameters
	//==========================

	// get VP model
	Variable* vp_var = NULL ;
	Grid* vp_grid = NULL ;
	Grid_1D_float *vp_1D_grid = NULL ;
	vp_var = pModel->get_parameter(VP) ;
	if (vp_var == NULL)
	{
		print_error("IN FEM_1D::initialize_mesh --> VP model not found");
		return(RTN_CODE_KO) ;
	}
	// get vp grid
	vp_grid = vp_var->get_grid() ;
	if (vp_grid == NULL)
	{
		print_error("IN FEM_1D::initialize_mesh --> VP grid not initialized");
		return(RTN_CODE_KO) ;
	}
	vp_1D_grid = dynamic_cast<Grid_1D_float*>(vp_grid) ;
	if (vp_1D_grid == NULL)
	{
		print_error("IN FEM_1D::initialize_mesh --> VP grid is not Grid_1D_float");
		return(RTN_CODE_KO) ;
	}

	// determine element region
	//=========================================

	// $$$ tmp, does not allow variable grid spacing including random spacing
	Myfloat el_sizez = (pModel->get_zcoord()->get_max() - pModel->get_zcoord()->get_min()) / nz ;

	Myfloat dz_model = pModel->get_dz() ;

	print_debug(ALL, MID_DEBUG, "determine element region") ;
	Myint nelem_no_region = 0 ;

	Myint izGhostBeg  = 0 ;
	Myint izLayerBeg  = izGhostBeg + nghost_zBeg ;
	Myint izMediumBeg = izLayerBeg + nlayer_zBeg ;
	Myint izMediumEnd = izMediumBeg + nz-1 ;
	Myint izLayerEnd  = izMediumEnd + nlayer_zEnd ;
	Myint izGhostEnd  = izLayerEnd + nghost_zEnd ;

	Myint iel = -1 ;
	for (Myint iz = 0; iz < nelem_z; iz++)
	{
		iel++ ;

		// determine region
		if ((iz >= izMediumBeg) && (iz<= izMediumEnd))
		{
			nelem_med++ ;
			pElement[iel].region = MEDIUM ;
		}
		else if ((iz >= izLayerBeg) && (iz<= izLayerEnd))
		{
			nelem_lay++ ;
			pElement[iel].region = LAYER ;
		}
		else if ((iz >= izGhostBeg) && (iz<= izGhostEnd))
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
	// check region distribution
	if (nelem != (nelem_ghost + nelem_lay + nelem_med))
	{
		print_error(" Error in FEM_1D::initialize_mesh, error in region distribution", nelem_no_region) ;
		return(RTN_CODE_KO) ;
	}

	// retrieve min velocity
	// this value will be used for p-adaptivity
	//=========================================
	print_debug(ALL, MID_DEBUG, "retrieve min velocity") ;

	// determine mesh origin coordinate
	// $$$ tmp, does not allow variable grid spacing including random spacing
	const Myfloat z_origin_mesh =  0.0 - el_sizez * (nlayer_zBeg + nghost_zBeg) ;

	iel = -1 ;
	for (Myint iz = 0; iz < nelem_z; iz++)
	{
		iel++ ;

		Myint iz_model = (z_origin_mesh + ((iz + 0) * el_sizez)) / dz_model ;
		iz_model = max(0,  iz_model) ;
		iz_model = min(nz, iz_model) ;
		Myfloat v1 = vp_1D_grid->pArray[iz_model] ;

		iz_model = (z_origin_mesh + ((iz + 1) * el_sizez)) / dz_model ;
		iz_model = max(0,  iz_model) ;
		iz_model = min(nz_model-1, iz_model) ;
		Myfloat v2 = vp_1D_grid->pArray[iz_model] ;

		iz_model = (z_origin_mesh + ((iz + 0.5) * el_sizez)) / dz_model ;
		iz_model = max(0,  iz_model) ;
		iz_model = min(nz_model-1, iz_model) ;
		Myfloat v0 = vp_1D_grid->pArray[iz_model] ;

		// retrieve vp_min
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
			vp_min = (v1 + v2) / 2.0 ;
			vp_max = vp_min ;
		}

		// PROP_LINEAR
		else if (prop == PROP_LINEAR)
		{
			// take minimum of values at 4 corners
			vp_min = min(v1,v2) ;

			// take maximum of values at 4 corners
			vp_max = max(v1,v2) ;
		}

		pElement[iel].vmax = vp_max ;
		if (eq_type == ACOUSTIC)
		{
			nelem_ac_iso++ ;
			pElement[iel].vmin = vp_min ;
		}
		else if (eq_type == AC_LOSSY)
		{
			nelem_ac_lossy++ ;
			pElement[iel].vmin = vp_min ;
		}
	}

	// check ac/el distribution
	if (nelem != (nelem_ac_iso + nelem_ac_lossy + nelem_el_iso))
	{
		print_error(" Error in FEM_1D::initialize_mesh, error in region ac/el distribution") ;
		return(RTN_CODE_KO) ;
	}

	// determine neighbour elements
	//=============================

	print_debug(ALL, MID_DEBUG, "determine neighbour elements") ;
	iel = -1 ;
	for (Myint iz = 0; iz < nelem_z; iz++)
	{
		iel++ ;

		if (iz == izGhostBeg)
		{
			pElement[iel].neigh[I_ZPREV] = NO_NEIGH ;
		}
		else
		{
			pElement[iel].neigh[I_ZPREV] = iel - 1 ;
		}

		if (iz == izGhostEnd)
		{
			pElement[iel].neigh[I_ZNEXT] = NO_NEIGH ;
		}
		else
		{
			pElement[iel].neigh[I_ZNEXT] = iel + 1 ;
		}
	}

	// initialize size of elements
	//============================

	print_debug(ALL, MID_DEBUG, "initialize size of elements") ;

	// get pointer to zcoord array
	if (pModel->get_zcoord() == NULL)
	{
		print_error(" Error in FEM_1D::initialize_mesh, pModel->get_zcoord() is NULL") ;
		return(RTN_CODE_KO) ;
	}
	Myfloat* zcoord = pModel->get_zcoord()->pArray ;
	if (zcoord == NULL)
	{
		print_error(" Error in FEM_1D::initialize_mesh, zcoord is NULL") ;
		return(RTN_CODE_KO) ;
	}

	// loop on all elements
	Myint ielem = -1 ;
	for (Myint iz = 0; iz < nelem_z; iz++)
	{
		ielem++ ;
		// case model type is GRID
		// elements are defined by 2 consecutive points in the model
		//pElement[ielem].size = zcoord[ielem+1] - zcoord[ielem] ;

		pElement[ielem].size = el_sizez ;
	}

	// initialize number of nodes per element
	//=======================================

	print_debug(ALL, MID_DEBUG, "initialize number of nodes per element") ;

	// check node distribution
	if ((node_distrib != NODE_DISTRIB_UNIFORM) &&
			(node_distrib != NODE_DISTRIB_MODEL) &&
			(node_distrib != NODE_DISTRIB_RANDOM))
	{
		print_error(" Error in FEM_1D::initialize_mesh, unsupported node_distrib", node_distrib) ;
		return(RTN_CODE_KO) ;
	}

	// check number of node per element
	Myint polynom_order = Singleton::Instance()->space_order ;
	if ( (polynom_order < 0) || (polynom_order > MAX_POLY_ORDER) )
	{
		print_error(" Error in FEM_1D::initialize_mesh, invalid polynomial order ", polynom_order) ;
		return(RTN_CODE_KO) ;
	}

	// special case
	// for CGM, one node element is not allowed
	if ((type == SCHEME_CGM) && (polynom_order == 0))
	{
		print_error(" Error in FEM_1D::initialize_mesh, one node element not allowed for CGM") ;
		return(RTN_CODE_KO) ;
	}

	// special case
	// for CGM, random number of node element is not yet supported
	if ((type == SCHEME_CGM) && (node_distrib == NODE_DISTRIB_RANDOM))
	{
		print_error(" Error in FEM_1D::initialize_mesh, random node distribution not yet supported for CGM") ;
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
					print_error(" Error in FEM_1D::initialize_mesh, invalid nnode_1D", pElement[ielem].nnode_1D) ;
					return(RTN_CODE_KO) ;
				}
			}
		}
		//===========================================================================
		// L O C A L   P - A D A P T I V I T Y
		//===========================================================================

		else if (node_distrib == NODE_DISTRIB_MODEL)
		{
			// determine nelem / lambda min
			Myfloat min_size = pElement[ielem].size   ;
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
		// in 1D, nnode = nnode_1D
		pElement[ielem].nnode = pElement[ielem].nnode_1D ;
	}

	if (node_distrib == NODE_DISTRIB_MODEL)
	{
		print_info(MASTER, " * p-adapt. mesh / fmax (Hz)", fmax) ;
	}

	// determine total number of nodes in the mesh
	//============================================

	print_debug(ALL, MID_DEBUG, "determine total number of nodes in the mesh") ;
	nnode = 0 ;
	Myint nnode_z ;
	for (Myint ielem = 0; ielem < nelem; ielem++)
	{
		// DGM case
		// nodes are not shared between elements
		if ((type == SCHEME_DGM) || (type == SCHEME_MGM))
		{
			nnode += pElement[ielem].nnode ;
		}

		// CGM case
		// nodes are shared between elements
		else if (type == SCHEME_CGM)
		{
			if (pElement[ielem].neigh[I_ZPREV] == NO_NEIGH)
			{
				// no previous element
				nnode += pElement[ielem].nnode ;
			}
			else
			{
				// one node is shared with previous element (-1)
				nnode += pElement[ielem].nnode - 1 ;
			}
		}
	}

	// allocate global node array
	//===========================

	print_debug(ALL, MID_DEBUG, "allocate global node array") ;
	if (nnode <= 0)
	{
		print_error("IN FEM_1D::initialize_mesh --> nnode <= 0");
		return(RTN_CODE_KO) ;
	}
	pNode = allocate_array<Node_1D_struct_type>(nnode) ;

	// check node type
	if ((node_type != NODE_TYPE_EQD) &&
			(node_type != NODE_TYPE_GLL))
	{
		print_error(" Error in FEM_1D::initialize_mesh, unsupported node_type", node_type) ;
		return(RTN_CODE_KO) ;
	}

	// initialize global node array
	// compute node coordinates
	//=============================

	print_debug(ALL, MID_DEBUG, "initialize global node array") ;

	Myint inode_global = 0 ;
	Myint inode_local  = 0 ;
	Myfloat z_origin_el = z_origin_mesh ;

	for (Myint ielem = 0; ielem < nelem; ielem++)
	{
		// store number of node along one direction
		nnode_z = pElement[ielem].nnode_1D ;

		// allocate local node array
		pElement[ielem].inode = allocate_array<Myint>(pElement[ielem].nnode) ;

		// get coordinates in reference element
		Myint ng = nnode_z ;
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

			// assign new node
			if (new_node_z)
			{
				// new node
				if ((inode_global > nnode) || (inode_global < 0))
				{
					print_error("IN FEM_1D::initialize_mesh --> invalid inode_global", inode_global);
					return(RTN_CODE_KO) ;
				}

				// assign global node index
				pElement[ielem].inode[inode_local] = inode_global ;

				// node coordinates
				pNode[inode_global].zcoord = z_origin_el + (zg[izn] + 1.0) * pElement[ielem].size / 2. ;

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
						// --> this is the last node
						Myint izn_neigh = pElement[ielem_neigh].nnode - 1 ;
						Myint inode_neigh = izn_neigh ;
						pElement[ielem].inode[inode_local] = pElement[ielem_neigh].inode[inode_neigh] ;
					}
				}
			}

			// increment local node index
			inode_local++ ;

		} // for (Myint izn = 0; izn < nnode_z ; izn++)

		// update element origin
		z_origin_el += pElement[ielem].size ;

	} // for (Myint ielem = 0; ielem < nelem; ielem++)

	// check all nodes have been assigned
	if (inode_global != nnode)
	{
		print_error(" Error in FEM_1D::initialize_mesh, inode_global != nnode", inode_global) ;
		return(RTN_CODE_KO) ;
	}

	// assign boundary type for each node
	//-----------------------------------
	for (Myint ielem = 0; ielem < nelem; ielem++)
	{
		for (Myint inode = 0; inode < pElement[ielem].nnode; inode++)
		{
			// $$$ tmp: does not allow mix of boundaries
			// boundary: NO_BOUNDARY, PML, RANDOM, SPG, FREESURF, RIGID, INTERNAL
			// region: NO_REGION, GHOST, LAYER, MEDIUM
			Myint inode_global = pElement[ielem].inode[inode] ;

			// free surface is mandatory at the edges of the mesh
			// if not, scheme is unstable
			if ((pElement[ielem].neigh[I_ZPREV] == NO_NEIGH) && (inode == 0))
			{
				pNode[inode_global].boundary = FREESURF ;
			}
			else if ((pElement[ielem].neigh[I_ZNEXT] == NO_NEIGH) && (inode == pElement[ielem].nnode - 1))
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
		}
	}

	// flux type
	for (Myint ielem = 0; ielem < nelem; ielem++)
	{
		for (Myint ineigh = 0; ineigh < 2; ineigh++)
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
	//----------------------------------------------------------------------------------

	// initialize pVec_vp_glob
	pVec_vp_glob = new Grid_1D_float(nnode, 0.) ;
	Myfloat *vp_glob = pVec_vp_glob->pArray ;

	iel = -1 ;
	for (Myint iz = 0; iz < nelem_z; iz++)
	{
		iel++ ;
		// determine vp
		//if (pElement[iel].region == MEDIUM)
		{
			// retrieve velocity at 1 corners of element
			// and take minimum value
			Myint iz_model = (z_origin_mesh + ((iz + 0) * el_sizez)) / dz_model ;
			iz_model = max(0,  iz_model) ;
			iz_model = min(nz_model-1, iz_model) ;
			Myfloat v1 = vp_1D_grid->pArray[iz_model] ;

			iz_model = (z_origin_mesh + ((iz + 1) * el_sizez)) / dz_model ;
			iz_model = max(0,  iz_model) ;
			iz_model = min(nz_model-1, iz_model) ;
			Myfloat v2 = vp_1D_grid->pArray[iz_model] ;

			// assign vp to each node
			for (Myint inode = 0; inode < pElement[iel].nnode ; inode++)
			{
				// PROP_CONST
				if (prop == PROP_CONST)
				{
					iz_model = (z_origin_mesh + ((iz + 0.5) * el_sizez)) / dz_model ;
					iz_model = max(0,  iz_model) ;
					iz_model = min(nz_model-1, iz_model) ;
					Myfloat v0 = vp_1D_grid->pArray[iz_model] ;
					vp_glob[ pElement[iel].inode[inode] ] = v0 ;
				}
				// PROP_AVER
				else if (prop == PROP_AVER)
				{
					vp_glob[ pElement[iel].inode[inode] ] = (v1 + v2) / 2.0 ;
				}
				// linear properties per element
				// PROP_LINEAR
				else if (prop == PROP_LINEAR)
				{
					Myfloat iz0 = z_origin_mesh + iz * el_sizez ;
					Myfloat alpha_z = (pNode[ pElement[iel].inode[inode] ].zcoord - iz0) / el_sizez ;
					vp_glob[ pElement[iel].inode[inode] ] = (1.0 - alpha_z) * v1 + alpha_z * v2 ;
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
		print_error(" Error in FEM_1D::initialize_mesh, min. vp < 0") ;
		return(RTN_CODE_KO) ;
	}
	if (maxval < 0.0)
	{
		print_error(" Error in FEM_1D::initialize_mesh, max. vp < 0") ;
		return(RTN_CODE_KO) ;
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

		// allocate array
		pSponge_coef_node = new Grid_1D_float(nnode, 0.) ;
		Myfloat* sponge_coef_node = pSponge_coef_node->pArray ;

		// compute sponge coef alpha
		// if coef is set to zero in xml, then constant coef is used in sponges
		Myfloat zLayerBeg_alpha = -1 ;
		if (get_boundary_coef(ZBEG) > 0) zLayerBeg_alpha = sqrt(-log(get_boundary_coef(ZBEG))) ;
		Myfloat zLayerEnd_alpha = -1 ;
		if (get_boundary_coef(ZEND) > 0) zLayerEnd_alpha = sqrt(-log(get_boundary_coef(ZEND))) ;

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
			}
		}
	}

	// save nodes coordinates on disk
	ofstream out_file ;
	out_file.open(NODE_COORD_OUT_FILE, ios::trunc | ios::out) ;
	for (Myint inode = 0; inode < nnode; inode++)
	{
		out_file << pNode[inode].zcoord << "\n" ;
	}
	out_file.close() ;

	// display mesh info
	Rtn_code rtn_code = mesh_info() ;
	if (rtn_code != RTN_CODE_OK) return (rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D::initialize_mesh");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FEM_1D::initialize(Model* pModel, Myfloat fmax)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D::initialize");

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
		pElement[ielem].vnz[I_ZPREV] = -1.0 ;

		pElement[ielem].vnz[I_ZNEXT] = +1.0 ;
	}

	// compute optimal time step
	Myfloat optimal_dt = compute_optimal_time_step() ;

	// set appropriate nt and dt
	if (set_nt_and_dt(optimal_dt) == RTN_CODE_KO) return(RTN_CODE_KO) ;

	// call parent initialization
	rtn_code = FEM::initialize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// save mesh in VTK format
	write_mesh_VTK() ;
	write_node_VTK() ;

	// initialize global vector and matrix
	//===============================================================

	// initialize pVec_k_glob
	pVec_k_glob = new Grid_1D_float(nnode, 0.) ;

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
		Myint const mat_idx_i = pElement[iel].nnode_1D - 1 ;
		Myfloat** const Mat_M_inv = pMat_M_inv[mat_idx_i]->pArray ;
		Myfloat cvol = Myfloat(2.0) / pElement[iel].size ;
		for (Myint ii = 0; ii < ni; ii++) {
			Myint jj = ii ;
			Myint kk = pElement[iel].inode[ii] ;
			Mat_M_inv_glob[kk] += 1.0 / (Mat_M_inv[ii][jj] * cvol) ;
		}
	}

	// compute inverse of Mat_M_inv_glob
	for (Myint inode = 0; inode < nnode; inode++)
	{
		Mat_M_inv_glob[inode] = 1.0 / Mat_M_inv_glob[inode] ;
		print_debug(ALL, MID_DEBUG, " Mat_M_inv_glob ", Mat_M_inv_glob[inode]) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D::initialize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FEM_1D::finalize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D::finalize");

	// deallocate global node array
	deallocate_array<Node_1D_struct_type>(pNode, nnode) ;

	// deallocate local node array
	for (Myint ielem = 0; ielem < nelem; ielem++)
	{
		deallocate_array<Myint>(pElement[ielem].inode, pElement[ielem].nnode) ;
	}

	// deallocate element array
	deallocate_array<Segment_struct_type>(pElement, nelem) ;

	// call parent finalize
	Rtn_code rtn_code = FEM::finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	// deallocate stiffness matrices
	for (Myint imat = 0; imat <= MAX_POLY_ORDER; imat++)
	{
		if (pMat_Dz[imat] != NULL) delete(pMat_Dz[imat]) ;
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D::finalize");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FEM_1D::compute_eigen_error(Myfloat* pr, Myint it, Acquisition* pAcquisition) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D::compute_eigen_error");

	// compute error only at output time steps
	Myint decim = round(dt_out / dt) ;
	if (it%decim != 0) return(RTN_CODE_OK) ;

	Myfloat ttime = dt * it ;

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

				if (abs(pr[inode]) > MAX_EIGEN_VAL)
				{
					print_error(" MAX_EIGEN_VAL HAS BEEN REACHED") ;
					return(RTN_CODE_KO) ;
				}

				Myfloat eigen_sol = -sin(M_PI*znode * eigen_nmode) * sin(M_PI*ttime * eigen_nmode) ;
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
			Myfloat eigen_sol = -sin(M_PI*znode * eigen_nmode) * sin(M_PI*ttime * eigen_nmode) ;

			Myfloat pr_rec = interpolate_variable(pr, ielem_rec[irec], eta_rec[irec]) ;

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

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D::compute_eigen_error");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat FEM_1D::compute_optimal_time_step(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D::compute_optimal_time_step");

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

		// P0 special case
		if (pElement[ielem].nnode_1D == 1)
		{
			Myfloat dist = pElement[ielem].size ;
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
					Myfloat dist = abs(pNode[inode1].zcoord - pNode[inode2].zcoord) ;

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
		Myfloat dt_i = min_dist_node_i / pElement[ielem].vmax ;
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

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D::compute_optimal_time_step");
	return(optim_dt) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FEM_1D::locate_pixel_in_mesh(Snapshot *pSnapshot)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D::locate_pixel_in_mesh");

	print_warning(" locate_pixel_in_mesh not yet implemented") ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D::locate_pixel_in_mesh");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FEM_1D::locate_src_and_rec_in_mesh(Acquisition *acquisition)
{
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D::locate_src_and_rec_in_mesh");

	// retrieve nodes for the source excitation
	//#########################################

	nnode_src = 0 ;
	Myfloat sum_src = 0.0 ;
	Myfloat* const Mat_M_inv_glob = pMat_M_inv_glob->pArray ;

	// case one, NO source smoothing
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
				Myfloat dist =  abs(acquisition->zsrc - pNode[pElement[ielem].inode[ii]].zcoord) ;
				if (dist <= min_dist)
				{
					min_dist = dist ;
					inode_src[0] = pElement[ielem].inode[ii] ;
					wnode_src[0] = Mat_M_inv_glob[inode_src[0]] ;
					sum_src = 1.0 ;
				}
			}
		}

		// $$$ tmp, need to check source is in mesh
		// if (min_dist > max_elem_size)
		// 	{
		// 	  print_warning(" Source is not in domain") ;
		// 	  inode_src[0] = NOT_FOUND ;
		// 	}
		// else
		{
			nnode_src = 1 ;
		}
	}

	// case one, WITH source GAUSSIAN
	// excitation on several node
	//================================
	else if (src_type == SRC_GAUSSIAN)

	{
		// 1st loop to get number of node
		Myfloat radius = src_sigma*3.0 ;

		for (Myint inode = 0; inode < nnode; inode++)
		{
			Myfloat dist =  abs(acquisition->zsrc - pNode[inode].zcoord) ;

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
				Myfloat dist =  abs(acquisition->zsrc - pNode[inode].zcoord) ;
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

	// case one, WITH source SINC
	// excitation on several node
	//================================
	else if (src_type == SRC_SINC)

	{
		Myfloat radius = src_sigma*3.0 ;
		// 1st loop to get number of node
		for (Myint inode = 0; inode < nnode; inode++)
		{
			Myfloat dist =  abs(acquisition->zsrc - pNode[inode].zcoord) ;
			if (dist <= radius)
			{
				nnode_src++ ;
			}
		}

		Bessik bessel = Bessik() ;
		Myfloat bb = 2.94 ;
		Myfloat b1 = bessel.in(0, bb) ;

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
				Myfloat dist =  abs(acquisition->zsrc - pNode[inode].zcoord) ;
				if (dist <= radius)
				{
					inode_src[nnode_src] = inode ;
					Myfloat term = bb * sqrt(1.0  - (dist*dist/(radius*radius)) ) ;
					Myfloat b2 = bessel.in(0, term) ;
					Myfloat kaiser = b2 / b1 ;
					if (dist > 0.001)
					{
						wnode_src[nnode_src] = kaiser * sin(PI*dist) / (PI*dist) ;
					}
					else
					{
						wnode_src[nnode_src] = 1.0 ;
					}
					sum_src += wnode_src[nnode_src] / Mat_M_inv_glob[inode_src[nnode_src]] ;
					cout << "coef. src." << wnode_src[nnode_src] << "\n" ;

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
	Myfloat znode_min, znode_max ;

	// --> exact location with interpolation
	//----------------------------------------

	print_debug(ALL, MID_DEBUG, "receiver located at exact position with interpolation") ;
	ielem_rec = allocate_array<Myint>(nrec) ;
	eta_rec   = allocate_array<Myfloat>(nrec) ;

	// define epsilon
	Myfloat eps_dist = min_dist_node / 1000.0 ;

	// loop on all receivers
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
			}
			// P0 element is special case
			else
			{
				znode_min = pNode[pElement[ielem].inode[0]].zcoord - pElement[ielem].size / 2.0 ;
				znode_max = pNode[pElement[ielem].inode[0]].zcoord + pElement[ielem].size / 2.0 ;
			}

			// check rec is within the extremums
			if ((acquisition->zrec[irec] >= (znode_min-eps_dist)) && (acquisition->zrec[irec] <= (znode_max+eps_dist)))
			{
				ielem_rec[irec] = ielem ;
				eta_rec[irec] = -1.0 + 2.0*(acquisition->zrec[irec] - znode_min)/pElement[ielem].size  ;
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
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D::locate_src_and_rec_in_mesh");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code FEM_1D::write_trace(Variable* pVar, Myint it)
{
	print_debug(ALL, FULL_DEBUG, "IN FEM_1D::write_trace");

	// check if ouput is required at the current time step
	Myint decim = round(dt_out / dt) ;
	string file_name = pVar->get_name() + TIME_REC_OUT_FILE ;

	if (it%decim != 0)
	{
		print_debug(ALL, FULL_DEBUG, "OUT FEM_1D::write_trace");
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
			u_tmp[irec] = interpolate_variable(time_u, ielem_rec[irec], eta_rec[irec]) ;
		}
	}

	// open, write and close
	ofstream pFile ;
	pFile.open(file_name.c_str(), ios::binary | ios::app | ios::out) ;
	assert(pFile.is_open());
	pFile.write((char*)u_tmp, nrec * sizeof(Myfloat)) ;
	pFile.close() ;

	print_debug(ALL, FULL_DEBUG, "OUT FEM_1D::write_trace");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat FEM_1D::interpolate_variable(Myfloat* pVal, Myint ielem, Myfloat eta_coord)
{
	print_debug(ALL, FULL_DEBUG, "IN FEM_1D::interpolate_variable");

	Myfloat return_val = 0.0 ;

	// get number of nodes in element
	Myint ni = pElement[ielem].nnode ;

	// P0 element
	if (ni == 1)
	{
		print_debug(ALL, FULL_DEBUG, "OUT FEM_1D::interpolate_variable");
		return(pVal[pElement[ielem].inode[0]]) ;
	}

	// build vector of nodes
	VecDoub *xnode = &(xnode_ref_elem_1D[ni-1]) ;

	// loop over the nodes of the element
	for (Myint in=0; in < ni; in++)
	{
		Myfloat poly_val = lagran(*xnode, in, eta_coord) ;
		return_val += poly_val * pVal[pElement[ielem].inode[in]] ;
	}

	print_debug(ALL, FULL_DEBUG, "OUT FEM_1D::interpolate_variable");
	return(return_val) ;
}

//-----------------------------------------------------------------------------------------
// compute mass matrix with quadrature rule such as Gauss-Legendre
// in the reference element with xi in [-1, 1]
//
// Mij = Int_elem_vol (phi_i x phi_j) dxi
//
// INPUT
// xg (VecDoub) = array of quadrature coordinates
// wg (VecDoub) = array of quadrature weights
// xnode (VecDoub) = array of node coordinates in element
//
// OUTPUT
// Mat_M (MatDoub) mass matrix
//-----------------------------------------------------------------------------------------
Rtn_code FEM_1D::compute_mass_matrix(VecDoub &xg, VecDoub& wg, VecDoub& xnode, MatDoub &Mat_M)
{   
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D::compute_mass_matrix");

	// # nodes in element
	Myint nn = xnode.size() ;
	if (nn == 0)
	{
		print_error(" Error in FEM_1D::compute_mass_matrix, xnode.size() = 0") ;
		return(RTN_CODE_KO) ;
	}

	// check size of matrix
	if (nn != Mat_M.ncols())
	{
		print_error(" Error in FEM_1D::compute_mass_matrix, nn != Mat_M.ncols()") ;
		return(RTN_CODE_KO) ;
	}
	if (nn != Mat_M.nrows())
	{
		print_error(" Error in FEM_1D::compute_mass_matrix, nn != Mat_M.nrows()") ;
		return(RTN_CODE_KO) ;
	}

	// # GL points
	Myint ng = xg.size() ;
	if (ng == 0)
	{
		print_error(" Error in FEM_1D::compute_mass_matrix, xg.size() = 0") ;
		return(RTN_CODE_KO) ;
	}

	// initialize mass matrix
	for (Myint ii=0; ii< nn; ii++)
	{
		for (Myint jj=0; jj< nn; jj++)
		{
			Mat_M[ii][jj] = 0.0 ;
		}
	}

	// loop on element nodes i
	for (Myint in=0; in<nn; in++)
	{
		// loop on element nodes j
		for (Myint jn=0; jn<nn; jn++)
		{
			// loop on quadrature points
			for (Myint ig=0; ig<ng; ig++)
			{
				// phi_in at xg[ig]
				Doub phi_in = lagran(xnode, in, xg[ig]) ;

				// phi_jn at xg[ig]
				Doub phi_jn = lagran(xnode, jn, xg[ig]) ;

				// update mass matrix
				Mat_M[in][jn] = Mat_M[in][jn] + wg[ig] * phi_in * phi_jn ;
			}
		}
	}
	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D::compute_mass_matrix");
	return(RTN_CODE_OK) ;
}

//-----------------------------------------------------------------------------------------
// compute stiffness matrix with quadrature rule such as Gauss-Legendre
// in the reference element with xi in [-1, 1]
//
// For 1ST ORDER:
// (Ki)jr = Int_elem_i_vol (d_ph_ir/dxi x phi_ij) dv
//
// For 2ND ORDER:
// (Ki)jr = Int_elem_i_vol (d_ph_ir/dxi x d_phi_ij/dxi) dv
//
// INPUT
// xg (float*) = array of quadrature coordinates
// wg (float*) = array of quadrature weights
// xnode (float*) = array of node coordinates in element
// axis (Axis_type) = axis for the derivative of basis function
//
// OUTPUT
// Mat_K (MatDoub) stiffness matrix
//-----------------------------------------------------------------------------------------
Rtn_code FEM_1D::compute_stiffness_matrix(VecDoub &xg, VecDoub& wg, VecDoub& xnode, MatDoub &Mat_K, Axis_type axis)
{   

	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D::compute_stiffness_matrix");

	// only Z_AXIS is allowed in 1D
	if (axis != Z_AXIS)
	{
		print_error(" Error in FEM_1D::compute_stiffness_matrix, axis != Z_AXIS", axis) ;
		return(RTN_CODE_KO) ;
	}

	// # nodes in element
	Myint nn = xnode.size() ;
	if (nn == 0)
	{
		print_error(" Error in FEM_1D::compute_stiffness_matrix, xnode.size() = 0") ;
		return(RTN_CODE_KO) ;
	}

	// check size of matrix
	if (nn != Mat_K.ncols())
	{
		print_error(" Error in FEM_1D::compute_stiffness_matrix, nn != Mat_K.ncols()") ;
		return(RTN_CODE_KO) ;
	}
	if (nn != Mat_K.nrows())
	{
		print_error(" Error in FEM_1D::compute_stiffness_matrix, nn != Mat_K.nrows()") ;
		return(RTN_CODE_KO) ;
	}

	// # GL points
	Myint ng = xg.size() ;
	if (ng == 0)
	{
		print_error(" Error in FEM_1D::compute_stiffness_matrix, xg.size() = 0") ;
		return(RTN_CODE_KO) ;
	}

	// initialize stiffness matrix
	for (Myint ii=0; ii< nn; ii++)
	{
		for (Myint jj=0; jj< nn; jj++)
		{
			Mat_K[ii][jj] = 0.0 ;
		}
	}

	// loop on element nodes i
	for (Myint in=0; in<nn; in++)
	{
		// loop on element nodes j
		for (Myint jn=0; jn<nn; jn++)
		{
			// loop on quadrature points
			for (Myint ig=0; ig<ng; ig++)
			{
				// d_phi_in at xg[ig]
				Doub d_phi_in = dlagran(xnode, in, xg[ig]) ;

				// phi_jn at xg[ig]
				Doub phi_jn ;
				if (eq_order == ORDER_1ST) phi_jn = lagran(xnode, jn, xg[ig]) ;
				if (eq_order == ORDER_2ND) phi_jn = dlagran(xnode, jn, xg[ig]) ;

				// update stiffness matrix
				Mat_K[in][jn] = Mat_K[in][jn] + wg[ig] * d_phi_in * phi_jn ;
			}
		}
	}
	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D::compute_stiffness_matrix");
	return(RTN_CODE_OK) ;
}

//-----------------------------------------------------------------------------------------
Rtn_code FEM_1D::write_mesh_VTK(void)
{   
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D::write_mesh_VTK");

	ofstream out_file ;
	out_file.open(MESH_VTK_OUT_FILE, ios::trunc | ios::out) ;
	out_file << "# vtk DataFile Version 2.0\n" ;
	out_file << "MESH FILE DJANGO\n" ;
	out_file << "ASCII\n" ;
	out_file << "DATASET POLYDATA\n" ;
	out_file << "POINTS " << nnode << " float\n" ;
	for (Myint inode = 0; inode < nnode; inode++)
	{
		out_file << "0.0 0.0 " << pNode[inode].zcoord << "\n" ;
	}
	out_file.close() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D::write_mesh_VTK");
	return(RTN_CODE_OK) ;
}

//-----------------------------------------------------------------------------------------
Rtn_code FEM_1D::write_node_VTK(void)
{   
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D::write_node_VTK");

	ofstream out_file ;
	out_file.open(NODE_VTK_OUT_FILE, ios::trunc | ios::out) ;
	out_file << "# vtk DataFile Version 2.0\n" ;
	out_file << "MESH FILE DJANGO\n" ;
	out_file << "ASCII\n" ;
	out_file << "DATASET POLYDATA\n" ;
	out_file << "POINTS " << nnode << " float\n" ;
	for (Myint inode = 0; inode < nnode; inode++)
	{
		out_file << "0.0 0.0 " << pNode[inode].zcoord << "\n" ;
	}
	out_file.close() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D::write_node_VTK");
	return(RTN_CODE_OK) ;
}

//-----------------------------------------------------------------------------------------
Rtn_code FEM_1D::mesh_info(void)
{   
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D::mesh_info");

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
			if (pElement[ielem].region == MEDIUM)
			{
				min_z_med = min(min_z_med, pNode[pElement[ielem].inode[ii]].zcoord) ;
				max_z_med = max(max_z_med, pNode[pElement[ielem].inode[ii]].zcoord) ;
			}
		}
		min_nnode = min(min_nnode, pElement[ielem].nnode) ;
		max_nnode = max(max_nnode, pElement[ielem].nnode) ;
	}

	// retrieve min and max element size
	Myfloat min_elem_size = FLT_MAX ;
	Myfloat max_elem_size = FLT_MIN ;
	for (Myint ielem = 0; ielem < nelem; ielem++)
	{
		min_elem_size = min(min_elem_size, pElement[ielem].size) ;
		max_elem_size = max(max_elem_size, pElement[ielem].size) ;
	}

	print_info(MASTER, "") ;
	print_info(MASTER, " min z in mesh (m)", min_z_mesh) ;
	print_info(MASTER, " max z in mesh (m)", max_z_mesh) ;
	print_info(MASTER, "") ;
	print_info(MASTER, " min z medium reg. (m)", min_z_med) ;
	print_info(MASTER, " max z medium reg. (m)", max_z_med) ;
	print_info(MASTER, "") ;
	print_info(MASTER, " min elem size (m)", min_elem_size) ;
	print_info(MASTER, " max elem size (m)", max_elem_size) ;

	if (min_elem_size <= 0.0)
	{
		print_error(" Error in FEM_1D::mesh_info, min_elem_size <= 0.0", min_elem_size) ;
		return(RTN_CODE_KO) ;
	}
	if (max_elem_size <= 0.0)
	{
		print_error(" Error in FEM_1D::mesh_info, max_elem_size <= 0.0", max_elem_size) ;
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
		for (Myint ii=0; ii<2; ii++)
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
		for (Myint ineigh = 0; ineigh < 2; ineigh++)
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

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D::mesh_info");
	return(RTN_CODE_OK) ;
}

//-----------------------------------------------------------------------------------------
Rtn_code FEM_1D::compute_element_matrices(void)
{   
	print_debug(ALL, LIGHT_DEBUG, "IN FEM_1D::compute_element_matrices");

	// open output file with matrix layout
	ofstream out_file ;
	out_file.open(FEM_MATRIX_OUT_FILE, ios::trunc | ios::out) ;

	//-----------------------------------------------------------------------------------------
	// loop on supported polynomials orders
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
			print_error(" Error in FEM_1D::compute_element_matrices, numerical integration not supported", node_integ) ;
			return(RTN_CODE_KO) ;
		}

		for (Myint ii=0; ii<xg.size(); ii++)
		{
			out_file << "xg[" << ii << "] " << xg[ii] << "\n" ;
			out_file << "wg[" << ii << "] " << wg[ii] << "\n" ;
		}

		//--------------------------------------
		// compute coordinates of element nodes
		//--------------------------------------
		VecDoub xnode(ni) ;

		//--------------------------------------
		// --> case equidistant node
		//--------------------------------------
		if (node_type == NODE_TYPE_EQD)
		{
			eqd(xnode) ;
		}
		//--------------------------------------
		// --> case GLL node
		//--------------------------------------
		else if (node_type == NODE_TYPE_GLL)
		{
			if (ni == 1)
			{
				xnode[0] = 0.0 ;
			}
			else
			{
				gll(xnode) ;
			}
		}
		else
		{
			print_error(" Error in FEM_1D::compute_element_matrices, node type not supported", node_integ) ;
			return(RTN_CODE_KO) ;
		}

		out_file << "\n*** Element nodes ***\n" ;
		for (Myint ii=0; ii<ni; ii++)
		{
			out_file << "xnode[" << ii << "] " << xnode[ii] << "\n" ;
		}

		//---------------------
		// compute mass matrix
		//---------------------
		MatDoub Mat_M2(ni, ni) ;
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
		pMat_M_inv[imat] = new Grid_2D_float(ni, ni, 1, 1) ;
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

		//----------------------------
		// compute stiffness matrix Dz
		//----------------------------
		MatDoub Mat_K2(ni, ni) ;
		Axis_type axis = Z_AXIS ;
		rtn_code = compute_stiffness_matrix(xg, wg, xnode, Mat_K2, axis) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		out_file << "\n*** Matrix Mat_Dz *** \n" ;
		pMat_Dz[imat] = new Grid_2D_float(ni, ni, 1, 1) ;
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

		//----------------------
		// compute flux matrices
		//----------------------
		rtn_code = compute_flux_matrices(imat, &out_file, ni, xg, wg, xnode) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	} // for (Myint imat=0; imat <= MAX_POLY_ORDER; imat++)

	print_info(MASTER, " Compute FEM matrices", "OK") ;
	print_info(MASTER, " Max polynomial order", MAX_POLY_ORDER) ;
	out_file.close() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT FEM_1D::compute_element_matrices");
	return(RTN_CODE_OK) ;
}

} // namespace django
