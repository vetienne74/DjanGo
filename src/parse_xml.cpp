//-------------------------------------------------------------------------------------------------------
//
// PARSE XML FILE
// AND CREATE APPROPRIATE OBJECTS
//
//-------------------------------------------------------------------------------------------------------

#include "parse_xml.h"

#include <iostream>

#include "FEM.h"
#include "output_report.h"
#include "program.h"
#include "program_fwi.h"
#include "program_guitar.h"
#include "program_modelling.h"
#include "program_rtm.h"
#include "singleton.h"
#include "scheme_factory.h"
#include "type_def.h"

namespace django {

//-------------------------------------------------------------------------------------------------------

void XMLCALL django_startElement(void *userData, const char *name, const char **attr)
{

	// display current element and its attributes
	print_debug(ALL, LIGHT_DEBUG, " django_startElement", name) ;
	for (Myint i = 0; attr[i]; i += 2)
		print_debug(ALL, LIGHT_DEBUG, attr[i], attr[i + 1]);


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// element django
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if (XML_DJANGO.compare(name) == 0) {

		Myfloat xml_version = 0.0 ;

		for (Myint i = 0; attr[i]; i += 2) {

			if (XML_DJANGO_VERSION.compare(attr[i]) == 0) {
				xml_version = atof(attr[i+1]) ;
			}
		}
		Singleton::Instance()->xml_version = xml_version ;
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// element program
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if (XML_PROG.compare(name) == 0) {

		Prog_type type ;
		for (Myint i = 0; attr[i]; i += 2) {

			//###################################################################################################
			// instantiate program
			//###################################################################################################
			if (XML_PROG_TYPE.compare(attr[i]) == 0) {
				if (XML_PROG_TYPE_MOD.compare(attr[i+1]) == 0)
				{
					type = MODELLING_PROG ;
					Singleton::Instance()->pProgram = new Program_Modelling(type) ;
				}
				else if (XML_PROG_TYPE_FWI.compare(attr[i+1]) == 0)
				{
					type = FWI_PROG ;
					Singleton::Instance()->pProgram = new Program_Fwi(type) ;
				}
				else if (XML_PROG_TYPE_RTM.compare(attr[i+1]) == 0)
				{
					type = RTM_PROG ;
					Singleton::Instance()->pProgram = new Program_Rtm(type) ;
				}
				else if (XML_PROG_TYPE_GUITAR.compare(attr[i+1]) == 0)
				{
					type = GUITAR_PROG ;
					Singleton::Instance()->pProgram = new Program_Guitar(type) ;
				}
			}
		}
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// element domain
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if (XML_DOMAIN.compare(name) == 0) {

		//###################################################################################################
		// instantiate domain
		//###################################################################################################
		Singleton::Instance()->pProgram->pDomain = new Domain() ;
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// element model
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if (XML_MODEL.compare(name) == 0) {

		Space_dim dim           = NO_DIM ;
		Model_type type         = NO_MODEL_TYPE ;

		for (Myint i = 0; attr[i]; i += 2) {

			// space dimension
			if (XML_MODEL_DIM.compare(attr[i]) == 0) {
				if (XML_MODEL_DIM_1D.compare(attr[i+1]) == 0)
				{
					dim = ONE ;
				}
				else if (XML_MODEL_DIM_2D.compare(attr[i+1]) == 0)
				{
					dim = TWO ;
				}
				else if (XML_MODEL_DIM_3D.compare(attr[i+1]) == 0)
				{
					dim = THREE ;
				}
			}

			// type
			if (XML_MODEL_TYPE.compare(attr[i]) == 0) {
				if (XML_MODEL_TYPE_GRID.compare(attr[i+1]) == 0)
				{
					type = GRID ;
				}
			}
		}

		//###################################################################################################
		// instantiate model
		//###################################################################################################
		if (Singleton::Instance()->pProgram != NULL)
		{
			if (Singleton::Instance()->pProgram->pDomain != NULL)
			{
				Singleton::Instance()->pProgram->pDomain->pModel = new Model(dim, type) ;
			}
		}
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// element scheme
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if (XML_SCHEME.compare(name) == 0) {

		//###################################################################################################
		// instantiate scheme factory
		//###################################################################################################
		if (Singleton::Instance()->pProgram != NULL)
		{
			if (Singleton::Instance()->pProgram->pDomain != NULL)
			{
				Singleton::Instance()->pProgram->pDomain->pScheme_Factory = new Scheme_Factory() ;
				Singleton::Instance()->pProgram->pDomain->pScheme_Factory->set_model(Singleton::Instance()->pProgram->pDomain->pModel) ;
			}
		}

		//###################################################################################################
		// instantiate scheme
		//###################################################################################################

		Scheme_method method = NO_SCHEME_METHOD ;
		Scheme_type   type   = NO_SCHEME_TYPE ;

		for (Myint i = 0; attr[i]; i += 2) {
			if (XML_SCHEME_METH.compare(attr[i]) == 0) {

				if (XML_SCHEME_METH_FDM.compare(attr[i+1]) == 0)
				{
					method = SCHEME_FDM ;
				}
				else if (XML_SCHEME_METH_FEM.compare(attr[i+1]) == 0)
				{
					method = SCHEME_FEM ;
				}
			}
			if (XML_SCHEME_TYPE.compare(attr[i]) == 0) {

				if (XML_SCHEME_TYPE_STG.compare(attr[i+1]) == 0)
				{
					type = SCHEME_STAGGERED ;
				}
				else if (XML_SCHEME_TYPE_CGM.compare(attr[i+1]) == 0)
				{
					type = SCHEME_CGM ;
				}
				else if (XML_SCHEME_TYPE_DGM.compare(attr[i+1]) == 0)
				{
					type = SCHEME_DGM ;
				}
				else if (XML_SCHEME_TYPE_MGM.compare(attr[i+1]) == 0)
				{
					type = SCHEME_MGM ;
				}
			}
		}
		Singleton::Instance()->method = method ;
		Singleton::Instance()->type   = type ;
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// element equation
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (XML_EQ.compare(name) == 0) {

		Eq_type  type  = NO_EQ_TYPE ;
		Eq_order order = NO_EQ_ORDER ;

		for (Myint i = 0; attr[i]; i += 2) {
			// order
			if (XML_EQ_ORD.compare(attr[i]) == 0) {
				if (XML_EQ_ORD_1.compare(attr[i+1]) == 0)
				{
					order = ORDER_1ST ;
				}
				else if (XML_EQ_ORD_2.compare(attr[i+1]) == 0)
				{
					order = ORDER_2ND ;
				}
			}

			// type
			if (XML_EQ_TYPE.compare(attr[i]) == 0) {
				if (XML_EQ_TYPE_AC.compare(attr[i+1]) == 0)
				{
					type = ACOUSTIC ;
				}
				else if (XML_EQ_TYPE_EL.compare(attr[i+1]) == 0)
				{
					type = ELASTIC ;
				}
				else if (XML_EQ_TYPE_AC_LOSSY.compare(attr[i+1]) == 0)
				{
					type = AC_LOSSY ;
				}
			}
		}

		Singleton::Instance()->eq_type  = type ;
		Singleton::Instance()->eq_order = order ;

	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// element accuracy
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (XML_ACCURACY.compare(name) == 0) {

		Myint space_order = 0 ;
		Myint time_order  = 0;
		Myfloat dt        = 0.0 ;
		Myfloat ratio_cfl = 1.0 ;

		for (Myint i = 0; attr[i]; i += 2) {

			// space order
			if (XML_ACCURACY_SPACE.compare(attr[i]) == 0) {
				space_order = atoi(attr[i+1]) ;
			}
			// time order
			if (XML_ACCURACY_TIME.compare(attr[i]) == 0) {
				time_order = atoi(attr[i+1]) ;
			}
			// dt
			if (XML_ACCURACY_TIME_DT.compare(attr[i]) == 0) {
				dt = atof(attr[i+1]) ;
			}
			// ratio cfl
			if (XML_ACCURACY_RATIO_CFL.compare(attr[i]) == 0) {
				ratio_cfl = atof(attr[i+1]) ;
			}
		}

		Singleton::Instance()->space_order = space_order ;
		Singleton::Instance()->time_order  = time_order ;
		Singleton::Instance()->dt          = dt ;
		Singleton::Instance()->ratio_cfl   = ratio_cfl ;

	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// element node
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (XML_NODE.compare(name) == 0) {

		Node_type type = NO_NODE_TYPE ;
		Node_distrib distrib = NO_NODE_DISTRIB ;
		Node_integ integ = NO_NODE_INTEG ;

		for (Myint i = 0; attr[i]; i += 2) {

			// node type
			if (XML_NODE_TYPE.compare(attr[i]) == 0) {
				if (XML_NODE_TYPE_EQD.compare(attr[i+1]) == 0)
				{
					type = NODE_TYPE_EQD ;
				}
				else if (XML_NODE_TYPE_GLL.compare(attr[i+1]) == 0)
				{
					type = NODE_TYPE_GLL ;
				}
			}

			// node distribution
			if (XML_NODE_DIST.compare(attr[i]) == 0) {
				if (XML_NODE_DIST_UNIFORM.compare(attr[i+1]) == 0)
				{
					distrib = NODE_DISTRIB_UNIFORM ;
				}
				else if (XML_NODE_DIST_RANDOM.compare(attr[i+1]) == 0)
				{
					distrib = NODE_DISTRIB_RANDOM ;
				}
				else if (XML_NODE_DIST_MODEL.compare(attr[i+1]) == 0)
				{
					distrib = NODE_DISTRIB_MODEL ;
				}
			}

			// numerical integration
			if (XML_NODE_INTEG.compare(attr[i]) == 0) {
				if (XML_NODE_INTEG_GL.compare(attr[i+1]) == 0)
				{
					integ = NODE_INTEG_GL ;
				}
				else if (XML_NODE_INTEG_GLL.compare(attr[i+1]) == 0)
				{
					integ = NODE_INTEG_GLL ;
				}
			}
		}

		Singleton::Instance()->node_type    = type ;
		Singleton::Instance()->node_distrib = distrib ;
		Singleton::Instance()->node_integ   = integ ;
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// element mesh
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (XML_MESH.compare(name) == 0) {

		Myint nelem_x  = 0 ;
		Myint nelem_z  = 0 ;
		Prop_type prop = PROP_CONST ;

		for (Myint i = 0; attr[i]; i += 2) {

			// number of element along x-axis
			if (XML_MESH_NELEM_X.compare(attr[i]) == 0)
			{
				nelem_x = atoi(attr[i+1]) ;
			}
			// number of element along z-axis
			else if (XML_MESH_NELEM_Z.compare(attr[i]) == 0)
			{
				nelem_z = atoi(attr[i+1]) ;
			}
			// properties within elements
			else if (XML_MESH_PROP.compare(attr[i]) == 0)
			{
				if (XML_MESH_PROP_CONST.compare(attr[i+1]) == 0)
				{
					prop = PROP_CONST ;
				}
				else if (XML_MESH_PROP_AVER.compare(attr[i+1]) == 0)
				{
					prop = PROP_AVER ;
				}
				else if (XML_MESH_PROP_LINEAR.compare(attr[i+1]) == 0)
				{
					prop = PROP_LINEAR ;
				}
			}
		}

		Singleton::Instance()->nelem_x_med = nelem_x ;
		Singleton::Instance()->nelem_z_med = nelem_z ;
		Singleton::Instance()->prop    = prop ;
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// element flux
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (XML_FLUX.compare(name) == 0) {

		Flux_type type = NO_FLUX_TYPE ;

		for (Myint i = 0; attr[i]; i += 2) {

			// flux type
			if (XML_FLUX_TYPE.compare(attr[i]) == 0) {
				if (XML_FLUX_CENTER.compare(attr[i+1]) == 0)
				{
					type = FLUX_TYPE_CENTERED ;
				}
				else if (XML_FLUX_UPWIND.compare(attr[i+1]) == 0)
				{
					type = FLUX_TYPE_UPWIND ;
				}
			}
		}

		Singleton::Instance()->flux_type = type ;

	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// element dynamic
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (XML_DYN.compare(name) == 0) {

		Front_type front_type = FRONT_STATIC ;

		for (Myint i = 0; attr[i]; i += 2) {

			// Front type
			if (XML_DYN_FRONT.compare(attr[i]) == 0) {
				if (XML_DYN_FRONT_STATIC.compare(attr[i+1]) == 0)
				{
					front_type = FRONT_STATIC ;
				}
				else if (XML_DYN_FRONT_VMAX.compare(attr[i+1]) == 0)
				{
					front_type = FRONT_DYN_VMAX ;
				}
			}
		}

		Singleton::Instance()->front_type = front_type ;
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// element adaptivity
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (XML_SCHEME_ADAPT.compare(name) == 0) {

		Myint   pmin         = 0 ;
		Myint   pmax         = MAX_POLY_ORDER ;
		Myfloat fmax0        = 0.0 ;
		Adapt_type adaptType = NO_ADAPT ;

		for (Myint i = 0; attr[i]; i += 2) {

			// max min allowed polynomial
			if (XML_SCHEME_ADAPT_PMIN.compare(attr[i]) == 0)
			{
				pmin = atoi(attr[i+1]) ;
			}
			// max allowed polynomial
			else if (XML_SCHEME_ADAPT_PMAX.compare(attr[i]) == 0)
			{
				pmax = atoi(attr[i+1]) ;
			}
			else if (XML_SCHEME_ADAPT_FMAX0.compare(attr[i]) == 0)
			{
				fmax0 = atof(attr[i+1]) ;
			}
			else if (XML_SCHEME_ADAPT_TYPE.compare(attr[i]) == 0)
			{
				if (XML_SCHEME_ADAPT_TYPE_NONE.compare(attr[i+1]) == 0)
				{
					adaptType = NO_ADAPT ;
				}
				else if (XML_SCHEME_ADAPT_TYPE_STATIC.compare(attr[i+1]) == 0)
				{
					adaptType = ADAPT_STATIC ;
				}
				else if (XML_SCHEME_ADAPT_TYPE_FREQTIME.compare(attr[i+1]) == 0)
				{
					adaptType = ADAPT_FREQTIME ;
				}
			}
		}

		Singleton::Instance()->adaptType = adaptType ;
		Singleton::Instance()->pmin      = pmin ;
		Singleton::Instance()->pmax      = pmax ;
		Singleton::Instance()->fmax0     = fmax0 ;
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// element freqTime
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (XML_SCHEME_ADAPT_FREQTIME.compare(name) == 0) {

		Myfloat tmin = 0.0 ;
		Myfloat fmax = 0.0 ;

		for (Myint i = 0; attr[i]; i += 2) {

			// tmin
			if (XML_SCHEME_ADAPT_FREQTIME_TMIN.compare(attr[i]) == 0)
			{
				tmin = atof(attr[i+1]) ;
			}
			// fmax
			else if (XML_SCHEME_ADAPT_FREQTIME_FMAX.compare(attr[i]) == 0)
			{
				fmax = atof(attr[i+1]) ;
			}
		}

		Singleton::Instance()->adaptTmin.push_back(tmin) ;
		Singleton::Instance()->adaptFmax.push_back(fmax) ;
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// element acquisition
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (XML_ACQUI.compare(name) == 0) {

		string acqui_file = "UNSPECIFIED" ;

		for (Myint i = 0; attr[i]; i += 2) {

			// file
			if (XML_ACQUI_FILE.compare(attr[i]) == 0) {
				acqui_file = attr[i + 1] ;
			}
		}

		if (Singleton::Instance()->pProgram != NULL)
		{
			if (Singleton::Instance()->pProgram->pDomain != NULL)
			{
				if (Singleton::Instance()->pProgram->pModelling != NULL)
				{
					Singleton::Instance()->pProgram->pModelling->set_acquisition(acqui_file) ;
				}
			}
		}
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// element time
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (XML_TIME.compare(name) == 0) {

		Myfloat dt = 0.0 ;

		for (Myint i = 0; attr[i]; i += 2) {

			// dt
			if (XML_TIME_DT.compare(attr[i]) == 0) {
				dt = atof(attr[i+1]) ;
			}
		}

		//Singleton::Instance()->pModelling->set_modelling_dt(dt) ;
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// boundary
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (XML_BOUND.compare(name) == 0) {

		Edge_type edge         = NO_EDGE ;
		Boundary_type boundary = NO_BOUNDARY ;
		Myint width            = 0 ;
		Myfloat coef           = 0.0 ;

		for (Myint i = 0; attr[i]; i += 2) {

			// edge
			if (XML_BOUND_EDGE.compare(attr[i]) == 0) {
				if (XML_BOUND_EDGE_ZBEG.compare(attr[i+1]) == 0)
				{
					edge = ZBEG ;
				}
				else if (XML_BOUND_EDGE_ZEND.compare(attr[i+1]) == 0)
				{
					edge = ZEND ;
				}
				else if (XML_BOUND_EDGE_XBEG.compare(attr[i+1]) == 0)
				{
					edge = XBEG ;
				}
				else if (XML_BOUND_EDGE_XEND.compare(attr[i+1]) == 0)
				{
					edge = XEND ;
				}
				else if (XML_BOUND_EDGE_YBEG.compare(attr[i+1]) == 0)
				{
					edge = YBEG ;
				}
				else if (XML_BOUND_EDGE_YEND.compare(attr[i+1]) == 0)
				{
					edge = YEND ;
				}
			}

			// type
			if (XML_BOUND_TYPE.compare(attr[i]) == 0) {
				if (XML_BOUND_TYPE_CPML.compare(attr[i+1]) == 0)
				{
					boundary = PML ;
				}
				else if (XML_BOUND_TYPE_SPG.compare(attr[i+1]) == 0)
				{
					boundary = SPG ;
				}
				else if (XML_BOUND_TYPE_SPG2.compare(attr[i+1]) == 0)
				{
					boundary = SPG2 ;
				}
				else if (XML_BOUND_TYPE_RAND.compare(attr[i+1]) == 0)
				{
					boundary = RANDOM ;
				}
				else if (XML_BOUND_TYPE_FS.compare(attr[i+1]) == 0)
				{
					boundary = FREESURF ;
				}
				else if (XML_BOUND_TYPE_RIG.compare(attr[i+1]) == 0)
				{
					boundary = RIGID ;
				}
			}

			// width
			if (XML_BOUND_WIDTH.compare(attr[i]) == 0) {
				width = atoi(attr[i+1]) ;
			}

			// coef
			if (XML_BOUND_COEF.compare(attr[i]) == 0) {
				coef = atof(attr[i+1]) ;
			}

		}

		Singleton::Instance()->create_boundary(edge, boundary, width, coef) ;
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// frequency
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (XML_FREQ.compare(name) == 0) {

		Myfloat fmin = 0.0 ;
		Myfloat fmax = 0.0 ;
		Myfloat df   = 0.0 ;

		for (Myint i = 0; attr[i]; i += 2) {

			// fmin
			if (XML_FREQ_FMIN.compare(attr[i]) == 0) {
				fmin = atof(attr[i+1]) ;
			}

			// fmax
			if (XML_FREQ_FMAX.compare(attr[i]) == 0) {
				fmax = atof(attr[i+1]) ;
			}

			// df
			if (XML_FREQ_DF.compare(attr[i]) == 0) {
				df = atof(attr[i+1]) ;
			}
		}
		if (Singleton::Instance()->pProgram != NULL)
		{
			if (Singleton::Instance()->pProgram->pModelling != NULL)
			{
				Singleton::Instance()->pProgram->pModelling->set_freq_param(fmin, fmax, df) ;
			}
		}
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// snapshot
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (XML_SNAPSHOT.compare(name) == 0) {

		Myfloat tmin ;
		Myfloat tmax ;
		Myfloat dt ;

		for (Myint i = 0; attr[i]; i += 2) {

			// tmin
			if (XML_SNAPSHOT_TMIN.compare(attr[i]) == 0) {
				tmin = atof(attr[i+1]) ;
			}

			// tmax
			if (XML_SNAPSHOT_TMAX.compare(attr[i]) == 0) {
				tmax = atof(attr[i+1]) ;
			}

			// dt
			if (XML_SNAPSHOT_DT.compare(attr[i]) == 0) {
				dt = atof(attr[i+1]) ;
			}
		}

		if (Singleton::Instance()->pProgram != NULL)
		{
			if (Singleton::Instance()->pProgram->pModelling != NULL)
			{
				Singleton::Instance()->pProgram->pModelling->set_snapshot_time_param(tmin, tmax, dt) ;
			}
		}

	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// pixel
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (XML_PIXEL.compare(name) == 0) {

		Myfloat xmin = NOT_SPECIFIED ;
		Myfloat xmax = NOT_SPECIFIED ;
		Myint nx     = NOT_SPECIFIED ;
		Myfloat zmin = NOT_SPECIFIED ;
		Myfloat zmax = NOT_SPECIFIED ;
		Myint nz     = NOT_SPECIFIED ;

		for (Myint i = 0; attr[i]; i += 2) {

			// xmin
			if (XML_PIXEL_XMIN.compare(attr[i]) == 0) {
				xmin = atof(attr[i+1]) ;
			}

			// xmax
			if (XML_PIXEL_XMAX.compare(attr[i]) == 0) {
				xmax = atof(attr[i+1]) ;
			}

			// nx
			if (XML_PIXEL_NX.compare(attr[i]) == 0) {
				nx = atoi(attr[i+1]) ;
			}

			// zmin
			if (XML_PIXEL_ZMIN.compare(attr[i]) == 0) {
				zmin = atof(attr[i+1]) ;
			}

			// zmax
			if (XML_PIXEL_ZMAX.compare(attr[i]) == 0) {
				zmax = atof(attr[i+1]) ;
			}

			// nz
			if (XML_PIXEL_NZ.compare(attr[i]) == 0) {
				nz = atoi(attr[i+1]) ;
			}
		}

		if (Singleton::Instance()->pProgram != NULL)
		{
			if (Singleton::Instance()->pProgram->pModelling != NULL)
			{
				Singleton::Instance()->pProgram->pModelling->set_snapshot_pixel_param(xmin, xmax, nx,
						zmin, zmax, nz) ;
			}
		}

	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// time
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (XML_TIME.compare(name) == 0) {

		Myfloat tmax      = 0.0 ;
		Myfloat dt        = 0.0 ;
		Myint rate        = 0 ;
		bool energy_flag  = false ;

		for (Myint i = 0; attr[i]; i += 2) {

			// tmax
			if (XML_TIME_TMAX.compare(attr[i]) == 0) {
				tmax = atof(attr[i+1]) ;
			}

			// dt
			if (XML_TIME_DT.compare(attr[i]) == 0) {
				dt = atof(attr[i+1]) ;
			}

			// rate
			if (XML_TIME_RATE.compare(attr[i]) == 0) {
				rate = atoi(attr[i+1]) ;
				dt = Myfloat(1.) / Myfloat(rate) ;
			}

			// energy
			if (XML_OUTPUT_ENERGY.compare(attr[i]) == 0) {
				if (XML_OUTPUT_ENERGY_ON.compare(attr[i+1]) == 0)
				{
					energy_flag = true ;
				}
				else if (XML_OUTPUT_ENERGY_OFF.compare(attr[i+1]) == 0)
				{
					energy_flag = false ;
				}
			}

			// frequency
			if (XML_OUTPUT_FREQ.compare(attr[i]) == 0) {
				if (XML_OUTPUT_FREQ_ON.compare(attr[i+1]) == 0)
				{
					//frequency_flag = YES ;
				}
				else if (XML_OUTPUT_FREQ_OFF.compare(attr[i+1]) == 0)
				{
					//frequency_flag = NO ;
				}
			}
		}

		if (Singleton::Instance()->pProgram != NULL)
		{
			if (Singleton::Instance()->pProgram->pModelling != NULL)
			{
				Singleton::Instance()->pProgram->pModelling->set_time_param(tmax, dt) ;
				Singleton::Instance()->pProgram->pModelling->set_energy_param(energy_flag) ;
			}
		}
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// source
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (XML_SOURCE.compare(name) == 0) {

		Myfloat freq        = 0 ;
		Myfloat src_t0      = 0 ;
		Src_func src_func   = NO_SRC_FUNC ;
		Src_type src_type   = NO_SRC_TYPE ;
		Src_stype src_stype = NO_SRC_STYPE ;
		Myfloat src_dt      = 0.0 ;
		string src_file     = UNSPECIFIED ;
		Myfloat src_sigma   = 0.0 ;

		for (Myint i = 0; attr[i]; i += 2) {

			// sigma
			if (XML_SOURCE_SIGMA.compare(attr[i]) == 0) {
				src_sigma = atof(attr[i+1]) ;
			}

			// freq
			if (XML_SOURCE_FREQ.compare(attr[i]) == 0) {
				freq = atof(attr[i+1]) ;
			}

			// function
			if (XML_SOURCE_FUNC.compare(attr[i]) == 0) {
				if (XML_SOURCE_FUNC_RIC.compare(attr[i+1]) == 0)
				{
					src_func = RICKER ;
				}
				else if (XML_SOURCE_FUNC_RIC_P.compare(attr[i+1]) == 0)
				{
					src_func = RICKER_P ;
				}
				else if (XML_SOURCE_FUNC_RIC_PP.compare(attr[i+1]) == 0)
				{
					src_func = RICKER_PP ;
				}
				else if (XML_SOURCE_FUNC_RIC_D.compare(attr[i+1]) == 0)
				{
					src_func = RICKER_D ;
				}
				else if (XML_SOURCE_FUNC_MONO.compare(attr[i+1]) == 0)
				{
					src_func = MONO_FREQ ;
				}
				else if (XML_SOURCE_FUNC_FILE.compare(attr[i+1]) == 0)
				{
					src_func = SRC_FILE ;
				}
			}

			// type
			if (XML_SOURCE_TYPE.compare(attr[i]) == 0) {
				if (XML_SOURCE_TYPE_POINT.compare(attr[i+1]) == 0)
				{
					src_type = SRC_POINT ;
				}
				else if (XML_SOURCE_TYPE_GAUSS.compare(attr[i+1]) == 0)
				{
					src_type = SRC_GAUSSIAN ;
				}
				else if (XML_SOURCE_TYPE_SINC.compare(attr[i+1]) == 0)
				{
					src_type = SRC_SINC ;
				}
			}

			// subtype
			if (XML_SOURCE_STYPE.compare(attr[i]) == 0) {
				if (XML_SOURCE_STYPE_EXP.compare(attr[i+1]) == 0)
				{
					src_stype = EXPLOSIVE ;
				}
				else if (XML_SOURCE_STYPE_FZ.compare(attr[i+1]) == 0)
				{
					src_stype = FORCE_Z ;
				}
			}

			// file
			if (XML_SOURCE_FUNC_NAME.compare(attr[i]) == 0) {
				src_file = attr[i + 1] ;
			}

			// dt
			if (XML_SOURCE_FUNC_DT.compare(attr[i]) == 0) {
				src_dt = atof(attr[i+1]) ;
			}

			// t0
			if (XML_SOURCE_FUNC_T0.compare(attr[i]) == 0) {
				src_t0 = atof(attr[i+1]) ;
			}

		}
		if (Singleton::Instance()->pProgram != NULL)
		{
			if (Singleton::Instance()->pProgram->pModelling != NULL)
			{
				Singleton::Instance()->pProgram->pModelling->set_src_param(freq, src_func,
						src_type, src_stype,
						src_file, src_dt,
						src_t0, src_sigma) ;
			}
		}
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// model size
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (XML_SIZE.compare(name) == 0) {

		Myint nx = 0 ;
		Myint ny = 0 ;
		Myint nz = 0 ;

		for (Myint i = 0; attr[i]; i += 2) {

			// nx
			if (XML_SIZE_NX.compare(attr[i]) == 0) {
				nx = atoi(attr[i+1]) ;
			}

			// ny
			else if (XML_SIZE_NY.compare(attr[i]) == 0) {
				ny = atoi(attr[i+1]) ;
			}

			// nz
			else if (XML_SIZE_NZ.compare(attr[i]) == 0) {
				nz = atoi(attr[i+1]) ;
			}

			if (Singleton::Instance()->pProgram != NULL)
			{
				if (Singleton::Instance()->pProgram->pDomain != NULL)
				{
					if (Singleton::Instance()->pProgram->pDomain->pModel != NULL)
					{
						Singleton::Instance()->pProgram->pDomain->pModel->set_size(nx, ny, nz) ;
					}
				}
			}
		}
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// sampling
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (XML_SAMPLING.compare(name) == 0) {

		Myfloat dx = 0.0 ;
		Myfloat dy = 0.0 ;
		Myfloat dz = 0.0 ;
		Myfloat dxrandom = 0.0 ;
		Myfloat dyrandom = 0.0 ;
		Myfloat dzrandom = 0.0 ;

		for (Myint i = 0; attr[i]; i += 2) {

			// dx
			if (XML_SAMPLING_DX.compare(attr[i]) == 0) {
				dx = atof(attr[i+1]) ;
			}

			// dy
			if (XML_SAMPLING_DY.compare(attr[i]) == 0) {
				dy = atof(attr[i+1]) ;
			}

			// dz
			if (XML_SAMPLING_DZ.compare(attr[i]) == 0) {
				dz = atof(attr[i+1]) ;
			}

			// dxrandom
			if (XML_SAMPLING_DX_RAND.compare(attr[i]) == 0) {
				dxrandom = atof(attr[i+1]) ;
			}

			// dyrandom
			if (XML_SAMPLING_DY_RAND.compare(attr[i]) == 0) {
				dyrandom = atof(attr[i+1]) ;
			}

			// dyrandom
			if (XML_SAMPLING_DZ_RAND.compare(attr[i]) == 0) {
				dzrandom = atof(attr[i+1]) ;
			}

			if (Singleton::Instance()->pProgram != NULL)
			{
				if (Singleton::Instance()->pProgram->pDomain != NULL)
				{
					if (Singleton::Instance()->pProgram->pDomain->pModel != NULL)
					{
						Singleton::Instance()->pProgram->pDomain->pModel->set_sampling(dx, dy, dz) ;
						Singleton::Instance()->pProgram->pDomain->pModel->set_sampling_random(dxrandom, dyrandom, dzrandom) ;
					}
				}
			}
		}
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// param
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (XML_PARAM.compare(name) == 0) {

		string file_name         = UNSPECIFIED ;
		string param_name        = UNSPECIFIED ;
		Var_type type            = NO_VAR_TYPE ;
		Myfloat val              = 0.0 ;
		Self_init_type self_init = NO_SELF_INIT ;

		for (Myint i = 0; attr[i]; i += 2) {

			// type
			if (XML_PARAM_TYPE.compare(attr[i]) == 0) {
				if (XML_PARAM_TYPE_VP.compare(attr[i+1]) == 0) {
					type = VP ;
				}
				else if (XML_PARAM_TYPE_VS.compare(attr[i+1]) == 0) {
					type = VS ;
				}
				if (XML_PARAM_TYPE_RHO.compare(attr[i+1]) == 0) {
					type = RHO ;
				}
				if (XML_PARAM_TYPE_LOSS1.compare(attr[i+1]) == 0) {
					type = LOSS1 ;
				}
				if (XML_PARAM_TYPE_LOSS2.compare(attr[i+1]) == 0) {
					type = LOSS2 ;
				}
				param_name = attr[i+1] ;
			}

			// file name
			if (XML_PARAM_FILE.compare(attr[i]) == 0) {
				file_name = attr[i+1] ;
				self_init = INIT_FROM_FILE ;
			}

			// file name
			if (XML_PARAM_CONST.compare(attr[i]) == 0) {
				val = atof(attr[i+1]) ;
				self_init = INIT_FROM_CONST ;
			}
		}

		//###################################################################################################
		// create parameter
		//###################################################################################################
		if (Singleton::Instance()->pProgram != NULL)
		{
			if (Singleton::Instance()->pProgram->pDomain != NULL)
			{
				if (Singleton::Instance()->pProgram->pDomain->pModel != NULL)
				{

					Variable* pVar = Singleton::Instance()->pProgram->pDomain->pModel->register_parameter(type, param_name) ;
					pVar->set_self_init_mode(self_init) ;
					if (self_init == INIT_FROM_CONST)
					{
						pVar->set_self_init_const(val) ;
					}
					else if (self_init == INIT_FROM_FILE)
					{
						pVar->set_self_init_file(file_name) ;
					}
				}
			}
		}
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// model
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (XML_INV.compare(name) == 0) {

		Myint niter = 0 ;
		Myint ntry  = 0 ;
		Myfloat init_try = 0.0 ;

		for (Myint i = 0; attr[i]; i += 2) {

			// niter
			if (XML_INV_NITER.compare(attr[i]) == 0) {
				niter = atoi(attr[i+1]) ;
			}
			// ntry
			if (XML_INV_NTRY.compare(attr[i]) == 0) {
				ntry = atoi(attr[i+1]) ;
			}
			// init_try
			if (XML_INV_INIT_TRY.compare(attr[i]) == 0) {
				init_try = atof(attr[i+1]) ;
			}
		}

		//Singleton::Instance()->pInversion->set_niter(niter) ;
		//Singleton::Instance()->pInversion->set_ntry(ntry) ;
		//Singleton::Instance()->pInversion->set_init_try(init_try) ;
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// domain
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (XML_DOMAIN.compare(name) == 0) {

		Myfloat min = 0.0 ;
		Myfloat max = 0.0 ;
		Myfloat delta = 0.0 ;
		Domain_type domain_type = NO_DOMAIN ;

		for (Myint i = 0; attr[i]; i += 2) {

			// type
			if (XML_DOMAIN_TYPE.compare(attr[i]) == 0) {
				if (XML_DOMAIN_TYPE_FREQ.compare(attr[i+1]) == 0) {
					domain_type = FREQ ;
				}
			}

			// min
			if (XML_DOMAIN_MIN.compare(attr[i]) == 0) {
				min = atof(attr[i+1]) ;
			}
			// max
			if (XML_DOMAIN_MAX.compare(attr[i]) == 0) {
				max = atof(attr[i+1]) ;
			}
			// delta
			if (XML_DOMAIN_DELTA.compare(attr[i]) == 0) {
				delta = atof(attr[i+1]) ;
			}
		}

		//if (domain_type == FREQ) Singleton::Instance()->pInversion->set_freq_inv(min, max, delta) ;
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// modelling
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (XML_MODELLING.compare(name) == 0) {

		// create a modelling object
		if (Singleton::Instance()->pProgram != NULL)
		{
			Singleton::Instance()->pProgram->pModelling = new Modelling() ;
		}

		for (Myint i = 0; attr[i]; i += 2) {

			// case
			if (XML_MODELLING_CASE.compare(attr[i]) == 0) {
				if (XML_MODELLING_CASE_STD.compare(attr[i+1]) == 0) {
					Singleton::Instance()->pProgram->pModelling->set_case(MODELLING_STD) ;
				}
				else if (XML_MODELLING_CASE_EIG.compare(attr[i+1]) == 0)
				{
					Singleton::Instance()->pProgram->pModelling->set_case(MODELLING_EIGEN) ;
				}
			}

			// parameter
			if (XML_MODELLING_PARAM.compare(attr[i]) == 0) {
				Myfloat param = atof(attr[i+1]) ;
				Singleton::Instance()->pProgram->pModelling->set_param(param) ;
			}
		}

	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// left_hand
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (XML_LEFT_HAND.compare(name) == 0) {

		string left_hand_file = "UNSPECIFIED" ;
		Myfloat left_hand_file_dt = 0.0 ;

		for (Myint i = 0; attr[i]; i += 2) {

			// file
			if (XML_LEFT_HAND_FILE.compare(attr[i]) == 0) {
				left_hand_file = attr[i + 1] ;
			}

			// dt
			if (XML_LEFT_HAND_DT.compare(attr[i]) == 0) {
				left_hand_file_dt = atof(attr[i + 1]) ;
			}
		}

		Program_Guitar *pProgram_Guitar = dynamic_cast<Program_Guitar*>(Singleton::Instance()->pProgram) ;
		if (pProgram_Guitar != NULL)
		{
			pProgram_Guitar->set_left_hand(left_hand_file, left_hand_file_dt) ;
		}
	}
}

//-------------------------------------------------------------------------------------------------------

void XMLCALL django_endElement(void *userData, const char *name)
{
	print_debug(ALL, LIGHT_DEBUG, "django_endElement", name) ;

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// element scheme
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if (XML_SCHEME.compare(name) == 0) {

		//###################################################################################################
		// instantiate scheme
		//###################################################################################################
		if (Singleton::Instance()->pProgram != NULL)
		{
			if (Singleton::Instance()->pProgram->pDomain != NULL)
			{
				Singleton::Instance()->pProgram->pDomain->pScheme =
						Singleton::Instance()->pProgram->pDomain->pScheme_Factory->create();
			}
		}
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// element model
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if (XML_MODEL.compare(name) == 0)
	{
		if (Singleton::Instance()->pProgram != NULL)
		{
			if (Singleton::Instance()->pProgram->pDomain != NULL)
			{
				if (Singleton::Instance()->pProgram->pDomain->pModel != NULL)
				{
					Singleton::Instance()->pProgram->pDomain->pModel->initialize() ;
				}
			}
		}
	}
}

} // namespace django
