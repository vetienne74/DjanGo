#ifndef DJANGO_TYPE_DEF_H_
#define DJANGO_TYPE_DEF_H_

#include <complex>

using namespace std;

namespace django {

//--------------------------------------------------
//
//          T Y P E   D E F I N I T I O N S
//
//--------------------------------------------------

// Basic types
#ifdef _DOUBLE_PRECISION_
typedef double Myfloat ;
typedef complex <double>  Mycomplex ;
#else
typedef float Myfloat ;
typedef complex <float>  Mycomplex ;
#endif

typedef int        Myint ;
typedef int        Myint32 ;
typedef long int   Myint64 ;
typedef float      Myfloat32 ;
typedef double     Myfloat64 ;

// Return codes used by all methods
enum Rtn_code {RTN_CODE_OK=0, RTN_CODE_KO=-1} ;

// Debug level
enum Debug_level {NO_DEBUG=-1, LIGHT_DEBUG, MID_DEBUG, FULL_DEBUG} ;

// Program mode
enum Prog_type {NO_PROG=-1, MODELLING_PROG, FWI_PROG, RTM_PROG, GUITAR_PROG} ;

// Modelling case
enum Modelling_case_type {MODELLING_STD, MODELLING_EIGEN} ;

// Space dimension
enum Space_dim {NO_DIM=-1, ONE=1, TWO=2, THREE=3} ;

// Numerical scheme (method)
enum Scheme_method {NO_SCHEME_METHOD, SCHEME_FDM, SCHEME_FEM} ;

// Numerical scheme (type)
enum Scheme_type {NO_SCHEME_TYPE, SCHEME_STAGGERED, SCHEME_CGM, SCHEME_DGM, SCHEME_MGM} ;

// Equation type
enum Eq_type {NO_EQ_TYPE, ACOUSTIC, ELASTIC, AC_LOSSY} ;

// Equation order
enum Eq_order {NO_EQ_ORDER, ORDER_1ST, ORDER_2ND} ;

// Edge
enum Edge_type {NO_EDGE, ZBEG, ZEND, XBEG, XEND, YBEG, YEND} ;

// Type of boundary
// PML      = Komatitsch's CPML
// SPG      = Cerjan's sponge
// SPG2     = Israeli's sponge
// FREESURF = Mirror condition
enum Boundary_type {NO_BOUNDARY, PML, RANDOM, SPG, SPG2, FREESURF, RIGID, INTERNAL} ;

// Acquisition type
enum Acqui_type {SEISMIC_FIXED=1, SEISMIC_STREAMER=2} ;

// For display in output report
enum Display_type {MASTER=0, ALL=1} ;

// Incident and adjoint wavefield
enum Wavefield_type {INCIDENT=0, ADJOINT=1} ;

// Gradient type
enum Grad_type {ADJOINT_METHOD=0, FD_FCOST=1} ;

// Fcost type
enum Fcost_type {DIFF_L2=0, DIFF_L1=1, LOG_L2=2, LOG_L2_PHASE=3,  LOG_L2_AMP=4} ;

// source time function
enum Src_func {NO_SRC_FUNC, RICKER_PP, RICKER_P, RICKER, RICKER_D, MONO_FREQ, SRC_FILE} ;

// source type
enum Src_type {NO_SRC_TYPE, SRC_POINT, SRC_GAUSSIAN, SRC_SINC} ;

// source type
enum Src_stype {NO_SRC_STYPE, EXPLOSIVE, FORCE_Z} ;

// self-initialisation type (for variables)
enum Self_init_type {NO_SELF_INIT, INIT_FROM_FILE, INIT_FROM_CONST} ;

// time or freq. domain
enum Domain_type {NO_DOMAIN, TIME, FREQ} ;

// dynamic or static front (expansible computing domain)
enum Front_type {FRONT_STATIC, FRONT_DYN_VMAX} ;

// model type
enum Model_type {NO_MODEL_TYPE, GRID} ;

// model sub-type
enum Model_sub_type {NO_MODEL_SUBTYPE, REGULAR, UNREGULAR} ;

// variable type
enum Var_type {NO_VAR_TYPE, VP, VS, RHO, LOSS1, LOSS2, INV_RHOVP2, INV_RHOVS2, INV_RHOVP2MINUSVS2,
	VX, VY, VZ, PR, PRN, PRC, PRP, SXX, SZZ, SXZ, TAU, TAUP, TAUPP} ;

// coef type
enum Coef_type {NO_COEF, COEF1, COEF2, COEF3} ;

// node type
enum Node_type {NO_NODE_TYPE, NODE_TYPE_EQD, NODE_TYPE_GLL} ;

// node distribution
enum Node_distrib {NO_NODE_DISTRIB, NODE_DISTRIB_UNIFORM, NODE_DISTRIB_RANDOM, NODE_DISTRIB_MODEL} ;

// node integration
enum Node_integ {NO_NODE_INTEG, NODE_INTEG_GL, NODE_INTEG_GLL} ;

// properties within element
enum Prop_type {NO_PROP_TYPE, PROP_CONST, PROP_AVER, PROP_LINEAR} ;

// flux type
enum Flux_type {NO_FLUX_TYPE, FLUX_TYPE_CENTERED, FLUX_TYPE_UPWIND} ;

// region type
enum Region_type {NO_REGION, GHOST, LAYER, MEDIUM} ;

// axis
enum Axis_type {Z_AXIS, X_AXIS, Y_AXIS} ;

// neighbours
enum Neigh_type {I_ZPREV, I_ZNEXT, I_XPREV, I_XNEXT, I_YPREV, I_YNEXT} ;

// adaptivity type
enum Adapt_type {NO_ADAPT, ADAPT_STATIC, ADAPT_FREQTIME} ;

} // namespace django

#endif
