#ifndef DJANGO_PARSE_XML_H_
#define DJANGO_PARSE_XML_H_

#include "expat.h"
#ifdef XML_LARGE_SIZE
#if defined(XML_USE_MSC_EXTENSIONS) && _MSC_VER < 1400
#define XML_FMT_INT_MOD "I64"
#else
#define XML_FMT_INT_MOD "ll"
#endif
#else
#define XML_FMT_INT_MOD "l"
#endif

#include "type_def.h"

namespace django {

void XMLCALL django_startElement(void *userData, const char *name, const char **attr) ;
void XMLCALL django_endElement(void *userData, const char *name) ;

//-------------------------------------------------------------------------------------------------------
// XML ELEMENTS AND ATTRIBUTES
//-------------------------------------------------------------------------------------------------------

// CLASS DJANGO
const string XML_DJANGO           = "django" ;
const string XML_DJANGO_VERSION   = "version" ;

// CLASS PROGRAM
const string XML_PROG               = "program" ;
const string XML_PROG_TYPE          = "type" ;
const string XML_PROG_TYPE_MOD      = "modelling" ;
const string XML_PROG_TYPE_RTM      = "rtm" ;
const string XML_PROG_TYPE_FWI      = "fwi" ;
const string XML_PROG_TYPE_GUITAR   = "guitar" ;

// CLASS SCHEME
const string XML_SCHEME             = "scheme" ;
const string XML_SCHEME_METH        = "method" ;
const string XML_SCHEME_METH_FDM    = "fdm" ;
const string XML_SCHEME_METH_FEM    = "fem" ;
const string XML_SCHEME_TYPE        = "type" ;
const string XML_SCHEME_TYPE_STG    = "staggered" ;
const string XML_SCHEME_TYPE_CGM    = "continuous" ;
const string XML_SCHEME_TYPE_DGM    = "discontinuous" ;
const string XML_SCHEME_TYPE_MGM    = "mixed" ;

const string XML_SCHEME_ADAPT               = "adaptivity" ;
const string XML_SCHEME_ADAPT_TYPE          = "type" ;
const string XML_SCHEME_ADAPT_TYPE_NONE     = "none" ;
const string XML_SCHEME_ADAPT_TYPE_STATIC   = "static" ;
const string XML_SCHEME_ADAPT_TYPE_FREQTIME = "freqTime" ;
const string XML_SCHEME_ADAPT_PMAX          = "pmax" ;
const string XML_SCHEME_ADAPT_PMIN          = "pmin" ;
const string XML_SCHEME_ADAPT_FMAX0         = "fmax0" ;

const string XML_SCHEME_ADAPT_FREQTIME      = "freqTime" ;
const string XML_SCHEME_ADAPT_FREQTIME_TMIN = "tmin" ;
const string XML_SCHEME_ADAPT_FREQTIME_FMAX = "fmax" ;

const string XML_ACCURACY           = "accuracy" ;
const string XML_ACCURACY_TIME      = "time" ;
const string XML_ACCURACY_TIME_O2   = "2" ;
const string XML_ACCURACY_TIME_DT   = "dt" ;
const string XML_ACCURACY_SPACE     = "space" ;
const string XML_ACCURACY_SPACE_O2  = "2" ;
const string XML_ACCURACY_SPACE_O4  = "4" ;
const string XML_ACCURACY_SPACE_O8  = "8" ;
const string XML_ACCURACY_SPACE_O12 = "12" ;
const string XML_ACCURACY_SPACE_O16 = "16" ;
const string XML_ACCURACY_RATIO_CFL = "ratio_cfl" ;

const string XML_NODE               = "node" ;
const string XML_NODE_TYPE          = "type" ;
const string XML_NODE_TYPE_EQD      = "equidistant" ;
const string XML_NODE_TYPE_GLL      = "GLL" ;
const string XML_NODE_DIST          = "distribution" ;
const string XML_NODE_DIST_UNIFORM  = "uniform" ;
const string XML_NODE_DIST_RANDOM   = "random" ;
const string XML_NODE_DIST_MODEL    = "model" ;
const string XML_NODE_INTEG         = "integration" ;
const string XML_NODE_INTEG_GL      = "GL" ;
const string XML_NODE_INTEG_GLL     = "GLL" ;

const string XML_FLUX               = "flux" ;
const string XML_FLUX_TYPE          = "type" ;
const string XML_FLUX_CENTER        = "centered" ;
const string XML_FLUX_UPWIND        = "upwind" ;

const string XML_DYN                = "dynamic" ;
const string XML_DYN_FRONT          = "front" ;
const string XML_DYN_FRONT_STATIC   = "static" ;
const string XML_DYN_FRONT_VMAX     = "vmax" ;

const string XML_EQ                 = "equation" ;
const string XML_EQ_ORD             = "order" ;
const string XML_EQ_ORD_1           = "1" ;
const string XML_EQ_ORD_2           = "2" ;
const string XML_EQ_TYPE            = "type" ;
const string XML_EQ_TYPE_AC         = "acoustic" ;
const string XML_EQ_TYPE_EL         = "elastic" ;
const string XML_EQ_TYPE_AC_LOSSY   = "acoustic_lossy" ;

const string XML_MESH               = "mesh" ;
const string XML_MESH_NELEM_X       = "nelx" ;
const string XML_MESH_NELEM_Z       = "nelz" ;
const string XML_MESH_PROP          = "properties" ;
const string XML_MESH_PROP_CONST    = "constant" ;
const string XML_MESH_PROP_AVER     = "average" ;
const string XML_MESH_PROP_LINEAR   = "linear" ;

// MODELLING CLASS
const string XML_TIME               = "time" ;
const string XML_TIME_TMAX          = "tmax" ;
const string XML_TIME_DT            = "dt" ;
const string XML_TIME_RATE          = "rate" ;

const string XML_FREQ               = "frequency" ;
const string XML_FREQ_FMIN          = "fmin" ;
const string XML_FREQ_FMAX          = "fmax" ;
const string XML_FREQ_DF            = "df" ;

const string XML_TIME_ORDER         = "order" ;
const string XML_TIME_ORDER_2       = "2" ;

const string XML_ACQUI              = "acquisition" ;
const string XML_ACQUI_FILE         = "file" ;

const string XML_BOUND              = "boundary" ;
const string XML_BOUND_EDGE         = "edge" ;
const string XML_BOUND_EDGE_ZBEG    = "zbeg" ;
const string XML_BOUND_EDGE_ZEND    = "zend" ;
const string XML_BOUND_EDGE_XBEG    = "xbeg" ;
const string XML_BOUND_EDGE_XEND    = "xend" ;
const string XML_BOUND_EDGE_YBEG    = "ybeg" ;
const string XML_BOUND_EDGE_YEND    = "yend" ;
const string XML_BOUND_TYPE         = "type" ;
const string XML_BOUND_TYPE_CPML    = "cpml" ;
const string XML_BOUND_TYPE_SPG     = "sponge" ;
const string XML_BOUND_TYPE_SPG2    = "sponge2" ;
const string XML_BOUND_TYPE_FS      = "freesurf" ;
const string XML_BOUND_TYPE_RAND    = "random" ;
const string XML_BOUND_TYPE_RIG     = "rigid" ;
const string XML_BOUND_WIDTH        = "width" ;
const string XML_BOUND_COEF         = "coef" ;

const string XML_OUTPUT_ENERGY      = "energy" ;
const string XML_OUTPUT_ENERGY_ON   = "on" ;
const string XML_OUTPUT_ENERGY_OFF  = "off" ;
const string XML_OUTPUT_FREQ        = "frequency" ;
const string XML_OUTPUT_FREQ_ON     = "on" ;
const string XML_OUTPUT_FREQ_OFF    = "off" ;

const string XML_SNAPSHOT           = "snapshot" ;
const string XML_SNAPSHOT_TMIN      = "tmin" ;
const string XML_SNAPSHOT_TMAX      = "tmax" ;
const string XML_SNAPSHOT_DT        = "dt" ;

const string XML_PIXEL              = "pixel" ;
const string XML_PIXEL_XMIN         = "xmin" ;
const string XML_PIXEL_XMAX         = "xmax" ;
const string XML_PIXEL_NX           = "nx" ;
const string XML_PIXEL_ZMIN         = "zmin" ;
const string XML_PIXEL_ZMAX         = "zmax" ;
const string XML_PIXEL_NZ           = "nz" ;

const string XML_SOURCE             = "source" ;
const string XML_SOURCE_FREQ        = "frequency" ;
const string XML_SOURCE_FUNC        = "function" ;
const string XML_SOURCE_FUNC_RIC    = "ricker" ;
const string XML_SOURCE_FUNC_RIC_P  = "rickerp" ;
const string XML_SOURCE_FUNC_RIC_PP = "rickerpp" ;
const string XML_SOURCE_FUNC_RIC_D  = "rickerd" ;
const string XML_SOURCE_FUNC_MONO   = "mono" ;
const string XML_SOURCE_FUNC_FILE   = "file" ;
const string XML_SOURCE_FUNC_NAME   = "name" ;
const string XML_SOURCE_FUNC_DT     = "dt" ;
const string XML_SOURCE_FUNC_T0     = "t0" ;
const string XML_SOURCE_TYPE        = "type" ;
const string XML_SOURCE_TYPE_POINT  = "point" ;
const string XML_SOURCE_TYPE_GAUSS  = "gaussian" ;
const string XML_SOURCE_TYPE_SINC   = "sinc" ;
const string XML_SOURCE_STYPE       = "subtype" ;
const string XML_SOURCE_STYPE_FZ    = "fz" ;
const string XML_SOURCE_STYPE_EXP   = "explosive" ;
const string XML_SOURCE_SIGMA       = "sigma" ;

// CLASS MODEL
const string XML_MODEL              = "model" ;
const string XML_MODEL_TYPE         = "type" ;
const string XML_MODEL_TYPE_GRID    = "grid" ;
const string XML_MODEL_DIM          = "dimension" ;
const string XML_MODEL_DIM_1D       = "1" ;
const string XML_MODEL_DIM_2D       = "2" ;
const string XML_MODEL_DIM_3D       = "3" ;

const string XML_SIZE               = "size" ;
const string XML_SIZE_NX            = "nx" ;
const string XML_SIZE_NY            = "ny" ;
const string XML_SIZE_NZ            = "nz" ;

const string XML_SAMPLING           = "sampling" ;
const string XML_SAMPLING_DX        = "dx" ;
const string XML_SAMPLING_DY        = "dy" ;
const string XML_SAMPLING_DZ        = "dz" ;
const string XML_SAMPLING_DX_RAND   = "dxrandom" ;
const string XML_SAMPLING_DY_RAND   = "dyrandom" ;
const string XML_SAMPLING_DZ_RAND   = "dzrandom" ;

const string XML_PARAM              = "parameter" ;
const string XML_PARAM_TYPE         = "type" ;
const string XML_PARAM_TYPE_VP      = "vp" ;
const string XML_PARAM_TYPE_VS      = "vs" ;
const string XML_PARAM_TYPE_RHO     = "rho" ;
const string XML_PARAM_TYPE_LOSS1   = "loss1" ;
const string XML_PARAM_TYPE_LOSS2   = "loss2" ;
const string XML_PARAM_FILE         = "file" ;
const string XML_PARAM_CONST        = "constant" ;

const string XML_INV                = "inversion" ;
const string XML_INV_NITER          = "niter" ;
const string XML_INV_NTRY           = "ntry" ;
const string XML_INV_INIT_TRY       = "init_try" ;

const string XML_DOMAIN             = "domain" ;
const string XML_DOMAIN_TYPE        = "type" ;
const string XML_DOMAIN_TYPE_FREQ   = "freq" ;
const string XML_DOMAIN_MIN         = "min" ;
const string XML_DOMAIN_MAX         = "max" ;
const string XML_DOMAIN_DELTA       = "delta" ;

const string XML_MODELLING          = "modelling" ;
const string XML_MODELLING_CASE     = "case" ;
const string XML_MODELLING_CASE_STD = "standard" ;
const string XML_MODELLING_CASE_EIG = "eigen" ;
const string XML_MODELLING_PARAM    = "parameter" ;

const string XML_LEFT_HAND          = "left_hand" ;
const string XML_LEFT_HAND_FILE     = "file" ;
const string XML_LEFT_HAND_DT       = "dt" ;

} // namespace django

#endif
