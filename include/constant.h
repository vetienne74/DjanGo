#ifndef DJANGO_CONSTANT_H_
#define DJANGO_CONSTANT_H_

#include <string>
#include <complex>

#include "type_def.h"

namespace django {

const Myint CURRENT_VERSION = 1 ;

//-------------------------------------------------------------------------------------------------------
// Useful numerical values
//-------------------------------------------------------------------------------------------------------

const Myfloat   PI         = 3.1415926535897 ;
const Mycomplex ZERO_CMPLX (0.,0.) ;
const Mycomplex I_CMPLX    (0.,1.) ;

#ifdef _DOUBLE_PRECISION_
const Myfloat MY_EPSILON    = 1.e-12 ;
#else
const Myfloat MY_EPSILON    = 1.e-5 ;
#endif

//-------------------------------------------------------------------------------------------------------
// Flags
//-------------------------------------------------------------------------------------------------------
const Myint  NOT_FOUND     = -1 ;
const Myint  NO_NEIGH      = -1 ;
const Myint  NOT_SPECIFIED = -1 ;
const string UNSPECIFIED = "UNSPECIFIED" ;

//-------------------------------------------------------------------------------------------------------
// Limits
//-------------------------------------------------------------------------------------------------------
const Myint MAX_BOUNDARY = 100 ;

// max allowed time since compilation in seconds
const time_t BINARY_VALIDITY_TIME = 30 * 24 * 60 * 60 ; // 30 days


//-------------------------------------------------------------------------------------------------------
// Miscellaneous
//-------------------------------------------------------------------------------------------------------
const Myfloat DEFAULT_AMP_SRC = 1.0 ;

//-------------------------------------------------------------------------------------------------------
// File names
//-------------------------------------------------------------------------------------------------------

// input files
//=============

const char GRAD_CONFIG_IN_FILE[]         = "gradient.config" ;
const char DATA_CONFIG_IN_FILE[]         = "data.config" ;

// output files
//==============

const char VX_TIME_REC_OUT_FILE[]        = "vx.time.rec.django.out.bin" ;
const char VY_TIME_REC_OUT_FILE[]        = "vy.time.rec.django.out.bin" ;
const char VZ_TIME_REC_OUT_FILE[]        = "vz.time.rec.django.out.bin" ;
const char PR_TIME_REC_OUT_FILE[]        = "pr.time.rec.django.out.bin" ;

const char TIME_SNAPSHOT_OUT_FILE[]      = ".time.snapshot.django.out.bin" ;
const char TIME_REC_OUT_FILE[]           = ".time.rec.django.out.bin" ;

const char VX_TIME_SNAPSHOT_OUT_FILE[]   = "vx.time.snapshot.django.out.bin" ;
const char VY_TIME_SNAPSHOT_OUT_FILE[]   = "vy.time.snapshot.django.out.bin" ;
const char VZ_TIME_SNAPSHOT_OUT_FILE[]   = "vz.time.snapshot.django.out.bin" ;
const char PR_TIME_SNAPSHOT_OUT_FILE[]   = "pr.time.snapshot.django.out.bin" ;

const char ENERGY_TIME_REC_OUT_FILE[]    = "energy.time.rec.django.out.ascii" ;

const char VX_FREQ_REC_OUT_FILE[]        = "vx.freq.rec.django.out.bin" ;
const char VY_FREQ_REC_OUT_FILE[]        = "vy.freq.rec.django.out.bin" ;
const char VZ_FREQ_REC_OUT_FILE[]        = "vz.freq.rec.django.out.bin" ;
const char PR_FREQ_REC_OUT_FILE[]        = "pr.freq.rec.django.out.bin" ;

const char PR_FREQ_ADJ_SRC_OUT_FILE[]    = "pr.freq.rec.adjoint.django.out.bin" ;

const char VX_FREQ_GRID_OUT_FILE[]       = "vx.freq.grid.django.out.bin" ;
const char VY_FREQ_GRID_OUT_FILE[]       = "vy.freq.grid.django.out.bin" ;
const char VZ_FREQ_GRID_OUT_FILE[]       = "vz.freq.grid.django.out.bin" ;
const char PR_FREQ_GRID_OUT_FILE[]       = "pr.freq.grid.django.out.bin" ;

const char VP_OUT_FILE[]                 = "vp.django.out.bin" ;

const char GRADIENT_VP_OUT_FILE[]        = "gradient.vp.django.out.bin" ;
const char GRAD_PRECOND_VP_OUT_FILE[]    = "gradient.precond.vp.django.out.bin" ;
const char PRECOND_VP_OUT_FILE[]         = "precond.vp.django.out.bin" ;

const char FCOST_OUT_FILE[]              = "fcost.django.out.ascii" ;

const char SRC_WAVELET_OUT_FILE[]        = "src.wavelet.django.out.bin" ;
const char ESTIM_SRC_WAVELET_OUT_FILE[]  = "estim.src.wavelet.django.out.bin" ;

const char PERF_OUT_FILE[]               = "perf.django.out.ascii" ;

const char EIGEN_ERROR_NODE_OUT_FILE[]   = "eigen.error.node.django.out.ascii" ;
const char EIGEN_ERROR_REC_OUT_FILE[]    = "eigen.error.rec.django.out.ascii" ;

const char FEM_MATRIX_OUT_FILE[]         = "fem.matrix.django.out.ascii" ;

const char GRID_POINT_OUT_FILE[]         = "grid.point.django.out.ascii" ;
const char NODE_COORD_OUT_FILE[]         = "node.coord.django.out.ascii" ;

const char MESH_VTK_OUT_FILE[]           = "mesh.django.out.vtk" ;
const char NODE_VTK_OUT_FILE[]           = "node.django.out.vtk" ;
const char REC_VTK_OUT_FILE[]            = "rec.django.out.vtk" ;
const char SRC_VTK_OUT_FILE[]            = "src.django.out.vtk" ;

const char MESH_ELEM_OUT_FILE[]          = "mesh.elem.django.out.bin" ;
const char MESH_VP_OUT_FILE[]            = "mesh.vp.django.out.bin" ;
const char MESH_VS_OUT_FILE[]            = "mesh.vs.django.out.bin" ;
const char MESH_RHO_OUT_FILE[]           = "mesh.rho.django.out.bin" ;
const char MESH_ORDER_OUT_FILE[]         = "mesh.order.django.out.bin" ;
const char MESH_NELPERLAMBDA_OUT_FILE[]  = "mesh.nelperlambda.django.out.bin" ;
const char MESH_TMIN_OUT_FILE[]          = "mesh.tmin.django.out.bin" ;

const char STAT_TIME_STEP_OUT_FILE[]     = "stat.time.step.django.out.ascii" ;

} // namespace django

#endif
