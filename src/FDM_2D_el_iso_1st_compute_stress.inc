//===========================================================================
// *** this is an include source file, it can not be compiled on its own ***
//===========================================================================

// for vectorization with Intel compiler

Grid_2D_float *vz_grid = (Grid_2D_float*) (Singleton::Instance()->get_variable(VZ))->get_grid() ;
Myfloat** const vz2 = vz_grid->pArray ;

Grid_2D_float *vx_grid = (Grid_2D_float*) (Singleton::Instance()->get_variable(VX))->get_grid() ;
Myfloat** const vx2 = vx_grid->pArray ;

Grid_2D_float *sxx_grid = (Grid_2D_float*) (Singleton::Instance()->get_variable(SXX))->get_grid() ;
Myfloat** const sxx2 = sxx_grid->pArray ;

Grid_2D_float *szz_grid = (Grid_2D_float*) (Singleton::Instance()->get_variable(SZZ))->get_grid() ;
Myfloat** const szz2 = szz_grid->pArray ;

Grid_2D_float *sxz_grid = (Grid_2D_float*) (Singleton::Instance()->get_variable(SXZ))->get_grid() ;
Myfloat** const sxz2 = sxz_grid->pArray ;

Myfloat ** const coef12 = coef1->pArray ;
Myfloat ** const coef22 = coef2->pArray ;
Myfloat ** const coef32 = coef3->pArray ;
Myfloat ** const coef42 = coef4->pArray ;
Myfloat ** const coef52 = coef5->pArray ;

Myint izMin = izBeg1 ;
Myint izMax = izEnd1 ;
Myint ixMin = ixBeg1 ;
Myint ixMax = ixEnd1 ;

//---------------
// medium points
//---------------

#pragma omp parallel for        
for (Myint ix = ixMin; ix < ixMax; ix++)
{
#pragma ivdep
	for (Myint iz = izMin; iz < izMax; iz++)
	{
		Myfloat d_vx_x = D_X(vx2, iz, ix) ;
		Myfloat d_vz_z = D_Z(vz2, iz, ix) ;
		Myfloat d_vx_z = D_Z(vx2, iz+1, ix) ;
		Myfloat d_vz_x = D_X(vz2, iz, ix+1) ;

		sxx2[ix][iz] += coef32[ix][iz] * d_vx_x + coef42[ix][iz] * d_vz_z ;
		szz2[ix][iz] += coef32[ix][iz] * d_vz_z + coef42[ix][iz] * d_vx_x ;
		sxz2[ix][iz] += coef52[ix][iz] * (d_vx_z + d_vz_x) ;
	}
}
nb_op_kernel += (izMax-izMin) * (ixMax-ixMin) * (4*nb_op_d + 11) ;
ncell       += (izMax-izMin) * (ixMax-ixMin) ;

//---------------
// CPML points
//---------------

// z- layer
if (get_boundary_type(ZBEG) == PML)
{
	// for vectorization with Intel compiler
	Myfloat ** const mem_vz_z    = mem_vz_zBeg->pArray ;
	Myfloat ** const mem_vx_z    = mem_vx_zBeg->pArray ;
	Myfloat ** const apml_z      = apml_zBeg->pArray ;
	Myfloat ** const bpml_z      = bpml_zBeg->pArray ;
	Myfloat ** const apml_half_z = apml_half_zBeg->pArray ;
	Myfloat ** const bpml_half_z = bpml_half_zBeg->pArray ;

	Myint iz1 = izBeg1 ;
	Myint iz2 = izBeg2 ;
	Myint ix1 = ixBeg1 ;
	Myint ix2 = ixEnd1 ;

#pragma omp parallel for
	for (Myint ix = ix1; ix < ix2; ix++)
	{
		Myint ipmlx = ix ;
#pragma ivdep       
		for (Myint iz = iz1; iz < iz2; iz++)
		{
			Myint ipmlz = iz ;
			Myfloat d_vz_z = D_Z(vz2, iz, ix) ;
			mem_vz_z[ipmlx][ipmlz] = bpml_z[ipmlx][ipmlz] * mem_vz_z[ipmlx][ipmlz] + apml_z[ipmlx][ipmlz] * d_vz_z ;
			sxx2[ix][iz] += coef42[ix][iz] *  mem_vz_z[ipmlx][ipmlz] ;
			szz2[ix][iz] += coef32[ix][iz] *  mem_vz_z[ipmlx][ipmlz] ;

			Myfloat d_vx_z = D_Z(vx2, iz+1, ix) ;
			mem_vx_z[ipmlx][ipmlz] = bpml_half_z[ipmlx][ipmlz] * mem_vx_z[ipmlx][ipmlz] + apml_half_z[ipmlx][ipmlz] * d_vx_z ;
			sxz2[ix][iz] += coef52[ix][iz] *  mem_vx_z[ipmlx][ipmlz] ;
		}
	}
	nb_op_bound += (iz2-iz1) * (ix2-ix1) * (2*nb_op_d + 12) ;
}

// z+ layer
if (get_boundary_type(ZEND) == PML)
{
	// for vectorization with Intel compiler
	Myfloat ** const mem_vz_z    = mem_vz_zEnd->pArray ;
	Myfloat ** const mem_vx_z    = mem_vx_zEnd->pArray ;
	Myfloat ** const apml_z      = apml_zEnd->pArray ;
	Myfloat ** const bpml_z      = bpml_zEnd->pArray ;
	Myfloat ** const apml_half_z = apml_half_zEnd->pArray ;
	Myfloat ** const bpml_half_z = bpml_half_zEnd->pArray ;

	Myint iz1 = izEnd2 ;
	Myint iz2 = izEnd1 ;
	Myint ix1 = ixBeg1 ;
	Myint ix2 = ixEnd1 ;
	Myint const offset_z = nlayer_zBeg + nz ;
#pragma omp parallel for
	for (Myint ix = ix1; ix < ix2; ix++)
	{
		Myint ipmlx = ix ;
#pragma ivdep       
		for (Myint iz = iz1; iz < iz2; iz++)
		{
			Myint ipmlz = iz - offset_z ;
			Myfloat d_vz_z = D_Z(vz2, iz, ix) ;
			mem_vz_z[ipmlx][ipmlz] = bpml_z[ipmlx][ipmlz] * mem_vz_z[ipmlx][ipmlz] + apml_z[ipmlx][ipmlz] * d_vz_z ;
			sxx2[ix][iz] += coef42[ix][iz] *  mem_vz_z[ipmlx][ipmlz] ;
			szz2[ix][iz] += coef32[ix][iz] *  mem_vz_z[ipmlx][ipmlz] ;

			Myfloat d_vx_z = D_Z(vx2, iz, ix) ;
			mem_vx_z[ipmlx][ipmlz] = bpml_half_z[ipmlx][ipmlz] * mem_vx_z[ipmlx][ipmlz] + apml_half_z[ipmlx][ipmlz] * d_vx_z ;
			sxz2[ix][iz-1] += coef52[ix][iz-1] *  mem_vx_z[ipmlx][ipmlz] ;
		}
	}
	nb_op_bound += (iz2-iz1) * (ix2-ix1) * (2*nb_op_d + 12) ;
}

// x- layer
if (get_boundary_type(XBEG) == PML)
{
	// for vectorization with Intel compiler
	Myfloat ** const mem_vz_x    = mem_vz_xBeg->pArray ;
	Myfloat ** const mem_vx_x    = mem_vx_xBeg->pArray ;
	Myfloat ** const apml_x      = apml_xBeg->pArray ;
	Myfloat ** const bpml_x      = bpml_xBeg->pArray ;
	Myfloat ** const apml_half_x = apml_half_xBeg->pArray ;
	Myfloat ** const bpml_half_x = bpml_half_xBeg->pArray ;

	Myint iz1 = izBeg1 ;
	Myint iz2 = izEnd1 ;
	Myint ix1 = ixBeg1 ;
	Myint ix2 = ixBeg2 ;

#pragma omp parallel for
	for (Myint ix = ix1; ix < ix2; ix++)
	{
		Myint ipmlx = ix ;
#pragma ivdep       
		for (Myint iz = iz1; iz < iz2; iz++)
		{
			Myint ipmlz = iz ;

			Myfloat d_vx_x = D_X(vx2, iz, ix) ;
			mem_vx_x[ipmlx][ipmlz] = bpml_x[ipmlx][ipmlz] * mem_vx_x[ipmlx][ipmlz] + apml_x[ipmlx][ipmlz] * d_vx_x ;
			sxx2[ix][iz] += coef32[ix][iz] *  mem_vx_x[ipmlx][ipmlz] ;
			szz2[ix][iz] += coef42[ix][iz] *  mem_vx_x[ipmlx][ipmlz] ;

			Myfloat d_vz_x = D_X(vz2, iz, ix+1) ;
			mem_vz_x[ipmlx][ipmlz] = bpml_half_x[ipmlx][ipmlz] * mem_vz_x[ipmlx][ipmlz] + apml_half_x[ipmlx][ipmlz] * d_vz_x ;
			sxz2[ix][iz] += coef52[ix][iz] *  mem_vz_x[ipmlx][ipmlz] ;
		}
	}
	nb_op_bound += (iz2-iz1) * (ix2-ix1) * (2*nb_op_d + 12) ;
}

// x+ layer
if (get_boundary_type(XEND) == PML)
{
	// for vectorization with Intel compiler
	Myfloat ** const mem_vz_x    = mem_vz_xEnd->pArray ;
	Myfloat ** const mem_vx_x    = mem_vx_xEnd->pArray ;
	Myfloat ** const apml_x      = apml_xEnd->pArray ;
	Myfloat ** const bpml_x      = bpml_xEnd->pArray ;
	Myfloat ** const apml_half_x = apml_half_xEnd->pArray ;
	Myfloat ** const bpml_half_x = bpml_half_xEnd->pArray ;

	Myint iz1 = izBeg1 ;
	Myint iz2 = izEnd1 ;
	Myint ix1 = ixEnd2 ;
	Myint ix2 = ixEnd1 ;
	Myint const offset_x = nlayer_xBeg + nx ;

#pragma omp parallel for
	for (Myint ix = ix1; ix < ix2; ix++)
	{
		Myint ipmlx = ix - offset_x ;
#pragma ivdep       
		for (Myint iz = iz1; iz < iz2; iz++)
		{
			Myint ipmlz = iz ;

			Myfloat d_vx_x = D_X(vx2, iz, ix) ;
			mem_vx_x[ipmlx][ipmlz] = bpml_x[ipmlx][ipmlz] * mem_vx_x[ipmlx][ipmlz] + apml_x[ipmlx][ipmlz] * d_vx_x ;
			sxx2[ix][iz] += coef32[ix][iz] *  mem_vx_x[ipmlx][ipmlz] ;
			szz2[ix][iz] += coef42[ix][iz] *  mem_vx_x[ipmlx][ipmlz] ;

			Myfloat d_vz_x = D_X(vz2, iz, ix) ;
			mem_vz_x[ipmlx][ipmlz] = bpml_half_x[ipmlx][ipmlz] * mem_vz_x[ipmlx][ipmlz] + apml_half_x[ipmlx][ipmlz] * d_vz_x ;
			sxz2[ix-1][iz] += coef52[ix-1][iz] *  mem_vz_x[ipmlx][ipmlz] ;
		}
	}
	nb_op_bound += (iz2-iz1) * (ix2-ix1) * (2*nb_op_d + 12) ;
}

//------------------------------------
// free surface with image method
// anti-symetry of stress wavefield
//------------------------------------

// // old
// if (get_boundary_type(ZBEG) == FREESURF)
//   {
//     Myint iz0 = izBeg1 ;
// #pragma omp parallel for        
//     for (Myint ix = ixMin; ix < ixMax; ix++)
//       {

// 	// szz and sxx at the free surface z=izBeg1
// 	szz2[ix][iz0] = 0. ;
// 	sxx2[ix][iz0] = 0. ;


// #pragma ivdep
// 	for (Myint ii = 1; ii <= lstencil; ii++)
// 	  {
// 	    // sxz  
// 	    sxz2[ix][iz0-ii] = -sxz2[ix][iz0+ii-1] ;

// 	    // szz 
// 	    szz2[ix][iz0-ii] = -szz2[ix][iz0+ii] ;	  
// 	  }
//       }
//   }

// Moczo p 220, W formulation
// surface on Vz

if (get_boundary_type(ZBEG) == FREESURF)
{
	Myint iz0 = izBeg1 ;
#pragma omp parallel for        
	for (Myint ix = ixMin; ix < ixMax; ix++)
	{
		// szx at surface = 0
		sxz2[ix][iz0] = 0. ;


#pragma ivdep
		for (Myint ii = 1; ii <= lstencil; ii++)
		{
			// anti-symetry of szz
			szz2[ix][iz0-ii+1] = -szz2[ix][iz0+ii] ;

			// anti-symetry of szx
			sxz2[ix][iz0-ii] = -sxz2[ix][iz0+ii] ;
		}
	}
}


