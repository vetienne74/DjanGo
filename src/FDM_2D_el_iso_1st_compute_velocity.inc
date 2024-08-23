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
		Myfloat d_sxx_x = D_X(sxx2, iz, ix+1) ;
		Myfloat d_szz_z = D_Z(szz2, iz+1, ix) ;
		Myfloat d_sxz_x = D_X(sxz2, iz, ix) ;
		Myfloat d_sxz_z = D_Z(sxz2, iz, ix) ;

		vx2[ix][iz] += coef12[ix][iz] * (d_sxx_x + d_sxz_z) ;
		vz2[ix][iz] += coef22[ix][iz] * (d_sxz_x + d_szz_z) ;
	}
}
nb_op_kernel += (izMax-izMin) * (ixMax-ixMin) * (4*nb_op_d + 6) ;

//---------------
// CPML points
//---------------

// z- layer
if (get_boundary_type(ZBEG) == PML)
{
	// for vectorization with Intel compiler
	Myfloat ** const mem_sxz_z   = mem_sxz_zBeg->pArray ;
	Myfloat ** const mem_szz_z   = mem_szz_zBeg->pArray ;
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

			Myfloat d_sxz_z = D_Z(sxz2, iz, ix) ;
			mem_sxz_z[ipmlx][ipmlz] = bpml_z[ipmlx][ipmlz] * mem_sxz_z[ipmlx][ipmlz] + apml_z[ipmlx][ipmlz] * d_sxz_z ;
			vx2[ix][iz] += coef12[ix][iz] *  mem_sxz_z[ipmlx][ipmlz] ;

			Myfloat d_szz_z = D_Z(szz2, iz+1, ix) ;
			mem_szz_z[ipmlx][ipmlz] = bpml_half_z[ipmlx][ipmlz] * mem_szz_z[ipmlx][ipmlz] + apml_half_z[ipmlx][ipmlz] * d_szz_z ;
			vz2[ix][iz] += coef22[ix][iz] *  mem_szz_z[ipmlx][ipmlz] ;
		}
	}
	nb_op_bound += (iz2-iz1) * (ix2-ix1) * (2*nb_op_d + 10) ;
}

// z+ layer
if (get_boundary_type(ZEND) == PML)
{
	// for vectorization with Intel compiler
	Myfloat ** const mem_sxz_z   = mem_sxz_zEnd->pArray ;
	Myfloat ** const mem_szz_z   = mem_szz_zEnd->pArray ;
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

			Myfloat d_sxz_z = D_Z(sxz2, iz, ix) ;
			mem_sxz_z[ipmlx][ipmlz] = bpml_z[ipmlx][ipmlz] * mem_sxz_z[ipmlx][ipmlz] + apml_z[ipmlx][ipmlz] * d_sxz_z ;
			vx2[ix][iz] += coef12[ix][iz] *  mem_sxz_z[ipmlx][ipmlz] ;

			Myfloat d_szz_z = D_Z(szz2, iz, ix) ;
			mem_szz_z[ipmlx][ipmlz] = bpml_half_z[ipmlx][ipmlz] * mem_szz_z[ipmlx][ipmlz] + apml_half_z[ipmlx][ipmlz] * d_szz_z ;
			vz2[ix][iz-1] += coef22[ix][iz-1] *  mem_szz_z[ipmlx][ipmlz] ;
		}
	}
	nb_op_bound += (iz2-iz1) * (ix2-ix1) * (2*nb_op_d + 10) ;
}

// x- layer
if (get_boundary_type(XBEG) == PML)
{
	// for vectorization with Intel compiler
	Myfloat ** const mem_sxx_x   = mem_sxx_xBeg->pArray ;
	Myfloat ** const mem_sxz_x   = mem_sxz_xBeg->pArray ;
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

			Myfloat d_sxx_x = D_X(sxx2, iz, ix+1) ;
			mem_sxx_x[ipmlx][ipmlz] = bpml_half_x[ipmlx][ipmlz] * mem_sxx_x[ipmlx][ipmlz] + apml_half_x[ipmlx][ipmlz] * d_sxx_x ;
			vx2[ix][iz] += coef12[ix][iz] *  mem_sxx_x[ipmlx][ipmlz] ;

			Myfloat d_sxz_x = D_X(sxz2, iz, ix) ;
			mem_sxz_x[ipmlx][ipmlz] = bpml_x[ipmlx][ipmlz] * mem_sxz_x[ipmlx][ipmlz] + apml_x[ipmlx][ipmlz] * d_sxz_x ;
			vz2[ix][iz] += coef22[ix][iz] *  mem_sxz_x[ipmlx][ipmlz] ;
		}
	}
	nb_op_bound += (iz2-iz1) * (ix2-ix1) * (2*nb_op_d + 10) ;
}

// x+ layer
if (get_boundary_type(XEND) == PML)
{
	// for vectorization with Intel compiler
	Myfloat ** const mem_sxx_x   = mem_sxx_xEnd->pArray ;
	Myfloat ** const mem_sxz_x   = mem_sxz_xEnd->pArray ;
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

			Myfloat d_sxx_x = D_X(sxx2, iz, ix) ;
			mem_sxx_x[ipmlx][ipmlz] = bpml_half_x[ipmlx][ipmlz] * mem_sxx_x[ipmlx][ipmlz] + apml_half_x[ipmlx][ipmlz] * d_sxx_x ;
			vx2[ix-1][iz] += coef12[ix-1][iz] *  mem_sxx_x[ipmlx][ipmlz] ;

			Myfloat d_sxz_x = D_X(sxz2, iz, ix) ;
			mem_sxz_x[ipmlx][ipmlz] = bpml_x[ipmlx][ipmlz] * mem_sxz_x[ipmlx][ipmlz] + apml_x[ipmlx][ipmlz] * d_sxz_x ;
			vz2[ix][iz] += coef22[ix][iz] *  mem_sxz_x[ipmlx][ipmlz] ;
		}
	}
	nb_op_bound += (iz2-iz1) * (ix2-ix1) * (2*nb_op_d + 10) ;
}

//-----------------------------------
// free surface with image method
// symetry of velocity wavefield
//----------------------------------

// // old
// if (get_boundary_type(ZBEG) == FREESURF)
//   {
//     Myint iz0 = izBeg1 ;
// #pragma omp parallel for        
//     for (Myint ix = ixMin; ix < ixMax; ix++)
//       {
// #pragma ivdep
// 	for (Myint ii = 1; ii <= lstencil; ii++)
// 	  {
// 	    vx2[ix][iz0-ii] = vx2[ix][iz0+ii] ;
// 	    vz2[ix][iz0-ii] = vz2[ix][iz0+ii-1] ;
// 	  }
//       }
//   }

// formulation W, Moczo p 220
// surface on Vz
if (get_boundary_type(ZBEG) == FREESURF)
{
	Myint iz0 = izBeg1 ;
	//#pragma omp parallel for
	for (Myint ix = ixMin; ix < ixMax; ix++)
	{
		// vx at -h/2
		// obtained from 2nd order approximation of
		// d(vx)/dz = - d(vz)/dx
		vx2[ix][iz0] = vx2[ix][iz0+1] + (vz2[ix+1][iz0] - vz2[ix][iz0]) * dz/dx ;
	}

	for (Myint ix = ixMin; ix < ixMax; ix++)
	{
		// vz at -h
		// obtained from 2nd order approximation of d(sxx[z=-h/2])/dt = -d(sxx[z=h/2])/dt

		// compute d(sxx[z=h/2])/dt
		Myfloat dszz_dt = (coef32[ix][iz0+1] * (vx2[ix][iz0+1]- vx2[ix-1][iz0+1]) / dx +
				coef42[ix][iz0+1] * (vz2[ix][iz0+1]- vz2[ix][iz0]) / dz) ;

		// get vz at -h
		Myfloat tmp = coef32[ix][iz0+1] * (vx2[ix][iz0]- vx2[ix-1][iz0]) / dx + dszz_dt ;
		vz2[ix][iz0-1] = -dz * tmp / coef42[ix][iz0+1] + vz2[ix][iz0] ;
	}
}
