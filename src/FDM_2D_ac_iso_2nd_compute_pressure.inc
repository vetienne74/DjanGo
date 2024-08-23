//===========================================================================
// *** this is an include source file, it can not be compiled on its own ***
//===========================================================================

// get pointers to variables
Grid_2D_float *prc_grid = (Grid_2D_float*) (Singleton::Instance()->get_variable(PRC))->get_grid() ;
Myfloat ** const prc = prc_grid->pArray ;

Grid_2D_float *prn_grid = (Grid_2D_float*) (Singleton::Instance()->get_variable(PRN))->get_grid() ;
Myfloat ** const prn = prn_grid->pArray ;

Myfloat ** const coef12 = coef1->pArray ;

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
		Myfloat d2_pr_x = D2_X(prc, iz, ix) ;
		Myfloat d2_pr_z = D2_Z(prc, iz, ix) ;

		prn[ix][iz] = coef12[ix][iz] * (d2_pr_x + d2_pr_z) - prn[ix][iz] + 2 * prc[ix][iz] ;
	}
}
nb_op_kernel += (izMax-izMin) * (ixMax-ixMin) * (2*nb_op_d2 + 5) ;
ncell       += (izMax-izMin) * (ixMax-ixMin) ;

//---------------
// CPML points
//---------------

// z- layer
if (get_boundary_type(ZBEG) == PML)
{
	// for vectorization with Intel compiler
	Myfloat ** const mem_vz_z    = mem_vz_zBeg->pArray ;
	Myfloat ** const mem_pr_z    = mem_pr_zBeg->pArray ;
	Myfloat ** const apml_z      = apml_zBeg->pArray ;
	Myfloat ** const bpml_z      = bpml_zBeg->pArray ;
	Myfloat ** const apml_half_z = apml_half_zBeg->pArray ;
	Myfloat ** const bpml_half_z = bpml_half_zBeg->pArray ;

	const Myfloat coef2 = dt / (TEMP_RHO_CONST * dz) ;
	const Myfloat invcoef2 = 1/coef2 ;

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
		Myfloat d_pr_z = D_Z(prc, iz+1, ix) ;
		mem_pr_z[ipmlx][ipmlz] = bpml_half_z[ipmlx][ipmlz] * mem_pr_z[ipmlx][ipmlz] + apml_half_z[ipmlx][ipmlz] * d_pr_z ;
	}
}
#pragma omp parallel for
for (Myint ix = ix1; ix < ix2; ix++)
{
	Myint ipmlx = ix ;
#pragma ivdep       
	for (Myint iz = iz1; iz < iz2; iz++) 
	{
		Myint ipmlz = iz ;
		Myfloat d_mem_pr_z = D_Z(mem_pr_z, ipmlz, ipmlx) ;
		Myfloat d2_pr_z = D2_Z(prc, iz, ix) ;
		mem_vz_z[ipmlx][ipmlz] = bpml_z[ipmlx][ipmlz] * mem_vz_z[ipmlx][ipmlz] + apml_z[ipmlx][ipmlz] * coef2 * (d2_pr_z + d_mem_pr_z) ;
		prn[ix][iz] += coef12[ix][iz] * (d_mem_pr_z + invcoef2 * mem_vz_z[ipmlx][ipmlz]) ;
	}
}
nb_op_bound += (iz2-iz1) * (ix2-ix1) * (2*nb_op_d + nb_op_d2 + 12) ;
}

// z+ layer
if (get_boundary_type(ZEND) == PML)
{
	// for vectorization with Intel compiler
	Myfloat ** const mem_vz_z    = mem_vz_zEnd->pArray ;
	Myfloat ** const mem_pr_z    = mem_pr_zEnd->pArray ;
	Myfloat ** const apml_z      = apml_zEnd->pArray ;
	Myfloat ** const bpml_z      = bpml_zEnd->pArray ;
	Myfloat ** const apml_half_z = apml_half_zEnd->pArray ;
	Myfloat ** const bpml_half_z = bpml_half_zEnd->pArray ;

	const Myfloat coef2 = dt / (TEMP_RHO_CONST * dz) ;
	const Myfloat invcoef2 = 1/coef2 ;

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
		Myint ipmlz = iz  - offset_z ;
		Myfloat d_pr_z = D_Z(prc, iz, ix) ;
		mem_pr_z[ipmlx][ipmlz] = bpml_half_z[ipmlx][ipmlz] * mem_pr_z[ipmlx][ipmlz] + apml_half_z[ipmlx][ipmlz] * d_pr_z ;
	}
}
#pragma omp parallel for
for (Myint ix = ix1; ix < ix2; ix++)
{
	Myint ipmlx = ix ;
#pragma ivdep       
	for (Myint iz = iz1; iz < iz2; iz++) 
	{
		Myint ipmlz = iz - offset_z ;
		Myfloat d_mem_pr_z = D_Z(mem_pr_z, ipmlz+1, ipmlx) ;
		Myfloat d2_pr_z = D2_Z(prc, iz, ix) ;
		mem_vz_z[ipmlx][ipmlz] = bpml_z[ipmlx][ipmlz] * mem_vz_z[ipmlx][ipmlz] + apml_z[ipmlx][ipmlz] * coef2 * (d2_pr_z + d_mem_pr_z) ;
		prn[ix][iz] += coef12[ix][iz] * (d_mem_pr_z + invcoef2 * mem_vz_z[ipmlx][ipmlz]) ;
	}
}
nb_op_bound += (iz2-iz1) * (ix2-ix1) * (2*nb_op_d + nb_op_d2 + 12) ;
}

// x- layer
if (get_boundary_type(XBEG) == PML)
{
	// for vectorization with Intel compiler
	Myfloat ** const mem_vx_x    = mem_vx_xBeg->pArray ;
	Myfloat ** const mem_pr_x    = mem_pr_xBeg->pArray ;
	Myfloat ** const apml_x      = apml_xBeg->pArray ;
	Myfloat ** const bpml_x      = bpml_xBeg->pArray ;
	Myfloat ** const apml_half_x = apml_half_xBeg->pArray ;
	Myfloat ** const bpml_half_x = bpml_half_xBeg->pArray ;

	const Myfloat coef2 = dt / (TEMP_RHO_CONST * dx) ;
	const Myfloat invcoef2 = 1/coef2 ;

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
		Myfloat d_pr_x = D_X(prc, iz, ix+1) ;
		mem_pr_x[ipmlx][ipmlz] = bpml_half_x[ipmlx][ipmlz] * mem_pr_x[ipmlx][ipmlz] + apml_half_x[ipmlx][ipmlz] * d_pr_x ;
	}
}
#pragma omp parallel for
for (Myint ix = ix1; ix < ix2; ix++)
{
	Myint ipmlx = ix ;
#pragma ivdep       
	for (Myint iz = iz1; iz < iz2; iz++) 
	{
		Myint ipmlz = iz ;
		Myfloat d_mem_pr_x = D_X(mem_pr_x, ipmlz, ipmlx) ;
		Myfloat d2_pr_x = D2_X(prc, iz, ix) ;
		mem_vx_x[ipmlx][ipmlz] = bpml_x[ipmlx][ipmlz] * mem_vx_x[ipmlx][ipmlz] + apml_x[ipmlx][ipmlz] * coef2 * (d2_pr_x + d_mem_pr_x) ;
		prn[ix][iz] += coef12[ix][iz] * (d_mem_pr_x + invcoef2 * mem_vx_x[ipmlx][ipmlz]) ;
	}
}
nb_op_bound += (iz2-iz1) * (ix2-ix1) * (2*nb_op_d + nb_op_d2 + 12) ;
}

// x+ layer
if (get_boundary_type(XEND) == PML)
{
	// for vectorization with Intel compiler
	Myfloat ** const mem_vx_x    = mem_vx_xEnd->pArray ;
	Myfloat ** const mem_pr_x    = mem_pr_xEnd->pArray ;
	Myfloat ** const apml_x      = apml_xEnd->pArray ;
	Myfloat ** const bpml_x      = bpml_xEnd->pArray ;
	Myfloat ** const apml_half_x = apml_half_xEnd->pArray ;
	Myfloat ** const bpml_half_x = bpml_half_xEnd->pArray ;

	const Myfloat coef2 = dt / (TEMP_RHO_CONST * dx) ;
	const Myfloat invcoef2 = 1/coef2 ;

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
		Myfloat d_pr_x = D_X(prc, iz, ix) ;
		mem_pr_x[ipmlx][ipmlz] = bpml_half_x[ipmlx][ipmlz] * mem_pr_x[ipmlx][ipmlz] + apml_half_x[ipmlx][ipmlz] * d_pr_x ;
	}
}
#pragma omp parallel for
for (Myint ix = ix1; ix < ix2; ix++)
{
	Myint ipmlx = ix - offset_x ;
#pragma ivdep       
	for (Myint iz = iz1; iz < iz2; iz++) 
	{
		Myint ipmlz = iz ;
		Myfloat d_mem_pr_x = D_X(mem_pr_x, ipmlz, ipmlx+1) ;
		Myfloat d2_pr_x = D2_X(prc, iz, ix) ;
		mem_vx_x[ipmlx][ipmlz] = bpml_x[ipmlx][ipmlz] * mem_vx_x[ipmlx][ipmlz] + apml_x[ipmlx][ipmlz] * coef2 * (d2_pr_x + d_mem_pr_x) ;
		prn[ix][iz] += coef12[ix][iz] * (d_mem_pr_x + invcoef2 * mem_vx_x[ipmlx][ipmlz]) ;
	}
}
nb_op_bound += (iz2-iz1) * (ix2-ix1) * (2*nb_op_d + nb_op_d2 + 12) ;
}

//--------------------------------
// free surface and rigid surface
//--------------------------------
{
	Myint ii_sign, pr_sign, x_axe, z_axe ;
	Myint ixmin, ixmax, izmin, izmax ;
	bool switch1, switch2 ;

	// loop on boundaries
	for (Myint ib=1; ib<=nBoundary; ib++)
	{
		switch1 = false ;
		switch2 = false ;

		//  set corect index
		if (pBoundary[ib-1]->get_edge() == ZBEG)
		{
			izmin   = izBeg2 ;
			izmax   = izBeg2 ;
			ixmin   = ixBeg ;
			ixmax   = ixEnd ;
			ii_sign = 1 ;
			x_axe   = 0 ;
			z_axe   = 1 ;
			switch1 = true ;
		}
		else if (pBoundary[ib-1]->get_edge() == ZEND)
		{
			izmin   = izEnd2-1 ;
			izmax   = izEnd2-1 ;
			ixmin   = ixBeg ;
			ixmax   = ixEnd ;
			ii_sign = -1 ;
			x_axe   = 0 ;
			z_axe   = 1 ;
			switch1 = true ;
		}
		else if (pBoundary[ib-1]->get_edge() == XBEG)
		{
			izmin   = izBeg ;
			izmax   = izEnd ;
			ixmin   = ixBeg2 ;
			ixmax   = ixBeg2 ;
			ii_sign = 1 ;
			x_axe   = 1 ;
			z_axe   = 0 ;
			switch1 = true ;
		}
		else if (pBoundary[ib-1]->get_edge() == XEND)
		{
			izmin   = izBeg ;
			izmax   = izEnd ;
			ixmin   = ixEnd2-1 ;
			ixmax   = ixEnd2-1 ;
			ii_sign = -1 ;
			x_axe   = 1 ;
			z_axe   = 0 ;
			switch1 = true ;
		}

		// free surface with image method
		// anti-symetry of pressure wavefield
		if (pBoundary[ib-1]->get_type() == FREESURF)
		{
			pr_sign = -1.0 ;
			switch2 = true ;
		}

		// free surface with image method
		// anti-symetry of pressure wavefield
		else if (pBoundary[ib-1]->get_type() == RIGID)
		{
			pr_sign = 1.0 ;
			switch2 = true ;
		}

		// apply boundary condition
		if (switch1 && switch2)
		{
#pragma omp parallel for        
			for (Myint ix = ixmin; ix <= ixmax; ix++)
			{
				for (Myint iz = izmin; iz <= izmax; iz++)
				{
					if (pBoundary[ib-1]->get_type() == FREESURF) prn[ix][iz] = 0. ;
#pragma ivdep
					for (Myint ii = 1; ii <= lstencil; ii++)
					{
						prn[ix-x_axe*(ii_sign*ii)][iz-z_axe*(ii_sign*ii)]
												   = pr_sign * prn[ix+x_axe*(ii_sign*ii)][iz+z_axe*(ii_sign*ii)] ;
					}
				}
			}
		}
	}
}
