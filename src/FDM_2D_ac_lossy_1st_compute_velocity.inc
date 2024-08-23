//===========================================================================
// *** this is an include source file, it can not be compiled on its own ***
//===========================================================================

// for vectorization with Intel compiler
Grid_2D_float *pr_grid = (Grid_2D_float*) (Singleton::Instance()->get_variable(PR))->get_grid() ;
Myfloat** const pr2 = pr_grid->pArray ;

Grid_2D_float *prp_grid = (Grid_2D_float*) (Singleton::Instance()->get_variable(PRP))->get_grid() ;
Myfloat** const prp2 = prp_grid->pArray ;

Grid_2D_float *vz_grid = (Grid_2D_float*) (Singleton::Instance()->get_variable(VZ))->get_grid() ;
Myfloat** const vz2 = vz_grid->pArray ;

Grid_2D_float *vx_grid = (Grid_2D_float*) (Singleton::Instance()->get_variable(VX))->get_grid() ;
Myfloat** const vx2 = vx_grid->pArray ;

Myfloat ** const coef12 = coef1->pArray ;

Myint izMin = izBeg1 ;
Myint izMax = izEnd1-1 ;
Myint ixMin = ixBeg1 ;
Myint ixMax = ixEnd1-1 ;

// loss terms
const Myfloat sigma_l = loss2->get_max() ;

//---------------
// medium points
//---------------

// if (layer_type == SPONGE)
//   {
//     Myfloat ** const sponge_coef_vz2 = sponge_coef_vz->pArray ;
//     Myfloat ** const sponge_coef_vx2 = sponge_coef_vx->pArray ;

// #pragma omp parallel for 
//     for (Myint ix = ixMin; ix < ixMax; ix++)
//       {
// #pragma ivdep
// 	for (Myint iz = izMin; iz < izMax; iz++)
// 	  {
// 	    Myfloat d_pr_x = D_X(pr2, iz, ix+1) ;
// 	    Myfloat d_pr_z = D_Z(pr2, iz+1, ix) ;

// 	    vx2[ix][iz] = sponge_coef_vx2[ix][iz] * (vx2[ix][iz] + (coef2 * d_pr_x)) ;
// 	    vz2[ix][iz] = sponge_coef_vz2[ix][iz] * (vz2[ix][iz] + (coef2 * d_pr_z)) ;
// 	  }
//       }
//     nb_op_kernel += (izMax-izMin) * (ixMax-ixMin) * (2*nb_op_d + 6) ;
//   }
//  else
//    {
#pragma omp parallel for 
for (Myint ix = ixMin; ix < ixMax; ix++)
{
#pragma ivdep
	for (Myint iz = izMin; iz < izMax; iz++)
	{
		Myfloat d_pr_x = D_X(pr2, iz, ix+1) ;
		Myfloat d_pr_z = D_Z(pr2, iz+1, ix) ;

		Myfloat d_prp_x = D_X(prp2, iz, ix+1) ;
		Myfloat d_prp_z = D_Z(prp2, iz+1, ix) ;

		Myfloat term1 = sigma_l * (d_pr_x - d_prp_x) / dt ;
		vx2[ix][iz] += coef2 * (d_pr_x + term1) ;

		Myfloat term2 = sigma_l * (d_pr_z - d_prp_z) / dt ;
		vz2[ix][iz] += coef2 * (d_pr_z + term2) ;
	}
}
nb_op_kernel += (izMax-izMin) * (ixMax-ixMin) * (4*nb_op_d + 12) ;

//--------------------------------
// sponge
//--------------------------------
if (is_there_boundary_type(SPG))
{
	Myfloat ** const sponge_coef_vz2 = sponge_coef_vz->pArray ;
	Myfloat ** const sponge_coef_vx2 = sponge_coef_vx->pArray ;
	for (Myint ix = ixMin; ix < ixMax; ix++)
	{
#pragma ivdep
		for (Myint iz = izMin; iz < izMax; iz++)
		{
			vx2[ix][iz] *= sponge_coef_vx2[ix][iz] ;
			vz2[ix][iz] *= sponge_coef_vz2[ix][iz] ;
		}
	}
}

//---------------
// CPML points
//---------------

// z- layer
if (get_boundary_type(ZBEG) == PML)
{
	// for vectorization with Intel compiler
	Myfloat ** const mem_pr_z = mem_pr_zBeg->pArray ;
	Myfloat ** const apml_z   = apml_half_zBeg->pArray ;
	Myfloat ** const bpml_z   = bpml_half_zBeg->pArray ;

	Myint iz1 = izBeg1 ;
	Myint iz2 = izBeg2 ;
	Myint ix1 = ixBeg1 ;
	Myint ix2 = ixEnd1-1 ;

#pragma omp parallel for
	for (Myint ix = ix1; ix < ix2; ix++)
	{
		Myint ipmlx = ix ;
#pragma ivdep       
		for (Myint iz = iz1; iz < iz2; iz++)
		{
			Myint ipmlz = iz ;
			Myfloat d_pr_z = D_Z(pr2, iz+1, ix) ;
			mem_pr_z[ipmlx][ipmlz] = bpml_z[ipmlx][ipmlz] * mem_pr_z[ipmlx][ipmlz] + apml_z[ipmlx][ipmlz] * d_pr_z ;
			vz2[ix][iz] += coef2 * mem_pr_z[ipmlx][ipmlz] ;
		}
	}
	nb_op_bound += (iz2-iz1) * (ix2-ix1) * (nb_op_d + 5) ;
}

// z+ layer
if (get_boundary_type(ZEND) == PML)
{
	// for vectorization with Intel compiler
	Myfloat ** const mem_pr_z = mem_pr_zEnd->pArray ;
	Myfloat ** const apml_z   = apml_half_zEnd->pArray ;
	Myfloat ** const bpml_z   = bpml_half_zEnd->pArray ;

	Myint iz1 = izEnd2 ;
	Myint iz2 = izEnd1 ;
	Myint ix1 = ixBeg1 ;
	Myint ix2 = ixEnd1-1 ;
	Myint const offset_z = nlayer_zBeg + nz ;
#pragma omp parallel for
	for (Myint ix = ix1; ix < ix2; ix++)
	{
		Myint ipmlx = ix ;
#pragma ivdep       
		for (Myint iz = iz1; iz < iz2; iz++)
		{
			Myint ipmlz = iz - offset_z ;
			Myfloat d_pr_z = D_Z(pr2, iz, ix) ;
			mem_pr_z[ipmlx][ipmlz] = bpml_z[ipmlx][ipmlz] * mem_pr_z[ipmlx][ipmlz] + apml_z[ipmlx][ipmlz] * d_pr_z ;
			vz2[ix][iz-1] += coef2 * mem_pr_z[ipmlx][ipmlz] ;
		}
	}
	nb_op_bound += (iz2-iz1) * (ix2-ix1) * (nb_op_d + 5) ;
}

// x- layer
if (get_boundary_type(XBEG) == PML)
{
	// for vectorization with Intel compiler
	Myfloat ** const mem_pr_x = mem_pr_xBeg->pArray ;
	Myfloat ** const apml_x   = apml_half_xBeg->pArray ;
	Myfloat ** const bpml_x   = bpml_half_xBeg->pArray ;

	Myint iz1 = izBeg1 ;
	Myint iz2 = izEnd1-1 ;
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
			Myfloat d_pr_x = D_X(pr2, iz, ix+1) ;
			mem_pr_x[ipmlx][ipmlz] = bpml_x[ipmlx][ipmlz] * mem_pr_x[ipmlx][ipmlz] + apml_x[ipmlx][ipmlz] * d_pr_x ;
			vx2[ix][iz] += coef2 * mem_pr_x[ipmlx][ipmlz] ;
		}
	}
	nb_op_bound += (iz2-iz1) * (ix2-ix1) * (nb_op_d + 5) ;
}

// x+ layer
if (get_boundary_type(XEND) == PML)
{
	// for vectorization with Intel compiler
	Myfloat ** const mem_pr_x = mem_pr_xEnd->pArray ;
	Myfloat ** const apml_x   = apml_half_xEnd->pArray ;
	Myfloat ** const bpml_x   = bpml_half_xEnd->pArray ;

	Myint iz1 = izBeg1 ;
	Myint iz2 = izEnd1-1 ;
	Myint ix1 = ixEnd2-1 ;
	Myint ix2 = ixEnd1-1 ;
	Myint const offset_x = nlayer_xBeg + nx -1;
#pragma omp parallel for
	for (Myint ix = ix1; ix < ix2; ix++)
	{
		Myint ipmlx = ix - offset_x ;
#pragma ivdep       
		for (Myint iz = iz1; iz < iz2; iz++)
		{
			Myint ipmlz = iz ;
			Myfloat d_pr_x = D_X(pr2, iz, ix+1) ;
			mem_pr_x[ipmlx][ipmlz] = bpml_x[ipmlx][ipmlz] * mem_pr_x[ipmlx][ipmlz] + apml_x[ipmlx][ipmlz] * d_pr_x ;
			vx2[ix][iz] += coef2 * mem_pr_x[ipmlx][ipmlz] ;
		}
	}
	nb_op_bound += (iz2-iz1) * (ix2-ix1) * (nb_op_d + 5) ;
}

//--------------------------------
// free surface and rigid surface
//--------------------------------

// loop on boundaries
Myint ii_sign, vel_sign, x_axe, z_axe ;
Myint ixmin, ixmax, izmin, izmax ;
Myint x_offset1, x_offset2, z_offset1, z_offset2 ;
bool switch1, switch2 ;
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
		ixmax   = ixEnd-1 ;
		ii_sign = 1 ;
		x_axe   = 0 ;
		z_axe   = 1 ;
		x_offset1 = 0 ;
		x_offset2 = 0 ;
		z_offset1 = 0 ;
		z_offset2 = -1 ;
		switch1 = true ;
	}
	else if (pBoundary[ib-1]->get_edge() == ZEND)
	{
		izmin   = izEnd2-1 ;
		izmax   = izEnd2-1 ;
		ixmin   = ixBeg ;
		ixmax   = ixEnd-1 ;
		ii_sign = -1 ;
		x_axe   = 0 ;
		z_axe   = 1 ;
		x_offset1 = 0 ;
		x_offset2 = 0 ;
		z_offset1 = -1 ;
		z_offset2 = 0 ;
		switch1 = true ;
	}
	else if (pBoundary[ib-1]->get_edge() == XBEG)
	{
		izmin   = izBeg ;
		izmax   = izEnd-1 ;
		ixmin   = ixBeg2 ;
		ixmax   = ixBeg2 ;
		ii_sign = 1 ;
		x_axe   = 1 ;
		z_axe   = 0 ;
		x_offset1 = 0 ;
		x_offset2 = -1 ;
		z_offset1 = 0 ;
		z_offset2 = 0 ;
		switch1 = true ;
	}
	else if (pBoundary[ib-1]->get_edge() == XEND)
	{
		izmin   = izBeg ;
		izmax   = izEnd-1 ;
		ixmin   = ixEnd2-1 ;
		ixmax   = ixEnd2-1 ;
		ii_sign = -1 ;
		x_axe   = 1 ;
		z_axe   = 0 ;
		x_offset1 = -1 ;
		x_offset2 = 0 ;
		z_offset1 = 0 ;
		z_offset2 = 0 ;
		switch1 = true ;
	}

	// free surface with image method
	// anti-symetry of pressure wavefield
	if (pBoundary[ib-1]->get_type() == FREESURF)
	{
		vel_sign = 1.0 ;
		switch2 = true ;
	}

	// free surface with image method
	// anti-symetry of pressure wavefield
	else if (pBoundary[ib-1]->get_type() == RIGID)
	{
		vel_sign = -1.0 ;
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
#pragma ivdep
				for (Myint ii = 1; ii <= lstencil; ii++)
				{
					vx2[ix-x_axe*(ii_sign*ii-x_offset1)][iz-z_axe*(ii_sign*ii-z_offset1)]
														 = vel_sign * vx2[ix+x_axe*(ii_sign*ii+x_offset2)][iz+z_axe*(ii_sign*ii+z_offset2)] ;

					vz2[ix-x_axe*(ii_sign*ii-x_offset1)][iz-z_axe*(ii_sign*ii-z_offset1)]
														 = vel_sign * vz2[ix+x_axe*(ii_sign*ii+x_offset2)][iz+z_axe*(ii_sign*ii+z_offset2)] ;
				}
			}
		}
	}
}

