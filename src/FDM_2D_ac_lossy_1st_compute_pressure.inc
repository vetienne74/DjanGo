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

Myfloat** const coef12 = coef1->pArray ;

Myint izMin = izBeg1 ;
Myint izMax = izEnd1 ;
Myint ixMin = ixBeg1 ;
Myint ixMax = ixEnd1 ;

// loss terms
const Myfloat sigma_0 = loss1->get_max() ;

//---------------
// medium points
//---------------

// if (layer_type == SPONGE)
//   {
//     Myfloat ** const sponge_coef_pr2 = sponge_coef_pr->pArray ; 
// #pragma omp parallel for        
//     for (Myint ix = ixMin; ix < ixMax; ix++)
//       {
// #pragma ivdep
// 	for (Myint iz = izMin; iz < izMax; iz++)
// 	  {
// 	    Myfloat d_vx_x = D_X(vx2, iz, ix) ; 
// 	    Myfloat d_vz_z = D_Z(vz2, iz, ix) ; 

// 	    pr2[ix][iz] = sponge_coef_pr2[ix][iz] * (pr2[ix][iz] + (coef12[ix][iz] * (d_vx_x + d_vz_z))) ;
// 	  }
//       }
//     nb_op_kernel += (izMax-izMin) * (ixMax-ixMin) * (2*nb_op_d + 4) ;
//   }
//  else
//    {
#pragma omp parallel for        
for (Myint ix = ixMin; ix < ixMax; ix++)
{
#pragma ivdep
	for (Myint iz = izMin; iz < izMax; iz++)
	{
		Myfloat term1 = Myfloat(1.0) + (coef12[ix][iz] * sigma_0 / Myfloat(2.0)) ;
		Myfloat term2 = Myfloat(1.0) - (coef12[ix][iz] * sigma_0 / Myfloat(2.0)) ;
		Myfloat d_vx_x = D_X(vx2, iz, ix) ;
		Myfloat d_vz_z = D_Z(vz2, iz, ix) ;
		pr2[ix][iz] = (coef12[ix][iz] * (d_vx_x + d_vz_z) + term2 * prp2[ix][iz]) / term1 ;
	}
}
nb_op_kernel += (izMax-izMin) * (ixMax-ixMin) * (2*nb_op_d + 11) ;
ncell       += (izMax-izMin) * (ixMax-ixMin) ;

//--------------------------------
// sponge
//--------------------------------
if (is_there_boundary_type(SPG))
{
	Myfloat ** const sponge_coef_pr2 = sponge_coef_pr->pArray ;
#pragma omp parallel for        
	for (Myint ix = ixMin; ix < ixMax; ix++)
	{
#pragma ivdep
		for (Myint iz = izMin; iz < izMax; iz++)
		{
			pr2[ix][iz] *= sponge_coef_pr2[ix][iz] ;
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
	Myfloat ** const mem_vz_z = mem_vz_zBeg->pArray ;
	Myfloat ** const apml_z   = apml_zBeg->pArray ;
	Myfloat ** const bpml_z   = bpml_zBeg->pArray ;

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
			pr2[ix][iz] += coef12[ix][iz] *  mem_vz_z[ipmlx][ipmlz] ;
		}
	}
	nb_op_bound += (iz2-iz1) * (ix2-ix1) * (nb_op_d + 5) ;
}

// z+ layer
if (get_boundary_type(ZEND) == PML)
{
	// for vectorization with Intel compiler
	Myfloat ** const mem_vz_z = mem_vz_zEnd->pArray ;
	Myfloat ** const apml_z   = apml_zEnd->pArray ;
	Myfloat ** const bpml_z   = bpml_zEnd->pArray ;

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
			pr2[ix][iz] += coef12[ix][iz] *  mem_vz_z[ipmlx][ipmlz] ;
		}
	}
	nb_op_bound += (iz2-iz1) * (ix2-ix1) * (nb_op_d + 5) ;
}

// x- layer
if (get_boundary_type(XBEG) == PML)
{
	// for vectorization with Intel compiler
	Myfloat ** const mem_vx_x = mem_vx_xBeg->pArray ;
	Myfloat ** const apml_x   = apml_xBeg->pArray ;
	Myfloat ** const bpml_x   = bpml_xBeg->pArray ;

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
			pr2[ix][iz] += coef12[ix][iz] *  mem_vx_x[ipmlx][ipmlz] ;
		}
	}
	nb_op_bound += (iz2-iz1) * (ix2-ix1) * (nb_op_d + 5) ;
}

// x+ layer
if (get_boundary_type(XEND) == PML)
{
	// for vectorization with Intel compiler
	Myfloat ** const mem_vx_x = mem_vx_xEnd->pArray ;
	Myfloat ** const apml_x   = apml_xEnd->pArray ;
	Myfloat ** const bpml_x   = bpml_xEnd->pArray ;

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
			pr2[ix][iz] += coef12[ix][iz] *  mem_vx_x[ipmlx][ipmlz] ;
		}
	}
	nb_op_bound += (iz2-iz1) * (ix2-ix1) * (nb_op_d + 5) ;
}

//--------------------------------
// free surface and rigid surface
//--------------------------------

// loop on boundaries
Myint ii_sign, pr_sign, x_axe, z_axe ;
Myint ixmin, ixmax, izmin, izmax ;
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
				if (pBoundary[ib-1]->get_type() == FREESURF) pr2[ix][iz] = 0. ;
#pragma ivdep
				for (Myint ii = 1; ii <= lstencil; ii++)
				{
					pr2[ix-x_axe*(ii_sign*ii)][iz-z_axe*(ii_sign*ii)]
											   = pr_sign * pr2[ix+x_axe*(ii_sign*ii)][iz+z_axe*(ii_sign*ii)] ;
				}
			}
		}
	}
}

