//===========================================================================
// *** this is an include source file, it can not be compiled on its own ***
//===========================================================================

// for vectorization with Intel compiler
Grid_1D_float *pr_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(PR))->get_grid() ;
Myfloat * const pr = pr_grid->pArray ;

Grid_1D_float *prp_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(PRP))->get_grid() ;
Myfloat * const prp = prp_grid->pArray ;

Grid_1D_float *vz_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(VZ))->get_grid() ;
Myfloat * const vz = vz_grid->pArray ;

Myfloat * const coef12 = coef1->pArray ;

Myint izMin = izBeg1 ;
Myint izMax = izEnd1-1 ;

// loss terms
const Myfloat sigma_l = loss2->get_max() ;

//---------------
// medium points
//---------------
#pragma omp parallel for 
#pragma ivdep
for (Myint iz = izMin; iz < izMax; iz++)
{
	Myfloat d_pr_z  = D_Z(pr,  iz+1) ;
	Myfloat d_prp_z = D_Z(prp, iz+1) ;
	Myfloat term1 = sigma_l * (d_pr_z - d_prp_z) / dt ;
	vz[iz] += coef2 * (d_pr_z + term1) ;
}
nb_op_kernel += (izMax-izMin) * (2*nb_op_d + 6) ;

//--------------------------------
// sponge
//--------------------------------
if (is_there_boundary_type(SPG))
{
	Myfloat * const sponge_coef_vz2 = sponge_coef_vz->pArray ;
#pragma omp parallel for 
#pragma ivdep
	for (Myint iz = izMin; iz < izMax; iz++)
	{
		vz[iz] *= sponge_coef_vz2[iz] ;
	}
	nb_op_kernel += (izMax-izMin) ;
}

//--------------------------------
// free surface and rigid surface
//--------------------------------

// loop on boundaries
Myint iz0, z_sign, z_offset1, z_offset2 ;
Myfloat vz_sign ;
bool switch1, switch2 ;
for (Myint ib=1; ib<=nBoundary; ib++)
{
	switch1 = false ;
	switch2 = false ;

	//  set corect index
	if (pBoundary[ib-1]->get_edge() == ZBEG)
	{
		iz0     = izBeg1 ;
		z_sign  = 1 ;
		z_offset1 = 0 ;
		z_offset2 = -1 ;
		switch1 = true ;
	}
	else if (pBoundary[ib-1]->get_edge() == ZEND)
	{
		iz0    = izEnd1-1 ;
		z_sign = -1 ;
		z_offset1 = -1 ;
		z_offset2 = 0 ;
		switch1 = true ;
	}

	// free surface with image method
	// symetry of velocity wavefield
	if (pBoundary[ib-1]->get_type() == FREESURF)
	{
		vz_sign = 1.0 ;
		switch2 = true ;
	}

	// free surface with image method
	// anti-symetry of pressure wavefield
	else if (pBoundary[ib-1]->get_type() == RIGID)
	{
		vz_sign = -1.0 ;
		switch2 = true ;
	}

	// apply boundary condition
	if (switch1 && switch2)
	{
#pragma ivdep
		for (Myint ii = 1; ii <= lstencil; ii++)
		{
			vz[iz0-(z_sign*ii)+z_offset1] = vz_sign * vz[iz0+(z_sign*ii)+z_offset2] ;
		}
	}
}

//---------------
// CPML points
//---------------

// layer z-
if (get_boundary_type(ZBEG) == PML)
{
	Myint iz1 = izMin ;
	Myint iz2 = izBeg2 ;

	// for vectorization with Intel compiler
	Myfloat * mem_pr_z    = mem_pr_zBeg->pArray ;
	Myfloat * apml_half_z = apml_half_zBeg->pArray ;
	Myfloat * bpml_half_z = bpml_half_zBeg->pArray ;
#pragma omp parallel for 
#pragma ivdep
	for (Myint iz = iz1; iz < iz2; iz++)
	{
		Myint ipml = iz ;
		Myfloat d_pr_z = D_Z(pr, iz+1) ;
		mem_pr_z[ipml] = bpml_half_z[ipml] * mem_pr_z[ipml] + apml_half_z[ipml] * d_pr_z ;
		vz[iz] += coef2 * mem_pr_z[ipml] ;
	}
	nb_op_bound += (iz2-iz1) * (nb_op_d + 5) ;
}

// layer z+
if (get_boundary_type(ZEND) == PML)
{
	Myint iz1 = izEnd2-1 ;
	Myint iz2 = izMax ;
	// for vectorization with Intel compiler
	Myfloat * mem_pr_z    = mem_pr_zEnd->pArray ;
	Myfloat * apml_half_z = apml_half_zEnd->pArray ;
	Myfloat * bpml_half_z = bpml_half_zEnd->pArray ;
	Myint const offset_z = nlayer_zBeg + nz - 1;
#pragma omp parallel for 
#pragma ivdep
	for (Myint iz = iz1; iz < iz2; iz++)
	{
		Myint ipml = iz - offset_z ;
		Myfloat d_pr_z = D_Z(pr, iz+1) ;
		mem_pr_z[ipml] = bpml_half_z[ipml] * mem_pr_z[ipml] + apml_half_z[ipml] * d_pr_z ;
		vz[iz] += coef2 * mem_pr_z[ipml] ;
	}
	nb_op_bound += (iz2-iz1) * (nb_op_d + 5) ;
}

