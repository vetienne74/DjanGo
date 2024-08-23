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
Myint izMax = izEnd1 ;

// loss terms
const Myfloat sigma_0 = loss1->get_max() ;

//---------------
// medium points
//---------------
#pragma omp parallel for 
#pragma ivdep
for (Myint iz = izMin; iz < izMax; iz++)
{
	Myfloat d_vz_z = D_Z(vz, iz) ;
	Myfloat term1 = Myfloat(1.0) + (coef12[iz] * sigma_0 / Myfloat(2.0)) ;
	Myfloat term2 = Myfloat(1.0) - (coef12[iz] * sigma_0 / Myfloat(2.0)) ;
	pr[iz] = (coef12[iz] * d_vz_z + term2 * prp[iz]) / term1 ;
}
nb_op_kernel += (izMax-izMin) * (nb_op_d + 10) ;
ncell       += (izMax-izMin) ;

//--------------------------------
// sponge
//--------------------------------
if (is_there_boundary_type(SPG))
{
	Myfloat * const sponge_coef_pr2 = sponge_coef_pr->pArray ;
#pragma omp parallel for 
#pragma ivdep
	for (Myint iz = izMin; iz < izMax; iz++)
	{
		pr[iz] *= sponge_coef_pr2[iz] ;
	}
	nb_op_kernel += (izMax-izMin) ;
}

//--------------------------------
// free surface and rigid surface
//--------------------------------

// loop on boundaries
Myint iz0, z_sign ;
Myfloat pr_sign ;
bool switch1, switch2 ;
for (Myint ib=1; ib<=nBoundary; ib++)
{
	switch1 = false ;
	switch2 = false ;

	//  set corect index
	if (pBoundary[ib-1]->get_edge() == ZBEG)
	{
		iz0     = izBeg2 ;
		z_sign  = 1 ;
		switch1 = true ;
	}
	else if (pBoundary[ib-1]->get_edge() == ZEND)
	{
		iz0    = izEnd2-1 ;
		z_sign = -1 ;
		switch1 = true ;
	}

	// free surface with image method
	// anti-symetry of pressure wavefield
	if (pBoundary[ib-1]->get_type() == FREESURF)
	{
		pr_sign = -1.0 ;
		pr[iz0] = 0.0 ;
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
#pragma ivdep
		for (Myint ii = 1; ii <= lstencil; ii++)
		{
			pr[iz0-(z_sign*ii)] = pr_sign * pr[iz0+(z_sign*ii)] ;
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
	Myfloat * mem_vz_z = mem_vz_zBeg->pArray ;
	Myfloat * apml_z   = apml_zBeg->pArray ;
	Myfloat * bpml_z   = bpml_zBeg->pArray ;
#pragma omp parallel for
#pragma ivdep
	for (Myint iz = iz1; iz < iz2; iz++)
	{
		Myint ipml = iz ;
		Myfloat d_vz_z = D_Z(vz, iz) ;
		mem_vz_z[ipml] = bpml_z[ipml] * mem_vz_z[ipml] + apml_z[ipml] * d_vz_z ;
		pr[iz] += coef12[iz] * mem_vz_z[ipml] ;
	}
	nb_op_bound += (iz2-iz1) * (nb_op_d + 5) ;
}

// layer z+
if (get_boundary_type(ZEND) == PML)
{
	Myint iz1 = izEnd2 ;
	Myint iz2 = izMax ;

	// for vectorization with Intel compiler
	Myfloat * mem_vz_z = mem_vz_zEnd->pArray ;
	Myfloat * apml_z   = apml_zEnd->pArray ;
	Myfloat * bpml_z   = bpml_zEnd->pArray ;
	Myint const offset_z = nlayer_zBeg + nz ;
#pragma omp parallel for
#pragma ivdep 
	for (Myint iz = izEnd2; iz < iz2; iz++)
	{
		Myint ipml = iz - offset_z ;
		Myfloat d_vz_z = D_Z(vz, iz) ;
		mem_vz_z[ipml] = bpml_z[ipml] * mem_vz_z[ipml] + apml_z[ipml] * d_vz_z ;
		pr[iz] += coef12[iz] * mem_vz_z[ipml] ;
	}
	nb_op_bound += (iz2-iz1) * (nb_op_d + 5) ;
}

