//===========================================================================
// *** this is an include source file, it can not be compiled on its own ***
//===========================================================================

// get pointers to variables
Grid_1D_float *prc_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(PRC))->get_grid() ;
Myfloat * const prc = prc_grid->pArray ;

Grid_1D_float *prn_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(PRN))->get_grid() ;
Myfloat * const prn = prn_grid->pArray ;

Grid_1D_float *prp_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(PRP))->get_grid() ;
Myfloat * const prp = prp_grid->pArray ;  

Myfloat * const coef12 = coef1->pArray ;

Myint izMin = izBeg1 ;
Myint izMax = izEnd1 ;

//---------------
// medium points
//---------------

// loss terms
const Myfloat sigma_0 = loss1->get_max() ;
const Myfloat sigma_l = loss2->get_max() ;

{  
	// update medium points
#pragma omp parallel for 
#pragma ivdep  
	for (Myint iz = izMin; iz < izMax; iz++)
	{
		Myfloat d2_prc_z = D2_Z(prc, iz) ;
		Myfloat term1 = d2_prc_z ;

		Myfloat term2   = sigma_0 / Myfloat(2.0) / dt * prp[iz] ;

		Myfloat d2_prp_z = D2_Z(prp, iz) ;
		Myfloat term3   = sigma_l / dt * (d2_prc_z - d2_prp_z) ;

		Myfloat term4 = Myfloat(1.0) + (coef12[iz] * sigma_0 / ((Myfloat(2.0) * dt))) ;

		prn[iz] = (coef12[iz] * (term1 + term2 + term3) + Myfloat(2.0) * prc[iz] - prp[iz])
				/ term4 ;

	}
	nb_op_kernel += (izMax-izMin) * (2*nb_op_d2 + 17) ;
	ncell       += (izMax-izMin) ;
}

//------------------------------------
// free surface with image method
// anti-symetry of pressure wavefield
//------------------------------------

if (get_boundary_type(ZBEG) == FREESURF)
{
	Myint iz0 = izBeg1 ;
	prn[iz0] = 0. ;
#pragma ivdep
	for (Myint ii = 1; ii <= lstencil; ii++)
	{
		prn[iz0-ii] = -prn[iz0+ii] ;
	}
}

if (get_boundary_type(ZEND) == FREESURF)
{
	Myint iz0 = izEnd1-1 ;
	prn[iz0] = 0. ;
#pragma ivdep
	for (Myint ii = 1; ii <= lstencil; ii++)
	{
		prn[iz0+ii] = -prn[iz0-ii] ;
	}
}

//------------------------------------
// rigid surface with image method
// symetry of pressure wavefield
//------------------------------------

if (get_boundary_type(ZBEG) == RIGID)
{
	Myint iz0 = izBeg1 ;
#pragma ivdep
	for (Myint ii = 1; ii <= lstencil; ii++)
	{
		prn[iz0-ii] = prn[iz0+ii] ;
	}
}

if (get_boundary_type(ZEND) == RIGID)
{
	Myint iz0 = izEnd1-1 ;
#pragma ivdep
	for (Myint ii = 1; ii <= lstencil; ii++)
	{
		prn[iz0+ii] = prn[iz0-ii] ;
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
	Myfloat * const mem_vz_z    = mem_vz_zBeg->pArray ;
	Myfloat * const mem_pr_z    = mem_pr_zBeg->pArray ;
	Myfloat * const apml_z      = apml_zBeg->pArray ;
	Myfloat * const bpml_z      = bpml_zBeg->pArray ;
	Myfloat * const apml_half_z = apml_half_zBeg->pArray ;
	Myfloat * const bpml_half_z = bpml_half_zBeg->pArray ;

	const Myfloat coef2 = dt / (TEMP_RHO_CONST * dz) ;
	const Myfloat invcoef2 = 1/coef2 ;

#pragma omp parallel for 
#pragma ivdep
	for (Myint iz = iz1; iz < iz2; iz++)
	{
		Myint ipml = iz ;
		Myfloat d_pr_z = D_Z(prc, iz+1) ;
		mem_pr_z[ipml] = bpml_half_z[ipml] * mem_pr_z[ipml] + apml_half_z[ipml] * d_pr_z ;
	}
#pragma omp parallel for 
#pragma ivdep
	for (Myint iz =  iz1; iz < iz2; iz++)
	{
		Myint ipml = iz ;
		Myfloat d_mem_pr_z = D_Z(mem_pr_z, ipml) ;
		Myfloat d2_pr_z = D2_Z(prc, iz) ;
		mem_vz_z[ipml] = bpml_z[ipml] * mem_vz_z[ipml] + apml_z[ipml] * coef2 * (d2_pr_z + d_mem_pr_z) ;
		prn[iz] += coef12[iz] * (d_mem_pr_z + invcoef2 * mem_vz_z[ipml]) ;
	}
	nb_op_bound += (iz2-iz1) * (2*nb_op_d + nb_op_d2 + 12) ;
}

// layer z+
if (get_boundary_type(ZEND) == PML)
{
	Myint iz1 = izEnd2 ;
	Myint iz2 = izMax ;

	// for vectorization with Intel compiler
	Myfloat * const mem_vz_z    = mem_vz_zEnd->pArray ;
	Myfloat * const mem_pr_z    = mem_pr_zEnd->pArray ;
	Myfloat * const apml_z      = apml_zEnd->pArray ;
	Myfloat * const bpml_z      = bpml_zEnd->pArray ;
	Myfloat * const apml_half_z = apml_half_zEnd->pArray ;
	Myfloat * const bpml_half_z = bpml_half_zEnd->pArray ;

	const Myfloat coef2 = dt / (TEMP_RHO_CONST * dz) ;
	const Myfloat invcoef2 = 1/coef2 ;
	Myint const offset_z = nlayer_zBeg + nz ;

#pragma omp parallel for 
#pragma ivdep
	for (Myint iz = iz1; iz < iz2; iz++)
	{
		Myint ipml = iz - offset_z ;
		Myfloat d_pr_z = D_Z(prc, iz) ;
		mem_pr_z[ipml] = bpml_half_z[ipml] * mem_pr_z[ipml] + apml_half_z[ipml] * d_pr_z ;
	}
#pragma omp parallel for 
#pragma ivdep
	for (Myint iz =  iz1; iz < iz2; iz++)
	{
		Myint ipml = iz - offset_z ;
		Myfloat d_mem_pr_z = D_Z(mem_pr_z, ipml+1) ;
		Myfloat d2_pr_z = D2_Z(prc, iz) ;
		mem_vz_z[ipml] = bpml_z[ipml] * mem_vz_z[ipml] + apml_z[ipml] * coef2 * (d2_pr_z + d_mem_pr_z) ;
		prn[iz] += coef12[iz] * (d_mem_pr_z + invcoef2 * mem_vz_z[ipml]) ;
	}
	nb_op_bound += (iz2-iz1) * (2*nb_op_d + nb_op_d2 + 12) ;
}

