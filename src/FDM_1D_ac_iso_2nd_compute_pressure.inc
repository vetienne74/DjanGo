//===========================================================================
// *** this is an include source file, it can not be compiled on its own ***
//===========================================================================

// get pointers to variables
Grid_1D_float *prc_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(PRC))->get_grid() ;
Myfloat * const prc = prc_grid->pArray ;

Grid_1D_float *prn_grid = (Grid_1D_float*) (Singleton::Instance()->get_variable(PRN))->get_grid() ;
Myfloat * const prn = prn_grid->pArray ;    

Myfloat * const coef12 = coef1->pArray ;
//Myfloat * const sponge_coef_pr_z2 = sponge_coef_pr_z->pArray ;

Myint izMin = izBeg1 ;
Myint izMax = izEnd1 ;

//---------------
// medium points
//---------------

// if (layer_type == SPONGE)
//   {
//     Myfloat * const sponge_coef_pr_z2 = sponge_coef_pr_z->pArray ;
// #pragma omp parallel for 
// #pragma ivdep
//     for (Myint iz = izMin; iz < izMax; iz++)
//       {
// 	Myfloat d2_pr_z = sponge_coef_pr_z2[iz] * D2_Z(prc, iz) ;
// 	prn[iz] = sponge_coef_pr_z2[iz] * ((coef12[iz] * d2_pr_z) - prn[iz] + 2 * prc[iz]) ;
//       }	
//     nb_op_kernel += (izMax-izMin) * (nb_op_d2 + 5) ;
//   }
//  else
{  
#pragma omp parallel for 
#pragma ivdep  
	for (Myint iz = izMin; iz < izMax; iz++)
	{
		Myfloat d2_pr_z = D2_Z(prc, iz) ;
		prn[iz] = (coef12[iz] * d2_pr_z) - prn[iz] + Myfloat(2.0) * prc[iz] ;
	}
	nb_op_kernel += (izMax-izMin) * (nb_op_d2 + 4) ;
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

