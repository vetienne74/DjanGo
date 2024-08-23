#ifndef DJANGO_VECMAT_H_
#define DJANGO_VECMAT_H_

void full2sparse(Grid_1D_int **Mat_i, Grid_1D_int **Mat_j, Grid_1D_float **Mat_c, MatDoub &Mat)
{
	// 1st step, count non-zero values in full matrix
	Int nn = 0 ;
	for (Int ii=0; ii< Mat.nrows(); ii++)
	{
		for (Int jj=0; jj< Mat.ncols(); jj++)
		{
			if (Mat[ii][jj] != 0.0) nn++ ;
		}
	}

	// allocate arrays for sparse matrix
	*Mat_i = new Grid_1D_int(nn, 1) ;
	*Mat_j = new Grid_1D_int(nn, 1) ;
	*Mat_c = new Grid_1D_float(nn, 1) ;

	// 2nd step, store non-zero values in sparse matrix
	Int in = 0 ;
	for (Int ii=0; ii< Mat.nrows(); ii++)
	{
		for (Int jj=0; jj< Mat.ncols(); jj++)
		{
			if (Mat[ii][jj] != 0.0) {
				(*Mat_i)->pArray[in] = ii ;
				(*Mat_j)->pArray[in] = jj ;
				(*Mat_c)->pArray[in] = Mat[ii][jj] ;
				in++ ;
			}
		}
	}
}
#endif
