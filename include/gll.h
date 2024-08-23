#ifndef DJANGO_GLL_H_
#define DJANGO_GLL_H_

void gll(VecDoub_O &x, VecDoub_O &w)
{  	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%
	//% follows Numerical Recipes stype
	//% code adapted from below
	//%
	//% lglnodes.m
	//%
	//% Computes the Legendre-Gauss-Lobatto nodes, weights and the LGL Vandermonde
	//% matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x). Useful for numerical
	//% integration and spectral methods.
	//%
	//% Reference on LGL nodes and weights:
	//%   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
	//%   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
	//%
	//% Written by Greg von Winckel - 04/17/2004
	//% Contact: gregvw@chtm.unm.edu
	//%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	// set precision
	const Doub EPS  = 1.0e-14;
	const Int MAXIT = 1000 ;

	// number of points should be >=2
	if (x.size() <= 1)
	{
		throw("number of point should be >=2 in gll");
	}

	// polynomial order
	Int N = x.size() - 1 ;

	//% Truncation + 1
	Int N1=N+1;

	//% Use the Chebyshev-Gauss-Lobatto nodes as the first guess
	for (Int ii=0; ii<x.size(); ii++)
	{
		x[ii] = cos(M_PI*Doub(ii)/Doub(N)) ;
	}

	//% The Legendre Vandermonde Matrix
	MatDoub P(N1, N1) ;

	//% Compute P_(N) using the recursion relation
	//% Compute its first and second derivatives and
	//% update x using the Newton-Raphson method.

	VecDoub xold(N1) ;
	for (Int ii=0; ii<x.size(); ii++)
	{
		xold[ii] = 2.0 ;
	}

	Doub max_x ;
	Doub eps = 1.e-14 ;
	max_x = 0 ;
	for (Int ii=0; ii<x.size(); ii++)
	{
		//max_x = max(max_x, abs(xold[ii] - x[ii])) ;
		max_x = (max_x > abs(xold[ii] - x[ii])) ? max_x : abs(xold[ii] - x[ii]) ;
	}

	Int MAX_ITER = 1000 ;
	Int niter = 0 ;
	while ((max_x > eps) && (niter < MAXIT))
	{
		for (Int ii=0; ii<x.size(); ii++)
		{
			xold[ii] = x[ii] ;
		}

		for (Int ii=0; ii<P.nrows(); ii++)
		{
			P[ii][0] = 1.0 ;
			P[ii][1] = x[ii] ;
		}

		for (Int k=2; k<=N; k++)
		{
			for (Int ii=0; ii<P.nrows(); ii++)
			{
				P[ii][k] = ( Doub(2*k-1) * x[ii] * P[ii][k-1] - Doub(k-1) * P[ii][k-2]) / Doub(k) ;
			}
		}

		for (Int ii=0; ii<x.size(); ii++)
		{
			x[ii] = xold[ii] -( (x[ii] * P[ii][N1-1]) - P[ii][N-1] ) / ( N1 * P[ii][N1-1] );
		}

		niter ++ ;

		max_x = 0 ;
		for (Int ii=0; ii<x.size(); ii++)
		{
			//max_x = max(max_x, abs(xold[ii] - x[ii])) ;
			max_x = (max_x > abs(xold[ii] - x[ii])) ? max_x : abs(xold[ii] - x[ii]) ;
		}
	}

	if (niter > MAXIT) throw("too many iterations in gll");

	for (Int ii=0; ii<P.nrows(); ii++)
	{
		w[ii] = 2. / (N*N1 * P[ii][N1-1] * P[ii][N1-1] );
	}

	// swap
	VecDoub tmp(N1) ;
	for (Int ii=0; ii<x.size(); ii++)
	{
		Int itmp = x.size()-ii-1 ;
		tmp[itmp] = w[ii] ;
	}
	for (Int ii=0; ii<x.size(); ii++)
	{
		w[ii] = tmp[ii] ;
	}
	for (Int ii=0; ii<x.size(); ii++)
	{
		Int itmp = x.size()-ii-1 ;
		tmp[itmp] = x[ii] ;
	}
	for (Int ii=0; ii<x.size(); ii++)
	{
		x[ii] = tmp[ii] ;
	}

}

void gll(VecDoub_O &x)
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// just get GLL point coordinates
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
{
	// number of points should be >=2
	if (x.size() <= 1)
	{
		throw("number of point should be >=2 in gll");
	}
	VecDoub dummy(x.size()) ;
	gll(x, dummy) ;
}

#endif
