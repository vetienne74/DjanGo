#ifndef DJANGO_LAGRANGE_H_
#define DJANGO_LAGRANGE_H_

inline Doub lagran(VecDoub_I &xnode, Int ia, Doub xi)
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Computes the Lagrange polynomial associated to node ia at coordinate xi in
// the reference element [-1 1]
// Nodes of the polynomials are given by input vector xnode
//
// Ref: Hughes, The finite element method, page 127, eq. 3.6.1
//
// Input
//   xnode nodes of the Lagrange polynomial
//   ia node index associated to the Lagrange polynomial
//   xi coordinates for the evaluation of the Lagrange polynomial
//
// Ouput
//   returns the value of the Lagrange polynomial
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
{
	// number of points should be >=1
	Int nn = xnode.size() ;
	if (nn <= 0)
	{
		throw("number of points should be >=1 in lagran");
	}

	// check node index ia
	if ((ia < 0) || (ia > nn))
	{
		throw("invalid node index ia in lagran");
	}

	// compute denominator of Lagrange polynomials
	Doub den_phi = 1.0 ;
	for (Int ib=0; ib<nn; ib++)
	{
		if (ib == ia) continue ;
		den_phi = den_phi * (xnode[ia] - xnode[ib]) ;
	}

	// compute phi(ia) at xg(ig)
	Doub phi = 1.0 ;
	for (Int ib=0; ib<nn; ib++)
	{
		if (ib == ia) continue ;
		phi = phi * (xi - xnode[ib]) ;
	}
	phi = phi / den_phi ;

	return(phi) ;
}

Doub dlagran(VecDoub_I &xnode, Int ia, Doub xi)
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Computes the DERIVATIVE of Lagrange polynomial associated to node ia
// at coordinate xi in // the reference element [-1 1]
// Nodes of the polynomials are given by input vector xnode
//
// Input
//   xnode nodes of the Lagrange polynomial
//   ia node index associated to the Lagrange polynomial
//   xi coordinates for the evaluation of the DERIVATIVE of Lagrange polynomial
//
// Ouput
//   returns the value of the DERIVATIVE of the Lagrange polynomial
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
{
	// number of points should be >=1
	Int nn = xnode.size() ;
	if (nn <= 0)
	{
		throw("number of points should be >=1 in dlagran");
	}

	// check node index ia
	if ((ia < 0) || (ia > nn))
	{
		throw("invalid node index ia in dlagran");
	}

	// compute denominator of Lagrange polynomials
	Doub den_phi = 1.0 ;
	for (Int ib=0; ib<nn; ib++)
	{
		if (ib == ia) continue ;
		den_phi = den_phi * (xnode[ia] - xnode[ib]) ;
	}

	// compute phi(ia) and d_phi(ia) at xi
	Doub phi = 1.0 ;
	Doub d_phi = 0.0 ;
	for (Int ib=0; ib<nn; ib++)
	{
		if (ib == ia) continue ;

		// derivative d_phi(ia)
		d_phi = phi + d_phi * (xi - xnode[ib]) ;
		phi = phi * (xi - xnode[ib]) ;
	}
	phi = phi / den_phi ;
	d_phi = d_phi / den_phi ;

	return(d_phi) ;
}
#endif
