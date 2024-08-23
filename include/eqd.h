#ifndef DJANGO_EQD_H_
#define DJANGO_EQD_H_

void eqd(VecDoub_O &x)
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// compute coordinates of equidistant points in ref element [-1 1]
// number of points equals the length of the input vector
// input/output VecDoub_O &x
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
{
  // number of points should be >=1
  Int ni = x.size() ;
  if (ni <= 0)
    {
      throw("number of point should be >=1 in eqd");
    }

  else if (ni == 1)
    {
      x[0] = 0.0 ;
    }
  else
    {
      Doub x2 = 1.0 ;
      Doub x1= -1.0 ;
      Doub hnode = (x2-x1) / (ni-1) ; 
      for (Int ii=0; ii<ni; ii++)
	{
	  x[ii] = x1 + ii*hnode ;            
	}
    }
}
#endif
