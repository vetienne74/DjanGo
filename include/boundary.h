#ifndef DJANGO_BOUNDARY_H_
#define DJANGO_BOUNDARY_H_

#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Boundary

{
 public:

  Boundary(Edge_type, Boundary_type, Myint, Myfloat) ;

  void info(void) ;

  Boundary_type get_type(void) ;
  
  Myint get_width(void) ;

  Myfloat get_coef(void) ;

  Edge_type get_edge(void) ;
  
  Rtn_code reset(void) ;

 protected:

  Edge_type     edge ;
  Boundary_type type ;
  Myint         width ;
  Myfloat       coef ;

} ;

} // namespace django

#endif
