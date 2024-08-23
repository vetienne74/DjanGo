#ifndef DJANGO_PROGRAM_GUITAR_H_
#define DJANGO_PROGRAM_GUITAR_H_

#include <string>

#include "domain.h"
#include "program.h"
#include "scheme.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class Program_Guitar : public Program

{
public:

	Program_Guitar(Prog_type) ;
	Rtn_code initialize(void) ;
	Rtn_code run(void) ;
	Rtn_code finalize(void) ;
	Rtn_code info(void) ;

	Rtn_code set_left_hand(string name, Myfloat dt) ;

protected:
	Model* pModel ;
	Data*  pData ;
	string left_hand_file ;
	Myfloat left_hand_file_dt ;
	Myfloat* left_hand_time_function ;
	Rtn_code get_left_hand_time_function(void) ;
	Rtn_code free_left_hand_time_function(void) ;

} ;

} // namespace django

#endif
