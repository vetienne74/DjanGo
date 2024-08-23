//-------------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//  WRITE INFORMATION IN PROGRAM OUTPUT REPORT
//
//-------------------------------------------------------------------------------------------------------

#include "output_report.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <string>

#include "global.h"
#include "singleton.h"
#include "type_def.h"
#include "version_django.h"

using namespace std;

namespace django {

static const char LINE_REPORT_T1[] = "================================================================================" ;
static const char LINE_REPORT_T2[] = "--------------------------------------------------------------------------------" ;
static const char LINE_REPORT_T3[] = "................................................................................" ;
static const char LINE_REPORT_T4[] = "XxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXx" ;
static const char LINE_REPORT_T5[] = "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" ;

static const char LOGO_DJANGO1[] = " o-/\\-o" ;
static const char LOGO_DJANGO2[] = " o-\\_\\-o" ;
static const char LOGO_DJANGO3[] = "     \\\\                     ___" ;
static const char LOGO_DJANGO4[] = "    /__\\    O  __          /___\\" ;
static const char LOGO_DJANGO5[] = "   //  \\\\   \\ / _\\  /\\  / //  __   /^\\" ;
static const char LOGO_DJANGO6[] = "  //    ||  /  /  \\ \\ \\ \\ ||  \\_\\ /   \\" ;
static const char LOGO_DJANGO7[] = " / \\____//  \\  \\__/ /  \\/ \\\\___// \\___/" ;
static const char LOGO_DJANGO8[] = " \\______/ \\_/              \\___/" ;

//-------------------------------------------------------------------------------------------------------

Rtn_code print_header_of_output_report(void)
{
	if (myid_world == 0)
	{
		cout << "\n" << LINE_REPORT_T1 << "\n" ;
		cout << "\n\t\t" << LOGO_DJANGO1 << "\n" ;
		cout << "\t\t" << LOGO_DJANGO2 << "\n" ;
		cout << "\t\t" << LOGO_DJANGO3 << "\n" ;
		cout << "\t\t" << LOGO_DJANGO4 << "\n" ;
		cout << "\t\t" << LOGO_DJANGO5 << "\n" ;
		cout << "\t\t" << LOGO_DJANGO6 << "\n" ;
		cout << "\t\t" << LOGO_DJANGO7 << "\n" ;
		cout << "\t\t" << LOGO_DJANGO8 << "\n" ;
		cout << "\n\t\t\t D j a n G o - ver " << CURRENT_VERSION << " (2024)\n\n" ;

		cout << " Git version " << DJANGO_GIT_COMMIT << "\n" ;
		cout << " " << DJANGO_GIT_AUTHOR << "\n" ;
		cout << " " << DJANGO_GIT_DATE << "\n" ;

		// single or double precision
#ifdef _DOUBLE_PRECISION_
		cout << " Computations in DOUBLE PRECISION\n" ;
#else
		cout << " Computations in SINGLE PRECISION\n" ;
#endif

		// check binary validity
		time_t current_time = time(nullptr) ;
		time_t compile_time = DJANGO_COMPILE_DATE ;
		time_t validity_time =  BINARY_VALIDITY_TIME ;
		time_t time_since_compile = current_time - compile_time ;

		//if (time_since_compile > validity_time)
		//{
		//	cout << "\n *** *** DjanGo binary validity expired! *** ***\n" ;
		//	cout << " current use (days) " << time_since_compile / 60 / 60 / 24 << " / max validity (days) "
		//			<< validity_time / 60 / 60 / 24 << "\n" ;
		//	cout << " contact vetienne@rocketmail.com to get a valid binary \n" ;
		//	cout << LINE_REPORT_T1 << "\n" ;
		//	return(RTN_CODE_KO) ;
		//}

		cout << LINE_REPORT_T1 ;
	}

	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code print_end_of_output_report(void)
{
	if (myid_world == 0)
	{
		cout << "\n" ;
		if (current_mem > 0)
		{
			print_info(MASTER, " BYTES STILL IN MEMORY:", current_mem) ;

			Singleton::Instance()->print_all_variable() ;
		}
		else if (current_mem < 0)
		{
			print_info(MASTER, " SOMETHING IS WRONG, BYTES IN MEMORY:", current_mem) ;

			Singleton::Instance()->print_all_variable() ;
		}
		else
		{
			print_info(MASTER, " MEMORY CLEAN OK") ;
		}

		if (max_mem < 1e+3) {
			print_info(MASTER, " MAX. MEMORY USED (B):", Myint(max_mem)) ;
		} else if (max_mem < 1e+6) {
			print_info(MASTER, " MAX. MEMORY USED (KB):", Myfloat(max_mem/1e+3)) ;
		} else if (max_mem < 1e+9) {
			print_info(MASTER, " MAX. MEMORY USED (MB):", Myfloat(max_mem/1e+6)) ;
		}	else {
			print_info(MASTER, " MAX. MEMORY USED (GB):", Myfloat(max_mem/1e+9)) ;
		}

		cout << "\n" ;
		cout << LINE_REPORT_T1 ;
		cout << "\n" ;
		cout << "\t\t\t DJANGO TERMINATED SUCCESSFULLY\n" ;
		cout << LINE_REPORT_T1 ;
		cout << "\n" ;
	}

	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code print_modelling_parameter(void)
{
	cout << LINE_REPORT_T2 ;
	cout << "\n\t\t\t\t MODELLING PARAMETERS\n" ;
	cout << LINE_REPORT_T3 ;

	cout << LINE_REPORT_T2 ;

	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

void print_debug(Display_type display_t, Debug_level debug_l, char* text)
{
	if (debug_l <= debug)
	{
		if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
		{
			cout << "* " << text << "\n" << flush ;
		}
	}
}

//-------------------------------------------------------------------------------------------------------

void print_debug(Display_type display_t, Debug_level debug_l, char* text, string text2)
{
	if (debug_l <= debug)
	{
		if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
		{
			cout << "* " << text << " " << text2 << "\n" << flush ;
		}
	}
}

//-------------------------------------------------------------------------------------------------------

void print_debug(Display_type display_t, Debug_level debug_l, char* text, Myint64 val)
{
	if (debug_l <= debug)
	{
		if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
		{
			cout << "* " << text << " " << val << "\n" << flush ;
		}
	}
}

//-------------------------------------------------------------------------------------------------------

void print_debug(Display_type display_t, Debug_level debug_l, char* text, Myint nb)
{
	if (debug_l <= debug)
	{
		if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
		{
			cout << "* " << text << " " << nb << "\n" << flush ;
		}
	}
}

//-------------------------------------------------------------------------------------------------------

void print_debug(Display_type display_t, Debug_level debug_l, char* text, Myfloat nb)
{
	if (debug_l <= debug)
	{
		if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
		{
			cout << "* " << text << " " << nb << "\n" << flush ;
		}
	}
}

//-------------------------------------------------------------------------------------------------------

void print_debug(Display_type display_t, Debug_level debug_l, const char* text, const char* text2)
{
	if (debug_l <= debug)
	{
		if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
		{
			cout << text << "\t" << text2 << "\n" << flush ;
		}
	}
}

//-------------------------------------------------------------------------------------------------------

void print_line1(void)
{
	if (myid_world == 0) cout << LINE_REPORT_T1 << "\n" << flush ;
}
void print_line2(void)
{
	if (myid_world == 0) cout << LINE_REPORT_T2 << "\n" << flush ;
}
void print_line3(void)
{
	if (myid_world == 0) cout << LINE_REPORT_T3 << "\n" << flush ;
}
void print_line4(void)
{
	if (myid_world == 0) cout << LINE_REPORT_T4 << "\n" << flush ;
}
void print_line5(void)
{
	if (myid_world == 0) cout << LINE_REPORT_T5 << "\n" << flush ;
}

//-------------------------------------------------------------------------------------------------------

void print_info(Display_type display_t, char* text)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void print_info(Display_type display_t, const char* text)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void print_info(Display_type display_t, char* text, Myint nb)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\t" << nb << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void print_info(Display_type display_t, const char* text, Myint nb)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\t" << nb << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void print_info(Display_type display_t, const char* text, Myint64 nb)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\t" << nb << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void print_info(Display_type display_t, char* text, Myfloat nb)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\t" << nb << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void print_info(Display_type display_t, char* text, Myfloat nb, Myfloat nb2)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\t" << nb << "\t" << nb2 << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void print_info(Display_type display_t, const char* text, Myfloat nb)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\t" << nb << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void print_info(Display_type display_t, char* text, char* text2)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\t" << text2 << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void print_info(Display_type display_t, const char* text, char* text2)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\t" << text2 << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void print_info(Display_type display_t, const char* text, const char* text2)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\t" << text2 << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void print_info(Display_type display_t, const char* text, string text2)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\t" << text2 << "\n" << flush ;
	}
}


//-------------------------------------------------------------------------------------------------------

void print_warning(char* text, Myint nb)
{
	cout << "\n" ;
	cout << LINE_REPORT_T3 ;
	cout << "\n\t\t\t W A R N I N G \n" ;
	cout << text << "\t" << nb << "\n\n" ;
	cout << LINE_REPORT_T3 ;
	cout << "\n" ;
}

//-------------------------------------------------------------------------------------------------------

void print_warning(char* text)
{
	cout << "\n" ;
	cout << LINE_REPORT_T3 ;
	cout << "\n\t\t\t W A R N I N G \n" ;
	cout << text << "\n\n" ;
	cout << LINE_REPORT_T3 ;
	cout << "\n" ;
}

//-------------------------------------------------------------------------------------------------------

void print_warning(const char* text)
{
	cout << "\n" ;
	cout << LINE_REPORT_T3 ;
	cout << "\n\t\t\t W A R N I N G \n" ;
	cout << text << "\n\n" ;
	cout << LINE_REPORT_T3 ;
	cout << "\n" ;
}

//-------------------------------------------------------------------------------------------------------

void print_error(string* text)
{
	cout << "\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n\t\t\t F A T A L    E R R O R\n\n" ;
	cout << *text << "\n\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n" ;
}

//-------------------------------------------------------------------------------------------------------

void print_error(char* text, Myint nb)
{
	cout << "\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n\t\t\t F A T A L    E R R O R\n\n" ;
	cout << text << "\t" << nb << "\n\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n" ;
}

//-------------------------------------------------------------------------------------------------------

void print_error(char* text, string text2)
{
	cout << "\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n\t\t\t F A T A L    E R R O R\n\n" ;
	cout << text << "\t" << text2 << "\n\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n" ;
}

//-------------------------------------------------------------------------------------------------------

void print_error(char* text, Myint nb1, Myint nb2)
{
	cout << "\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n\t\t\t F A T A L    E R R O R\n\n" ;
	cout << text << "\t" << nb1 << "\t" << nb2 << "\n\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n" ;
}

//-------------------------------------------------------------------------------------------------------

void print_error(char* text)
{
	cout << "\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n\t\t\t F A T A L    E R R O R\n\n" ;
	cout << text << "\n\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n" ;
}

//-------------------------------------------------------------------------------------------------------

void print_error(char* text1, char* text2)
{
	cout << "\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n\t\t\t F A T A L    E R R O R\n\n" ;
	cout << text1 << "\t" << text2 << "\n\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n" ;
}

//-------------------------------------------------------------------------------------------------------

void print_error(const char* text)
{
	cout << "\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n\t\t\t F A T A L    E R R O R\n\n" ;
	cout << text << "\n\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n" ;
}

} // namespace django
