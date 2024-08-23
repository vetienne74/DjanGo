//-------------------------------------------------------------------------------------------------------
//                                                                                                       
// DESCRIPTION                                                                                           
//  PARSE COMMAND LINE ARGUMENT                                                      
//                                                                                                       
//-------------------------------------------------------------------------------------------------------

#include "parse_argument.h"

#include <cstdlib>
#include <iostream>
#include <string>

#include "global.h"
#include "output_report.h"
#include "singleton.h"
#include "type_def.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------

Rtn_code parse_argument(int argc, char* argv[])
{
	print_debug(ALL, LIGHT_DEBUG, "IN parse_argument");
	print_info(MASTER, "\n") ;

	// check for -version or -help or -h
	Myint ii = 0 ;
	while (ii < argc)
	{
		if ((string(argv[ii]) == "-version") ||  (string(argv[ii]) == "-v"))
		{
			return(RTN_CODE_KO);
		}
		else if ((argc == 1) || (string(argv[ii]) == "-help") || (string(argv[ii]) == "-h"))
		{
			print_info(MASTER, " list of command line parameters:") ;
			print_info(MASTER, " -cbx [N]       = block size along x-axis (grid points)") ;
			print_info(MASTER, " -cby [N]       = block size along y-axis (grid points)") ;
			print_info(MASTER, " -cbz [N]       = block size along z-axis (grid points)") ;
			print_info(MASTER, " -debug [OPT]   = debug trace OPT=none/light/mid/full") ;
			print_info(MASTER, " -dryrun        = dry run") ;
			print_info(MASTER, " -help or -h    = list of command line parameters") ;
			print_info(MASTER, " -iobufsize [N] = io buffer size (bytes)") ;
			print_info(MASTER, " -nsubx [N]     = number of subdomains along x-axis") ;
			print_info(MASTER, " -nsuby [N]     = number of subdomains along y-axis") ;
			print_info(MASTER, " -nsubz [N]     = number of subdomains along z-axis") ;
			print_info(MASTER, " -version or -v = print version information") ;
			print_info(MASTER, " -xml [NAME]    = NAME=path/name of xml config file") ;
			return(RTN_CODE_KO);
		}
		ii++ ;
	}

	print_line2() ;
	print_info(MASTER, " COMMAND LINE ARGUMENTS\n") ;

	if (argc == 0)
	{
		print_info(MASTER, " No input arguments") ;
	}
	else
	{
		print_info(MASTER, " Executable\t ", argv[0]) ;
		ii = 1 ;

		// loop on all arguments
		while (ii < argc)
		{

			if (string(argv[ii]) == "-cbx")
			{
				ii++ ;
				if (ii >= argc)
				{
					print_error(" parameter is needed after -cbx") ;
					return(RTN_CODE_KO) ;
				}
				Singleton::Instance()->cbx = atoi(argv[ii]);
				print_info(MASTER, " Cache block size X ", Singleton::Instance()->cbx) ;
			}

			else if (string(argv[ii]) == "-cby")
			{
				ii++ ;
				if (ii >= argc)
				{
					print_error(" parameter is needed after -cby") ;
					return(RTN_CODE_KO) ;
				}
				Singleton::Instance()->cby = atoi(argv[ii]);
				print_info(MASTER, " Cache block size Y ", Singleton::Instance()->cby) ;
			}

			else if (string(argv[ii]) == "-cbz")
			{
				ii++ ;
				if (ii >= argc)
				{
					print_error(" parameter is needed after -cbz") ;
					return(RTN_CODE_KO) ;
				}
				Singleton::Instance()->cbz = atoi(argv[ii]);
				print_info(MASTER, " Cache block size Z ", Singleton::Instance()->cbz) ;
			}

			else if (string(argv[ii]) == "-nsubx")
			{
				ii++ ;
				if (ii >= argc)
				{
					print_error(" parameter is needed after -nsubx") ;
					return(RTN_CODE_KO) ;
				}
				Singleton::Instance()->nsubx = atoi(argv[ii]);
				print_info(MASTER, " No. subdomains in X ", Singleton::Instance()->nsubx) ;
			}

			else if (string(argv[ii]) == "-nsuby")
			{
				ii++ ;
				if (ii >= argc)
				{
					print_error(" parameter is needed after -nsuby") ;
					return(RTN_CODE_KO) ;
				}
				Singleton::Instance()->nsuby = atoi(argv[ii]);
				print_info(MASTER, " No. subdomains in Y ", Singleton::Instance()->nsuby) ;
			}

			else if (string(argv[ii]) == "-nsubz")
			{
				ii++ ;
				if (ii >= argc)
				{
					print_error(" parameter is needed after -nsubz") ;
					return(RTN_CODE_KO) ;
				}
				Singleton::Instance()->nsubz = atoi(argv[ii]);
				print_info(MASTER, " No. subdomains in Z ", Singleton::Instance()->nsubz) ;
			}

			else if (string(argv[ii]) == "-iobufsize")
			{
				ii++ ;
				if (ii >= argc)
				{
					print_error(" parameter is needed after -iobufsize") ;
					return(RTN_CODE_KO) ;
				}
				Singleton::Instance()->iobufsize = atoi(argv[ii]);
				print_info(MASTER, " IO buffer (bytes) ", Singleton::Instance()->iobufsize) ;
			}

			else if (string(argv[ii]) == "-dryrun")
			{
				print_info(MASTER, " *** DRY RUN *** ") ;
				Singleton::Instance()->dryrun = true ;
			}

			else if (string(argv[ii]) == "-xml")
			{
				ii++ ;
				if (ii >= argc)
				{
					print_error(" parameter is needed after -xml") ;
					return(RTN_CODE_KO) ;
				}
				Singleton::Instance()->xml_config_file = argv[ii] ;
				print_info(MASTER, " XML config file", Singleton::Instance()->xml_config_file.c_str()) ;

			}

			else if (string(argv[ii]) == "-debug")
			{
				ii++ ;
				if (ii >= argc)
				{
					print_error(" parameter is needed after -debug") ;
					return(RTN_CODE_KO) ;
				}
				if (string(argv[ii]) == "none")
				{
					debug = NO_DEBUG ;
					print_info(MASTER, " Debug info level \tNONE") ;
				}
				else if (string(argv[ii]) == "light")
				{
					debug = LIGHT_DEBUG ;
					print_info(MASTER, " Debug info level \tLIGHT") ;
				}
				else if (string(argv[ii]) == "mid")
				{
					debug = MID_DEBUG ;
					print_info(MASTER, " Debug info level \tMEDIUM") ;
				}
				else if (string(argv[ii]) == "full")
				{
					debug = FULL_DEBUG ;
					print_info(MASTER, " Debug info level \tFULL") ;
				}
				else
				{
					print_error(" Invalid debug level", argv[ii]) ;
					return(RTN_CODE_KO) ;
				}
			}

			else
			{
				print_error(" Unknown argument ", argv[ii]) ;
				return(RTN_CODE_KO);
			}
			ii++ ;
		}
	}

	print_debug(ALL, LIGHT_DEBUG, "OUT parse_argument");
	return(RTN_CODE_OK);
}  

} // namespace django
