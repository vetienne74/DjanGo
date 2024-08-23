//-------------------------------------------------------------------------------------------------------
//
// CLASS DOMAIN
//
//-------------------------------------------------------------------------------------------------------

#include "domain.h"

#include <fstream>
#include <cassert>

#include "output_report.h"
#include "parse_xml.h"
#include "singleton.h"
#include "type_def.h"

namespace django {

//-------------------------------------------------------------------------------------------------------

Domain::Domain(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Domain::Domain");
	print_debug(ALL, LIGHT_DEBUG, "OUT Domain::Domain");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Domain::initialize(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Domain::initialize");

	//-----------------------------------------
	// parse XML configuration file with EXPAT
	//-----------------------------------------
	char buf[BUFSIZ];
	XML_Parser parser = XML_ParserCreate(NULL);
	int done;
	int depth = 0;

	// init parser
	XML_SetUserData(parser, &depth);
	XML_SetElementHandler(parser, django_startElement, django_endElement);

	// open xml config file
	ifstream config_file(Singleton::Instance()->xml_config_file.c_str()) ;
	assert(config_file.is_open());

	// parse file byte after byte up to the end
	do {
		int len = 1 ;
		config_file.read (buf,len);
		done = config_file.eof() ;
		if (XML_Parse(parser, buf, len, done) == XML_STATUS_ERROR) {
			fprintf(stderr,
					"%s at line %" XML_FMT_INT_MOD "u\n",
					XML_ErrorString(XML_GetErrorCode(parser)),
					XML_GetCurrentLineNumber(parser));
			print_error(" Error while parsing XML configuration file") ;
			return(RTN_CODE_KO) ;
		}
	} while ( !done ) ;

	// free parser
	XML_ParserFree(parser);
	config_file.close() ;

	Rtn_code rtn = this->info() ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Domain::initialize");
	return(rtn) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Domain::info(void)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Domain::info");

	print_info(MASTER, "") ;
	print_line2() ;
	print_info(MASTER, " DOMAIN PARAMETERS") ;

	// model
	if (pModel == NULL)
	{
		print_error(" Model non initialized") ;
		return(RTN_CODE_KO) ;
	}
	Rtn_code rtn = pModel->info() ;
	if (rtn != RTN_CODE_OK) return(rtn) ;

	// scheme
	if (pScheme == NULL)
	{
		print_error(" Scheme non initialized") ;
		pScheme_Factory->info() ;
		return(RTN_CODE_KO) ;
	}
	rtn = pScheme->info() ;
	if (rtn != RTN_CODE_OK) return(rtn) ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Domain::info");
	return(RTN_CODE_OK) ;
}

} // namespace django
