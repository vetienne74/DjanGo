#ifndef DJANGO_OUTPUT_REPORT_H_
#define DJANGO_OUTPUT_REPORT_H_

#include <string>

#include "type_def.h"

namespace django {

Rtn_code print_header_of_output_report(void) ;
Rtn_code print_end_of_output_report(void) ;
Rtn_code print_modelling_parameter(void) ;

void print_line1(void) ;
void print_line2(void) ;
void print_line3(void) ;
void print_line4(void) ;
void print_line5(void) ;

void print_info(Display_type, char*) ;
void print_info(Display_type, const char*) ;
void print_info(Display_type, char*, Myint) ;
void print_info(Display_type, const char*, Myint) ;
void print_info(Display_type, const char*, Myint64) ;
void print_info(Display_type, char*, Myfloat) ;
void print_info(Display_type, char*, Myfloat, Myfloat) ;
void print_info(Display_type, const char*, Myfloat) ;
void print_info(Display_type, char*, char*) ;
void print_info(Display_type, const char*, char*) ;
void print_info(Display_type, const char*, string) ;
void print_info(Display_type, const char*, const char*) ;

void print_debug(Display_type, Debug_level, char*) ;
void print_debug(Display_type, Debug_level, char*, string) ;
void print_debug(Display_type, Debug_level, char*, Myint) ;
void print_debug(Display_type, Debug_level, char*, Myint64) ;
void print_debug(Display_type, Debug_level, char*, Myfloat) ;
void print_debug(Display_type, Debug_level, const char*, const char*) ;

void print_warning(string* text) ;
void print_warning(char* text) ;
void print_warning(char* text, Myint) ;
void print_warning(const char* text) ;

void print_error(string* text) ;
void print_error(char* text) ;
void print_error(char* text, Myint) ;
void print_error(char* text, string) ;
void print_error(char* text, Myint, Myint) ;
void print_error(char* text1, char* text2) ;
void print_error(const char* text) ;

} // namespace django

#endif
