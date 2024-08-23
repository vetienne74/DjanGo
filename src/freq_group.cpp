//-------------------------------------------------------------------------------------------------------
//
// FREQUENCY GROUPS
//
//-------------------------------------------------------------------------------------------------------

#include "freq_group.h"

#include <cassert>
#include <fstream>

#include "allocate_array.h"
#include "output_report.h"
#include "type_def.h"

using namespace std;

namespace django {

//-------------------------------------------------------------------------------------------------------

Freq_group::~Freq_group(void)
{
	if (nb_freq > 0) deallocate_array<Myfloat>(pFreq_list, nb_freq) ;
}

//-------------------------------------------------------------------------------------------------------

Freq_group::Freq_group(Myfloat min, Myfloat max, Myfloat delta)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Freq_group::Freq_group");

	print_info(MASTER, "") ;
	print_line2() ;
	print_info(MASTER, " FREQUENCY GROUPS DEFINITION") ;
	print_info(MASTER, "") ;

	// nb freq
	if (delta == 0.0)
	{
		nb_freq = 1 ;
	}
	else
	{
		nb_freq = floor((max - min) / delta) + 1 ;
	}

	if (nb_freq >= 0)
	{
		print_info(MASTER, " Nb frequency:\t", nb_freq) ;
	}
	else
	{
		print_error(" Nb frequeny < 0 in config file") ;
	}

	// list of freq
	if (nb_freq >= 1)
	{
		pFreq_list = allocate_array<Myfloat>(nb_freq) ;
		for (int ifreq = 0; ifreq < nb_freq; ifreq++)
		{
			pFreq_list[ifreq] = min + ifreq * delta ;
			print_info(MASTER, " Freq (Hz):\t", pFreq_list[ifreq]) ;
		}
	}

	print_line2() ;

	print_debug(ALL, LIGHT_DEBUG, "IN Freq_group::Freq_group");
}

} // namespace django
