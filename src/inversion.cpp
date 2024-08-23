//-------------------------------------------------------------------------------------------------------
//
// INVERSION
//
//-------------------------------------------------------------------------------------------------------

#include "inversion.h"

#include <cassert>
#include <fstream>
#include <iostream>

#include "mpi.h"

#include "acquisition.h"
#include "data.h"
#include "freq_group.h"
#include "gradient.h"
#include "grid.h"
#include "model.h"
#include "modelling.h"
#include "output_report.h"
#include "parse_xml.h"
#include "singleton.h"
#include "snapshot.h"

using namespace std;

namespace django {

// constructor

Inversion::Inversion(Modelling* pModelling_in)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Inversion::Inversion");
	pModelling = pModelling_in ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Inversion::Inversion");
}

// setters
void Inversion::set_niter(Myint niter_in)
{
	niter = niter_in ;
}
void Inversion::set_ntry(Myint ntry_in)
{
	ntry = ntry_in ;
}
void Inversion::set_init_try(Myfloat init_try_in)
{
	init_try = init_try_in ;
}
void Inversion::set_freq_inv(Myfloat min, Myfloat max, Myfloat delta)
{
	pFreq_group = new Freq_group(min, max, delta) ;
}

//-------------------------------------------------------------------------------------------------------
void Inversion::open_config_file(const char* file_name) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Inversion::open_config_file");
	config_file.open(file_name) ;
	assert(config_file.is_open());
	print_debug(ALL, LIGHT_DEBUG, "OUT Inversion::open_config_file");

}

//-------------------------------------------------------------------------------------------------------
void Inversion::close_config_file(void) 
{
	print_debug(ALL, LIGHT_DEBUG, "IN Inversion::close_config_file");
	config_file.close() ;
	print_debug(ALL, LIGHT_DEBUG, "OUT Inversion::close_config_file");

}

//-------------------------------------------------------------------------------------------------------

Rtn_code Inversion::read_config()
{
	print_debug(ALL, LIGHT_DEBUG, "IN Inversion::read_config");

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
			printf( "%s at line %" XML_FMT_INT_MOD "u\n",
					XML_ErrorString(XML_GetErrorCode(parser)),
					XML_GetCurrentLineNumber(parser) );
			print_error(" Error while parsing XML configuration file") ;
			return(RTN_CODE_KO) ;
		}
	} while ( !done ) ;

	// free parser
	XML_ParserFree(parser);
	config_file.close() ;

	// check parameters

	print_info(MASTER, "") ;
	print_line2() ;
	print_info(MASTER, " INVERSION PARAMETERS") ;
	print_info(MASTER, "");

	// no. iterations
	if (niter > 0)
	{
		print_info(MASTER, " No. iteration:\t", niter) ;
	}
	else
	{
		print_error(" Invalid no. iterations") ;
		return(RTN_CODE_KO) ;
	}

	// no. try for line search
	if (ntry >= 0)
	{
		print_info(MASTER, " No. try line search:", ntry) ;
	}
	else
	{
		print_error(" Invalid no. try for line search") ;
		return(RTN_CODE_KO) ;
	}

	// no. try for line search
	if (init_try > 0.)
	{
		print_info(MASTER, " Initial step:\t", init_try) ;
	}
	else
	{
		print_error(" Invalid initial step for line search") ;
		return(RTN_CODE_KO) ;
	}

	print_line2() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Inversion::read_config");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Inversion::initialize()
{
	print_debug(ALL, LIGHT_DEBUG, "IN Inversion::initialize");

	//-------------------------------------------------------------------------------------------------------
	// read inversion configuration file
	//-------------------------------------------------------------------------------------------------------
	Rtn_code rtn_code = (*this).read_config() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_debug(ALL, LIGHT_DEBUG, "Out Inversion::initialize");
	return(RTN_CODE_OK) ;
}


//=======================================================================================================
//
// FULL WAVEFORM INVERSION
//
//=======================================================================================================

Rtn_code Inversion::perform_fwi(Acquisition* pAcquisition, Data* pData, Model* pModel)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Inversion::perform_fwi");

	//-------------------------------------------------------------------------------------------------------
	// initialize inversion
	//-------------------------------------------------------------------------------------------------------
	Rtn_code rtn_code = this->initialize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	//-------------------------------------------------------------------------------------------------------
	// initialize acquisition
	//-------------------------------------------------------------------------------------------------------
	rtn_code = pAcquisition->initialize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	//-------------------------------------------------------------------------------------------------------
	// initialize model
	//-------------------------------------------------------------------------------------------------------
	rtn_code = pModel->initialize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	rtn_code = pModel->read() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	//-------------------------------------------------------------------------------------------------------
	// instanciate gradient
	//-------------------------------------------------------------------------------------------------------
	pGradient = new Gradient(pAcquisition, pModel, pFreq_group) ;

	//-------------------------------------------------------------------------------------------------------
	// initialize modelling
	//-------------------------------------------------------------------------------------------------------
	rtn_code = pModelling->initialize(pModel) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
	pModelling->set_freq_list_to_extract(pFreq_group) ;

	//-------------------------------------------------------------------------------------------------------
	// loop over iterations
	//-------------------------------------------------------------------------------------------------------

	// open fcost
	ofstream fcost_file ;
	fcost_file.open(FCOST_OUT_FILE) ;

	double t0 = MPI_Wtime() ;
	Myfloat alpha_try, alpha_optim, alpha1, alpha2 ;
	Myfloat fcost0, fcost_new, fcost_alpha1, fcost_alpha2 ;

	for (Myint iter = 1; iter <= niter; iter++)
	{
		print_info(MASTER, "") ;
		print_line5() ;
		print_info(MASTER, " ITERATION", iter) ;
		print_info(MASTER, "");

		//-------------------------------------------------------------------------------------------------------
		// compute gradient
		//-------------------------------------------------------------------------------------------------------

		rtn_code = pGradient->compute(pAcquisition, pData, pModelling, pModel, &fcost0) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		print_info(MASTER, " FCOST CURRENT MODEL", fcost0) ;
		fcost_file << iter << "\t" << 0 << "\t" << 0. << "\t" << fcost0 << "\n" ;

		//-------------------------------------------------------------------------------------------------------
		// line search
		//-------------------------------------------------------------------------------------------------------

		alpha1       = 0. ;
		alpha2       = 0. ;
		fcost_alpha1 = 0. ;
		fcost_alpha2 = 0. ;

		if (iter == 1)
		{
			alpha_try = -init_try ;
		}
		else
		{
			alpha_try = alpha_optim ;
		}

		//Model_2D_ac_iso *pModel1 = dynamic_cast<Model_2D_ac_iso*>(pModel) ;
		//Model_2D_ac_iso *pModel2 = dynamic_cast<Model_2D_ac_iso*>(pGradient->pGrad_precond) ;
		Model* pModel1 = pModel ;
		Model* pModel2 = pGradient->pGrad_precond ;

		for (Myint itry = 1; itry <= ntry; itry++)
		{

			print_info(MASTER, "") ;
			print_line2() ;
			print_info(MASTER, " TRY", itry) ;
			print_info(MASTER, " ALPHA ", alpha_try) ;
			print_info(MASTER, "");

			// update model
			//Model_2D_ac_iso *pNew_model = new Model_2D_ac_iso(*pModel1, *pModel2, alpha_try) ;
			Model* pNew_model = new Model() ;
			pNew_model->info() ;

			// compute new cost function
			rtn_code = compute_fcost_freq_domain(pAcquisition, pData, pNew_model, &fcost_new) ;
			if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

			print_info(MASTER, " FCOST NEW MODEL", fcost_new) ;
			fcost_file << iter << "\t" << itry << "\t" << alpha_try << "\t" << fcost_new << "\n" ;

			// initialize alpha1 and alpha2
			if (itry == 1)
			{
				alpha1       = alpha_try ;
				fcost_alpha1 = fcost_new ;

				alpha2       = alpha_try ;
				fcost_alpha2 = fcost_new ;
			}

			//------------------------
			// case 1
			// fcost_alpha1 < fcost0
			//------------------------

			if (fcost_alpha1 < fcost0)
			{

				if (fcost_new > fcost_alpha1)
				{
					alpha2       = alpha_try ;
					fcost_alpha2 = fcost_new ;
					//delete((Model_2D_ac_iso*) pNew_model) ; // delete new model
					delete(pNew_model) ;
					break ;
				}
				else if (fcost_new <= fcost_alpha1)
				{
					alpha1       = alpha_try ;
					fcost_alpha1 = fcost_new ;
					alpha_try    = 2. * alpha1 ;
					//delete((Model_2D_ac_iso*) pNew_model) ; // delete new model
					delete(pNew_model) ;
				}
			}

			//------------------------
			// case 2
			// fcost_alpha1 > fcost0
			//------------------------

			if (fcost_alpha1 > fcost0)
			{

				if (fcost_new >= fcost_alpha1)
				{
					alpha1       = alpha_try ;
					fcost_alpha1 = fcost_new ;
					alpha_try    = 0.5 * alpha1 ;
					//delete((Model_2D_ac_iso*) pNew_model) ; // delete new model
					delete(pNew_model) ;
				}

				else if (fcost_new < fcost_alpha1)
				{
					alpha1       = alpha_try ;
					fcost_alpha1 = fcost_new ;
					//delete((Model_2D_ac_iso*) pNew_model) ; // delete new model
					delete(pNew_model) ;

					if (fcost_alpha1 > fcost0)
					{
						alpha_try    = 0.5 * alpha1 ;
					}
					else
					{
						break ;
					}
				}
			}

		} // for (Myint itry = 1; itry <= ntry; itry++)

		// determine optimal step length and new model
		// parabola fitting to be done

		cout << "alpha1       " << alpha1 << "\n" ;
		cout << "fcost_alpha1 " << fcost_alpha1 << "\n" ;
		cout << "alpha2       " << alpha2 << "\n" ;
		cout << "fcost_alpha2 " << fcost_alpha2 << "\n" ;

		alpha_optim = alpha1 ;
		//Model_2D_ac_iso *pNew_model = new Model_2D_ac_iso(*pModel1, *pModel2, alpha_optim) ;
		//Model* pNew_model = new Model(pModel1, pModel2, alpha_optim) ;
		Model* pNew_model = new Model() ;
		pNew_model->write(iter) ;
		pNew_model->info() ;
		delete(pModel1) ;     // delete old model
		pModel = pNew_model ;

		if (fcost_alpha1 > fcost0)
		{
			cout << " *** STOP ITERATION *** \n" ;
			break ;
		}

	} // for (Myint iter = 1; iter <= niter; iter++)

	rtn_code = pData->finalize(pAcquisition, FWI_PROG) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	//-------------------------------------------------------------------------------------------------------
	// delete objects
	//-------------------------------------------------------------------------------------------------------
	delete(pGradient) ;
	//delete((Model_2D_ac_iso*) pModel) ;
	delete(pModel) ;
	delete(pFreq_group) ;

	double t3 = MPI_Wtime() ;
	double total_time = t3 - t0 ;
	if (total_time < 60.0)
	{
		print_info(MASTER, " TOTAL TIME (sec):", (float) (total_time)) ;
	}
	else if (total_time < 3600.0)
	{
		print_info(MASTER, " TOTAL TIME (min):", (float) (total_time/60.)) ;
	}
	else
	{
		print_info(MASTER, " TOTAL TIME (hrs):", (float) (total_time/3600.)) ;
	}
	print_line2() ;

	// close fcost
	fcost_file.close() ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Inversion::perform_fwi");
	return(RTN_CODE_OK) ;
}

//=======================================================================================================
//
// COMPUTE FCOST IN FREQUENCY DOMAIN
//
//=======================================================================================================


Rtn_code Inversion::compute_fcost_freq_domain(Acquisition* pAcquisition, 
		Data* pData,
		Model* pModel,
		Myfloat* fcost)
{
	print_debug(ALL, LIGHT_DEBUG, "IN Inversion::compute_fcost_freq_domain");

	// return code
	Rtn_code rtn_code ;

	print_info(MASTER, "") ;
	print_line2() ;
	print_info(MASTER, " START MODELLING") ;
	print_info(MASTER, "") ;

	//-------------------------------------------------------------------------------------------------------
	// initialize modelling
	//-------------------------------------------------------------------------------------------------------
	rtn_code = pModelling->initialize(pModel) ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	double t0 = MPI_Wtime() ;

	// instantiate snapshot
	Snapshot* pSnapshot ;
	pSnapshot = NULL ;

	//-------------------------------------------------------------------------------------------------------
	// start loop on shots until all have been computed
	//-------------------------------------------------------------------------------------------------------

	pAcquisition->move_to_first_shot() ;
	*fcost = 0 ;

	while ( pAcquisition->remain_shot_to_compute() )

	{

		//-------------------------------------------------------------------------------------------------------
		// retrieve src and rec positions the current shot
		//-------------------------------------------------------------------------------------------------------

		rtn_code = pAcquisition->get_position_for_current_shot() ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		//-------------------------------------------------------------------------------------------------------
		// reset data
		//-------------------------------------------------------------------------------------------------------

		rtn_code = pData->reset(pAcquisition, FWI_PROG) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		//-------------------------------------------------------------------------------------------------------
		// compute INCIDENT wavefield
		//-------------------------------------------------------------------------------------------------------

		print_info(ALL, " COMPUTE INCIDENT WAVEFIELD FOR SOURCE:", pAcquisition->current_src) ;
		double t1 = MPI_Wtime() ;

		rtn_code = pModelling->pScheme->solve_current_shot(pAcquisition, pData, INCIDENT, pSnapshot) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

		//-------------------------------------------------------------------------------------------------------
		// evaluate obj function
		//-------------------------------------------------------------------------------------------------------

		Myfloat fcost_tmp = 0 ;
		rtn_code = pGradient->evaluate_fcost(pAcquisition, pData, pModelling, &fcost_tmp) ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;
		*fcost += fcost_tmp ;

		// increment shot counter

		rtn_code = pAcquisition->move_to_next_shot() ;
		if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	}

	//-------------------------------------------------------------------------------------------------------
	// end loop on shots until all have been computed
	//-------------------------------------------------------------------------------------------------------

	//rtn_code = this->deallocate_grid_for_freq_wavefield(pIncident) ;
	//if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	//-------------------------------------------------------------------------------------------------------
	// finalize modelling
	//-------------------------------------------------------------------------------------------------------

	rtn_code = pModelling->finalize() ;
	if (rtn_code != RTN_CODE_OK) return(rtn_code) ;

	print_info(MASTER, "") ;
	print_info(MASTER, " END MODELLING") ;
	print_info(MASTER, "") ;

	print_debug(ALL, LIGHT_DEBUG, "OUT Inversion::compute_fcost_freq_domain");
	return(RTN_CODE_OK) ;
}

} // namespace django



