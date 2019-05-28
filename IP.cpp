/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//																										  //	
//								NGUYEN DUY LINH -1511755												  //	
//								LVTN:IP for SMTI														  //	
//																										  //
///////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "IP.hpp"
#include <sstream>
#include "scip_exception.hpp"
#include "scip/pub_message.h"
#include <math.h> 
#include <cmath>  
using namespace std;
using namespace scipexamples;


/////////////////////////////////////////////////

////////////////1_NOBIN_1STA_NOMERGED/////////////////////////
/////////////////////////////////////////////////

/* constructor */
//scipexamples::IP::IP(Allocation& allo)
//	: _scip(0), n1(allo.nbChildren), n2(allo.nbFamilies), _cons()
//{
//
//	// initialize scip
//	SCIP_CALL_EXC(SCIPcreate(&_scip));
//
//	// load default plugins linke separators, heuristics, etc.
//	SCIP_CALL_EXC(SCIPincludeDefaultPlugins(_scip));
//
//	// disable scip output to stdout
//	SCIPmessagehdlrSetQuiet(SCIPgetMessagehdlr(_scip), TRUE);
//
//	// create an empty problem
//	SCIP_CALL_EXC(SCIPcreateProb(_scip, "SMTI", NULL, NULL, NULL, NULL, NULL, NULL, NULL));
//
//	// set the objective sense to maximize, default is minimize
//	SCIP_CALL_EXC(SCIPsetObjsense(_scip, SCIP_OBJSENSE_MAXIMIZE));
//
//	
//
//	 //create a binary variable xij
//	ostringstream namebuf;
//	_vars.resize(allo.nbChildren);
//
//	for (size_t i = 0; i < allo.nbChildren; ++i)
//	{
//		_vars[i].resize(allo.children[i].nbPref);
//
//		for (size_t j = 0; j < allo.children[i].nbPref; ++j)
//		{
//			_vars[i][j].resize(allo.children[i].preferences[j].size());
//			for (size_t k = 0; k < allo.children[i].preferences[j].size(); ++k)
//			{
//				SCIP_VAR* var;
//				namebuf.str("");
//				namebuf << "x#" << i << "#" << j << "#" << k;
//
//				// create the SCIP_VAR object
//				SCIP_CALL_EXC(SCIPcreateVar(_scip, &var, namebuf.str().c_str(), 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));
//
//				// add the SCIP_VAR object to the scip problem
//				SCIP_CALL_EXC(SCIPaddVar(_scip, var));
//
//				// storing the SCIP_VAR pointer for later access
//				_vars[i][j][k] = var;
//			}
//		}
//	}
//	
//	// create constraints
//	
//	// each child is matched with at most one family
//	for (size_t i = 0; i < allo.nbChildren; ++i)
//	{
//		SCIP_CONS* cons;
//		namebuf.str("");
//		namebuf << "child_" << i;
//		SCIP_CALL_EXC(SCIPcreateConsLinear(_scip, &cons, namebuf.str().c_str(), 0, NULL, NULL, 0.0, 1.0,
//			TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
//
//		for (size_t j = 0; j < allo.children[i].nbPref; ++j)
//		{
//			for (size_t k = 0; k < allo.children[i].preferences[j].size(); ++k)
//			{
//				SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _vars[i][j][k], 1.0));
//			}
//		}
//
//		SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
//		_cons.push_back(cons);
//	}
//
//	// each family is matched with at most one child
//	for (size_t j = 0; j < allo.nbFamilies; ++j)
//	{
//		SCIP_CONS* cons;
//		namebuf.str("");
//		namebuf << "Family_" << j;
//		SCIP_CALL_EXC(SCIPcreateConsLinear(_scip, &cons, namebuf.str().c_str(), 0, NULL, NULL, 0.0, 1.0,
//			TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
//
//		for (size_t i = 0; i < allo.families[j].nbPref; ++i)
//		{
//			for (size_t k = 0; k < allo.families[j].preferences[i].size(); ++k)
//			{
//				SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _vars[allo.families[j].preferences[i][k]][allo.families[j].ranks[i][k]][allo.families[j].positions[i][k]], 1.0));
//			}
//		}
//
//		SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
//		_cons.push_back(cons);
//	}
//	
//	//constraints stable
//	for (size_t i = 0; i < allo.nbChildren; ++i)
//	{
//		for (size_t j = 0; j < allo.children[i].nbPref; ++j)
//		{
//			for (size_t k = 0; k < allo.children[i].preferences[j].size(); ++k)
//			{
//				SCIP_CONS* cons;
//				namebuf.str("");
//				namebuf << "Stable_" << i << "_" << j << "_" << k;
//				SCIP_CALL_EXC(SCIPcreateConsLinear(_scip, &cons, namebuf.str().c_str(), 0, NULL, NULL, 1.0, 2.0,
//					TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
//
//
//				for (size_t l = 0; l < j + 1; ++l)
//				{
//					
//					for (size_t m = 0; m < allo.children[i].preferences[l].size(); ++m)
//					{						
//						SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _vars[i][l][m], 1.0));
//					}
//				}
//
//				for (size_t l = 0; l < allo.children[i].ranks[j][k]+1; ++l)
//				{
//					for (size_t m = 0; m < allo.families[allo.children[i].preferences[j][k]].preferences[l].size(); ++m)
//					{
//						SCIP_CALL_EXC(SCIPaddCoefLinear(
//							_scip,
//							cons,
//							_vars[
//								allo.families[allo.children[i].preferences[j][k]].preferences[l][m]
//							][
//								allo.families[allo.children[i].preferences[j][k]].ranks[l][m]
//							][
//								allo.families[allo.children[i].preferences[j][k]].positions[l][m]
//							],
//							1.0));
//					}
//				}
//
//
//				SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
//				_cons.push_back(cons);
//			}
//		}
//	}	
//	
//}

 /* display the solution */
void scipexamples::IP::disp(Allocation& allo, std::ostream& out)
{
	// get the best found solution from scip
	 SCIP_SOL* sol = SCIPgetBestSol(_scip);
	//SCIP_SOL** sol = SCIPgetSols(_scip);
	
	out << "solution for Stable matching" << endl << endl;

	// when SCIP did not succeed then sol is NULL
	if (sol == NULL)
	{
		out << "no solution found" << endl;
		return;
	}


	allo.assignmentByChild.resize(allo.nbChildren, -1);
	allo.assignmentByFamily.resize(allo.nbFamilies, -1);

	for (int i = 0; i < allo.nbChildren; i++) {
		for (int j = 0; j < allo.children[i].nbPref; j++) {
			for (int k = 0; k < allo.children[i].preferences[j].size(); k++) {
				if (SCIPgetSolVal(_scip, sol, _vars[i][j][k]) > 0.5) {
					int idxFam = allo.children[i].preferences[j][k];
					allo.assignmentByChild[i] = idxFam;
					allo.assignmentByFamily[idxFam] = i;
				}
			}
		}
	}


	for (size_t i = 0; i < n1; ++i)
	{

		for (size_t j = 0; j < allo.children[i].nbPref; ++j)
		{
			for (size_t k = 0; k < allo.children[i].preferences[j].size(); ++k)
			{
				if (SCIPgetSolVal(_scip, sol, _vars[i][j][k]) > 0.5)
				{
					out << "Child_" << i << " with Family" << allo.children[i].preferences[j][k] << endl;
				}

			}
		}
	}


}



/* destructor */
//scipexamples::IP::~IP(void)
//{
//	std::vector<std::vector<std::vector<size_t>>> sizeAllo(_vars.size());
//	for (size_t i = 0; i < _vars.size(); ++i) {
//		sizeAllo[i].resize(_vars[i].size());
//		for (size_t j = 0; j < _vars[i].size(); ++j) {
//			sizeAllo[i][j].resize(_vars[i][j].size());
//		}
//	}
//	
//   // since the SCIPcreateVar captures all variables, we have to release them now
//   for( size_t i = 0; i < sizeAllo.size() ; ++i )
//   {
//      for ( size_t j = 0; j < sizeAllo[i].size(); ++j )
//		  for (size_t k = 0; k < sizeAllo[i][j].size(); ++k) {
//			  SCIP_CALL_EXC(SCIPreleaseVar(_scip, &_vars[i][j][k])); /*lint !e1551 !e1546*/
//		  }
//   }
//   _vars.clear(); /*lint !e1551*/
//
//   // the same for the constraints
//   for( vector<SCIP_CONS *>::size_type i = 0; i < _cons.size(); ++i ) /*lint !e1551*/
//      SCIP_CALL_EXC( SCIPreleaseCons(_scip, &_cons[i]) ); /*lint !e1551 !e1546*/
//   _cons.clear(); /*lint !e1551*/
//
//   // after releasing all vars and cons we can free the scip problem
//   // remember this has allways to be the last call to scip
//   SCIP_CALL_EXC( SCIPfree(&_scip) ); /*lint !e1551 !e1546*/
//}


/* destructor */
scipexamples::IP::~IP(void)
{
	std::vector<std::vector<std::vector<size_t>>> sizeAllo(_vars.size());
	for (size_t i = 0; i < _vars.size(); ++i) {
		sizeAllo[i].resize(_vars[i].size());
		for (size_t j = 0; j < _vars[i].size(); ++j) {
			sizeAllo[i][j].resize(_vars[i][j].size());
		}
	}
	std::vector<std::vector<size_t>> sizeYF(_varYF.size());
	for (size_t i = 0; i < _varYF.size(); ++i) {
		sizeYF[i].resize(_varYF[i].size());

	}

	// since the SCIPcreateVar captures all variables, we have to release them now
	for (size_t i = 0; i < sizeAllo.size(); ++i)
	{
		for (size_t j = 0; j < sizeAllo[i].size(); ++j)
			for (size_t k = 0; k < sizeAllo[i][j].size(); ++k) {
				SCIP_CALL_EXC(SCIPreleaseVar(_scip, &_vars[i][j][k])); /*lint !e1551 !e1546*/
			}
	}
	_vars.clear(); /*lint !e1551*/

	for (size_t i = 0; i < sizeAllo.size(); ++i)
	{
		for (size_t j = 0; j < sizeAllo[i].size(); ++j)

			SCIP_CALL_EXC(SCIPreleaseVar(_scip, &_varYC[i][j])); /*lint !e1551 !e1546*/

	}
	_varYC.clear();

	for (size_t i = 0; i < sizeYF.size(); ++i)
	{
		for (size_t j = 0; j < sizeYF[i].size(); ++j)

			SCIP_CALL_EXC(SCIPreleaseVar(_scip, &_varYF[i][j])); /*lint !e1551 !e1546*/

	}
	_varYF.clear();
	// the same for the constraints
	for (vector<SCIP_CONS*>::size_type i = 0; i < _cons.size(); ++i) /*lint !e1551*/
		SCIP_CALL_EXC(SCIPreleaseCons(_scip, &_cons[i])); /*lint !e1551 !e1546*/
	_cons.clear(); /*lint !e1551*/

	// after releasing all vars and cons we can free the scip problem
	// remember this has allways to be the last call to scip
	SCIP_CALL_EXC(SCIPfree(&_scip)); /*lint !e1551 !e1546*/
}




/* solve the stable matching problem */
void IP::solve(void)
{
   // this tells scip to start the solution process
   SCIP_CALL_EXC( SCIPsolve(_scip) );
}

/////////////////////////////////////////////////

////////////////2_YESBIN_1STA_NOMERGED/////////////////////////
/////////////////////////////////////////////////


//scipexamples::IP::IP(Allocation& allo)
//	: _scip(0), n1(allo.nbChildren), n2(allo.nbFamilies), _cons()
//{
//
//	// initialize scip
//	SCIP_CALL_EXC(SCIPcreate(&_scip));
//
//	// load default plugins linke separators, heuristics, etc.
//	SCIP_CALL_EXC(SCIPincludeDefaultPlugins(_scip));
//
//	// disable scip output to stdout
//	SCIPmessagehdlrSetQuiet(SCIPgetMessagehdlr(_scip), TRUE);
//
//	// create an empty problem
//	SCIP_CALL_EXC(SCIPcreateProb(_scip, "SMTI", NULL, NULL, NULL, NULL, NULL, NULL, NULL));
//
//	// set the objective sense to maximize, default is minimize
//	SCIP_CALL_EXC(SCIPsetObjsense(_scip, SCIP_OBJSENSE_MAXIMIZE));
//
//
//
//	//create a binary variable xij, yC
//	ostringstream namebuf;
//	_vars.resize(allo.nbChildren);
//	_varYC.resize(allo.nbChildren);
//	for (size_t i = 0; i < allo.nbChildren; ++i)
//	{
//		_vars[i].resize(allo.children[i].nbPref);
//		_varYC[i].resize(allo.children[i].nbPref);
//		for (size_t j = 0; j < allo.children[i].nbPref; ++j)
//		{
//
//			SCIP_VAR* varYC;
//			namebuf.str("");
//			namebuf << "yC#" << i << "#" << j;
//
//			// create the SCIP_VAR object
//			SCIP_CALL_EXC(SCIPcreateVar(_scip, &varYC, namebuf.str().c_str(), 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));
//
//			// add the SCIP_VAR object to the scip problem
//			SCIP_CALL_EXC(SCIPaddVar(_scip, varYC));
//
//			// storing the SCIP_VAR pointer for later access
//			_varYC[i][j] = varYC;
//			
//			//////////////////////////////////////////////////////////////
//			//////////////////////////////////////////////////////////////
//
//			_vars[i][j].resize(allo.children[i].preferences[j].size());
//			for (size_t k = 0; k < allo.children[i].preferences[j].size(); ++k)
//			{
//				SCIP_VAR* var;
//				namebuf.str("");
//				namebuf << "x#" << i << "#" << j << "#" << k;
//
//				// create the SCIP_VAR object
//				SCIP_CALL_EXC(SCIPcreateVar(_scip, &var, namebuf.str().c_str(), 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));
//
//				// add the SCIP_VAR object to the scip problem
//				SCIP_CALL_EXC(SCIPaddVar(_scip, var));
//
//				// storing the SCIP_VAR pointer for later access
//				_vars[i][j][k] = var;
//			}
//		}
//	}
//
//	//create a binary variable xij, yC, yF
//	_varYF.resize(allo.nbFamilies);
//	for (size_t i = 0; i < allo.nbFamilies; ++i)
//	{
//		_varYF[i].resize(allo.families[i].nbPref);
//		for (size_t j = 0; j < allo.families[i].nbPref; ++j)
//		{
//
//			SCIP_VAR* varYF;
//			namebuf.str("");
//			namebuf << "yF#" << i << "#" << j;
//
//			// create the SCIP_VAR object
//			SCIP_CALL_EXC(SCIPcreateVar(_scip, &varYF, namebuf.str().c_str(), 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));
//
//			// add the SCIP_VAR object to the scip problem
//			SCIP_CALL_EXC(SCIPaddVar(_scip, varYF));
//
//			// storing the SCIP_VAR pointer for later access
//			_varYF[i][j] = varYF;
//		}
//	}
//	// create constraints
//	
//
//	// init YC
//	for (size_t i = 0; i < allo.nbChildren; ++i)
//	{
//		/*SCIP_CONS* cons;
//		namebuf.str("");
//		namebuf << "child_" << i;
//		SCIP_CALL_EXC(SCIPcreateConsLinear(_scip, &cons, namebuf.str().c_str(), 0, NULL, NULL, 0.0, 0.0,
//			TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));*/
//
//		for (size_t j = 0; j < allo.children[i].nbPref; ++j)
//		{
//			SCIP_CONS* cons;
//			namebuf.str("");
//			namebuf << "child_" << i<<"_rank_"<<j;
//			SCIP_CALL_EXC(SCIPcreateConsLinear(_scip, &cons, namebuf.str().c_str(), 0, NULL, NULL, 0.0, 0.0,
//				TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
//			for (size_t k = 0; k < allo.children[i].preferences[j].size(); ++k)
//			{
//				SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _vars[i][j][k], 1.0));
//			}
//			if (j > 0) {
//				SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _varYC[i][j-1], 1.0));
//			}
//			SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _varYC[i][j], -1.0));
//			SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
//			_cons.push_back(cons);
//		}
//
//		/*SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
//		_cons.push_back(cons);*/
//	}
//
//	//init yF
//	for (size_t j = 0; j < allo.nbFamilies; ++j)
//	{
//		/*SCIP_CONS* cons;
//		namebuf.str("");
//		namebuf << "child_" << i;
//		SCIP_CALL_EXC(SCIPcreateConsLinear(_scip, &cons, namebuf.str().c_str(), 0, NULL, NULL, 0.0, 0.0,
//			TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));*/
//
//		for (size_t i = 0; i < allo.families[j].nbPref; ++i)
//		{
//			SCIP_CONS* cons;
//			namebuf.str("");
//			namebuf << "Fam_" << j << "_rank_" << i;
//			SCIP_CALL_EXC(SCIPcreateConsLinear(_scip, &cons, namebuf.str().c_str(), 0, NULL, NULL, 0.0, 0.0,
//				TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
//			for (size_t k = 0; k < allo.families[j].preferences[i].size(); ++k)
//			{
//				SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _vars[allo.families[j].preferences[i][k]][allo.families[j].ranks[i][k]][allo.families[j].positions[i][k]], 1.0));
//			}
//			if (i > 0) {
//				SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _varYF[j][i - 1], 1.0));
//			}
//			SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _varYF[j][i], -1.0));
//			SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
//			_cons.push_back(cons);
//		}
//
//		/*SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
//		_cons.push_back(cons);*/
//	}
//
//	//stable constraint
//	for (size_t i = 0; i < allo.nbChildren; ++i)
//	{
//		for (size_t j = 0; j < allo.children[i].nbPref; ++j)
//		{
//			for (size_t k = 0; k < allo.children[i].preferences[j].size(); ++k)
//			{
//				SCIP_CONS* cons;
//				namebuf.str("");
//				namebuf << "Stable_" << i << "_" << j << "_" << k;
//				SCIP_CALL_EXC(SCIPcreateConsLinear(_scip, &cons, namebuf.str().c_str(), 0, NULL, NULL, 1.0, 2.0,
//					TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
//
//
//				SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _varYC[i][j], 1.0));
//
//				SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _varYF[allo.children[i].preferences[j][k]][allo.children[i].ranks[j][k]], 1.0));
//
//				SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
//				_cons.push_back(cons);
//			}
//		}
//	}
//
//
//	// diag col up
//
//}



/////////////////////////////////////////////////

////////////////3_NOBIN_1STA_YESMERGED/////////////////////////
/////////////////////////////////////////////////




//scipexamples::IP::IP(Allocation& allo)
//	: _scip(0), n1(allo.nbChildren), n2(allo.nbFamilies), _cons()
//{
//
//	// initialize scip
//	SCIP_CALL_EXC(SCIPcreate(&_scip));
//
//	// load default plugins linke separators, heuristics, etc.
//	SCIP_CALL_EXC(SCIPincludeDefaultPlugins(_scip));
//
//	// disable scip output to stdout
//	SCIPmessagehdlrSetQuiet(SCIPgetMessagehdlr(_scip), TRUE);
//
//	// create an empty problem
//	SCIP_CALL_EXC(SCIPcreateProb(_scip, "SMTI", NULL, NULL, NULL, NULL, NULL, NULL, NULL));
//
//	// set the objective sense to maximize, default is minimize
//	SCIP_CALL_EXC(SCIPsetObjsense(_scip, SCIP_OBJSENSE_MAXIMIZE));
//
//	
//
//	 //create a binary variable xij
//	ostringstream namebuf;
//	_vars.resize(allo.nbChildren);
//
//	for (size_t i = 0; i < allo.nbChildren; ++i)
//	{
//		_vars[i].resize(allo.children[i].nbPref);
//
//		for (size_t j = 0; j < allo.children[i].nbPref; ++j)
//		{
//			_vars[i][j].resize(allo.children[i].preferences[j].size());
//			for (size_t k = 0; k < allo.children[i].preferences[j].size(); ++k)
//			{
//				SCIP_VAR* var;
//				namebuf.str("");
//				namebuf << "x#" << i << "#" << j << "#" << k;
//
//				// create the SCIP_VAR object
//				SCIP_CALL_EXC(SCIPcreateVar(_scip, &var, namebuf.str().c_str(), 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));
//
//				// add the SCIP_VAR object to the scip problem
//				SCIP_CALL_EXC(SCIPaddVar(_scip, var));
//
//				// storing the SCIP_VAR pointer for later access
//				_vars[i][j][k] = var;
//			}
//		}
//	}
//	
//	// create constraints
//	
//	//each child is matched with at most one family
//	for (size_t i = 0; i < allo.nbChildren; ++i)
//	{
//		SCIP_CONS* cons;
//		namebuf.str("");
//		namebuf << "child_" << i;
//		SCIP_CALL_EXC(SCIPcreateConsLinear(_scip, &cons, namebuf.str().c_str(), 0, NULL, NULL, 0.0, 1.0,
//			TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
//
//		for (size_t j = 0; j < allo.children[i].nbPref; ++j)
//		{
//			for (size_t k = 0; k < allo.children[i].preferences[j].size(); ++k)
//			{
//				SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _vars[i][j][k], 1.0));
//			}
//		}
//
//		SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
//		_cons.push_back(cons);
//	}
//
//	//each family is matched with at most one child
//	for (size_t j = 0; j < allo.nbFamilies; ++j)
//	{
//		SCIP_CONS* cons;
//		namebuf.str("");
//		namebuf << "Family_" << j;
//		SCIP_CALL_EXC(SCIPcreateConsLinear(_scip, &cons, namebuf.str().c_str(), 0, NULL, NULL, 0.0, 1.0,
//			TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
//
//		for (size_t i = 0; i < allo.families[j].nbPref; ++i)
//		{
//			for (size_t k = 0; k < allo.families[j].preferences[i].size(); ++k)
//			{
//				SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _vars[allo.families[j].preferences[i][k]][allo.families[j].ranks[i][k]][allo.families[j].positions[i][k]], 1.0));
//			}
//		}
//
//		SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
//		_cons.push_back(cons);
//	}
//	
//	//constraints stable
//	for (size_t i = 0; i < allo.nbChildren; ++i)
//	{
//		for (size_t j = 0; j < allo.children[i].nbPref; ++j)
//		{
//
//			SCIP_CONS* cons;
//			namebuf.str("");
//			namebuf << "Stable_" << i << "_rank_" << j ;
//			SCIP_CALL_EXC(SCIPcreateConsLinear(_scip, &cons, namebuf.str().c_str(), 0, NULL, NULL, double(allo.children[i].preferences[j].size()), double(2*allo.children[i].preferences[j].size()),
//				TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
//
//
//			for (size_t l = 0; l < j + 1; ++l)
//			{
//
//				for (size_t m = 0; m < allo.children[i].preferences[l].size(); ++m)
//				{
//					SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _vars[i][l][m], double(allo.children[i].preferences[j].size())));
//				}
//			}
//			for (size_t k = 0; k < allo.children[i].preferences[j].size(); ++k) {
//				for (size_t l = 0; l < allo.children[i].ranks[j][k] + 1; ++l)
//				{
//					for (size_t m = 0; m < allo.families[allo.children[i].preferences[j][k]].preferences[l].size(); ++m)
//					{
//						SCIP_CALL_EXC(SCIPaddCoefLinear(
//							_scip,
//							cons,
//							_vars[
//								allo.families[allo.children[i].preferences[j][k]].preferences[l][m]
//							][
//								allo.families[allo.children[i].preferences[j][k]].ranks[l][m]
//							][
//								allo.families[allo.children[i].preferences[j][k]].positions[l][m]
//							],
//									1.0));
//					}
//
//
//
//					
//				}
//			}
//			SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
//			_cons.push_back(cons);
//		}
//	}
//
//	
//	
//	
//}




/////////////////////////////////////////////////

////////////////4_YESBIN_1STA_YESMERGED/////////////////////////
/////////////////////////////////////////////////




//scipexamples::IP::IP(Allocation& allo)
//	: _scip(0), n1(allo.nbChildren), n2(allo.nbFamilies), _cons()
//{
//
//	// initialize scip
//	SCIP_CALL_EXC(SCIPcreate(&_scip));
//
//	// load default plugins linke separators, heuristics, etc.
//	SCIP_CALL_EXC(SCIPincludeDefaultPlugins(_scip));
//
//	// disable scip output to stdout
//	SCIPmessagehdlrSetQuiet(SCIPgetMessagehdlr(_scip), TRUE);
//
//	// create an empty problem
//	SCIP_CALL_EXC(SCIPcreateProb(_scip, "SMTI", NULL, NULL, NULL, NULL, NULL, NULL, NULL));
//
//	// set the objective sense to maximize, default is minimize
//	SCIP_CALL_EXC(SCIPsetObjsense(_scip, SCIP_OBJSENSE_MAXIMIZE));
//
//
//
//	//create a binary variable xij, yC
//	ostringstream namebuf;
//	_vars.resize(allo.nbChildren);
//	_varYC.resize(allo.nbChildren);
//	for (size_t i = 0; i < allo.nbChildren; ++i)
//	{
//		_vars[i].resize(allo.children[i].nbPref);
//		_varYC[i].resize(allo.children[i].nbPref);
//		for (size_t j = 0; j < allo.children[i].nbPref; ++j)
//		{
//
//			SCIP_VAR* varYC;
//			namebuf.str("");
//			namebuf << "yC#" << i << "#" << j;
//
//			// create the SCIP_VAR object
//			SCIP_CALL_EXC(SCIPcreateVar(_scip, &varYC, namebuf.str().c_str(), 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));
//
//			// add the SCIP_VAR object to the scip problem
//			SCIP_CALL_EXC(SCIPaddVar(_scip, varYC));
//
//			// storing the SCIP_VAR pointer for later access
//			_varYC[i][j] = varYC;
//
//			//////////////////////////////////////////////////////////////
//			//////////////////////////////////////////////////////////////
//
//			_vars[i][j].resize(allo.children[i].preferences[j].size());
//			for (size_t k = 0; k < allo.children[i].preferences[j].size(); ++k)
//			{
//				SCIP_VAR* var;
//				namebuf.str("");
//				namebuf << "x#" << i << "#" << j << "#" << k;
//
//				// create the SCIP_VAR object
//				SCIP_CALL_EXC(SCIPcreateVar(_scip, &var, namebuf.str().c_str(), 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));
//
//				// add the SCIP_VAR object to the scip problem
//				SCIP_CALL_EXC(SCIPaddVar(_scip, var));
//
//				// storing the SCIP_VAR pointer for later access
//				_vars[i][j][k] = var;
//			}
//		}
//	}
//
//	//create a binary variable xij, yC, yF
//	_varYF.resize(allo.nbFamilies);
//	for (size_t i = 0; i < allo.nbFamilies; ++i)
//	{
//		_varYF[i].resize(allo.families[i].nbPref);
//		for (size_t j = 0; j < allo.families[i].nbPref; ++j)
//		{
//
//			SCIP_VAR* varYF;
//			namebuf.str("");
//			namebuf << "yF#" << i << "#" << j;
//
//			// create the SCIP_VAR object
//			SCIP_CALL_EXC(SCIPcreateVar(_scip, &varYF, namebuf.str().c_str(), 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));
//
//			// add the SCIP_VAR object to the scip problem
//			SCIP_CALL_EXC(SCIPaddVar(_scip, varYF));
//
//			// storing the SCIP_VAR pointer for later access
//			_varYF[i][j] = varYF;
//		}
//	}
//	// create constraints
//
//
//	// init YC
//	for (size_t i = 0; i < allo.nbChildren; ++i)
//	{
//		/*SCIP_CONS* cons;
//		namebuf.str("");
//		namebuf << "child_" << i;
//		SCIP_CALL_EXC(SCIPcreateConsLinear(_scip, &cons, namebuf.str().c_str(), 0, NULL, NULL, 0.0, 0.0,
//			TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));*/
//
//		for (size_t j = 0; j < allo.children[i].nbPref; ++j)
//		{
//			SCIP_CONS* cons;
//			namebuf.str("");
//			namebuf << "child_" << i << "_rank_" << j;
//			SCIP_CALL_EXC(SCIPcreateConsLinear(_scip, &cons, namebuf.str().c_str(), 0, NULL, NULL, 0.0, 0.0,
//				TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
//			for (size_t k = 0; k < allo.children[i].preferences[j].size(); ++k)
//			{
//				SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _vars[i][j][k], 1.0));
//			}
//			if (j > 0) {
//				SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _varYC[i][j - 1], 1.0));
//			}
//			SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _varYC[i][j], -1.0));
//			SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
//			_cons.push_back(cons);
//		}
//
//		/*SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
//		_cons.push_back(cons);*/
//	}
//
//	//init yF
//	for (size_t j = 0; j < allo.nbFamilies; ++j)
//	{
//		/*SCIP_CONS* cons;
//		namebuf.str("");
//		namebuf << "child_" << i;
//		SCIP_CALL_EXC(SCIPcreateConsLinear(_scip, &cons, namebuf.str().c_str(), 0, NULL, NULL, 0.0, 0.0,
//			TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));*/
//
//		for (size_t i = 0; i < allo.families[j].nbPref; ++i)
//		{
//			SCIP_CONS* cons;
//			namebuf.str("");
//			namebuf << "Fam_" << j << "_rank_" << i;
//			SCIP_CALL_EXC(SCIPcreateConsLinear(_scip, &cons, namebuf.str().c_str(), 0, NULL, NULL, 0.0, 0.0,
//				TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
//			for (size_t k = 0; k < allo.families[j].preferences[i].size(); ++k)
//			{
//				SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _vars[allo.families[j].preferences[i][k]][allo.families[j].ranks[i][k]][allo.families[j].positions[i][k]], 1.0));
//			}
//			if (i > 0) {
//				SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _varYF[j][i - 1], 1.0));
//			}
//			SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _varYF[j][i], -1.0));
//			SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
//			_cons.push_back(cons);
//		}
//
//		/*SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
//		_cons.push_back(cons);*/
//	}
//
//	//stable constraint
//	for (size_t i = 0; i < allo.nbChildren; ++i)
//	{
//		for (size_t j = 0; j < allo.children[i].nbPref; ++j)
//		{
//
//			SCIP_CONS* cons;
//			namebuf.str("");
//			namebuf << "Stable_" << i << "_ranks_" << j ;
//			SCIP_CALL_EXC(SCIPcreateConsLinear(_scip, &cons, namebuf.str().c_str(), 0, NULL, NULL, double(allo.children[i].preferences[j].size()), double(2 * allo.children[i].preferences[j].size()),
//				TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
//
//
//			SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _varYC[i][j], double(allo.children[i].preferences[j].size())));
//			for (size_t k = 0; k < allo.children[i].preferences[j].size(); ++k)
//			{
//				SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _varYF[allo.children[i].preferences[j][k]][allo.children[i].ranks[j][k]], 1.0));
//			}
//
//			SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
//			_cons.push_back(cons);
//
//		}
//	}
//
//
//	// diag col up
//
//}
/* destructor */



/////////////////////////////////////////////////

////////////////7_NOBIN_2STA_YESMERGED/////////////////////////
/////////////////////////////////////////////////


scipexamples::IP::IP(Allocation& allo)
	: _scip(0), n1(allo.nbChildren), n2(allo.nbFamilies), _cons()
{

	// initialize scip
	SCIP_CALL_EXC(SCIPcreate(&_scip));

	// load default plugins linke separators, heuristics, etc.
	SCIP_CALL_EXC(SCIPincludeDefaultPlugins(_scip));

	// disable scip output to stdout
	SCIPmessagehdlrSetQuiet(SCIPgetMessagehdlr(_scip), TRUE);

	// create an empty problem
	SCIP_CALL_EXC(SCIPcreateProb(_scip, "SMTI", NULL, NULL, NULL, NULL, NULL, NULL, NULL));

	// set the objective sense to maximize, default is minimize
	SCIP_CALL_EXC(SCIPsetObjsense(_scip, SCIP_OBJSENSE_MAXIMIZE));

	

	 //create a binary variable xij
	ostringstream namebuf;
	_vars.resize(allo.nbChildren);

	for (size_t i = 0; i < allo.nbChildren; ++i)
	{
		_vars[i].resize(allo.children[i].nbPref);

		for (size_t j = 0; j < allo.children[i].nbPref; ++j)
		{
			_vars[i][j].resize(allo.children[i].preferences[j].size());
			for (size_t k = 0; k < allo.children[i].preferences[j].size(); ++k)
			{
				SCIP_VAR* var;
				namebuf.str("");
				namebuf << "x#" << i << "#" << j << "#" << k;

				// create the SCIP_VAR object
				SCIP_CALL_EXC(SCIPcreateVar(_scip, &var, namebuf.str().c_str(), 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));

				// add the SCIP_VAR object to the scip problem
				SCIP_CALL_EXC(SCIPaddVar(_scip, var));

				// storing the SCIP_VAR pointer for later access
				_vars[i][j][k] = var;
			}
		}
	}
	
	// create constraints
	
	//each child is matched with at most one family
	for (size_t i = 0; i < allo.nbChildren; ++i)
	{
		SCIP_CONS* cons;
		namebuf.str("");
		namebuf << "child_" << i;
		SCIP_CALL_EXC(SCIPcreateConsLinear(_scip, &cons, namebuf.str().c_str(), 0, NULL, NULL, 0.0, 1.0,
			TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));

		for (size_t j = 0; j < allo.children[i].nbPref; ++j)
		{
			for (size_t k = 0; k < allo.children[i].preferences[j].size(); ++k)
			{
				SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _vars[i][j][k], 1.0));
			}
		}

		SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
		_cons.push_back(cons);
	}

	//each family is matched with at most one child
	for (size_t j = 0; j < allo.nbFamilies; ++j)
	{
		SCIP_CONS* cons;
		namebuf.str("");
		namebuf << "Family_" << j;
		SCIP_CALL_EXC(SCIPcreateConsLinear(_scip, &cons, namebuf.str().c_str(), 0, NULL, NULL, 0.0, 1.0,
			TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));

		for (size_t i = 0; i < allo.families[j].nbPref; ++i)
		{
			for (size_t k = 0; k < allo.families[j].preferences[i].size(); ++k)
			{
				SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _vars[allo.families[j].preferences[i][k]][allo.families[j].ranks[i][k]][allo.families[j].positions[i][k]], 1.0));
			}
		}

		SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
		_cons.push_back(cons);
	}
	
	//constraints stable
	for (size_t i = 0; i < allo.nbChildren; ++i)
	{
		for (size_t j = 0; j < allo.children[i].nbPref; ++j)
		{

			SCIP_CONS* cons;
			namebuf.str("");
			namebuf << "Stable_" << i << "_rank_" << j ;
			SCIP_CALL_EXC(SCIPcreateConsLinear(_scip, &cons, namebuf.str().c_str(), 0, NULL, NULL, double(allo.children[i].preferences[j].size()), double(2*allo.children[i].preferences[j].size()),
				TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));


			for (size_t l = 0; l < j + 1; ++l)
			{

				for (size_t m = 0; m < allo.children[i].preferences[l].size(); ++m)
				{
					SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _vars[i][l][m], double(allo.children[i].preferences[j].size())));
				}
			}
			for (size_t k = 0; k < allo.children[i].preferences[j].size(); ++k) {
				for (size_t l = 0; l < allo.children[i].ranks[j][k] + 1; ++l)
				{
					for (size_t m = 0; m < allo.families[allo.children[i].preferences[j][k]].preferences[l].size(); ++m)
					{
						SCIP_CALL_EXC(SCIPaddCoefLinear(
							_scip,
							cons,
							_vars[
								allo.families[allo.children[i].preferences[j][k]].preferences[l][m]
							][
								allo.families[allo.children[i].preferences[j][k]].ranks[l][m]
							][
								allo.families[allo.children[i].preferences[j][k]].positions[l][m]
							],
									1.0));
					}



					
				}
			}
			SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
			_cons.push_back(cons);
		}
	}

	//constraints stable
	for (size_t j = 0; j < allo.nbFamilies; ++j)
	{
		for (size_t i = 0; i < allo.families[j].nbPref; ++i)
		{

			SCIP_CONS* cons;
			namebuf.str("");
			namebuf << "Stable_" << i << "_rank_" << j;
			SCIP_CALL_EXC(SCIPcreateConsLinear(_scip, &cons, namebuf.str().c_str(), 0, NULL, NULL, double(allo.families[j].preferences[i].size()), double(2 * allo.families[j].preferences[i].size()),
				TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));


			for (size_t l = 0; l < i + 1; ++l)
			{

				for (size_t m = 0; m < allo.families[j].preferences[l].size(); ++m)
				{
					SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _vars[allo.families[j].preferences[l][m]][allo.families[j].ranks[l][m]][allo.families[j].positions[l][m]], double(allo.families[j].preferences[i].size())));
				}
			}
			for (size_t k = 0; k < allo.families[j].preferences[i].size(); ++k) {
				for (size_t l = 0; l < allo.families[j].ranks[i][k] + 1; ++l)
				{
					for (size_t m = 0; m < allo.children[allo.families[j].preferences[i][k]].preferences[l].size(); ++m)
					{
						SCIP_CALL_EXC(SCIPaddCoefLinear(
							_scip,
							cons,
							_vars[
								allo.families[j].preferences[l][m]
							][
								allo.children[].ranks[l][m]
							][
								allo.children[allo.families[j].preferences[i][k]].positions[l][m]
							],
									1.0));
					}




				}
			}
			SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
			_cons.push_back(cons);
		}
	}

	
	
	
}