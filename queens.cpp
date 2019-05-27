/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the examples to                     */
/*                         An introduction to SCIP                           */
/*                                                                           */
/*    Copyright (C) 2007 Cornelius Schwarz                                   */
/*                                                                           */
/*                  2007 University of Bayreuth                              */
/*                                                                           */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**
 * @file queens.cpp
 * @author Cornelius Schwarz
 * @brief n-queens examlple implementation
 */

#include "queens.hpp"
#include <sstream>
#include "scip_exception.hpp"
#include "scip/pub_message.h"
#include <math.h> 
#include <cmath>  
using namespace std;
using namespace scipexamples;

/* constructor */
scipexamples::QueensSolver::QueensSolver(Allocation& allo)
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

	

	 //create a binary variable for every field (i,j) on the chess board
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
	// one queen per row

	// now we create the constraints for the diagonals
	// there is only one queen allowed in each diagonal, but there do not need to be one. Therefore we add a <= 1 constraint
	// in this problem case we can set the lhs to zero instead of -SCIPinfinity
	// diag col down

	// diag row down
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
	

	for (size_t i = 0; i < allo.nbChildren; ++i)
	{
		for (size_t j = 0; j < allo.children[i].nbPref; ++j)
		{
			for (size_t k = 0; k < allo.children[i].preferences[j].size(); ++k)
			{
				SCIP_CONS* cons;
				namebuf.str("");
				namebuf << "Stable_" << i << "_" << j << "_" << k;
				SCIP_CALL_EXC(SCIPcreateConsLinear(_scip, &cons, namebuf.str().c_str(), 0, NULL, NULL, 1.0, 2.0,
					TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));


				for (size_t l = 0; l < j + 1; ++l)
				{
					
					for (size_t m = 0; m < allo.children[i].preferences[l].size(); ++m)
					{						
						SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _vars[i][l][m], 1.0));
					}
				}

				for (size_t l = 0; l < allo.children[i].ranks[j][k]+1; ++l)
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


				SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
				_cons.push_back(cons);
			}
		}
	}

	
	// diag col up
	
}

 /* display the solution */
void scipexamples::QueensSolver::disp(Allocation& allo, std::ostream& out)
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


	/*for (size_t i = 0; i < n1; ++i)
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
	}*/


}



/* destructor */
scipexamples::QueensSolver::~QueensSolver(void)
{
	std::vector<std::vector<std::vector<size_t>>> sizeAllo(_vars.size());
	for (size_t i = 0; i < _vars.size(); ++i) {
		sizeAllo[i].resize(_vars[i].size());
		for (size_t j = 0; j < _vars[i].size(); ++j) {
			sizeAllo[i][j].resize(_vars[i][j].size());
		}
	}
	
   // since the SCIPcreateVar captures all variables, we have to release them now
   for( size_t i = 0; i < sizeAllo.size() ; ++i )
   {
      for ( size_t j = 0; j < sizeAllo[i].size(); ++j )
		  for (size_t k = 0; k < sizeAllo[i][j].size(); ++k) {
			  SCIP_CALL_EXC(SCIPreleaseVar(_scip, &_vars[i][j][k])); /*lint !e1551 !e1546*/
		  }
   }
   _vars.clear(); /*lint !e1551*/

   // the same for the constraints
   for( vector<SCIP_CONS *>::size_type i = 0; i < _cons.size(); ++i ) /*lint !e1551*/
      SCIP_CALL_EXC( SCIPreleaseCons(_scip, &_cons[i]) ); /*lint !e1551 !e1546*/
   _cons.clear(); /*lint !e1551*/

   // after releasing all vars and cons we can free the scip problem
   // remember this has allways to be the last call to scip
   SCIP_CALL_EXC( SCIPfree(&_scip) ); /*lint !e1551 !e1546*/
}

/* solve the n-queens problem */
void QueensSolver::solve(void)
{
   // this tells scip to start the solution process
   SCIP_CALL_EXC( SCIPsolve(_scip) );
}
