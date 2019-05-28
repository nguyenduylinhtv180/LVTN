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

/**@file queens_main.cpp
 * @brief main file for the queens example
 * @author Cornelius Schwarz
 */

#include "IP.hpp"
#include <cstdlib>
#include <iostream>
#include "scip_exception.hpp"



using namespace std;
using namespace scipexamples;


/** main function for queens example */
int
main(
     int args,
     char ** argv
     )
{
   cout << "********************************************" << endl;
   cout << "* n-queens solver based on SCIP            *" << endl;
   cout << "*                                          *" << endl;
   cout << "* (c) Cornelius Schwarz (2007)             *" << endl;
   cout << "********************************************" << endl << endl;

   //if (args < 2)
   //{
   //   cerr << "call " << argv[0] << " <number of queens>" << endl;
   //   exit(1);
   //}

   //// get the number of queens for commandline
   //size_t n = abs(atoi(argv[1])); /*lint !e732*/

	// local variables
   Allocation allo;
   //string filein = argv[2];
   string filein = "Text.txt";
   string path = "";
   string pathAndFileout = "Output.txt";

   // functions
   try
   {
	   allo.load(path, filein);
	   allo.reduction();

	   

	  

	   // initialize the queens solver
	   IP solver(allo);

	   // solve the queens problem
	   solver.solve();

	   // display the solution on stdout
	   solver.disp(allo);
	   allo.printSol();
	   allo.checkSolution();
	   cout << "aaaa" << endl;
	   allo.printInfo(pathAndFileout);

   }
   catch (const SCIPException & exc)
   {
	   cout << "BBBBB" << endl;
      cerr << exc.what() << endl;
      exit((int) exc.getRetcode());
   }




   system("pause");
   return EXIT_SUCCESS;



   //size_t n;
   //cout << "nhap vao n" << endl;
   //cin >> n;
   //try
   //{
   //   // initialize the queens solver
   //   QueensSolver solver(n);

   //   // solve the queens problem
   //   solver.solve();

   //   // display the solution on stdout
   //   solver.disp();

   //} catch(const SCIPException& exc)
   //{
   //   cerr << exc.what() << endl;
   //   exit((int) exc.getRetcode());
   //}
   //system("pause");
   //return EXIT_SUCCESS;
}
