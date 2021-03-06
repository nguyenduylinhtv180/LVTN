/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the examples to                     */
/*                         An introduction to SCIP                           */
/*                                                                           */
/*    Nguyen Duy Linh                                   */
/*                                                                           */
/*                  SMTI                              */
/*                                                                           */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  queens.hpp
 * @brief  n-queens example
 * @author Cornelius Schwarz
 */

#ifndef IP_H
#define IP_H

#include <vector>
#include <iostream>

#include <scip/scip.h>
#include <scip/scipdefplugins.h>

#include <math.h> 
#include <cmath>  


#include "E:\stable-matching\StableMatchingCodes-master\SMTI\1_NOBIN_1STA_NOMERGED\Allocation.h" 


namespace scipexamples
{
   /**@class QueensSolver
    * @brief solver class for the n-queens problem
    *
    *  this class implements a solver for the n-queens problem as an mip model, which will be solved using SCIP
    */
   class IP
   {
   private:

      /** @brief pointer to scip structure
       *
       *  SCIP organizes all the problem informations by itself, we can access them by the SCIP * pointer
       */
      SCIP * _scip;

      /** @brief number of queens  */
      size_t n1;

	  size_t n2;

      /** @brief one binary variable for each field (i,j) on the chess bord
       *
       * To access variable information (objective value, bounds,
       * etc.) use the SCIP_VAR * pointer. Since we want to know the
       * value of each variable in the solution, we have to store
       * these pointers.
       */

	  std::vector<std::vector<std::vector<SCIP_VAR*>> > _vars;

	  ///////2_YESBIN_1STA_NOMERGED//////////////

	  std::vector<std::vector<SCIP_VAR*> > _varYC;
	  std::vector<std::vector<SCIP_VAR*> > _varYF;


      /** @brief constraints for rows, cols and diags of the chess board
       *
       * To access constraint information (right hand side, left hand
       * side, dual values, etc.) use the SCIP_CONS * pointer. For the
       * n-queens problem we do not really need to store them but we
       * do for illustration.
       */
      std::vector<SCIP_CONS *> _cons;
	  

   public:
      /** @brief constructs the BP model for the n-queens problem
       *
       * the constructor builds a BP model in scip for the n-queens problem
       * @param[in] n the number of queens
       */
      IP(Allocation& allo);

      /** @brief destructor this is the place to release the SCIP_VAR
       * and SCIP_CONS pointers and to free the SCIP pointer
       * afterwards
       */
      ~IP();

      void solve(void); ///< solves the queens problem using SCIPsolve

      /** @brief display the solution
       *
       * a simplex ASCII output function to display the solution of the n - queens problem
       * @param[in,out] out ostream class for output(default cout)
       */
      void disp(Allocation& allo, std::ostream& out = std::cout);
   };
}
#endif
