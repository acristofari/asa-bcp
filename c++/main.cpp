// -------------------------------------------------------------------------
//
// This file is part of ASA-BCP, which is a solver for bound-constrained
// optimization problems of the following form:
//
//                                 min f(x)
//                           s.t. l <= x <= u
//
// with f(x) twice continuously differentiable.
//
// This is a driver for running ASA-BCP on user-defined problems.
// See the file 'README.txt' to know how to run the program.
//
// -------------------------------------------------------------------------
//
// Reference paper:
//
// A. Cristofari, M. De Santis, S. Lucidi, F. Rinaldi (2017). A Two-Stage
// Active-Set Algorithm for Bound-Constrained Optimization. Journal of
// Optimization Theory and Applications, 172(2), 369-401.
//
// -------------------------------------------------------------------------
//
// Authors:
// Andrea Cristofari (e-mail: andrea.cristofari@unipd.it)
// Marianna De Santis (e-mail: mdesantis@diag.uniroma1.it)
// Stefano Lucidi (e-mail: lucidi@diag.uniroma1.it)
// Francesco Rinaldi (e-mail: rinaldi@math.unipd.it)
//
// Last update of this file:
// December 18th, 2020
//
// Copyright 2017-2020 Andrea Cristofari, Marianna De Santis,
// Stefano Lucidi, Francesco Rinaldi.
//
// Licensing:
// This file is part of ASA-BCP.
// ASA-BCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ASA-BCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with ASA-BCP. If not, see <http://www.gnu.org/licenses/>.
//
// -------------------------------------------------------------------------


#include "problem.h"
#include "asa_bcp.h"
#include <iostream>
#include <time.h>
#include <vector>
#include <algorithm>

int main() {

    unsigned short int status;
    float elap_time;
    std::ofstream file_stats,file_opt_sol;
    clock_t t_start,t_end;

    std::ios_base::sync_with_stdio(false);
    std::cout.precision(15);

    // build an object of type Problem
    Problem p(status);

    // check if an error occurred when calling the Problem constructor
    // (something went wrong if 'status' > 0)
    if (status > 0) {
        std::cout << "\nerror when calling the Problem constructor\n";
        return 1;
    }

    // build an object of type Asa_bcp
    Asa_bcp alg(status,&p);
    
    // check if an error occurred when calling the Asa_bcp constructor
    // (something went wrong if 'status' > 0)
    if (status > 0) {
        std::cout << "\nerror when calling the Asa_bcp constructor\n";
        return 1;
    }

    // ------------------------------------------------------------------------------------------------------------
    // *** EXAMPLE OF HOW TO CHANGE ASA-BCP PARAMETERS ***
    // (see the description of Asa_bcp in the file 'asa_bcp.cpp' to know which parameters can be changed and their
    // default values)
    //
    // Instead of creating the object 'alg' by the above instruction 'Asa_bcp alg(status,&p);', do the following:
    //
    // (1) create an object of structure type asa_bcp_options (see its declaration in file 'asa_bcp.h'), e.g.,:
    //
    //       asa_bcp_options opts;
    //
    // (2) assign new values to (some of) the members of the new structure object, e.g.,:
    //
    //       opts.verbosity = 0;
    //
    // (3) pass the address of the structure object as third input argument when calling the Asa_bcp constructor, e.g.,:
    //
    //       Asa_bcp alg(status,&p,&opts);
    // ------------------------------------------------------------------------------------------------------------

    // call the solver
    t_start = clock();
    alg.solve();
    t_end = clock();

    // compute the elapsed time
    elap_time = std::max((float)(t_end-t_start)/CLOCKS_PER_SEC,(float)0e0);

    // write statistics to the screen and to file 'statistics.txt'
    file_stats.open("statistics.txt",std::ios::trunc);
    if (alg.get_flag() < 0) {
        file_stats << "infeasible problem\n";
        file_stats.close();
        return 0;
    }
    std::cout.precision(5);
    std::cout.setf(std::ios::scientific,std::ios::floatfield);
    file_stats.precision(5);
    file_stats.setf(std::ios::scientific,std::ios::floatfield);
    std::cout << "************************************************"
              << "\n\nAlgorithm: ASA-BCP"
              << "\n\nnumber of variables = " << *p.dim
              << "\n\nf = " << alg.get_f()
              << "\n\nsup-norm of the projected gradient = " << alg.get_sup_norm_proj_g()
              << "\nnumber of iterations = " << alg.get_it()
              << "\nnumber of function evaluations = " << alg.get_n_f()
              << "\nnumber of gradient evaluations = " << alg.get_n_g()
              << "\nnumber of Hessian-vector products = " << alg.get_n_hd() 
              << "\nnumber inner cg iterations = " << alg.get_inner_it()
              << "\nexit flag = " << alg.get_flag()
              << "\nelapsed time (s) = " << elap_time
              << "\n\n************************************************\n";
    file_stats << "************************************************"
               << "\n\nAlgorithm: ASA-BCP"
               << "\n\nnumber of variables = " << *p.dim
               << "\n\nf = " << alg.get_f()
               << "\n\nsup-norm of the projected gradient = " << alg.get_sup_norm_proj_g()
               << "\nnumber of iterations = " << alg.get_it()
               << "\nnumber of function evaluations = " << alg.get_n_f()
               << "\nnumber of gradient evaluations = " << alg.get_n_g()
               << "\nnumber of Hessian-vector products = " << alg.get_n_hd() 
               << "\nnumber inner cg iterations = " << alg.get_inner_it()
               << "\nexit flag = " << alg.get_flag()
               << "\nelapsed time (s) = " << elap_time
               << "\n\n************************************************\n";
    file_stats.close();

    // write the solution found by the algorithm to file 'opt_sol.txt'
    file_opt_sol.open("opt_sol.txt",std::ios::trunc);
    file_opt_sol.precision(5);
    file_opt_sol.setf(std::ios::scientific,std::ios::floatfield);
    const double *x_ptr = &(alg.get_x()[0]);
    for (unsigned int i=0; i<*p.dim; i++) {
        file_opt_sol << *(x_ptr+i) << "\n";
    }
    file_opt_sol.close();

    return 0;
}