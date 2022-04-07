// -------------------------------------------------------------------------
// 
// This file is part of ASA-BCP, which is a solver for bound-constrained
// optimization problems of the following form:
// 
//                                min f(x)
//                           s.t. l <= x <= u
// 
// with given vectors l, u and where f(x) is a twice continuously
// differentiable function.
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
// April 7th, 2022
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
// Copyright 2017-2022 Andrea Cristofari, Marianna De Santis,
// Stefano Lucidi, Francesco Rinaldi.
//
// -------------------------------------------------------------------------


#include "problem.h"
#include "asa_bcp.h"
#include <iostream>
#include <vector>
#include <algorithm>

// In this file, it is shown how to call ASA-BCP to solve a user-defined problem

int main() {

    unsigned short int status;
    std::ofstream file_stats,file_opt_sol;

    std::ios_base::sync_with_stdio(false);
    std::cout.precision(15);

    // (1) Build an object of type Problem
    Problem p(status);

    // (2) Check if an error occurred when calling the Problem constructor
    if (status > 0) {
        std::cout << "\nSomething went wrong when calling the Problem constructor.\n";
        return 1;
    }

    // (3) Build an object of type Asa_bcp
    Asa_bcp alg(status,&p);
    
    // (4) Check if an error occurred when calling the Asa_bcp constructor
    if (status > 0) {
        std::cout << "\nSomething went wrong when calling the Asa_bcp constructor.\n";
        return 1;
    }

    // ------------------------------------------------------------------------------------------------------------
    // *** EXAMPLE OF HOW TO CHANGE ASA-BCP PARAMETERS ***
    // (see the description of Asa_bcp in the file 'usage.txt' to know which parameters can be changed and their
    // default values)
    //
    // Instead of creating the above object 'alg', do the following:
    //
    // - create an object of structure type asa_bcp_options (see its declaration in file 'asa_bcp.h'), e.g.,
    //
    //     asa_bcp_options opts;
    //
    // - assign new values to (some of) the members of the new structure object, e.g.,
    //
    //     opts.verbosity = 0;
    //
    // - pass the address of the structure object as third input argument of the Asa_bcp constructor, e.g.,
    //
    //     Asa_bcp alg(status,&p,&opts);
    // ------------------------------------------------------------------------------------------------------------

    // (5) Call the solver
    alg.solve();

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