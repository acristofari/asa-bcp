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
// June 5th, 2020
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


#include "asa_bcp.h"
#include "problem.h"
#include <numeric>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <limits>


// constructor
//-------------------------------------------------------------------------------------
Asa_bcp::Asa_bcp(unsigned short int& status, Problem* p, const asa_bcp_options* opts){

    unsigned short int prob_status = 0;

    if (opts == NULL) {
        opts = new asa_bcp_options();
    }

    eps_opt = opts->eps_opt;
    min_gd = opts->min_gd;
    min_norm_proj_d = opts->min_norm_proj_d;
    min_stepsize = opts->min_stepsize;
    max_it = (unsigned int) opts->max_it;
    max_n_f = (unsigned int) opts->max_n_f;
    max_n_g = (unsigned int) opts->max_n_g;
    max_n_hd = (unsigned int) opts->max_n_hd;
    min_f = opts->min_f;
    m = (unsigned int) opts->m;
    z = (unsigned int) opts->z;
    hd_exact = opts->hd_exact;
    verbosity = (short unsigned int) opts->verbosity;

    if (eps_opt < 0e0) {
        std::cout << "error when calling asa_bcp: 'eps_opt' must be a number greater than or equal to 0\n";
        exit(-1);
    }
    if (min_gd < 0e0) {
        std::cout << "error when calling asa_bcp: 'min_gd' must be a number greater than or equal to 0\n";
        exit(-1);
    }
    if (min_norm_proj_d < 0e0) {
        std::cout << "error when calling asa_bcp: 'min_norm_proj_d' must be a number greater than or equal to 0\n";
        exit(-1);
    }
    if (min_stepsize < 0e0) {
        std::cout << "error when calling asa_bcp: 'min_stepsize' must be a number greater than or equal to 0\n";
        exit(-1);
    }
    if (opts->max_it < 0) {
        std::cout << "error when calling asa_bcp: 'max_it' must be a number greater than or equal to 0\n";
        exit(-1);
    }
    if (opts->max_n_f < 1) {
        std::cout << "error when calling asa_bcp: 'max_n_f' must be a number greater than or equal to 1\n";
        exit(-1);
    }
    if (opts->max_n_g < 1) {
        std::cout << "error when calling asa_bcp: 'max_n_g' must be a number greater than or equal to 1\n";
        exit(-1);
    }
    if (opts->max_n_hd < 0) {
        std::cout << "error when calling asa_bcp: 'max_n_hd' must be a number greater than or equal to 0\n";
        exit(-1);
    }
    if (opts->m < 1) {
        std::cout << "error when calling asa_bcp: 'm' must be a number greater than or equal to 1\n";
        exit(-1);
    }
    if (opts->z < 1) {
        std::cout << "error when calling asa_bcp: 'z' must be a number greater than or equal to 1\n";
        exit(-1);
    }
    if (opts->verbosity < 0) {
        std::cout << "error when calling asa_bcp: 'verbosity' must be a number between 0 and 2\n";
        exit(-1);
    }

    status = 0;
    prob = p;

    // initialize dimension, bounds and starting point
    n = *prob->dim; // we have checked that n > 0 when calling the constructor of Problem
    x.resize(n);
    l.resize(n);
    u.resize(n);
    prob->bounds(prob_status,l,u);
    if (prob_status > 0) {
        std::cout << "error with the lower and upper bounds\n";
        status = 1;
    } else {
        prob->starting_point(prob_status,x);
        if (prob_status > 0) {
            std::cout << "error with the starting point\n";
            status = 2;
        }
    }

}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
void Asa_bcp::solve() {

    double delta_dir,delta_act,sq_norm_d_act,beta_dir,tmp;

    //           'n' = number of variables
    //      'n_true' = 'n' - number of fixed variables
    //   'n_non_act' = number of estimated non-active variables
    // 'n_new_act_l' = number of variables x(i) that are estimated active at the lower bound
    //                 and such that x(i) > l(i)
    // 'n_new_act_u' = number of variables x(i) that are estimated active at the upper bound
    //                 and such that x(i) < u(i)

    // index vectors:
    // 'ind_non_act' : in the first 'n_non_act' positions there are the indices of the estimated non-active variables,
    //                 and in the last 'n'-'n_true' positions there are the indices of the fixed variables
    // 'ind_new_act' : in the first 'n_new_act_l' positions there are the indices of the variables x(i) that are
    //                 estimated active at the lower bound and such that x(i) > l(i),
    //                 and in the last 'n_new_act_u' positions there are the indices of the variables x(i) that are
    //                 estimated active at the upper bound and such that x(i) < u(i)

    l_min = -1e20; // no i-th lower bound if l(i) <= l_min
    u_max = 1e20; // no i-th upper bound if u(i) >= u_max

    if (verbosity > 0) {
        std::cout.precision(5);
        std::cout.setf(std::ios::scientific,std::ios::floatfield);
        file_output.open("iteration_history.txt",std::ios::trunc);
        file_output.precision(5);
        file_output.setf(std::ios::scientific,std::ios::floatfield);
    }

    // check problem feasibility, project the starting point onto the box
    // and identify fixed variabled
    ind_non_act.resize(n);
    n_true = n;
    delta_act = 0e0;
    for (unsigned int i=n-1; i>0; i--) {
        if (l[i] < u[i]) {
            if (l[i] > l_min) {
                x[i] = std::max(l[i],x[i]);
                if (u[i] < u_max) {
                    x[i] = std::min(x[i],u[i]);
                    delta_act += (u[i]-l[i])*(u[i]-l[i]);
                }
            }
        } else if (l[i] > u[i]) {
            x.assign(n,nan(""));
            f = nan("");
            sup_norm_proj_g = nan("");
            it = n_f = n_g = n_hd = it_cg_tot = 0;
            flag = -1;
            std::cout << "infeasible problem\n";
            if (verbosity > 0) {
                file_output << "infeasible problem\n";
            }
            clear_vectors();
            return;
        } else {
            n_true--;
            ind_non_act[n_true] = i;
            x[i] = l[i];
        }
    }
    if (l[0] < u[0]) {
        if (l[0] > l_min) {
            x[0] = std::max(l[0],x[0]);
            if (u[0] < u_max) {
                x[0] = std::min(x[0],u[0]);
                delta_act += (u[0]-l[0])*(u[0]-l[0]);
            }
        }
    } else if (l[0] > u[0]) {
        x.assign(n,nan(""));
        f = nan("");
        sup_norm_proj_g = nan("");
        it = n_f = n_g = n_hd = it_cg_tot = 0;
        flag = -1;
        std::cout << "infeasible problem\n";
        if (verbosity > 0) {
            file_output << "infeasible problem\n";
        }
        clear_vectors();
        return;
    } else {
        n_true--;
        ind_non_act[n_true] = 0;
        x[0] = l[0];
    }

    // the point obtained by setting the estimated active variables to the bounds
    // is accepted in case of sufficient decrease in the objective function
    // if the distance between the two points is less than or equal to 'delta_act'
    delta_act = std::min(1e30,std::max(1e3,sqrt(delta_act)));

    // the unit stepsize is accepted without evaluating the objective function
    // if the distance between the two points is less than or equal to 'delta_dir'
    delta_dir = delta0_dir = 1e3;
    beta_dir = 9e-1; // reduction factor of 'delta_dir' (must be >=0 and <1)

    // allocate vectors
    g.resize(n);
    hp.resize(n);
    if (hd_exact) {
        g1.resize(0); // g1 is needed only if H(x)*d is approximated
    } else {
        g1.resize(n);
    }
    ind_new_act.resize(n_true);
    lambda.resize(n_true);
    mu.resize(n_true);
    d.assign(n,0e0);
    g_q.assign(n,0e0);
    v.resize(n);

    // initialize counters
    it = k = it_nm = it_cg_tot = n_f = n_g = n_hd = 0;

    if (verbosity > 0) {
        std::cout << "\nnumber of variables = " << n << " (" << n-n_true << " fixed)";
        file_output << "number of variables = " << n << " (" << n-n_true << " fixed)";
    }
    
    // first objective function evaluation
    prob->funct(status,x,f);
    if (status > 0) {
        f_best = nan("");
        sup_norm_proj_g_best = nan("");
        flag = 9;
        err_obj();
        goto end_asa_bcp;
    }
    n_f = 1;
    f_computed = true;
    
    x_best = x;
    f_best = f_w = f;
    w.assign(m,std::numeric_limits<double>::lowest());

    // first gradient evaluation
    prob->grad(status,x,g);
    if (status > 0) {
        sup_norm_proj_g = nan("");
        flag = 10;
        err_obj();
        goto end_asa_bcp;
    }
    
    n_g = 1;
    compute_sup_norm_proj_g();
    
    g_best = g;

    gd = std::numeric_limits<double>::lowest();

    n_hd = 0;

    // compute the sup-norm of the projected gradient and initalize 'eps_act_set'
    sup_norm_proj_g = 0e0;
    eps_act_set = 0e0;
    j = 0;
    for (unsigned int i=0; i<n;) {
        ind_found = false;
        while (!ind_found) {
            if (j == n-n_true) {
                h = n;
                ind_found = true;
            } else {
                if (i == ind_non_act[j+n_true]) {
                    i++;
                    j++;
                } else {
                    h = ind_non_act[j+n_true];
                    ind_found = true;
                }
            }
        }
        while (i < h) {
            if (std::max(l_min,x[i]-g[i])<l[i]) {
                tmp = x[i] - l[i];
            } else if (std::min(u_max,x[i]-g[i])>u[i]) {
                tmp = x[i] - u[i];
            } else {
                tmp = g[i];
            }
            sup_norm_proj_g = std::max(sup_norm_proj_g,std::abs(tmp));
            eps_act_set += tmp*tmp;
            i++;
        }
    }
    // at this point, 'eps_act_set' is the squared norm of the projected gradient
    if (eps_act_set > 0) {
        eps_act_set = std::min(1e-6,std::pow(eps_act_set,-(3e0/2e0)));
    } else {
        eps_act_set = 1e-6;
    }

    sup_norm_proj_g_best = sup_norm_proj_g;

    z_nm = std::min(n_true,z);

    act_phase = false;
    gd_exit = false;
    dir_exit = false;
    stepsize_exit = false;
    min_f_exit = (f<=min_f);
    checkpoint = true;
    is_restarted = false;
    is_first_linesearch = true;

    flag = 0;
    
    //-----------------
    // START MAIN LOOP
    //-----------------

    while (!converged()) {
              
        if (!is_restarted) {

            //---------------------------------------------------------------
            //     MINIMIZATION STEP OVER THE ESTIMATED ACTIVE VARIABLES
            //---------------------------------------------------------------

            if (!act_phase) {
                // active-set estimate
                act_phase = true;
                compute_multipliers();
                estimate_active_set();
            }

            if (verbosity > 1) {
                std::cout << "\n\n--- iteration details ---\n\n"
                          << "number of estimated active variables not at the bounds = "
                          << n_new_act_l+n_new_act_u;
                file_output << "\n\n--- iteration details ---\n\n"
                            << "number of estimated active variables not at the bounds = "
                            << n_new_act_l+n_new_act_u;
            }

            if (act_phase) {

                if (!f_computed) {
                    prob->funct(status,x,f);
                    if (status > 0) {
                        flag = 9;
                        err_obj();
                        goto end_asa_bcp;
                    }
                    n_f++;
                    f_computed = true;
                    if (verbosity > 1) {
                        std::cout << "\nfunction control, f = " << f;
                        file_output << "\nfunction control, f = " << f;
                    }
                }

                f_decreased = false;
                v = x;
                while (act_phase && !f_decreased) {

                    sq_norm_d_act = 0e0;

                    // set x(i) = l(i), for all i in A_l(x)
                    for (unsigned int i=0; i<n_new_act_l; i++) {
                        j = ind_new_act[i];
                        sq_norm_d_act += (x[j]-l[j])*(x[j]-l[j]);
                        x[j] = l[j];
                    }

                    // set x(i) = u(i), for all i in A_u(x)
                    for (unsigned int i=n_true-n_new_act_u; i<n_true; i++) {
                        j = ind_new_act[i];
                        sq_norm_d_act += (x[j]-u[j])*(x[j]-u[j]);
                        x[j] = u[j];
                    }

                    // compute the objective function
                    prob->funct(status,x,tmp);
                    if (status > 0) {
                        flag = 9;
                        err_obj();
                        goto end_asa_bcp;
                    }
                    n_f++;

                    // check if the objective function is sufficiently decreased
                    if (tmp <= f-(sq_norm_d_act/(2e0*eps_act_set))) {
                        f_decreased = true;
                        if (sqrt(sq_norm_d_act) <= delta_act) { // new point accepted
                            delta_act *= beta_dir;
                        } else if (tmp > f_w) { // point not accepted
                            if (verbosity > 1) {
                                std::cout << "\npoint not accepted (f = " << tmp
                                          << ")\n\nrestart from the best point";
                                file_output << "\npoint not accepted (f = " << tmp
                                            << ")\n\nrestart from the best point";
                            }
                            restart();
                        }
                        if (!is_restarted) { // N.B. 'f_computed' remains true
                            f = tmp;
                            // compute the gradient and the sup-norm of the projected gradient
                            prob->grad(status,x,g);
                            if (status > 0) {
                                flag = 10;
                                err_obj();
                                goto end_asa_bcp;
                            }
                            n_g++;
                            compute_sup_norm_proj_g();
                            if (f < f_best) {
                                if (f <= min_f) {
                                    min_f_exit = true;
                                } else {
                                    x_best = x;
                                    f_best = f;
                                    g_best = g;
                                    sup_norm_proj_g_best = sup_norm_proj_g;
                                }
                            }
                            if (verbosity > 1) {
                                std::cout << "\npoint accepted (f = " << f << ")";
                                file_output << "\npoint accepted (f = " << f << ")";
                            }
                        }
                    } else { // set x to the previous value and reduce 'eps_act_set'
                             // to carry out a new active-set estimate
                        for (unsigned int i=0; i<n_new_act_l; i++) {
                            x[ind_new_act[i]] = v[ind_new_act[i]];
                        }
                        for (unsigned int i=n_true-n_new_act_u; i<n_true; i++) {
                            x[ind_new_act[i]] = v[ind_new_act[i]];
                        }
                        if (n_f < max_n_f) {
                            eps_act_set *= 1e-1;
                            estimate_active_set();
                            if (verbosity > 1) {
                                std::cout << "\npoint not accepted (f = " << tmp << ")\nreducing epsilon"
                                          << "\nnumber of estimated active variables not at the bounds = "
                                          << n_new_act_l+n_new_act_u;
                                file_output << "\npoint not accepted (f = " << tmp << ")\nreducing epsilon"
                                            << "\nnumber of estimated active variables not at the bounds = "
                                            << n_new_act_l+n_new_act_u;
                            }
                        } else {
                            act_phase = false;
                        }
                    }

                }

            }

        } else {

            is_restarted = false;
            if (verbosity > 1) {
                std::cout << "\n\n--- iteration details ---";
                file_output << "\n\n--- iteration details ---";
            }
        }
        
        //---------------------------------------------------------------
        //   MINIMIZATION STEP OVER THE ESTIMATED NON-ACTIVE VARIABLES
        //---------------------------------------------------------------

        if (!is_restarted && sup_norm_proj_g>eps_opt && n_f<max_n_f && n_g<max_n_g && !min_f_exit) {

            if (act_phase) {
                // active-set estimate
                act_phase = false;
                compute_multipliers();
                estimate_active_set();
            }

            if (verbosity > 1) {
                file_output << "\n\nnumber of estimated non-active variables = " << n_non_act;
                std::cout << "\n\nnumber of estimated non-active variables = " << n_non_act;
            }

            if (!act_phase) {

                // compute the norm of the gradient with respect to the esimated non-active variables
                norm_g_non_act = 0e0;
                for (unsigned int i=0; i<n_non_act; i++) {
                    norm_g_non_act += g[ind_non_act[i]]*g[ind_non_act[i]];
                }
                norm_g_non_act = sqrt(norm_g_non_act);

                if (norm_g_non_act > min_norm_proj_d) {

                    if (checkpoint) {
                        // N.B. 'f_computed' = true at this point
                        if (f < f_best) {
                            if (f <= min_f) {
                                min_f_exit = true;
                            } else {
                                x_best = x;
                                f_best = f;
                                g_best = g;
                                sup_norm_proj_g_best = sup_norm_proj_g;
                            }
                        }
                        if (!min_f_exit) {
                            update_w();
                            checkpoint = false;
                        }
                    }

                    // function control
                    if (it_nm >= z_nm) {
                        z_nm = std::min(z_nm+n_true,z);
                        if (!f_computed) {
                            prob->funct(status,x,f);
                            if (status > 0) {
                                flag = 9;
                                err_obj();
                                goto end_asa_bcp;
                            }
                            n_f++;
                        }
                        if (f >= f_w) {
                            if (verbosity > 1) {
                                std::cout << "\nfunction control not satisfied (f = " << f
                                          << ")\n\nrestart from the best point";
                                file_output << "\nfunction control not satisfied (f = " << f
                                            << ")\n\nrestart from the best point";
                            }
                            restart();
                        } else {
                            if (verbosity > 1) {
                                std::cout << "\nfunction control satisfied (f = " << f << ")";
                                file_output << "\nfunction control satisfied (f = " << f << ")";
                            }
                            if (f < f_best) {
                                if (f <= min_f) {
                                    min_f_exit = true;
                                } else {
                                    x_best = x;
                                    f_best = f;
                                    g_best = g;
                                    sup_norm_proj_g_best = sup_norm_proj_g;
                                }
                            }
                            if (!min_f_exit) {
                                update_w();
                                f_computed = true;
                            }
                        }
                    }

                    if (!is_restarted && !min_f_exit) {

                        max_inner_it = hd_exact ? std::min(2*n_non_act,max_n_hd-it_cg_tot) : std::min(2*n_non_act,max_n_g-it_cg_tot);

                        if (max_inner_it > 1) {

                            // compute the search direction
                            trunc_newton_dir();
                            it_cg_tot += it_cg;

                            if (flag==10 || flag==11) {
                                err_obj();
                                goto end_asa_bcp;
                            }

                            if (verbosity > 1) {
                                std::cout << "\nnumber of inner cg iterations = " << it_cg;
                                file_output << "\nnumber of inner cg iterations = " << it_cg;
                                if (warn_sing) {
                                    std::cout << "\nHessian matrix probably singular";
                                    file_output << "\nHessian matrix probably singular";
                                } else if (warn_noposdef) {
                                    std::cout << "\nHessian matrix probably not positive definite";
                                    file_output << "\nHessian matrix probably not positive definite";
                                } else if (warn_small) {
                                    std::cout << "\nnorm of the inner conjugate direction too small";
                                    file_output << "\nnorm of the inner conjugate direction too small";
                                } else if (warn_conjfail) {
                                    std::cout << "\nconjugacy failure when computing the search direction";
                                    file_output << "\nconjugacy failure when computing the search direction";
                                } else if (warn_maxcgit) {
                                    std::cout << "\nconjugate gradient method not converged";
                                    file_output << "\nconjugate gradient method not converged";
                                }
                                if (warn_grad) {
					                std::cout << "\nanti-gradient used as search direction";
                                    file_output << "\nanti-gradient used as search direction";
                                }
                                if (is_d_neg_curv) {
                                    std::cout << "\nnegative curvature direction";
                                    file_output << "\nnegative curvature direction";
                                }
                                std::cout << "\ndirectional derivative = " << gd;
                                file_output << "\ndirectional derivative = " << gd;
                            }

                            // check if the directional derivative is sufficiently negative
                            if (gd < -min_gd) {

                                // check if line search must be performed
                                if (warn_noposdef || is_first_linesearch) {
                                    ls = true;
                                } else {
                                    ls = false;
                                }

                                // try accepting the unit stepsize or prepare for the line search
                                if (!ls) {

                                    // set v = x and x = p[x + d]
                                    norm_proj_d = 0e0;
                                    for (unsigned int i=0; i<n_non_act; i++) {
                                        j = ind_non_act[i];
                                        v[j] = x[j];
                                        x[j] += d[j];
                                        if (std::max(l_min,x[j])<l[j]) {
                                            x[j] = l[j];
                                        } else if (std::min(u_max,x[j])>u[j]) {
                                            x[j] = u[j];
                                        }
                                        norm_proj_d += (x[j]-v[j])*(x[j]-v[j]);
                                    }
                                    norm_proj_d = sqrt(norm_proj_d);

                                    if (verbosity > 1) {
                                        std::cout << "\nnorm of the projected direction = " << norm_proj_d;
                                        file_output << "\nnorm of the projected direction = " << norm_proj_d;
                                    }

                                    //  check if the norm of the projected direction is sufficiently large
                                    if (norm_proj_d > min_norm_proj_d) {
                                        if (norm_proj_d <= delta_dir) { // unit stepsize accepted
                                            delta_dir *= beta_dir;
                                            // compute the gradient and the sup-norm of the projected gradient
                                            prob->grad(status,x,g);
                                            if (status > 0) {
                                                flag = 10;
                                                err_obj();
                                                goto end_asa_bcp;
                                            } else {
                                                n_g++;
                                                compute_sup_norm_proj_g();
                                                f_computed = false;
                                                it_nm++;
                                                if (verbosity > 1) {
                                                    std::cout << "\nunit stepsize accepted without computing f";
                                                    file_output << "\nunit stepsize accepted without computing f";
                                                }
                                            }
                                        } else if (!f_computed) { // check the objective function
                                            // restore x
                                            for (unsigned int i=0; i<n_non_act; i++) {
                                                j = ind_non_act[i];
                                                tmp = v[j];
                                                v[j] = x[j];
                                                x[j] = tmp;
                                            }
                                            prob->funct(status,x,f); // evaluate f(x)
                                            if (status > 0) {
                                                flag = 9;
                                                err_obj();
                                                goto end_asa_bcp;
                                            }
                                            n_f++;
                                            if (f < f_w) { // objective function decreased -> line search
                                                if (verbosity > 1) {
                                                    std::cout << "\nfunction control satisfied (f = " << f << ")";
                                                    file_output << "\nfunction control satisfied (f = " << f << ")";
                                                }
                                                if (f < f_best) {
                                                    if (f <= min_f) {
                                                        min_f_exit = true;
                                                    } else {
                                                        x_best = x;
                                                        f_best = f;
                                                        g_best = g;
                                                        sup_norm_proj_g_best = sup_norm_proj_g;
                                                    }
                                                }
                                                if (!min_f_exit) {
                                                    update_w();
                                                    f_computed = true;
                                                    ls = true;
                                                    // set again v = x and x = p[x + d]
                                                    for (unsigned int i=0; i<n_non_act; i++) {
                                                        j = ind_non_act[i];
                                                        tmp = v[j];
                                                        v[j] = x[j];
                                                        x[j] = tmp;
                                                    }
                                                }
                                            } else { // objective function not decreased -> restart
                                                if (verbosity > 1) {
                                                    std::cout << "\nfunction control not satisfied (f = " << f
                                                              << ")\n\nrestart from the best point";
                                                    file_output << "\nfunction control not satisfied (f = " << f
                                                                << ")\n\nrestart from the best point";
                                                }
                                                restart();
                                            }
                                        } else {
                                            ls = true;
                                        }
                                    } else {
                                        // restore x
                                        for (unsigned int i=0; i<n_non_act; i++) {
                                            x[ind_non_act[i]] = v[ind_non_act[i]];
                                        }
                                        dir_exit = true;
                                    }

                                } else if (!f_computed) {

                                    prob->funct(status,x,f); // evaluate f(x)
                                    if (status > 0) {
                                        flag = 9;
                                        err_obj();
                                        goto end_asa_bcp;
                                    }
                                    n_f++;
                                    if (f < f_w) { // objective function decreased -> line search
                                        // set v = x and x = p[x + d]
                                        f_computed = true;
                                        norm_proj_d = 0e0;
                                        for (unsigned int i=0; i<n_non_act; i++) {
                                            j = ind_non_act[i];
                                            v[j] = x[j];
                                            x[j] += d[j];
                                            if (std::max(l_min,x[j])<l[j]) {
                                                x[j] = l[j];
                                            } else if (std::min(u_max,x[j])>u[j]) {
                                                x[j] = u[j];
                                            }
                                            norm_proj_d += (x[j]-v[j])*(x[j]-v[j]);
                                        }
                                        norm_proj_d = sqrt(norm_proj_d);
                                        if (verbosity > 1) {
                                            std::cout << "\nnorm of the projected direction = " << norm_proj_d;
                                            file_output << "\nnorm of the projected direction = " << norm_proj_d;
                                        }
                                        if (norm_proj_d > min_norm_proj_d) {
                                            if (verbosity > 1) {
                                                std::cout << "\nfunction control satisfied (f = " << f << ")";
                                                file_output << "\nfunction control satisfied (f = " << f << ")";
                                            }
                                            if (f < f_best) {
                                                if (f <= min_f) {
                                                    for (unsigned int i=0; i<n_non_act; i++) {
                                                        x[ind_non_act[i]] = v[ind_non_act[i]];
                                                    }
                                                    min_f_exit = true;
                                                    ls = false; // to skip the line search and exit the while loop
                                                } else {
                                                    x_best = x;
                                                    f_best = f;
                                                    g_best = g;
                                                    sup_norm_proj_g_best = sup_norm_proj_g;
                                                }
                                            }
                                            if (!min_f_exit) {
                                                update_w();
                                            }
                                        } else {
                                            // restore x
                                            for (unsigned int i=0; i<n_non_act; i++) {
                                                x[ind_non_act[i]] = v[ind_non_act[i]];
                                            }
                                            dir_exit = true;
                                            ls = false; // to skip the line search and exit the while loop
                                        }
                                    } else { // objective function not decreased -> restart
                                             // (first check if the norm of the projected direction is sufficiently large)
                                        norm_proj_d = 0e0;
                                        for (unsigned int i=0; i<n_non_act; i++) {
                                            j = ind_non_act[i];
                                            tmp = x[j] + d[j];
                                            if (std::max(l_min,tmp)<l[j]) {
                                                norm_proj_d += (x[j]-l[j])*(x[j]-l[j]);
                                            } else if (std::min(u_max,tmp)>u[j]) {
                                                norm_proj_d += (u[j]-x[j])*(u[j]-x[j]);
                                            }
                                            norm_proj_d += (tmp-x[j])*(tmp-x[j]);
                                        }
                                        norm_proj_d = sqrt(norm_proj_d);
                                        if (verbosity > 1) {
                                            std::cout << "\nnorm of the projected direction = " << norm_proj_d;
                                            file_output << "\nnorm of the projected direction = " << norm_proj_d;
                                        }
                                        if (norm_proj_d > min_norm_proj_d) {
                                            if (verbosity > 1) {
                                                std::cout << "\npoint not accepted (f = " << f
                                                          << ")\n\nrestart from the best point";
                                                file_output << "\npoint not accepted (f = " << f
                                                            << ")\n\nrestart from the best point";
                                            }
                                            restart();
                                        } else {
                                            dir_exit = true;
                                            f_computed = true;
                                            ls = false; // to skip the line search and exit the while loop
                                        }                            
                                    }

                                } else {

                                    // set v = x and x = p[x + d]
                                    norm_proj_d = 0e0;
                                    for (unsigned int i=0; i<n_non_act; i++) {
                                        j = ind_non_act[i];
                                        v[j] = x[j];
                                        x[j] += d[j];
                                        if (std::max(l_min,x[j])<l[j]) {
                                            x[j] = l[j];
                                        } else if (std::min(u_max,x[j])>u[j]) {
                                            x[j] = u[j];
                                        }
                                        norm_proj_d += (x[j]-v[j])*(x[j]-v[j]);
                                    }
                                    norm_proj_d = sqrt(norm_proj_d);
                                    if (verbosity > 1) {
                                        std::cout << "\nnorm of the projected direction = " << norm_proj_d;
                                        file_output << "\nnorm of the projected direction = " << norm_proj_d;
                                    }
                                    if (norm_proj_d <= min_norm_proj_d) {
                                        // restore x
                                        for (unsigned int i=0; i<n_non_act; i++) {
                                            x[ind_non_act[i]] = v[ind_non_act[i]];
                                        }
                                        dir_exit = true;
                                        ls = false; // to skip the line search and exit the while loop
                                    }

                                }

                                // line search
                                if (ls && !is_restarted) {
                                    if (verbosity >= 2) {
                                        std::cout << "\nline search";
                                        file_output << "\nline search";
                                    }
                                    if (is_first_linesearch) {
                                        f_newton_first = f;
                                        delta_f0 = 0e0;
                                    }
                                    linesearch(); // N.B. 'f_computed' = true before and after line search
                                    if (flag == 9) {
                                        err_obj();
                                        goto end_asa_bcp;
                                    }
                                    if (flag != 5) {
                                        if (!stepsize_exit) { // check if the line search succeeded
                                            // compute the gradient and the sup-norm of the projected gradient
                                            prob->grad(status,x,g);
                                            if (status > 0) {
                                                flag = 10;
                                                err_obj();
                                                goto end_asa_bcp;
                                            }
                                            n_g++;
                                            compute_sup_norm_proj_g();
                                            checkpoint = true;
                                            if (is_first_linesearch) {
                                                delta_f0 = f_newton_first - f;
                                                delta_dir = delta0_dir*stepsize*norm_proj_d;
                                                is_first_linesearch = false;
                                            }
                                            if (verbosity > 1) {
                                                std::cout << "\nstepsize = " << stepsize;
                                                file_output << "\nstepsize = " << stepsize;
                                            }
                                        } else {
                                            // restore x and f(x)
                                            for (unsigned int i=0; i<n_non_act; i++) {
                                                x[ind_non_act[i]] = v[ind_non_act[i]];
                                            }
                                            f = fv;
                                        }
                                    } else {
                                        // restore x and f(x)
                                        for (unsigned int i=0; i<n_non_act; i++) {
                                            x[ind_non_act[i]] = v[ind_non_act[i]];
                                        }
                                        f = fv;
                                    }
                                }
                            } else {
                                gd_exit = true;
                            }
                        }
                    }
                } else {
                    dir_exit = true;
                }

            } else {

                gd = std::numeric_limits<double>::lowest();
                norm_proj_d = std::numeric_limits<double>::max();
                stepsize = std::numeric_limits<double>::max();
                it_nm++;

            }

        } else {

            gd = std::numeric_limits<double>::lowest();
            norm_proj_d = std::numeric_limits<double>::max();
            stepsize = std::numeric_limits<double>::max();

        }

        it++;
        k++;

    }

    end_asa_bcp:

    if (verbosity > 0) {
        std::cout << "\n\nWARNING: using 'verbosity = 0' may be faster\n\n";
        file_output << "\n\nWARNING: using 'verbosity = 0' may be faster\n";
        file_output.close();
    }

    clear_vectors();

}
//-------------------------------------------------------------------------------------




// other functions
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
bool Asa_bcp::converged() {

    if ( (sup_norm_proj_g > eps_opt) && (!gd_exit) && (!dir_exit) && (!stepsize_exit)
         && (it<max_it) && (n_f<max_n_f) && (n_g<max_n_g) && (n_hd<max_n_hd) && (!min_f_exit) ) {
        if (verbosity > 0) {
            main_prints();
        }
        return false;
    } else {
        if (!f_computed) {
            prob->funct(status,x,f);
            if (status > 0) {
                x = x_best;
                f = f_best;
                sup_norm_proj_g = sup_norm_proj_g_best;
                flag = 9;
                std::cout << "\n\nerror when computing the objective function, algorithm stopped";
                if (verbosity > 0) {
                    file_output << "\n\nerror when computing the objective function, algorithm stopped";
                }
                return true;
            }
            n_f++;
            f_computed = true;
            if (f <= min_f) {
                flag = 8;
                if (verbosity > 0) {
                    std::cout << "\n\nobjective value below the minimum threshold, algorithm stopped";
                    file_output << "\n\nobjective value below the minimum threshold, algorithm stopped";
                }
                return true;
            }
        }
        if (sup_norm_proj_g <= eps_opt) {
            if (verbosity > 0) {
                main_prints();
                std::cout << "\n\n================================================================================="
                          << "\noptimality condition satisfied: sup-norm of the projected gradient <= " << eps_opt
                          << "\n=================================================================================";
                file_output << "\n\n================================================================================="
                            << "\noptimality condition satisfied: sup-norm of the projected gradient <= " << eps_opt
                            << "\n=================================================================================";
            }
            if (f <= f_best) {
                flag = 0;
                return true;
            } else {
                return check_stop(true);
            }
        } else if (gd_exit) {
            if (verbosity > 0) {
                std::cout << "\n\ndirectional derivative not sufficiently negative";
                file_output << "\n\ndirectional derivative not sufficiently negative";
            }
            // active-set estimate
            act_phase = true;
            compute_multipliers();
            estimate_active_set();                
            if (act_phase) {
                 gd_exit = false;
                 return check_stop(false);
            } else {
                if (f <= f_best) {
                    it--;
                    flag = 1;
                    if (verbosity > 1) {
                        std::cout << ", algorithm stopped";
                        file_output << ", algorithm stopped";
                    }
                    return true;
                } else {
                    return check_stop(true);
                }
            }
        } else if (dir_exit) {
            if (verbosity > 0) {
                std::cout << "\n\nnorm of the projected direction too small";
                file_output << "\n\nnorm of the projected direction too small";
            }
            // active-set estimate
            act_phase = true;
            compute_multipliers();
            estimate_active_set();                
            if (act_phase) {
                 dir_exit = false;
                 return check_stop(false);
            } else {
                if (f <= f_best) {
                    it--;
                    flag = 2;
                    if (verbosity > 1) {
                        std::cout << ", algorithm stopped";
                        file_output << ", algorithm stopped";
                    }
                    return true;
                } else {
                    return check_stop(true);
                }
            }
        } else if (stepsize_exit) {
            if (verbosity > 0) {
                std::cout << "\n\nstepsize too small";
                file_output << "\n\nstepsize too small";
            }
            // active-set estimate
            act_phase = true;
            compute_multipliers();
            estimate_active_set();                
            if (act_phase) {
                 stepsize_exit = false;
                 return check_stop(false);
            } else {
                if (f <= f_best) {
                    it--;
                    flag = 3;
                    if (verbosity > 1) {
                        std::cout << ", algorithm stopped";
                        file_output << ", algorithm stopped";
                    }
                    return true;
                } else {
                    return check_stop(true);
                }
            }
        } else if (min_f_exit) {
            flag = 8;
            if (verbosity > 0) {
                std::cout << "\n\nobjective value below the minimum threshold, algorithm stopped";
                file_output << "\n\nobjective value below the minimum threshold, algorithm stopped";
            }
            return true;
        } else {
            if (f > f_best) {
                x = x_best;
                f = f_best;
                sup_norm_proj_g = sup_norm_proj_g_best;
            }
            if (it >= max_it) {
                flag = 4;
                if (verbosity > 0) {
                    main_prints();
                    std::cout << "\n\ntoo many iterations, algorithm stopped";
                    file_output << "\n\ntoo many iterations, algorithm stopped";
                }
            } else if (n_f >= max_n_f) {
                flag = 5;
                if (verbosity > 0) {
                    main_prints();
                    std::cout << "\n\ntoo many function evaluations, algorithm stopped";
                    file_output << "\n\ntoo many function evaluations, algorithm stopped";
                }
            } else if (n_g >= max_n_g) {
                flag = 6;
                if (verbosity > 0) {
                    main_prints();
                    std::cout << "\n\ntoo many gradient evaluations, algorithm stopped";
                    file_output << "\n\ntoo many gradient evaluations, algorithm stopped";
                }
            } else { // n_hd >= max_n_hd
                flag = 7;
                if (verbosity > 0) {
                    main_prints();
                    std::cout << "\n\ntoo many Hessian-vector products, algorithm stopped";
                    file_output << "\n\ntoo many Hessian-vector products, algorithm stopped";
                }
            }
            return true;
        }
    }
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
bool Asa_bcp::check_stop(bool rst) { // it returns true if the algorithm must stop, false otherwise
    if ( (it<max_it) && (n_f<max_n_f) && (n_g<max_n_g) && (n_hd<max_n_hd) ) {
        if (rst) {
            restart();
            if (verbosity > 0) {
                if (verbosity > 1) {
                    std::cout << "\n\nrestart from the best point";
                    file_output << "\n\nrestart from the best point";
                }
                main_prints();
            }
        }
        return false;
    } else {        
        if (rst || f>f_best) {
            x = x_best;
            f = f_best;
            sup_norm_proj_g = sup_norm_proj_g_best;
        }
        if (it >= max_it) {
            flag = 4;
            if (verbosity > 0) {
                std::cout << "\ntoo many iterations, algorithm stopped";
                file_output << "\ntoo many iterations, algorithm stopped";
            }
        } else if (n_f >= max_n_f) {
            flag = 5;
            if (verbosity > 0) {
                std::cout << "\ntoo many function evaluations, algorithm stopped";
                file_output << "\ntoo many function evaluations, algorithm stopped";
            }
        } else if (n_g >= max_n_g) {
            flag = 6;
            if (verbosity > 0) {
                std::cout << "\ntoo many gradient evaluations, algorithm stopped";
                file_output << "\ntoo many gradient evaluations, algorithm stopped";
            }
        } else { // n_hd >= max_n_hd
            flag = 7;
            if (verbosity > 0) {
                std::cout << "\ntoo many Hessian-vector products, algorithm stopped";
                file_output << "\ntoo many Hessian-vector products, algorithm stopped";
            }
        }
        return true;
    }
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
void Asa_bcp::main_prints() {
    std::cout << "\n\n--------------------------------------------------------\n\n"
              << "iteration " << it
              << "\n\nbest f = " << std::min(f_best,f)
              << "\nsup-norm of the projected gradient at the current point = " << sup_norm_proj_g
              << "\nnumber of function evaluations = " << n_f
              << "\nnumber of gradient evaluations = " << n_g
              << "\nnumber of Hessian-vector products = " << n_hd
              << "\nnumber of inner cg iterations = " << it_cg_tot;
    file_output << "\n\n--------------------------------------------------------\n\n"
                << "iteration " << it
                << "\n\nbest f = " << std::min(f_best,f)
                << "\nsup-norm of the projected gradient at the current point = " << sup_norm_proj_g
                << "\nnumber of function evaluations = " << n_f
                << "\nnumber of gradient evaluations = " << n_g
                << "\nnumber of Hessian-vector products = " << n_hd
                << "\nnumber of inner cg iterations = " << it_cg_tot;
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
void Asa_bcp::compute_sup_norm_proj_g() { // compute the sup-norm of (x - p[x-g(x)])
    double proj_g;

    sup_norm_proj_g = 0e0;
    j = 0;
    for (unsigned int i=0; i<n;) {
        ind_found = false;
        while (!ind_found) {
            if (j == n-n_true) {
                h = n;
                ind_found = true;
            } else {
                if (i == ind_non_act[j+n_true]) {
                    i++;
                    j++;
                } else {
                    h = ind_non_act[j+n_true];
                    ind_found = true;
                }
            }
        }
        while (i < h) {
            if (std::max(l_min,x[i]-g[i])<l[i]) {
                proj_g = x[i] - l[i];
            } else if (std::min(u_max,x[i]-g[i])>u[i]) {
                proj_g = x[i] - u[i];
            } else {
                proj_g = g[i];
            }
            sup_norm_proj_g = std::max(sup_norm_proj_g, std::abs(proj_g));
            i++;
        }
    }
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
void Asa_bcp::compute_multipliers() {
    double mult1,mult2;
    unsigned int z;

    j = z = 0;
    for (unsigned int i=0; i<n;) {
        ind_found = false;
        while (!ind_found) {
            if (j == n-n_true) {
                h = n;
                ind_found = true;
            } else {
                if (i == ind_non_act[j+n_true]) {
                    i++;
                    j++;
                } else {
                    h = ind_non_act[j+n_true];
                    ind_found = true;
                }
            }
        }
        while (i < h) {
            if (l[i] > l_min) {
                mult1 = (u[i]-x[i])*(u[i]-x[i]);
                mult2 = (x[i]-l[i])*(x[i]-l[i]);
                lambda[z] = g[i]*mult1/(mult1+mult2);
                if (u[i] < u_max) {
                    mu[z] = -g[i]*mult2/(mult1+mult2);
                }
            } else if (u[i] < u_max) {
                mu[z] = -g[i]*(x[i]-l[i])*(x[i]-l[i])/((u[i]-x[i])*(u[i]-x[i])+(x[i]-l[i])*(x[i]-l[i]));
            }
            i++;
            z++;
        }
    }
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
void Asa_bcp::estimate_active_set() {
    unsigned int count,z,t;
    bool est;

    est = true;
    count = 0; // to avoid cycling
    while (est && count<2) {
        if (act_phase) {
            n_new_act_l = n_new_act_u = 0;
            j = z = 0;
            t = n_true;
            for (unsigned int i=0; i<n;) {
                ind_found = false;
                while (!ind_found) {
                    if (j == n-n_true) {
                        h = n;
                        ind_found = true;
                    } else {
                        if (i == ind_non_act[j+n_true]) {
                            i++;
                            j++;
                        } else {
                            h = ind_non_act[j+n_true];
                            ind_found = true;
                        }
                    }
                }
                while (i < h) {
                    if ( (x[i]<=l[i]+eps_act_set*lambda[z]) && (g[i]>0e0) && (x[i]>l[i]) && (l[i]>l_min) ) {
                        ind_new_act[n_new_act_l] = i;
                        n_new_act_l++;
                    } else if ( (x[i]>=u[i]-eps_act_set*mu[z]) && (g[i]<0e0) && ((x[i]<u[i])) && (u[i]<u_max) ) {
                        t--;
                        ind_new_act[t] = i;
                    }
                    i++;
                    z++;
                }
            }
            n_new_act_u = n_true - t;
            act_phase = (n_new_act_l+n_new_act_u>0);
            est = !act_phase;
            count++;
        } else {
            n_non_act = 0;
            j = z = 0;
            t = n_true;
            for (unsigned int i=0; i<n;) {
                ind_found = false;
                while (!ind_found) {
                    if (j == n-n_true) {
                        h = n;
                        ind_found = true;
                    } else {
                        if (i == ind_non_act[j+n_true]) {
                            i++;
                            j++;
                        } else {
                            h = ind_non_act[j+n_true];
                            ind_found = true;
                        }
                    }
                }
                while (i < h) {
                    if ( !((x[i]<=l[i]+eps_act_set*lambda[z]) && (g[i]>0e0) && (l[i]>l_min)) &&
                         !((x[i]>=u[i]-eps_act_set*mu[z]) && (g[i]<0e0) && (u[i]<u_max)) ) {
                        ind_non_act[n_non_act] = i;
                        n_non_act++;
                    }
                    i++;
                    z++;
                }
            }
            act_phase = (n_non_act==0);
            est = act_phase;
            count++;
        }
    }
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
void Asa_bcp::update_w() { // update the vector of reference values
    f_w = f;
    for (unsigned int i=m-1; i>0; i--) {
        w[i] = w[i-1];
        f_w = std::max(f_w,w[i]);
    }
    w[0] = f;
    it_nm = 0;
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
void Asa_bcp::restart() {
    x = x_best;
    f = f_w = f_best;
    g = g_best;
    sup_norm_proj_g = sup_norm_proj_g_best;
    m = (m+4)/5;
    z_nm = std::min(n_true,z);
    w.assign(m,f_w);
    k = it_nm = 0;
    delta0_dir *= 1e-1;
    f_computed = true;
    checkpoint = false;
    is_restarted = true;
    is_first_linesearch = true;
    act_phase = true;
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
void Asa_bcp::trunc_newton_dir() {
    const double eps_cg_dir = std::min(1e0,pow(min_norm_proj_d,3e0/2e0)); // 'eps_cg_dir' should be <= 'min_norm_proj_d',
                                                                          // so that the norm of the first conjugate
                                                                          // direction is sufficiently large
    const double eps_cg_curv = 1e-9;

    double alpha,beta,curv,norm_p,gp,gqp,tmp,eps_cg_tr,norm_d;;
    std::vector<double> p(n,0e0),hd(n,0e0);
    bool trunc_exit;
    
    // In this function, 'p' are the (normalized) conjugate directions in the subspace
    // of the estimated non-active variables and 'd' is the resulting search direction

    warn_grad = warn_small = warn_noposdef = warn_conjfail = warn_maxcgit = is_d_neg_curv = false;
    it_cg = 1;
    
    // compute the first 'p' and normalize it
    for (unsigned int i=0; i<n_non_act; i++) {
        p[ind_non_act[i]] = -g[ind_non_act[i]]/norm_g_non_act;
    }

    // compute H(x)*p
    if (hd_exact) {
        prob->hd_prod(status,false,x,p,hp);
        if (status > 0) {
            flag = 11;
            return;
        }
        n_hd++;
    } else {
        approximate_hd(p);
        if (status > 0) {
            flag = 10;
            return;
        }
    }

    // compute curvature
    curv = 0e0;
    for (unsigned int i=0; i<n_non_act; i++) {
        curv += p[ind_non_act[i]]*hp[ind_non_act[i]];
    }

    // check curvature
    if (curv > eps_cg_curv) {

        // update 'd'
        for (unsigned int i=0; i<n_non_act; i++) {
            d[ind_non_act[i]] = -g[ind_non_act[i]]/curv;
        }

        // first truncated-Newton termination test
        if (curv <= 1e0) {
            eps_cg_tr = std::min(1e0,1e-1+exp(-1e-3*(double)k));
        } else {
            eps_cg_tr = std::min(1e0,1e-1+exp(-1e-3*(double)k))*std::min(1e0,norm_g_non_act/curv);
        }
        if (eps_cg_tr < 1e0) {

            gd = -norm_g_non_act*norm_g_non_act/curv;
            q = 5e-1*gd;

            tmp = norm_g_non_act/curv;
            for (unsigned int i=0; i<n_non_act; i++) {
                j = ind_non_act[i];
                hd[j] = hp[j]*tmp;
                g_q[j] = g[j] + hd[j];
            }

            // check if the maximum number of inner iterations has been reached
            if (max_inner_it > 1) {

                // update 'p'
                beta = 0e0;
                for (unsigned int i=0; i<n_non_act; i++) {
                    beta += g_q[ind_non_act[i]]*hp[ind_non_act[i]];
                }
                beta /= curv;
                gp = gqp = norm_p = 0e0;
                for (unsigned int i=0; i<n_non_act; i++) {
                    j = ind_non_act[i];
                    p[j] = -g_q[j] + beta*p[j];
                    gp += g[j]*p[j];
                    gqp += g_q[j]*p[j];
                    norm_p += p[j]*p[j];
                }
                
                // check conjugacy
                if((abs(gqp-gp)>abs(gqp)+1e-6) || ((gp*gqp<0e0) && (abs(gqp-gp)>1e-6*(abs(gqp)+1e-9)))) {
                    warn_conjfail = true;
                }

                trunc_exit = false;

                // repeat for iteration 2,3,...
                while (!trunc_exit && !warn_conjfail && !warn_noposdef && !warn_small && !warn_maxcgit) {

                    it_cg++;

                    norm_p = 0e0;
                    for (unsigned int i=0; i<n_non_act; i++) {
                        norm_p += p[ind_non_act[i]]*p[ind_non_act[i]];
                    }
                    norm_p = sqrt(norm_p);

                    // check if the norm of 'p' is sufficiently large
                    if (norm_p > eps_cg_dir) {

                        // normalize 'p'
                        for (unsigned int i=0; i<n_non_act; i++) {
                            p[ind_non_act[i]] = p[ind_non_act[i]]/norm_p;
                        }

                        // compute H(x)*p
                        if (hd_exact) {
                            prob->hd_prod(status,true,x,p,hp);
                            if (status > 0) {
                                flag = 11;
                                return;
                            }
                            n_hd++;
                        } else {
                            approximate_hd(p);
                            if (status > 0) {
                                flag = 10;
                                return;
                            }
                        }

                        // compute curvature
                        curv = 0e0;
                        for (unsigned int i=0; i<n_non_act; i++) {
                            curv += p[ind_non_act[i]]*hp[ind_non_act[i]];
                        }

                        // check curvature
                        if (curv > eps_cg_curv) {

                            // update 'd'
                            alpha = 0e0;
                            for (unsigned int i=0; i<n_non_act; i++) {
                                alpha -= p[ind_non_act[i]]*g_q[ind_non_act[i]];
                            }
                            alpha /= curv;
                            for (unsigned int i=0; i<n_non_act; i++) {
                                j = ind_non_act[i];
                                d[j] += alpha*p[j];
                                hd[j] += alpha*hp[j];
                                g_q[j] += alpha*hp[j];
                            }

                            // truncated-Newton termination test
                            norm_d = 0e0;
                            for (unsigned int i=0; i<n_non_act; i++) {
                                norm_d += d[ind_non_act[i]]*d[ind_non_act[i]];
                            }
                            norm_d = sqrt(norm_d);
                            if (norm_d >= norm_g_non_act) {
                                eps_cg_tr = std::min(1e0,1e-1+exp(-1e-3*(double)k));
                            } else {
                                eps_cg_tr = std::min(1e0,1e-1+exp(-1e-3*(double)k))*std::min(1e0,norm_d);
                            }
                            q_old = q;
                            gd_old = gd;
                            q = gd = 0e0;
                            for (unsigned int i=0; i<n_non_act; i++) {
                                j = ind_non_act[i];
                                gd += g[j]*d[j];
                                q += g_q[j]*d[j];
                            }
                            q = (q+gd)/2e0;
                            if (!(abs((q-q_old-3e0*(gd-gd_old)/2e0) / (q-3e0*gd/2e0))*(double)it_cg <= eps_cg_tr)) {

                                // check if the maximum number of inner iterations has been reached
                                if (it_cg < max_inner_it) {

                                    // update 'p'
                                    beta = 0e0;
                                    for (unsigned int i=0; i<n_non_act; i++) {
                                        beta += g_q[ind_non_act[i]]*hp[ind_non_act[i]];
                                    }
                                    beta /= curv;
                                    gp = gqp = 0e0;
                                    for (unsigned int i=0; i<n_non_act; i++) {
                                        j = ind_non_act[i];
                                        p[j] = -g_q[j] + beta*p[j];
                                        gp += g[j]*p[j];
                                        gqp += g_q[j]*p[j];
                                    }

                                    // check conjugacy
                                    if((abs(gqp-gp)>abs(gqp)+1e-6) || ((gp*gqp<0e0) && (abs(gqp-gp)>1e-6*(abs(gqp)+1e-9)))) {
                                        warn_conjfail = true;
                                    }

                                } else {
                                    warn_maxcgit = true;
                                }
                            } else {
                                trunc_exit = true;
                            }

                        } else {

                            // Hessian matrix not (sufficiently) positive definite
                            warn_noposdef = true;
                            if (curv < -eps_cg_curv) {
                                alpha = 0e0;
                                for (unsigned int i=0; i<n_non_act; i++) {
                                    alpha -= p[ind_non_act[i]]*g_q[ind_non_act[i]];
                                }
                                alpha /= curv;
                                double sq_norm_d = 0e0;
                                gd = 0e0;
                                curv = 0e0;
                                for (unsigned int i=0; i<n_non_act; i++) {
                                    j = ind_non_act[i];
                                    d[j] -= alpha*p[j];
                                    hd[j] -= alpha*hp[j];
                                    sq_norm_d += d[j]*d[j];
                                    gd += g[j]*d[j];
                                    curv += d[j]*hd[j];
                                }
                                if (curv < eps_cg_curv*sq_norm_d) {
                                    is_d_neg_curv = (curv < -eps_cg_curv*sq_norm_d);
                                    warn_sing = !is_d_neg_curv;
                                }
                            } else {
                                warn_sing = true;
                            }

                        }

                    } else {
                        warn_small = true;
                    }
                }
            } else {
                warn_maxcgit = true;
            }
        } else {
            gd = -norm_g_non_act*norm_g_non_act/curv;
        }

    } else {
        
        // Hessian matrix not (sufficiently) positive definite
        warn_noposdef = true;

        if (curv < -eps_cg_curv) {
            gd = 0e0;
            for (unsigned int i=0; i<n_non_act; i++) {
                d[ind_non_act[i]] = g[ind_non_act[i]]/curv;
            }
            gd = norm_g_non_act*norm_g_non_act/curv;
            is_d_neg_curv = true;
        } else {
            for (unsigned int i=0; i<n_non_act; i++) {
                d[ind_non_act[i]] = -g[ind_non_act[i]];
            }
            gd = -norm_g_non_act*norm_g_non_act;
            warn_sing = true;
            warn_grad = true;
        }

    }

}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
void Asa_bcp::linesearch() {
    const double gamma = 1e-6;
    const double delta = 5e-1;
    const double f_dec = 1e6;
    double f_ref,aa,a1;
    bool red_param_computed;
    
    // N.B. 'v' is the current point, 'x' is the point obtained by using
    // the search direction with unit stepsize (and projecting onto the box)
    // and 'f' is f(v)
    
    red_param_computed = false;

    fv = f;
    prob->funct(status,x,f);
    if (status > 0) {
        flag = 9;
        return;
    }
    n_f++;

    f_ref = ((f-f_newton_first>=f_dec*delta_f0) || warn_noposdef) ? fv : f_w;
    
    stepsize = 1e0;

    while (f > f_ref+gamma*stepsize*gd) {

        if (n_f >= max_n_f) {
            flag = 5;
            break;
        }

        // update the stepsize
        if ((f-fv)/std::min(1e12,std::max(1e-12,-gamma*stepsize*gd)) < 1e10) {
            stepsize *= delta;
        } else {
            if (!red_param_computed) {
                a1 = 1e0*std::max(1e0,sqrt(std::inner_product(v.begin(),v.end(),v.begin(),0e0)))/std::max(1e-15,norm_proj_d);
                red_param_computed = true;
            }
            a1 = std::min(a1,stepsize*delta*delta);
            aa = 2e-12*stepsize;
            stepsize = std::max(aa,a1) > min_stepsize ? std::max(aa,a1) : delta*stepsize;
            
        }

        if (stepsize <= min_stepsize) {
            stepsize_exit = true;
            break;
        }
        
        // compute a new point
        for (unsigned int i=0; i<n_non_act; i++) {
            j = ind_non_act[i];
            x[j] = v[j] + stepsize*d[j];
            if (std::max(l_min,x[j])<l[j]) {
                x[j] = l[j];
            } else if (std::min(u_max,x[j])>u[j]) {
                x[j] = u[j];
            }
        }
        prob->funct(status,x,f);
        if (status > 0) {
            flag = 9;
            return;
        }
        n_f++;

    }
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
void Asa_bcp::approximate_hd(const std::vector<double>& p) { // approximate H(x)*p in the truncated-Newton method
    const double eps_approx = 1e-6;

    for (unsigned int i=0; i<n_non_act; i++) {
        j = ind_non_act[i];
        v[j] = x[j];
        x[j] += eps_approx*p[j];
    }
    prob->grad(status,x,g1);
    if (status > 0) {
        flag = 10;
        return;
    }
    n_g++;
    for (unsigned int i=0; i<n_non_act; i++) {
        hp[ind_non_act[i]] = (g1[ind_non_act[i]]-g[ind_non_act[i]])/eps_approx;
    }
    // restore x
    for (unsigned int i=0; i<n_non_act; i++) {
        x[ind_non_act[i]] = v[ind_non_act[i]];
    }
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
void Asa_bcp::err_obj() {
    x = x_best;
    f = f_best;
    sup_norm_proj_g = sup_norm_proj_g_best;
    if (flag == 9) {
        std::cout << "\n\nerror when computing the objective function, algorithm stopped";
        if (verbosity > 0) {
            file_output << "\n\nerror when computing the objective function, algorithm stopped";
        }
    } else if (flag == 10) {
        std::cout << "\n\nerror when computing the gradient of the objective function, algorithm stopped";
        if (verbosity > 0) {
            file_output << "\n\nerror when computing the gradient of the objective function, algorithm stopped";
        } 
    } else { // flag == 11
        std::cout << "\n\nerror when computing the Hessian-vector product, algorithm stopped";
        if (verbosity > 0) {
            file_output << "\n\nerror when computing the Hessian-vector product, algorithm stopped";
        }
    }
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
void Asa_bcp::clear_vectors() {
    // clear vectors to free memory (only 'x' is maintained)
    g.clear();
    g.shrink_to_fit();
    v.clear();
    v.shrink_to_fit();
    w.clear();
    w.shrink_to_fit();
    l.clear();
    u.shrink_to_fit();
    x_best.clear();
    x_best.shrink_to_fit();
    g_best.clear();
    g_best.shrink_to_fit();
    d.clear();
    d.shrink_to_fit();
    g_q.clear();
    g_q.shrink_to_fit();
    hp.clear();
    hp.shrink_to_fit();
    lambda.clear();
    lambda.shrink_to_fit();
    mu.clear();
    mu.shrink_to_fit();
    g1.clear();
    g1.shrink_to_fit();
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
const std::vector<double>& Asa_bcp::get_x() {
    return x;
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
const double Asa_bcp::get_f() {
    return f;
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
const double Asa_bcp::get_sup_norm_proj_g() {
    return sup_norm_proj_g;
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
const unsigned int Asa_bcp::get_it() {
    return it;
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
const unsigned int Asa_bcp::get_inner_it() {
    return it_cg_tot;
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
const unsigned int Asa_bcp::get_n_f() {
    return n_f;
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
const unsigned int Asa_bcp::get_n_g() {
    return n_g;
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
const unsigned int Asa_bcp::get_n_hd() {
    return n_hd;
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
const short int Asa_bcp::get_flag() {
    return flag;
}
//-------------------------------------------------------------------------------------