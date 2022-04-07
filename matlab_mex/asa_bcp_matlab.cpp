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


#include "mex.h"
#include <string>
#include "problem.h"
#include "../c++/asa_bcp.cpp"

#ifndef mxGetDoubles
#define mxGetDoubles mxGetPr
typedef double mxDouble;
#endif

// This is a MEX file for Matlab.
// See the file 'README.md' to know how to build the MEX file and run the
// program.

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    unsigned int n;
    unsigned short int status;
    const double *cdbl_ptr;
    bool obj_set,grad_set,hd_prod_set;
    const char *cchar_ptr;
    mxDouble *mdbl_ptr;
    mxArray *tmp_mxArray;

    // check the number of inupts and outputs
    if (nrhs<4) {
        mexErrMsgTxt("At least four inputs are required.\n");
    }
    if (nrhs>5) {
        mexErrMsgTxt("At most five inputs are required.\n");
    }
    if (nlhs<1) {
        mexErrMsgTxt("At least one input is required.\n");
    }
    if (nlhs>3) {
        mexErrMsgTxt("At most three outputs are required.\n");
    }

    // check inputs
    if (!mxIsStruct(prhs[0]) || mxGetNumberOfElements(prhs[0])>1) {
        mexErrMsgTxt("The first input must be a structure of function handle elements.");
    }
    if (mxGetNumberOfDimensions(prhs[1])>2 || mxGetN(prhs[1])!=1 || !mxIsDouble(prhs[1]) ||
        mxIsComplex(prhs[1]) || mxIsSparse(prhs[1])) {
        mexErrMsgTxt("The second input must be a real column vector.");
    }
    if (mxGetNumberOfDimensions(prhs[2])>2 || mxGetN(prhs[2])!=1 || !mxIsDouble(prhs[2]) ||
        mxIsComplex(prhs[2]) || mxIsSparse(prhs[2])) {
        mexErrMsgTxt("The third input must be a real column vector.");
    }
    if (mxGetNumberOfDimensions(prhs[3])>2 || mxGetN(prhs[3])!=1 || !mxIsDouble(prhs[3]) ||
        mxIsComplex(prhs[3]) || mxIsSparse(prhs[3])) {
        mexErrMsgTxt("The fourth input must be a real column vector.");
    }
    obj_set = grad_set = hd_prod_set = false;
    for (int i=0; i<mxGetNumberOfFields(prhs[0]); i++) {
        tmp_mxArray = mxGetFieldByNumber(prhs[0],0,i);
        cchar_ptr = mxGetFieldNameByNumber(prhs[0],i);
        if (std::string(cchar_ptr).compare(std::string("funct")) == 0) {
            if(!mxIsClass(tmp_mxArray,"function_handle")) {
                mexErrMsgTxt("In the structure passed as first input, 'funct' must be a function handle.");
            }
            obj_set = true;
        } else if (std::string(cchar_ptr).compare(std::string("grad")) == 0) {
            if(!mxIsClass(tmp_mxArray,"function_handle")) {
                mexErrMsgTxt("In the structure passed as first input, 'grad' must be a function handle.");
            }
            grad_set = true;
        } else if (std::string(cchar_ptr).compare(std::string("hd_prod")) == 0) {
            if(!mxIsClass(tmp_mxArray,"function_handle")) {
                mexErrMsgTxt("In the structure passed as first input, 'hd_prod' must be a function handle.");
            }
            hd_prod_set = true;
        } else {
            mexErrMsgTxt("Not valid field name in the structure passed as first input.");
        }
    }

    // get the problem dimension
    n = (unsigned int) mxGetM(prhs[1]);
    if (n != mxGetM(prhs[2])) {
        mexErrMsgTxt("Lower bound dimension and problem dimension must agree.");
    }
    if (n !=mxGetM(prhs[3])) {
        mexErrMsgTxt("Upper bound dimension and problem dimension must agree.");
    }

   // build an object of type Problem
    Problem p;
    *(p.dim) = n; // set dimension
    p.set_obj(prhs[0]); // set the objective function, the gradient and (if present) the Hessian-vector product
    p.set_starting_point(prhs[1]); // set starting point
    p.set_bounds(prhs[2],prhs[3]); // set bounds

    // get options
    asa_bcp_options opts;
    if (nrhs > 4) {
        if (!mxIsStruct(prhs[4]) || mxGetNumberOfElements(prhs[4])>1) {
            mexErrMsgTxt("The fifth input (which is optional) must be a structure.");
        }
        for (int i=0; i<mxGetNumberOfFields(prhs[4]); i++) {
            tmp_mxArray = mxGetFieldByNumber(prhs[4],0,i);
            cchar_ptr = mxGetFieldNameByNumber(prhs[4],i);
            if (std::string(cchar_ptr).compare(std::string("eps_opt")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("In the options, 'eps_opt' must be a non-negative number.");
                }
                opts.eps_opt = *mxGetDoubles(tmp_mxArray);
            } else if (std::string(cchar_ptr).compare(std::string("max_it")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("In the options, 'max_it' must be a number greater than or equal to 1.");
                }
                opts.max_it = (int)floor(*mxGetDoubles(tmp_mxArray));
            } else if (std::string(cchar_ptr).compare(std::string("max_n_f")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("In the options, 'max_n_f' must be a number greater than or equal to 1.");
                }
                opts.max_n_f = (int)floor(*mxGetDoubles(tmp_mxArray));
            } else if (std::string(cchar_ptr).compare(std::string("max_n_g")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("In the options, 'max_n_g' must bea number greater than or equal to 1.");
                }
                opts.max_n_g = (int)floor(*mxGetDoubles(tmp_mxArray));
            } else if (std::string(cchar_ptr).compare(std::string("max_n_hd")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("In the options, 'max_n_hd' must bea number greater than or equal to 0.");
                }
                opts.max_n_hd = (int)floor(*mxGetDoubles(tmp_mxArray));
            } else if (std::string(cchar_ptr).compare(std::string("min_f")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("In the options, 'min_f' must be a real number.");
                }
                opts.min_f = *mxGetDoubles(tmp_mxArray);
            } else if (std::string(cchar_ptr).compare(std::string("min_gd")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("In the options, 'min_gd' must be a non-negative number.");
                }
                opts.min_gd = *mxGetDoubles(tmp_mxArray);
            } else if (std::string(cchar_ptr).compare(std::string("min_norm_proj_d")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("In the options, 'min_norm_proj_d' must be a non-negative number.");
                }
                opts.min_norm_proj_d = *mxGetDoubles(tmp_mxArray);
            } else if (std::string(cchar_ptr).compare(std::string("min_stepsize")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("In the options, 'min_stepsize' must be a non-negative number.");
                }
                opts.min_stepsize = *mxGetDoubles(tmp_mxArray);
            } else if (std::string(cchar_ptr).compare(std::string("ls_memory")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("In the options, 'ls_memory' must be a number greater than or equal to 1.");
                }
                opts.ls_memory = (int)floor(*mxGetDoubles(tmp_mxArray));
            } else if (std::string(cchar_ptr).compare(std::string("z")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("In the options, 'z' must be a non-negative number.");
                }
                opts.z = (int)floor(*mxGetDoubles(tmp_mxArray));
            } else if (std::string(cchar_ptr).compare(std::string("hd_exact")) == 0) {
                if (!mxIsLogicalScalar(tmp_mxArray)) {
                    mexErrMsgTxt("In the options, 'hd_exact' must be a logical.");
                }
                opts.hd_exact = *mxGetLogicals(tmp_mxArray);
            } else if (std::string(cchar_ptr).compare(std::string("verbosity")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("In the options, 'verbosity' must be a number between 0 and 2.");
                }
                opts.verbosity = (int)floor(*mxGetDoubles(tmp_mxArray));
            } else {
                mexErrMsgTxt("Not valid field name in the structure of options.");
            }
        }
        if (!(hd_prod_set || !opts.hd_exact)) {
            mexErrMsgTxt("the Hessian-vector product must be specified (or set 'hd_exact' to false for approximating Hessian-vector products).");
        }
    }

    // build an object of type Asa_bcp
    Asa_bcp alg(status,&p,&opts);

    // check if an error occurred when calling the Asa_bcp constructor (something went wrong if 'status' > 0)
    if (status > 0) {
        mexErrMsgTxt("Something went wrong when calling the Asa_bcp constructor.");
    }

    // call the solver
    alg.solve();

    // set outputs
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    mdbl_ptr = mxGetDoubles(plhs[0]);
    cdbl_ptr = &(alg.get_x()[0]);
    for (unsigned int i=0; i<n; i++) {
        mdbl_ptr[i] = cdbl_ptr[i];
    }
    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleScalar(alg.get_f());
        if (nlhs > 2) {
            const char* asa_bcp_info_field_names[] = {"sup_norm_proj_g","it","inner_it","n_f","n_g","n_hd","flag"};
            mwSize asa_bcp_info_dims[2] = {1,1};
            plhs[2] = mxCreateStructArray(2,asa_bcp_info_dims,7,asa_bcp_info_field_names);
            mxSetField(plhs[2],0,asa_bcp_info_field_names[0],mxCreateDoubleScalar(alg.get_sup_norm_proj_g()));
            mxSetField(plhs[2],0,asa_bcp_info_field_names[1],mxCreateDoubleScalar((double)alg.get_it()));
            mxSetField(plhs[2],0,asa_bcp_info_field_names[2],mxCreateDoubleScalar((double)alg.get_inner_it()));
            mxSetField(plhs[2],0,asa_bcp_info_field_names[3],mxCreateDoubleScalar((double)alg.get_n_f()));
            mxSetField(plhs[2],0,asa_bcp_info_field_names[4],mxCreateDoubleScalar((double)alg.get_n_g()));
            mxSetField(plhs[2],0,asa_bcp_info_field_names[5],mxCreateDoubleScalar((double)alg.get_n_hd()));
            mxSetField(plhs[2],0,asa_bcp_info_field_names[6],mxCreateDoubleScalar((double)alg.get_flag()));
        }
    }
}

void Problem::set_obj(const mxArray* obj) {
    funct_fhandle = mxGetField(obj,0,"funct");
    grad_fhandle = mxGetField(obj,0,"grad");
    if (mxGetNumberOfFields(obj) == 3) {
        hd_prod_fhandle = mxGetField(obj,0,"hd_prod");
    }
}

void Problem::set_starting_point(const mxArray* x0_in){
    mxDouble* x0_ptr = mxGetDoubles(x0_in);
    x0.resize(n);
    for (unsigned int i=0; i<n; i++){
        x0[i] = x0_ptr[i];
    }
}

void Problem::set_bounds(const mxArray* l_in, const mxArray* u_in){
    mxDouble *l_ptr,*u_ptr;
    l_ptr = mxGetDoubles(l_in);
    u_ptr = mxGetDoubles(u_in);
    l.resize(n);
    u.resize(n);
    for (unsigned int i=0; i<n; i++){
        l[i] = l_ptr[i];
        u[i] = u_ptr[i];
    }
}

void Problem::funct(short unsigned int& status, const std::vector<double>& x, double& f) {
    mxDouble *x_ptr;
    mxArray *lhs[2],*rhs[2];

    rhs[0] = funct_fhandle;
    rhs[1] = mxCreateDoubleMatrix(n,1,mxREAL);
    x_ptr = mxGetDoubles(rhs[1]);
    for (unsigned int i=0; i<n; i++) {
        x_ptr[i] = x[i];
    }
    status = (short unsigned int) mexCallMATLAB(2,lhs,2,rhs,"feval");
    if (status == 0) {
        status = (short unsigned int) *mxGetDoubles(lhs[1]);
    }
    f = status == 0 ? *mxGetDoubles(lhs[0]) : nan("");
    mxDestroyArray(lhs[0]);
    mxDestroyArray(lhs[1]);
}

void Problem::grad(short unsigned int& status, const std::vector<double>& x, std::vector<double>& g) {
    mxDouble *xg_ptr;
    mxArray *lhs[2],*rhs[2];

    rhs[0] = grad_fhandle;
    rhs[1] = mxCreateDoubleMatrix(n,1,mxREAL);
    xg_ptr = mxGetDoubles(rhs[1]);
    for (unsigned int i=0; i<n; i++) {
        xg_ptr[i] = x[i];
    }
    status = (short unsigned int) mexCallMATLAB(2,lhs,2,rhs,"feval");
    if (status == 0) {
        status = (short unsigned int) *mxGetDoubles(lhs[1]);
    }
    if (status == 0) {
        g.resize(n);
        xg_ptr = mxGetDoubles(lhs[0]);
        for (unsigned int i=0; i<n; i++) {
            g[i] = xg_ptr[i];
        }
    } else {
        g.assign(n,nan(""));
    }
    mxDestroyArray(lhs[0]);
    mxDestroyArray(lhs[1]);
}

void Problem::hd_prod(short unsigned int& status, bool goth, const std::vector<double>& x, const std::vector<double>& d, std::vector<double>& hd) {
    mxDouble *xhd_ptr,*d_ptr;
    mxArray *lhs[2],*rhs[4];

    rhs[0] = hd_prod_fhandle;
    rhs[1] = mxCreateLogicalScalar(goth);
    rhs[2] = mxCreateDoubleMatrix(n,1,mxREAL);
    rhs[3] = mxCreateDoubleMatrix(n,1,mxREAL);
    xhd_ptr = mxGetDoubles(rhs[2]);
    d_ptr = mxGetDoubles(rhs[3]);
    for (unsigned int i=0; i<n; i++) {
        xhd_ptr[i] = x[i];
        d_ptr[i] = d[i];
    }
    status = (short unsigned int) mexCallMATLAB(2,lhs,4,rhs,"feval");
    if (status == 0) {
        status = (short unsigned int) *mxGetDoubles(lhs[1]);
    }
    if (status == 0) {
        hd.resize(n);
        xhd_ptr = mxGetDoubles(lhs[0]);
        for (unsigned int i=0; i<n; i++) {
            hd[i] = xhd_ptr[i];
        }
    } else {
        hd.assign(n,nan(""));
    }
    mxDestroyArray(lhs[0]);
    mxDestroyArray(lhs[1]);
}

void Problem::starting_point(short unsigned int& status, std::vector<double>& x) {
    x = x0;
    status = 0;
}

void Problem::bounds(short unsigned int& status, std::vector<double>& lb, std::vector<double>& ub) {
    lb = l;
    ub = u;
    status = 0;
}