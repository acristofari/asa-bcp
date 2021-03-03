// -------------------------------------------------------------------------
// 
// This file is part of ASA-BCP, which is a solver for bound-constrained
// optimization problems of the following form:
// 
//                                min f(x)
//                           s.t. l <= x <= u
// 
// where f(x) is a twice continuously differentiable.
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
// March 3rd, 2021
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
// Copyright 2017-2021 Andrea Cristofari, Marianna De Santis,
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
// See the file 'README.txt' to know how to build the MEX file and run the
// program.

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    unsigned int n;
    unsigned short int status;
    bool obj_set,grad_set,hd_prod_set;
    const char *tmp_char;
    mxDouble *x_ptr_out;
    mxArray *tmp_mxArray;

    // check the number of inupts and outputs
    if (nrhs<4 || nrhs>5) {
        mexErrMsgTxt("error when calling asa_bcp: the number of input arguments must be either 4 or 5.\n");
    }
    if (nlhs<1 || nlhs>3) {
        mexErrMsgTxt("error when calling asa_bcp: the number of output arguments must be between 1 and 3.\n");
    }

    // check inputs
    if (!mxIsStruct(prhs[0]) || mxGetNumberOfElements(prhs[0])>1) {
        mexErrMsgTxt("error when calling asa_bcp: the first input argument must be a structure of function handle elements.");
    }
    if (mxIsScalar(prhs[1]) || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
        mxGetNumberOfDimensions(prhs[1])>2 || mxIsSparse(prhs[1]) || mxGetN(prhs[1])!=1) {
        mexErrMsgTxt("error when calling asa_bcp: the second input argument must be a full column vector of real numbers.");
    }
    if (mxIsScalar(prhs[2]) || !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
        mxGetNumberOfDimensions(prhs[2])>2 || mxIsSparse(prhs[2]) || mxGetN(prhs[2])!=1) {
        mexErrMsgTxt("error when calling asa_bcp: the third input argument must be a full column vector of real numbers.");
    }
    if (mxIsScalar(prhs[3]) || !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
        mxGetNumberOfDimensions(prhs[3])>2 || mxIsSparse(prhs[3]) || mxGetN(prhs[3])!=1) {
        mexErrMsgTxt("error when calling asa_bcp: the fourth input argument must be a full column vector of real numbers.");
    }
    obj_set = grad_set = hd_prod_set = false;
    for (int i=0; i<mxGetNumberOfFields(prhs[0]); i++) {
        tmp_mxArray = mxGetFieldByNumber(prhs[0],0,i);
        tmp_char = mxGetFieldNameByNumber(prhs[0],i);
        if (std::string(tmp_char).compare(std::string("funct")) == 0) {
            if(!mxIsClass(tmp_mxArray,"function_handle")) {
                mexErrMsgTxt("error when calling asa_bcp: the objective function must be a function handle.");
            }
            obj_set = true;
        } else if (std::string(tmp_char).compare(std::string("grad")) == 0) {
            if(!mxIsClass(tmp_mxArray,"function_handle")) {
                mexErrMsgTxt("error when calling asa_bcp: the gradient of the objective function must be a function handle.");
            }
            grad_set = true;
        } else if (std::string(tmp_char).compare(std::string("hd_prod")) == 0) {
            if(!mxIsClass(tmp_mxArray,"function_handle")) {
                mexErrMsgTxt("error when calling asa_bcp: the Hessian-vector product must be a function handle.");
            }
            hd_prod_set = true;
        } else {
            mexErrMsgTxt("error when calling asa_bcp: not valid field name present in the structure of the first input argument.");
        }
    }

    // get the problem dimension
    n = (unsigned int) mxGetM(prhs[1]);
    if (n!=mxGetM(prhs[2]) || n!=mxGetM(prhs[3])) {
        mexErrMsgTxt("error when calling asa_bcp: the dimension of the starting point and the dimension of the variable bounds must agree.");
    }

   // build an object of type Problem
    Problem p;
    *(p.dim) = n; // set dimension
    p.set_obj(prhs[0]); // set the objective function, the gradient and (if present) the Hessian-vector product
    p.set_starting_point(prhs[1]); // set starting point
    p.set_bounds(prhs[2],prhs[3]); // set bounds

    // set asa-bcp options
    asa_bcp_options opts;
    if (nrhs > 4) {
        if (!mxIsStruct(prhs[4]) || mxGetNumberOfElements(prhs[4])>1) {
            mexErrMsgTxt("error when calling asa_bcp: the fifth input argument (which is optional) must be a structure.");
        }
        for (int i=0; i<mxGetNumberOfFields(prhs[4]); i++) {
            tmp_mxArray = mxGetFieldByNumber(prhs[4],0,i);
            tmp_char = mxGetFieldNameByNumber(prhs[4],i);
            if (std::string(tmp_char).compare(std::string("eps_opt")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("error when calling asa_bcp: 'eps_opt' must be a number greater than or equal to 0.");
                }
                opts.eps_opt = *mxGetDoubles(tmp_mxArray);
            } else if (std::string(tmp_char).compare(std::string("min_gd")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("error when calling asa_bcp: 'min_gd' must be a number greater than or equal to 0.");
                }
                opts.min_gd = *mxGetDoubles(tmp_mxArray);
            } else if (std::string(tmp_char).compare(std::string("min_norm_proj_d")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("error when calling asa_bcp: 'min_norm_proj_d' must be a number greater than or equal to 0.");
                }
                opts.min_norm_proj_d = *mxGetDoubles(tmp_mxArray);
            } else if (std::string(tmp_char).compare(std::string("min_stepsize")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("error when calling asa_bcp: 'min_stepsize' must be a number greater than or equal to 0.");
                }
                opts.min_stepsize = *mxGetDoubles(tmp_mxArray);
            } else if (std::string(tmp_char).compare(std::string("max_it")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("error when calling asa_bcp: 'max_it' must be a number greater than or equal to 0.");
                }
                opts.max_it = (int)floor(*mxGetDoubles(tmp_mxArray));
            } else if (std::string(tmp_char).compare(std::string("max_n_f")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("error when calling asa_bcp: 'max_n_f' must be a number greater than or equal to 1.");
                }
                opts.max_n_f = (int)floor(*mxGetDoubles(tmp_mxArray));
            } else if (std::string(tmp_char).compare(std::string("max_n_g")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("error when calling asa_bcp: 'max_n_g' must be a number greater than or equal to 1.");
                }
                opts.max_n_g = (int)floor(*mxGetDoubles(tmp_mxArray));
            } else if (std::string(tmp_char).compare(std::string("max_n_hd")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("error when calling asa_bcp: 'max_n_hd' must be a number greater than or equal to 0.");
                }
                opts.max_n_hd = (int)floor(*mxGetDoubles(tmp_mxArray));
            } else if (std::string(tmp_char).compare(std::string("min_f")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("error when calling asa_bcp: 'min_f' must be a real number.");
                }
                opts.min_f = *mxGetDoubles(tmp_mxArray);
            } else if (std::string(tmp_char).compare(std::string("m")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("error when calling asa_bcp: 'm' must be a number greater than or equal to 1.");
                }
                opts.m = (int)floor(*mxGetDoubles(tmp_mxArray));
            } else if (std::string(tmp_char).compare(std::string("z")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("error when calling asa_bcp: 'z' must be a number greater than or equal to 1.");
                }
                opts.z = (int)floor(*mxGetDoubles(tmp_mxArray));
            } else if (std::string(tmp_char).compare(std::string("hd_exact")) == 0) {
                if (!mxIsLogicalScalar(tmp_mxArray)) {
                    mexErrMsgTxt("error when calling asa_bcp: 'hd_exact' must be a logical.");
                }
                opts.hd_exact = *mxGetLogicals(tmp_mxArray);
            } else if (std::string(tmp_char).compare(std::string("verbosity")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("error when calling asa_bcp: 'verbosity' must be a number between 0 and 2.");
                }
                opts.verbosity = (int)floor(*mxGetDoubles(tmp_mxArray));
            } else {
                mexErrMsgTxt("error when calling asa_bcp: not valid field name in the structure of the fifth input argument (which is optional).");
            }
        }
        if (!(hd_prod_set || !opts.hd_exact)) {
            mexErrMsgTxt("error when calling asa_bcp: the Hessian-vector product must be specified (or set 'hd_exact' to false for approximating Hessian-vector products).");
        }
    }

    // build an object of type Asa_bcp
    Asa_bcp alg(status,&p,&opts);

    // check if an error occurred when calling the Asa_bcp constructor (something went wrong if 'status' > 0)
    if (status > 0) {
        mexErrMsgTxt("error when calling the Asa_bcp constructor");
    }

    // call the solver
    alg.solve();

    // assign value to outputs
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    x_ptr_out = mxGetDoubles(plhs[0]);
    const double *x_ptr = &(alg.get_x()[0]);
    for (unsigned int i=0; i<n; i++) {
        *(x_ptr_out+i) = *(x_ptr+i);
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
        x0[i] = *(x0_ptr+i);
    }
}

void Problem::set_bounds(const mxArray* l_in, const mxArray* u_in){
    mxDouble *l_ptr,*u_ptr;
    l_ptr = mxGetDoubles(l_in);
    u_ptr = mxGetDoubles(u_in);
    l.resize(n);
    u.resize(n);
    for (unsigned int i=0; i<n; i++){
        l[i] = *(l_ptr+i);
        u[i] = *(u_ptr+i);
    }
}

void Problem::funct(short unsigned int& status, const std::vector<double>& x, double& f) {
    mxDouble *x_ptr;
    mxArray *lhs[2],*rhs[2];

    rhs[0] = funct_fhandle;
    rhs[1] = mxCreateDoubleMatrix(n,1,mxREAL);
    x_ptr = mxGetDoubles(rhs[1]);
    for (unsigned int i=0; i<n; i++) {
        *(x_ptr+i) = x[i];
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
        *(xg_ptr+i) = x[i];
    }
    status = (short unsigned int) mexCallMATLAB(2,lhs,2,rhs,"feval");
    if (status == 0) {
        status = (short unsigned int) *mxGetDoubles(lhs[1]);
    }
    if (status == 0) {
        g.resize(n);
        xg_ptr = mxGetDoubles(lhs[0]);
        for (unsigned int i=0; i<n; i++) {
            g[i] = *(xg_ptr+i);
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
        *(xhd_ptr+i) = x[i];
        *(d_ptr+i) = d[i];
    }
    status = (short unsigned int) mexCallMATLAB(2,lhs,4,rhs,"feval");
    if (status == 0) {
        status = (short unsigned int) *mxGetDoubles(lhs[1]);
    }
    if (status == 0) {
        hd.resize(n);
        xhd_ptr = mxGetDoubles(lhs[0]);
        for (unsigned int i=0; i<n; i++) {
            hd[i] = *(xhd_ptr+i);
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