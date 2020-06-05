#include "problem.h"

// constructor
//
// Required input/output arguments:
//   status (out) is 0 if no error occurred, >0 otherwise
//
//-------------------------------------------------------------------------------------
Problem::Problem(unsigned short int& status) {
    n = 1000; // problem dimension
    status = n > 0 ? 0 : 1;
}
//-------------------------------------------------------------------------------------


// ******************************************************
//            GENERALIZED ROSENBROCK FUNCTION
// ******************************************************


// OBJECTIVE FUNCTION
//
// Required input/output arguments:
//   status (out) is 0 if no error occurred, >0 otherwise
//   x (in) is the point where the objective function will be computed
//   f (out) is the value of objective function at x
//
//-------------------------------------------------------------------------------------
void Problem::funct(short unsigned int& status, const std::vector<double>& x, double& f) {
    const double c = 1e2;
    double t;

    f = 0e0;
    for (unsigned int i=0; i<n-1; i++) {
        t = x[i]*x[i];
        f += (x[i]-1e0)*(x[i]-1e0) + c*(x[i+1]-t)*(x[i+1]-t);
    }

    status = 0;
}
//-------------------------------------------------------------------------------------


// GRADIENT OF THE OBJECTIVE FUNCTION
//
// Required input/output arguments:
//   status (out) is 0 if no error occurred, >0 otherwise
//   x (in) is the point where the gradient of the objective function will be computed
//   g (out) is the gradient of the objective function at x
//
//-------------------------------------------------------------------------------------
void Problem::grad(short unsigned int& status, const std::vector<double>& x, std::vector<double>& g) {
    const double c = 1e2;

    g[0] = 2e0*(x[0]-1e0) + 4e0*c*x[0]*((x[0]*x[0])-x[1]);
    for (unsigned int i=1; i<n-1; i++) {
        g[i] = 2e0*(x[i]-1e0) + 4e0*c*x[i]*((x[i]*x[i])-x[i+1]) - 2e0*c*(x[i-1]*x[i-1]-x[i]);
    }
    g[n-1] = -2e0*c*(x[n-2]*x[n-2]-x[n-1]);

    status = 0;
}
//-------------------------------------------------------------------------------------


// HESSIAN-VECTOR PRODUCT
//
// Required input/output arguments:
//   status (out) is 0 if no error occurred, >0 otherwise
//   goth (in) is true if x is the same point used in the previous call, false otherwise
//             (it can be used to store values of the Hessian matrix)
//   x (in) is the point where the Hessian-vector product will be computed to be multiplied with d
//   d (in) is a vector whose product with the Hessian matrix is required
//   hd (out) is the Hessian-vector product
//
// N.B. This function can be ignored (in the sense that it can return any dummy value)
// if Hessian-vector products are approximated (i.e., if 'hd_exact' is equal to false
// in function 'asa_bcp')
//
//-------------------------------------------------------------------------------------
void Problem::hd_prod(short unsigned int& status, const bool goth, const std::vector<double>& x,
    const std::vector<double>& d, std::vector<double>& hd) {
    const double c = 1e2;
    static std::vector<double> hd1(n-1,0e0),hd2(n-1,0e0);

    if (!goth) {
        hd2[0] = -4e0*c*x[0];
        hd2[1] = -4e0*c*x[1];
        hd1[0] = 2e0 + 12e0*c*x[0]*x[0] + hd2[1];
        for (unsigned int i=1; i<n-2; i++) {
            hd2[i+1] = -4e0*c*x[i+1];
            hd1[i] = 2e0*(c+1e0) + 12e0*c*x[i]*x[i] + hd2[i+1];
        }
        hd1[n-2] = 2e0*(c+1e0) + 12e0*c*x[n-2]*x[n-2] - 4e0*c*x[n-1];
    }

    hd.assign(n,0e0);
    hd[0] = hd1[0]*d[0] + hd2[0]*d[1];
    hd[1] = hd2[0]*d[0];
    for (unsigned int i=1; i<n-1; i++) {
        hd[i] += hd1[i]*d[i] + hd2[i]*d[i+1];
        hd[i+1] += hd2[i]*d[i];
    }
    hd[n-1] += 2e0*c*d[n-1];

    status = 0;
}
//-------------------------------------------------------------------------------------


// BOUNDS
//
// Required input/output arguments:
//   status (out) is 0 if no error occurred, >0 otherwise
//   l (out) is the vector of lower bound for the variables
//   u (out) is the vector of upper bound for the variables
//
// N.B. There is no i-th lower bound if l(i) <= -1e20 and there is no i-th upper bound
// if u(i) >= 1e20.
//
//-------------------------------------------------------------------------------------
void Problem::bounds(short unsigned int& status, std::vector<double>& l, std::vector<double>& u) {
    for (unsigned int i=0; i<n-1; i+=2) {
        l[i] = -15e-1;
        l[i+1] = 5e-1;
        u[i] = 5e-1;
        u[i+1] = 2e0;
    }
    status = 0;
}
//-------------------------------------------------------------------------------------


// STARTING POINT
//
// Required input/output arguments:
//   status (out) is 0 if no error occurred, >0 otherwise
//   x (out) is the starting point for the algorithm
//
//-------------------------------------------------------------------------------------
void Problem::starting_point(short unsigned int& status, std::vector<double>& x){
    for (unsigned int i=0; i<n-1; i+=2) {
        x[i] = -12e-1;
        x[i+1] = 1e0;
    }
    status = 0;
}
//-------------------------------------------------------------------------------------