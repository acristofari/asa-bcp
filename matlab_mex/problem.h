#ifndef __PROBLEM_H_INCLUDED__
#define __PROBLEM_H_INCLUDED__
#include "mex.h"
#include <vector>

class Problem {
private:
    unsigned int n;
    std::vector<double> x0,l,u;
    mxArray *funct_fhandle,*grad_fhandle,*hd_prod_fhandle;
public:
    unsigned int* const dim = &n;
    void funct(short unsigned int&, const std::vector<double>&, double&);
    void grad(short unsigned int&, const std::vector<double>&, std::vector<double>&);
    void hd_prod(short unsigned int&, const bool, const std::vector<double>&, const std::vector<double>&, std::vector<double>&);
    void bounds(short unsigned int&, std::vector<double>&, std::vector<double>&);
    void starting_point(short unsigned int&, std::vector<double>&);
    void set_obj(const mxArray*);
    void set_starting_point(const mxArray*);
    void set_bounds(const mxArray*,const mxArray*);
};

#endif