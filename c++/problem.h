#ifndef __PROBLEM_H_INCLUDED__
#define __PROBLEM_H_INCLUDED__
#include <vector>

class Problem {
private:
    unsigned int n;
public:
    Problem(unsigned short int&);
    unsigned int* const dim = &n;
    void funct(short unsigned int&, const std::vector<double>&, double&);
    void grad(short unsigned int&, const std::vector<double>&, std::vector<double>&);
    void hd_prod(short unsigned int&, const bool, const std::vector<double>&, const std::vector<double>&, std::vector<double>&);
    void bounds(short unsigned int&, std::vector<double>&, std::vector<double>&);
    void starting_point(short unsigned int&, std::vector<double>&);
};

#endif