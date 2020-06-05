ASA-BCP is a solver for bound-constrained optimization problems of the
following form:

                                min f(x)
                          s.t. l <= x <= u

with f(x) twice continuously differentiable.

--------------------------------------------------------------------------

Reference paper:

A. Cristofari, M. De Santis, S. Lucidi, F. Rinaldi (2017). A Two-Stage
Active-Set Algorithm for Bound-Constrained Optimization. Journal of
Optimization Theory and Applications, 172(2), 369-401.

--------------------------------------------------------------------------

Authors:
Andrea Cristofari (e-mail: andrea.cristofari@unipd.it)
Marianna De Santis (e-mail: mdesantis@diag.uniroma1.it)
Stefano Lucidi (e-mail: lucidi@diag.uniroma1.it)
Francesco Rinaldi (e-mail: rinaldi@math.unipd.it)

--------------------------------------------------------------------------


How to use ASA-BCP (C++ implementation)
==========================================================================

1 - This directory should contain the following files:

    - 'asa_bcp.cpp',
    - 'asa_bcp.h',
    - 'main.cpp',
    - 'problem.cpp',
    - 'problem.h',
    - 'README.txt',
    - 'syntax.txt'.

2 - See the file 'syntax.txt' to know how to set the problem, change
    algorithm parameters and get output values.

3 - See the files 'main.cpp' and 'problem.cpp' for an example. To run the
    example, first create the executable file by compiling 'asa_bcp.cpp',
    'main.cpp' and 'problem.cpp', and then run the executable file. The
    solution found by the algorithm will be reported in the file
    'opt_sol.txt' and the final statistics will be reported in the file
    'statistics.txt'.
