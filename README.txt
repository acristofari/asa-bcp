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


How to use ASA-BCP
==========================================================================

This directory should contain the files 'COPYING.txt' and 'README.txt',
plus the following subdirectories:

- 'c++', with a C++ implementation of ASA-BCP;
- 'fortran', with a Fortran implementation of ASA-BCP and an interface to
  CUTEst;
- 'matlab', with a Matlab implementation of ASA-BCP (a first way to run
  ASA-BCP via Matlab);
- 'matlab_mex', with a Matlab interface to the C++ implementation of
  ASA-BCP (a second way to run ASA-BCP via Matlab).

Each of these subdirectories contains a file 'README.txt' where it is
explained how to use the program.
