# ASA-BCP

_Active-Set Algorithm for Box-Constrained Problems_ (ASA-BCP) is a solver for bound-constrained
optimization problems of the following form:

    min f(x)
    l <= x <= u

where _f(x)_ is a twice continuously differentiable function.

ASA-BCP combines an active-set strategy with a truncated-Newton search direction and a non-monotone line search.

## Reference paper

[A. Cristofari, M. De Santis, S. Lucidi, F. Rinaldi (2017). _A Two-Stage Active-Set Algorithm for Bound-Constrained Optimization._
Journal of Optimization Theory and Applications, 172(2), 369-401.](https://link.springer.com/article/10.1007/s10957-016-1024-9)

## Authors

* Andrea Cristofari (e-mail: [andrea.cristofari@unipd.it](mailto:andrea.cristofari@unipd.it))
* Marianna De Santis (e-mail: [mdesantis@diag.uniroma1.it](mailto:mdesantis@diag.uniroma1.it))
* Stefano Lucidi (e-mail: [lucidi@diag.uniroma1.it](mailto:lucidi@diag.uniroma1.it))
* Francesco Rinaldi (e-mail: [rinaldi@math.unipd.it](mailto:rinaldi@math.unipd.it))

## How to use ASA-BCP

This directory should contain the files `COPYING.txt` and `README.md`,
plus the following subdirectories:

* `c++`, with a C++ implementation of ASA-BCP;
* `fortran`, with a Fortran implementation of ASA-BCP and an interface to
  CUTEst;
* `matlab`, with a Matlab implementation of ASA-BCP (a first way to run
  ASA-BCP via Matlab);
* `matlab_mex`, with a Matlab interface to the C++ implementation of
  ASA-BCP (a second way to run ASA-BCP via Matlab).

Each of these subdirectories contains a file `README.md` where it is explained how to use the program.
