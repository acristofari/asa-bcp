# Active-Set Algorithm for Box-Constrained Problems (ASA-BCP)

_Active-Set Algorithm for Box-Constrained Problems_ (ASA-BCP) is a solver for bound-constrained
optimization problems of the following form:

         min f(x)
    s.t. l <= x <= u

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

## Licensing

ASA-BCP is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
ASA-BCP is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with ASA-BCP. If not, see <http://www.gnu.org/licenses/>.

Copyright 2017-2021 Andrea Cristofari, Marianna De Santis,
Stefano Lucidi, Francesco Rinaldi.

## How to use ASA-BCP (via C++ source MEX file for Matlab)

1. This directory should contain the following files:
    * `asa_bcp_matlab.cpp`,
    * `make.m`,
    * `problem.h`,
    * `README.md`.

    You will also need the files `asa_bcp.h` and `asa_bcp.cpp`, located in
    `../c++`, and the files `main.m` and `problem.m`, located in
    `../matlab`.

2. Run `make.m` to build the MEX file.

3. See the file `../matlab/syntax.txt` to know how to call ASA-BCP,
    change algorithm parameters and get output values.

4. Copy the files `main.m` and `problem.m` from `../matlab` to the current
   directory for running an example. After copying the files, to run the 
   example just run `main.m`. The solution found by the algorithm will be
   reported in the file `opt_sol.txt` and the final statistics will be
   reported in the file `statistics.txt`.
