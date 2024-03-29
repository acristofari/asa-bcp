# Active-Set Algorithm for Box-Constrained Problems (ASA-BCP)

_Active-Set Algorithm for Box-Constrained Problems_ (ASA-BCP) is a solver for bound-constrained
optimization problems of the following form:

         min f(x)
    s.t. l <= x <= u

with given vectors _l, u_ and where _f(x)_ is a twice continuously differentiable function.

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

Copyright 2017-2022 Andrea Cristofari, Marianna De Santis,
Stefano Lucidi, Francesco Rinaldi.

## How to use ASA-BCP (C++ implementation)

1. This directory should contain the following files:
    * `asa_bcp.cpp`,
    * `asa_bcp.h`,
    * `main.cpp`,
    * `problem.cpp`,
    * `problem.h`,
    * `README.md`,
    * `syntax.txt`.

2. See the file `usage.txt` to know how to call ASA-BCP in C++, change
   algorithm parameters and get output values.

3. See the files `main.cpp` and `problem.cpp` for an example.
   To run the example, first create the executable file by compiling `asa_bcp.cpp`,
   `main.cpp` and `problem.cpp`, then run the executable file.
   The solution found by the algorithm will be reported in the file
   `opt_sol.txt` and the final statistics will be reported in the file
   `statistics.txt`.