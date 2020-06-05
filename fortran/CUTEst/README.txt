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


How to run ASA-BCP on CUTEst problems (for Linux system)
==========================================================================

Before continuing, make sure you have downloaded CUTEst with the related
packages to your computer and have installed all correctly.

1 - This directory should contain the following files:

    - 'asa_bcp',
    - 'asa_bcp_main.f90',
    - 'makemaster',
    - 'README.txt'.

    You will also need the file 'asa_bcp.f90', located in the parent
    directory.

2 - Copy the file 'asa_bcp' to '[CUTEst_path]/cutest/packages/defaults'.

3 - In '[CUTEst_path]/cutest/src', create the directory 'asa_bcp'.

4 - Copy the files 'asa_bcp.f90' (from the parent directory),
    'asa_bcp_main.f90' and 'makemaster' to '[CUTEst_path]/cutest/src/asa_bcp'.

5 - In '[CUTEst_path]/cutest/src/asa_bcp', compile both the .f90 files to
    create the object files 'asa_bcp.o' and 'asa_bcp_main.o'.

6 - Now you can run ASA-BCP on CUTEst problems. For instance, on the
    command prompt you may type

    runcutest -p asa_bcp -D problem_name

    where 'problem_name' is the CUTEst problem you want to solve.

See the file 'syntax.txt' in the parent directory and the file
'asa_bcp_main.f90' to know how to set the problem, change algorithm
parameters and get output values.
