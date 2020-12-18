! -------------------------------------------------------------------------
!
! This file is part of ASA-BCP, which is a solver for bound-constrained
! optimization problems of the following form:
!
!                                 min f(x)
!                           s.t. l <= x <= u
!
! with f(x) twice continuously differentiable.
!
! This is a driver for running ASA-BCP on user-defined problems.
! See the file 'README.txt' to know how to run the program.
!
! -------------------------------------------------------------------------
!
! Reference paper:
!
! A. Cristofari, M. De Santis, S. Lucidi, F. Rinaldi (2017). A Two-Stage
! Active-Set Algorithm for Bound-Constrained Optimization. Journal of
! Optimization Theory and Applications, 172(2), 369-401.
!
! -------------------------------------------------------------------------
!
! Authors:
! Andrea Cristofari (e-mail: andrea.cristofari@unipd.it)
! Marianna De Santis (e-mail: mdesantis@diag.uniroma1.it)
! Stefano Lucidi (e-mail: lucidi@diag.uniroma1.it)
! Francesco Rinaldi (e-mail: rinaldi@math.unipd.it)
!
! Last update of this file:
! December 18th, 2020
!
! Copyright 2017-2020 Andrea Cristofari, Marianna De Santis,
! Stefano Lucidi, Francesco Rinaldi.
!
! Licensing:
! This file is part of ASA-BCP.
! ASA-BCP is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! ASA-BCP is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with ASA-BCP. If not, see <http://www.gnu.org/licenses/>.
!
! -------------------------------------------------------------------------


program main

    use asa_bcp_opts

    implicit none
    
    integer :: status,n,it,n_f,n_g,n_hd,inner_it,flag
    double precision, allocatable :: x(:),l(:),u(:)
    double precision :: f,sup_norm_proj_g
    real :: time_start,time_end,elap_time
    type(asa_bcp_options) :: opts
    
    ! get problem dimension
    call prob_dim(status,n)
    
    ! check if an error occurred with the problem dimension
    ! (something went wrong if 'status' > 0)
    if (status > 0) then
        write(*,*) "error with the problem dimension"
        stop
    endif
    
    ! allocate vectors
    allocate(x(n),l(n),u(n))
    
    ! get the lower and upper bounds
    call bounds(status,n,l,u)
    ! check if an error occurred with the lower and upper bounds
    ! (something went wrong if 'status' > 0)
    if (status > 0) then
        write(*,*) "error with the lower and upper bounds"
        deallocate(x,l,u)
        stop
    endif
    
    ! get the starting point
    call starting_point(status,n,x)
    
    ! check if an error occurred with the starting point
    ! (something went wrong if 'status' > 0)
    if (status > 0) then
        write(*,*) "error with the starting point"
        deallocate(x,l,u)
        stop
    endif
    
    ! ------------------------------------------------------------------------------------
    ! *** EXAMPLE OF HOW TO CHANGE ASA-BCP PARAMETERS ***
    ! (see the description of asa_bcp in the file 'asa_bcp.f90' to know which parameters can
    ! be changed and their default values)
    !
    ! Assign new values to (some of) the members of the object 'opts' of derived data type
    ! 'asa_bcp_options' (see its declaration in the module 'asa_bcp_opts'), e.g.,:
    !
    !   opts%verbosity = 0
    ! ------------------------------------------------------------------------------------
    
    ! call the solver
    call cpu_time(time_start)
    call asa_bcp(n,x,f,l,u,opts,flag)
    call cpu_time(time_end)
    
    if (flag.eq.-1) then ! infeasible problem
        deallocate(x,l,u)
        stop
    endif
    
    ! compute the elapsed time
    elap_time = amax1(time_end-time_start,0.e0)
    
    ! get statistics
    call get_sup_norm_proj_g(sup_norm_proj_g)
    call get_it(it)
    call get_n_f(n_f)
    call get_n_g(n_g)
    call get_n_hd(n_hd)
    call get_inner_it(inner_it)
    
    ! write statistics to the screen and to file 'statistics.txt'
    open(11,file='statistics.txt')
    write(*,200) n,f,sup_norm_proj_g,it,n_f,n_g,n_hd,inner_it,flag,elap_time
    write(11,200) n,f,sup_norm_proj_g,it,n_f,n_g,n_hd,inner_it,flag,elap_time
    close(11)
    
    ! write the solution found by the algorithm to file 'opt_sol.txt'
    open(12,file='opt_sol.txt')
    write(12,'(d13.6)') x
    close(12)
    
    ! free allocated memory
    deallocate(x,l,u)
    
    ! non-executable statement
200 format(1x, 50('*') //, &
    ' Algorithm: ASA-BCP', //, &
    ' number of variables = ', i0, //, &
    ' f = ', d13.6, /, &
    ' sup-norm of the projected gradient = ', d13.6, /, &
    ' number of iterations = ', i0, /, &
    ' number of function evaluations = ', i0, /, &
    ' number of gradient evaluations = ', i0, /, &
    ' number of Hessian-vector products = ', i0, /, &
    ' number of inner cg iterations = ', i0, /, &
    ' exit flag = ', i1, /, &
    ' elapsed time (s) = ', e12.5, //, 1x, 50('*'))

end program