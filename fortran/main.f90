! -------------------------------------------------------------------------
! 
! This file is part of ASA-BCP, which is a solver for bound-constrained
! optimization problems of the following form:
! 
!                                min f(x)
!                           s.t. l <= x <= u
! 
! with given vectors l, u and where f(x) is a twice continuously
! differentiable function.
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
! April 7th, 2022
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
! Copyright 2017-2022 Andrea Cristofari, Marianna De Santis,
! Stefano Lucidi, Francesco Rinaldi.
! 
! -------------------------------------------------------------------------


program main

    use asa_bcp_opts

    implicit none
    
    integer :: status,n,it,n_f,n_g,n_hd,inner_it,flag
    double precision, allocatable :: x(:),l(:),u(:)
    double precision :: f,sup_norm_proj_g
    type(asa_bcp_options) :: opts
    
    ! In this file, it is shown how to call ASA-BCP to solve a user-defined problem
    
    ! (1) Get problem dimension
    call prob_dim(status,n)
    
    ! (2) Check if an error occurred with the problem dimension
    !     (something went wrong if 'status' > 0)
    if (status > 0) then
        write(*,*) "Something went wrong with the problem dimension."
        stop
    endif
    
    ! (3) Allocate vectors
    allocate(x(n),l(n),u(n))
    
    ! (4) Get lower and upper bounds
    call bounds(status,n,l,u)
    
    ! (5) Check if an error occurred with the lower and upper bounds
    !     (something went wrong if 'status' > 0)
    if (status > 0) then
        write(*,*) "Something went wrong with the lower and upper bounds."
        deallocate(x,l,u)
        stop
    endif
    
    ! (6) Get the starting point
    call starting_point(status,n,x)
    
    ! (7) Check if an error occurred with the starting point
    !     (something went wrong if 'status' > 0)
    if (status > 0) then
        write(*,*) "Something went wrong with the starting point."
        deallocate(x,l,u)
        stop
    endif
    
    ! ------------------------------------------------------------------------------------
    ! *** EXAMPLE OF HOW TO CHANGE ASA-BCP PARAMETERS ***
    ! (see the description of asa_bcp in the file 'usage.txt' to know which parameters can
    ! be changed and their default values)
    !
    ! Assign new values to (some of) the members of the object 'opts' of derived data type
    ! 'asa_bcp_options' (see its declaration in the module 'asa_bcp_opts'), e.g.,
    !
    !   opts%verbosity = 0
    ! ------------------------------------------------------------------------------------
    
    ! (8) call ASA-BCP
    call asa_bcp(n,x,f,l,u,opts,flag)
        
    if (flag.eq.-1) then ! infeasible problem
        deallocate(x,l,u)
        stop
    endif
        
    ! get statistics
    call get_sup_norm_proj_g(sup_norm_proj_g)
    call get_it(it)
    call get_n_f(n_f)
    call get_n_g(n_g)
    call get_n_hd(n_hd)
    call get_inner_it(inner_it)
    
    ! write statistics to the screen and to file 'statistics.txt'
    open(11,file='statistics.txt')
    write(*,200) n,f,sup_norm_proj_g,it,n_f,n_g,n_hd,inner_it,flag
    write(11,200) n,f,sup_norm_proj_g,it,n_f,n_g,n_hd,inner_it,flag
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
    ' exit flag = ', i1, //, 1x, 50('*'))

end program