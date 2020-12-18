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
! This is a driver for running ASA-BCP on CUTEst problems.
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


program asa_bcp_main

    use asa_bcp_opts

    implicit none
    integer :: i,n,m,status
    integer :: flag,it,n_f,n_g,n_hd,inner_it
    integer, parameter :: input=47,io_buffer=11,out=6
    double precision :: f,sup_norm_proj_g
    double precision, dimension(:), allocatable :: x,lb,ub
    real :: time_start,time_end,elap_time
    character(len=10) :: pname
    type(asa_bcp_options) :: opts
    
    ! open problem file
    open(input,file='OUTSDIF.d',form='FORMATTED',status='OLD')
    rewind input
    
    ! check if the problem has only bound constraints
    call CUTEST_cdimen(status,input,n,m)
    if (status.ne.0) goto 910
    if (m.gt.0) then
        write(6,'(a)') 'ERROR: ASA-BCP can only handle bound constraints'
        stop
    elseif (m.lt.0) then
        write(6,'(a)') 'ERROR: reading OUTSDIF.d'
        stop
    endif
    
    ! set up SIF data from the problem file
    allocate(x(n),lb(n),ub(n))
    call CUTEST_usetup(status,input,out,io_buffer,n,x,lb,ub)
    if (status.ne.0) then
        goto 910
    endif
        
    ! get problem name
    call CUTEST_pname(status,input,pname)
    if (status.ne.0) then
        goto 910
    endif
    
    ! call the solver
    opts%verbosity = 0 ! to suppress verbosity
    call cpu_time(time_start)
    call asa_bcp(n,x,f,lb,ub,opts,flag)
    call cpu_time(time_end)
    
    ! compute the elapsed time
    elap_time = amax1(time_end-time_start,0.e0)
    
    ! get statistics
    call get_sup_norm_proj_g(sup_norm_proj_g)
    call get_it(it)
    call get_n_f(n_f)
    call get_n_g(n_g)
    call get_n_hd(n_hd)
    call get_inner_it(inner_it)
    
    ! write statistics to the screen and to the file 'statistics.txt'
    open(11,file='statistics.txt')
    write(*,200) pname,n,f,sup_norm_proj_g,it,n_f,n_g,n_hd,inner_it,flag,elap_time
    write(11,200) pname,n,f,sup_norm_proj_g,it,n_f,n_g,n_hd,inner_it,flag,elap_time
    close(11)
    
    ! write the solution found by the algorithm to the file 'opt_sol.txt'
    open(12,file='opt_sol.txt')
    write(12,'(d13.6)') x
    close(12)
        
    ! close problem file
    close(input)
    
    ! free allocated memory
    deallocate(x,lb,ub)
    
    ! exit
    call CUTEST_uterminate(status)
    stop
    
    910 continue
    write(out,"('CUTEst error, status = ',i0,', stopping')") status
    stop
    
! non-executable statement
200 format(1x, 50('*') //, &
    ' Algorithm: ASA-BCP', //, &
    ' Problem name: ' a10, /, &
    ' number of variables = ', i0, //, &
    ' f = ', d12.5, /, &
    ' sup-norm of the projected gradient = ', d12.5, /, &
    ' number of iterations = ', i0, /, &
    ' number of function evaluations = ', i0, /, &
    ' number of gradient evaluations = ', i0, /, &
    ' number of Hessian-vector products = ', i0, /, &
    ' number of inner cg iterations = ', i0, /, &
    ' exit flag = ', i1, /, &
    ' elapsed time (s) = ', e12.5, //, 1x, 50('*'))

end program

subroutine funct(status,n,x,f)
    implicit none
    
    integer, intent(in) :: n
    integer, intent(out) :: status
    double precision, intent(in) :: x(n)
    double precision, intent(out) :: f
    
    call CUTEST_ufn(status,n,x,f)
    return
end subroutine
    
subroutine grad(status,n,x,g)
    implicit none
    
    integer, intent(in) :: n
    integer, intent(out) :: status
    double precision, intent(in) :: x(n)
    double precision, intent(out) :: g(n)
    
    call CUTEST_ugr(status,n,x,g)
    return
end subroutine

subroutine hd_prod(status,n,goth,x,d,hd)
    implicit none
    
    integer, intent(in) :: n
    integer, intent(out) :: status
    double precision, intent(in) :: x(n),d(n)
    double precision, intent(out) :: hd(n)
    logical, intent(in) :: goth
    
    call CUTEST_uhprod(status,n,goth,x,d,hd)
    return
end subroutine    