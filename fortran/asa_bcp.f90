! -------------------------------------------------------------------------
!
! This file is part of ASA-BCP, which is a solver for bound-constrained
! optimization problems of the following form:
!
!                                min f(x)
!                           s.t. l <= x <= u
!
! where f(x) is a twice continuously differentiable.
!
! -------------------------------------------------------------------------
!
! Reference paper:
!
! A. Cristofari, M. De Santis, S. Lucidi, F. Rinaldi (2017). A Two-Stage
! Active-Set Algorithm for Bound-Constrained Optimization. Journal of
! Optimization Theory and Applications, 172(2), 369-401.
!
!-------------------------------------------------------------------------
!
! Authors:
! Andrea Cristofari (e-mail: andrea.cristofari@unipd.it)
! Marianna De Santis (e-mail: mdesantis@diag.uniroma1.it)
! Stefano Lucidi (e-mail: lucidi@diag.uniroma1.it)
! Francesco Rinaldi (e-mail: rinaldi@math.unipd.it)
!
! Last update of this file:
! January 31st, 2022
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


! modules
!--------------------------------------------

module asa_bcp_opts
    
    implicit none

    type :: asa_bcp_options

        ! *** DO NOT CHANGE DEFINITIONS IN THIS DERIVED DATA TYPE ***
        !
        ! to change ASA-BCP parameters, in the main program create an
        ! object of derived data type asa_bcp_options and assign new
        ! values to (some of) its members
        ! (see file 'asa_bcp_main.f90' for an example).

        ! ====================================
        ! DEFAULT VALUES OF ASA-BCP PARAMETERS
        ! ====================================

        ! PARAMETERS FOR TERMINATION (see the description of ASA_BCP above)
        double precision :: eps_opt = 1.d-5
        double precision :: min_gd = 1.d-15
        double precision :: min_norm_proj_d = 1.d-9
        double precision :: min_stepsize = 1.d-20
        double precision :: min_f = -1.d90
        integer :: max_it = 1000000
        integer :: max_n_f = 1000000
        integer :: max_n_g = 1000000
        integer :: max_n_hd = 1000000

        ! OTHER ALGORITHM PARAMETERS (see the description of ASA_BCP above)
        integer :: m = 100
        integer :: z = 20
        logical :: hd_exact = .true.
        integer :: verbosity = 1

    end type

end module

!-----------

module asa_bcp_active_set
    
    implicit none

    integer :: n_act_new,n_non_act,n_fix
    double precision, allocatable :: lambda(:),mu(:)
    logical, allocatable :: ind_act_l(:),ind_act_u(:),ind_non_act(:),ind_fix(:)

end module

!-----------

module asa_bcp_gradients
    
    implicit none

    double precision :: sup_norm_proj_g
    double precision, allocatable :: g(:)

end module

!-----------

module asa_bcp_dir
    
    implicit none

    double precision, allocatable :: d(:)

end module

!-----------

module asa_bcp_trunc_dir
    
    implicit none
    
    integer :: max_it_cg
    double precision :: gd,curv_g_non_act
    double precision :: eps_cg_dir,eps_cg_curv,eps_approx
    logical :: warn_grad,warn_small,warn_noposdef,warn_sing,warn_conjfail,warn_maxcgit,is_d_neg_curv

end module

!-----------

module asa_bcp_line_search
    
    implicit none
    
    double precision :: stepsize,f_newton_first,delta_f0,norm_proj_d
    double precision, allocatable :: v(:)

end module

!-----------

module asa_bcp_counters
    
    implicit none

    integer :: it,it_cg_tot,k,it_cg,n_f,n_g,n_hd

end module

!-----------

module asa_bcp_nm_ref
    
    implicit none

    integer :: mj
    double precision :: w_max
    double precision, allocatable :: w(:)

end module

!-----------

module asa_bcp_max_bounds
    
    implicit none

    double precision, parameter :: l_min = -1.d20 ! no i-th lower bound if l(i) <= l_min
    double precision, parameter :: u_max = 1.d20 ! no i-th upper bound if u(i) >= u_max

end module

!--------------------------------------------


! ASA-BCP solver

! Required input/output arguments:
!   n (in) is the problem dimension
!   x (in/out) is the starting point of the algorithm as input and the final solution as output
!   f (out) is the objective value at the final solution
!   l (in) is the lower bound on the variables
!   u (in) is the upper bound on the variables
!   opts (in) contains the algorithm parameters
!   flag (out) is the exit condition of the algorithm

subroutine asa_bcp(n,x,f,l,u,opts,flag)
    
    use asa_bcp_opts
    use asa_bcp_active_set
    use asa_bcp_gradients
    use asa_bcp_dir
    use asa_bcp_trunc_dir
    use asa_bcp_line_search
    use asa_bcp_counters
    use asa_bcp_nm_ref
    use asa_bcp_max_bounds
    
    implicit none
    
    integer, intent(in) :: n
    double precision, intent(in) :: l(n),u(n)
    type(asa_bcp_options), intent(in) :: opts
    double precision, intent(inout) :: x(n)
    integer, intent(out) :: flag
    double precision, intent(out) :: f
    
    integer :: max_it,max_n_f,max_n_g,max_n_hd,m,z,verbosity
    integer :: n_true,i,lj,z_nm,status
    double precision :: eps_opt,min_gd,min_norm_proj_d,min_stepsize,min_f
    double precision :: f_best,sup_norm_proj_g_best,fv,sq_d_act,proj_g,delta0_dir
    double precision :: delta_dir,delta_act,beta_dir,eps_act_set,x_best(n),g_best(n)
    logical :: hd_exact,act_phase,f_computed,is_first_linesearch,f_decreased,checkpoint,restart
    logical :: gd_exit,dir_exit,stepsize_exit
    
    eps_opt = opts%eps_opt
    min_gd = opts%min_gd
    min_norm_proj_d = opts%min_norm_proj_d
    min_stepsize = opts%min_stepsize
    min_f = opts%min_f
    max_it = opts%max_it
    max_n_f = opts%max_n_f
    max_n_g = opts%max_n_g
    max_n_hd = opts%max_n_hd
    m = opts%m
    z = opts%z
    hd_exact = opts%hd_exact
    verbosity = opts%verbosity

    if (eps_opt.lt.0.d0) then
        write(*,'(a)') "error when calling asa_bcp: 'eps_opt' must be a number greater than or equal to 0"
        stop
    endif
    if (min_gd.lt.0.d0) then
        write(*,'(a)') "error when calling asa_bcp: 'min_gd' must be a number greater than or equal to 0"
        stop
    endif
    if (min_norm_proj_d.lt.0.d0) then
        write(*,'(a)') "error when calling asa_bcp: 'min_norm_proj_d' must be a number greater than or equal to 0"
        stop
    endif
    if (min_stepsize.lt.0.d0) then
        write(*,'(a)') "error when calling asa_bcp: 'min_stepsize' must be a number greater than or equal to 0"
        stop
    endif
    if (max_it.lt.0) then
        write(*,'(a)') "error when calling asa_bcp: 'max_it' must be a number greater than or equal to 0"
        stop
    endif
    if (max_n_f.lt.1) then
        write(*,'(a)') "error when calling asa_bcp: 'max_n_f' must be a number greater than or equal to 1"
        stop
    endif
    if (max_n_g.lt.1) then
        write(*,'(a)') "error when calling asa_bcp: 'max_n_g' must be a number greater than or equal to 1"
        stop
    endif
    if (max_n_hd.lt.0) then
        write(*,'(a)') "error when calling asa_bcp: 'max_n_hd' must be a number greater than or equal to 0"
        stop
    endif
    if (m.lt.1) then
        write(*,'(a)') "error when calling asa_bcp: 'm' must be a number greater than or equal to 1"
        stop
    endif
    if (z.lt.1) then
        write(*,'(a)') "error when calling asa_bcp: 'z' must be a number greater than or equal to 1"
        stop
    endif
    if (verbosity.lt.0) then
        write(*,'(a)') "error when calling asa_bcp: 'verbosity' must be a number between 0 and 2"
        stop
    endif
    
    ! allocate vectors
    allocate(ind_act_l(n),ind_act_u(n),ind_non_act(n),ind_fix(n),lambda(n),mu(n),g(n),d(n),v(n),w(m))

    flag = -1
    
    lambda = 0.d0
    mu = 0.d0
    
    if (verbosity.gt.0) then
        open(10, file = 'iteration_history.txt')
    endif
    
    ! check problem feasibility, project the starting point onto the box,
    ! identify fixed variabled, initialize 'delta_act' and 'delta0_dir'
    n_fix = 0
    delta_act = 0.d0
    do i = 1,n
        if (l(i).lt.u(i)) then
            ind_fix(i) = .false.
            if (l(i) > l_min) then
                x(i) = dmax1(l(i),x(i))
                if (u(i) < u_max) then
                    x(i) = dmin1(x(i),u(i))
                    delta_act = delta_act + (u(i)-l(i))*(u(i)-l(i))
                endif
            endif
        elseif (l(i).gt.u(i)) then
            x = 0.d0
            f = 0.d0
            sup_norm_proj_g = 0.d0
            it = 0
            n_f = 0
            n_g = 0
            n_hd = 0
            it_cg_tot = 0
            write(*,*) 'infeasible problem'
            if (verbosity.gt.0) then
                write(*,*) 'infeasible problem'
            endif
            deallocate(ind_act_l,ind_act_u,ind_non_act,ind_fix,lambda,mu,g,d,v,w) ! free allocated memory
            return
        else
            x(i) = l(i)
            ind_fix(i) = .true.
            n_fix = n_fix + 1
        endif
    enddo
    
    ! the point obtained by setting the estimated active variables to the bounds
    ! is accepted in case of sufficient decrease in the objective function
    ! if the distance between the two points is less than or equal to 'delta_act'
    delta_act = dmin1(1.d30,dmax1(1.d3,dsqrt(delta_act)))
    
    ! the unit stepsize is accepted without evaluating the objective function
    ! if the distance between the two points is less than or equal to 'delta_dir's
    delta0_dir = 1.d3
    delta_dir = delta0_dir
    beta_dir = 9.d-1 ! reduction factor of 'delta_dir' (must be >=0 and <1)
    
    n_true = n - n_fix
    
    if (verbosity.gt.0) then
        write(*,'(a,i0,a,i0,a)') ' number of variables = ', n, ' (', n_fix, ' fixed)'
        write(10,'(a,i0,a,i0,a)') ' number of variables = ', n, ' (', n_fix, ' fixed)'
    endif
    
    ! first objective function evaluation
    call funct(status,n,x,f)
    if (status.ne.0) then
        f = huge(1.d0)
        sup_norm_proj_g = huge(1.d0)
        it = 0
        n_f = 0
        n_g = 0
        n_hd = 0
        it_cg_tot = 0
        flag = 9
        write(*,'(/,a)') ' error when computing the objective function, algorithm stopped'
        if (verbosity > 0) then
            write(10,'(/,a)') ' error when computing the objective function, algorithm stopped'
        endif
        goto 19
    endif
    n_f = 1
    f_computed = .true.
    
    x_best = x
    f_best = f
    v = x
    w = -huge(1.d0)
    w_max = f
    
    ! first gradient evaluation
    call grad(status,n,x,g)
    if (status.ne.0) then
        sup_norm_proj_g = huge(1.d0)
        n_g = 0
        it = 0
        n_hd = 0
        it_cg_tot = 0
        flag = 10
        write(*,'(/,a)') ' error when computing the gradient of objective function, algorithm stopped'
        if (verbosity > 0) then
            write(10,'(/,a)') ' error when computing the gradient of objective function, algorithm stopped'
        endif
        goto 19
    endif
    n_g = 1
        
    g_best = g
    
    gd = -huge(1.d0)
    
    n_hd = 0
    
    ! compute the sup-norm of the projected gradient and initalize 'eps_act_set'
    sup_norm_proj_g = 0.d0
    eps_act_set = 0.d0
    do i = 1,n
        if (.not.ind_fix(i)) then
            if (dmax1(l_min,x(i)-g(i)).lt.l(i)) then
                proj_g = x(i) - l(i)
            elseif (dmin1(u_max,x(i)-g(i)).gt.u(i)) then
                proj_g = u(i) - x(i)
            else
                proj_g = g(i)
            endif
            eps_act_set = eps_act_set + proj_g*proj_g
            sup_norm_proj_g = dmax1(sup_norm_proj_g,dabs(proj_g))
        endif
    enddo
    ! at this point, 'eps_act_set' is the squared norm of the projected gradient
    if (eps_act_set.gt.0.d0) then
        eps_act_set = dmin1(1.d-6,eps_act_set**(-3.d0/2.d0))
    else
        eps_act_set = 1.d-6
    endif
    
    sup_norm_proj_g_best = sup_norm_proj_g
    
    z_nm = min0(n_true,z)
    
    ! parameters for the truncated-Newton method
    eps_cg_dir = dmin1(1.d0,min_norm_proj_d**(3.d0/2.d0)) ! 'eps_cg_dir' should be <= 'min_norm_proj_d',
                                                          ! so that the norm of the first conjugate direction
                                                          ! in subroutine 'trunc_dir' is sufficiently large
    eps_cg_curv = 1.d-9
    eps_approx = 1.d-6
    
    ! initialize counters
    lj = -1
    mj = -1
    k = -1
    it = -1
    it_cg_tot = 0
    it_cg = 0
    
    act_phase = .false.
    gd_exit = .false.
    dir_exit = .false.
    stepsize_exit = .false.
    checkpoint = .true.
    restart = .false.
    is_first_linesearch = .true.
    
    if (f.le.min_f) then
        it = 0
        n_hd = 0
        it_cg_tot = 0
        if (verbosity.gt.0) then
            write(*,'(/,a)') ' objective value below the minimum threshold, algorithm stopped'
            write(10,'(/,a)') ' objective value below the minimum threshold, algorithm stopped'
        endif
        flag = 8
        goto 19
    endif
    
    !-------------------
    !  START MAIN LOOP
    !-------------------
    
    do
        
        k = k + 1
        it = it + 1
        
        if ((.not.act_phase).and.(it.gt.0)) then
            
            ! compute the gradient and the sup-norm of the projected gradient
            call grad(status,n,x,g)
            if (status.ne.0) then
                flag = 10
                goto 22
            endif
            n_g = n_g + 1
            sup_norm_proj_g = 0.d0
            do i = 1,n                    
                if (.not.ind_fix(i)) then
                    if (dmax1(l_min,x(i)-g(i)).lt.l(i)) then
                        sup_norm_proj_g = dmax1(sup_norm_proj_g,dabs(x(i)-l(i)))
                    elseif (dmin1(u_max,x(i)-g(i)).gt.u(i)) then
                        sup_norm_proj_g = dmax1(sup_norm_proj_g,dabs(u(i)-x(i)))
                    else
                        sup_norm_proj_g = dmax1(sup_norm_proj_g,dabs(g(i)))
                    endif
                endif
            enddo
            
        endif
        
        ! main prints
        if (verbosity.gt.0) then
            write(*,31) it,dmin1(f_best,f),sup_norm_proj_g,n_f,n_g,n_hd,it_cg_tot
            write(10,31) it,dmin1(f_best,f),sup_norm_proj_g,n_f,n_g,n_hd,it_cg_tot
        endif
        
        if ((sup_norm_proj_g.le.eps_opt).or.(it.ge.max_it).or.(n_f.ge.max_n_f).or.(n_g.ge.max_n_g)) then
            goto 20
        endif
        
        !---------------------------------------------------------------
        !     MINIMIZATION STEP OVER THE ESTIMATED ACTIVE VARIABLES
        !---------------------------------------------------------------
        
        if (restart) then
            
            ! active-set estimate
            act_phase = .false.
            call multipliers(n,x,l,u)
            call estimate_active_set(n,x,l,u,eps_act_set,act_phase)
            
            restart = .false.
            
            if (verbosity.gt.1) then
                write(*,'(/,a)') ' --- iteration details ---'
                write(10,'(/,a)') ' --- iteration details ---'
            endif
            
            goto 200
            
        endif
        
        if (.not.act_phase) then
            ! active-set estimate
            act_phase = .true.
            call multipliers(n,x,l,u)
            call estimate_active_set(n,x,l,u,eps_act_set,act_phase)
        endif
        
        if (verbosity.gt.1) then
            write(*,'(/a,//,a,i0)') ' --- iteration details ---', &
                                    ' number of estimated active variables not at the bounds =  ', n_act_new
            write(10,'(/a,//,a,i0)') ' --- iteration details ---', &
                                     ' number of estimated active variables not at the bounds =  ', n_act_new
        endif
        
        if (act_phase) then
            
            v = x
            f_decreased = .false.
            
            if (.not.f_computed) then
                call funct(status,n,v,f)
                if (status.ne.0) then
                    flag = 9
                    goto 22
                endif
                n_f = n_f + 1
                f_computed = .true.
                if (verbosity.gt.1) then
                    write(*,'(a,d13.6,a)') ' function control, f = ', f, ')'
                    write(10,'(a,d13.6,a)') ' function control, f = ', f, ')'
                endif
            endif
            
            do while (act_phase.and.(.not.f_decreased))
                
                ! set the estimated active variables to the bounds
                sq_d_act = 0.d0
                do i = 1,n
                    if ( (ind_act_l(i)).and.(x(i).gt.l(i)) ) then
                        sq_d_act = sq_d_act + (x(i)-l(i))*(x(i)-l(i))
                        x(i) = l(i)
                    elseif ( (ind_act_u(i)).and.(x(i).lt.u(i)) ) then
                        sq_d_act = sq_d_act + (x(i)-u(i))*(x(i)-u(i))
                        x(i) = u(i)
                    endif
                enddo
                
                call funct(status,n,x,fv)
                if (status.ne.0) then
                    flag = 9
                    goto 22
                endif
                n_f = n_f + 1
                
                ! check if the objective function is sufficiently decreased
                if (fv.le.f-(sq_d_act/(2d0*eps_act_set))) then
                    
                    f_decreased = .true.
                    
                    if (dsqrt(sq_d_act).le.delta_act) then
                        
                        delta_act = beta_dir*delta_act
                        
                    elseif (fv.gt.w_max) then ! point not accepted
                        
                        if (n_f.ge.max_n_f) then
                            goto 20
                        endif
                        
                        ! restart
                        if (verbosity.gt.1) then
                            write(*,'(a,d13.6,a,//,a)') ' point not accepted (f = ', fv, ')', &
                                                        ' restart from the best point'
                            write(10,'(a,d13.6,a,//,a)') ' point not accepted (f = ', fv, ')', &
                                                         ' restart from the best point'
                        endif
                        x = x_best
                        f = f_best
                        g = g_best
                        sup_norm_proj_g_best = sup_norm_proj_g
                        delta0_dir = 1.d-1*delta0_dir
                        w = -huge(1.d0)
                        m = m/5 + 1
                        k = -1
                        lj = 0
                        mj = 0
                        z_nm = min0(n_true,z)
                        w(m) = f
                        w_max = f
                        is_first_linesearch = .true.
                        restart = .true.
                        checkpoint = .false.
                        act_phase = .true.
                        
                        cycle
                        
                    endif   
                    
                    ! the new point has been accepted
                    f = fv
                    
                    ! compute the gradient and the sup-norm of the projected gradient
                    call grad(status,n,x,g)
                    if (status.ne.0) then
                        flag = 10
                        goto 22
                    endif
                    n_g = n_g + 1
                    sup_norm_proj_g = 0.d0
                    do i = 1,n                    
                        if (.not.ind_fix(i)) then
                            if (dmax1(l_min,x(i)-g(i)).lt.l(i)) then
                                sup_norm_proj_g = dmax1(sup_norm_proj_g,dabs(x(i)-l(i)))
                            elseif (dmin1(u_max,x(i)-g(i)).gt.u(i)) then
                                sup_norm_proj_g = dmax1(sup_norm_proj_g,dabs(x(i)-u(i)))
                            else
                                sup_norm_proj_g = dmax1(sup_norm_proj_g,dabs(g(i)))
                            endif
                        endif
                    enddo
                    
                    if (f.lt.f_best) then
                        if (f.le.min_f) then
                            goto 21
                        endif
                        f_best = f
                        x_best = x
                        g_best = g
                        sup_norm_proj_g_best = sup_norm_proj_g
                    endif
                    
                    if (verbosity.gt.1) then
                        write(*,'(a,d13.6,a)') ' point accepted (f = ', f, ')'
                        write(10,'(a,d13.6,a)') ' point accepted (f = ', f, ')'
                    endif
                    
                    if ((sup_norm_proj_g.le.eps_opt).or.(n_f.ge.max_n_f).or.(n_g.ge.max_n_g)) then
                        goto 20
                    endif
                    
                else ! set x to the previous value and reduce 'eps_act_set'
                     ! to carry out a new active-set estimate
                    
                    x = v
                    
                    if (n_f.ge.max_n_f) then
                        if (verbosity.gt.1) then
                            write(*,'(a,d13.6,a,/,a,/,a,i0)') ' point not accepted (f = ', fv, ')'
                            write(10,'(a,d13.6,a,/,a,/,a,i0)') ' point not accepted (f = ', fv, ')'
                            goto 20
                        endif
                    endif
                    
                    eps_act_set = eps_act_set*1.d-1
                    
                    call estimate_active_set(n,x,l,u,eps_act_set,act_phase)
                    
                    if (verbosity.gt.1) then
                        write(*,'(a,d13.6,a,/,a,/,a,i0)') &
                              ' point not accepted (f = ', fv, ')', &
                              ' reducing epsilon', &
                              ' number of estimated active variables not at the bounds =  ', n_act_new
                        write(10,'(a,d13.6,a,/,a,/,a,i0)') &
                              ' point not accepted (f = ', fv, ')', &
                              ' reducing epsilon', &
                              ' number of estimated active variables not at the bounds =  ', n_act_new
                    endif
                    
                endif
                
            enddo
            
        endif
        
        !---------------------------------------------------------------
        !   MINIMIZATION STEP OVER THE ESTIMATED NON-ACTIVE VARIABLES
        !---------------------------------------------------------------
        
        if (act_phase) then
            ! active-set estimate
            act_phase = .false.
            call multipliers(n,x,l,u)
            call estimate_active_set(n,x,l,u,eps_act_set,act_phase)
        endif
        
200     if (verbosity.gt.1) then
            write(*,'(/,a,i0)') ' number of estimated non-active variables = ', n_non_act
            write(10,'(/,a,i0)') ' number of estimated non-active variables = ', n_non_act
        endif
        
        if (.not.act_phase) then
            
            if (checkpoint) then
                if (f.lt.f_best) then
                    if (f.le.min_f) then
                        goto 21
                    endif
                    x_best = x
                    f_best = f
                    g_best = g
                    sup_norm_proj_g_best = sup_norm_proj_g
                endif
                call update_w(m,f)
                lj = k
                checkpoint = .false.
            endif
            
            ! function control
            if (k.ge.(lj+z_nm)) then
                
                z_nm = min0(z_nm+n_true,z)
                
                if (.not.f_computed) then
                    call funct(status,n,x,f)
                    if (status.ne.0) then
                        flag = 9
                        goto 22
                    endif
                    n_f = n_f + 1
                    f_computed = .true.
                endif
                
                if (f.gt.w_max) then
                    
                    if (n_f.ge.max_n_f) then
                        goto 20
                    endif
                    
                    ! restart
                    if (verbosity.gt.1) then
                        write(*,'(a,d13.6,a,//,a)') ' function control not satisfied (f = ', f, ')', &
                                                    ' restart from the best point'
                        write(10,'(a,d13.6,a,//,a)') ' function control not satisfied (f = ', f, ')', &
                                                     ' restart from the best point'
                    endif
                    x = x_best
                    f = f_best
                    g = g_best
                    sup_norm_proj_g = sup_norm_proj_g_best
                    delta0_dir = 1.d-1*delta0_dir
                    w = -huge(1.d0)
                    m = m/5 + 1
                    k = -1
                    lj = 0
                    mj = 0
                    z_nm = min0(n_true,z)
                    w(m) = f
                    w_max = f
                    is_first_linesearch = .true.
                    restart = .true.
                    checkpoint = .false.
                    act_phase = .true.
                    
                    cycle
                    
                else
                    
                    if (verbosity.gt.1) then
                        write(*,'(a,d13.6,a)') ' function control satisfied (f = ', f, ')'
                        write(10,'(a,d13.6,a)') ' function control satisfied (f = ', f, ')'
                    endif
                    if (f.lt.f_best) then
                        if (f.le.min_f) then
                            goto 21
                        endif
                        x_best = x
                        f_best = f
                        g_best = g
                        sup_norm_proj_g_best = sup_norm_proj_g
                    endif
                    call update_w(m,f)
                    lj = k
                    
                endif
                
            endif
            
            ! compute the search direction
            if (hd_exact) then
                max_it_cg = min0(2*n_non_act,max_n_hd-it_cg_tot)
            else
                max_it_cg = min0(2*n_non_act,max_n_g-it_cg_tot)
            endif
            if (max_it_cg.lt.1) then
                goto 20
            endif
            call dir(n_non_act,n,ind_non_act,x,hd_exact,flag)
            it_cg_tot = it_cg_tot + it_cg
            if ((flag.eq.10).or.(flag.eq.11)) then
                goto 22
            endif
            
            if (verbosity.gt.1) then
                write(*,'(a,i0)') ' number of inner cg iterations = ', it_cg
                write(10,'(a,i0)') ' number of inner cg iterations = ', it_cg
				if (warn_sing) then
                    write(*,*) 'Hessian matrix probably singular'
                    write(10,*) 'Hessian matrix probably singular'
				elseif (warn_noposdef) then
                    write(*,*) 'Hessian matrix probably not positive definite'
                    write(10,*) 'Hessian matrix probably not positive definite'
                elseif (warn_small) then
                    write(*,*) 'norm of the inner conjugate direction too small'
                    write(10,*) 'norm of the inner conjugate direction too small'
                elseif (warn_conjfail) then
                    write(*,*) 'conjugacy failure when computing the search direction'
                    write(10,*) 'conjugacy failure when computing the search direction'
                elseif (warn_maxcgit) then
                    write(*,*) 'conjugate gradient method not converged'
                    write(10,*) 'conjugate gradient method not converged'
                endif
				if (warn_grad) then
					write(*,*) 'anti-gradient used as search direction'
                    write(10,*) 'anti-gradient used as search direction'
                endif
                if (is_d_neg_curv) then
                    write(*,*) 'negative curvature direction'
                    write(10,*) 'negative curvature direction'
                endif
                write(*,'(a,d13.6)') ' directional derivative = ', gd
                write(10,'(a,d13.6)') ' directional derivative = ', gd
            endif
            
            if (gd.ge.-min_gd) then
                gd_exit = .true.
                goto 20
            endif
            
            ! set v = p(x + d)
            do i = 1,n                    
                if (.not.ind_fix(i)) then
                    v(i) = x(i) + d(i)
                    if (dmax1(l_min,v(i)).lt.l(i)) then
                        v(i) = l(i)
                    elseif (dmin1(u_max,v(i)).gt.u(i)) then
                        v(i) = u(i)
                    endif
                endif
            enddo
            
            norm_proj_d = dsqrt(dot_product(x-v,x-v))
            
            if (verbosity.gt.1) then 
                write(*,'(a,d13.6)') ' norm of the projected direction = ', norm_proj_d
                write(10,'(a,d13.6)') ' norm of the projected direction = ', norm_proj_d
            endif
            
            ! check if the norm of the projected direction is sufficiently large
            if (norm_proj_d.le.min_norm_proj_d) then
                dir_exit = .true.
                goto 20
            endif    
            
            ! try accepting the unit stepsize
            if ((norm_proj_d.le.delta_dir).and.(.not.(warn_noposdef.or.is_first_linesearch))) then
                
                if (verbosity.gt.1) then
                    write(*,*) 'stepsize = 1 accepted without computing f'
                    write(10,*) 'stepsize = 1 accepted without computing f'
                endif
                
                delta_dir = beta_dir*delta_dir
                
                x = v
                f_computed = .false.
                
                cycle
                
            endif
            
            if (n_f.ge.max_n_f) then
                goto 20
            endif
            
            ! unit stepsize not accepted -> line search
            if (k.ne.lj) then
                
                if (.not.f_computed) then
                    call funct(status,n,x,f)
                    if (status.ne.0) then
                        flag = 9
                        goto 22
                    endif
                    n_f = n_f + 1
                    f_computed = .true.
                endif
                
                ! check the objective function
                if (f.gt.w_max) then
                    
                    if (n_f.ge.max_n_f) then
                        goto 20
                    endif
                    
                    ! restart
                    if (verbosity.gt.1) then
                        write(*,'(a,d13.6,a,//,a)') ' function control not satisfied (f = ', f, ')', &
                                                    ' restart from the best point'
                        write(10,'(a,d13.6,a,//,a)') ' function control not satisfied (f = ', f, ')', &
                                                     ' restart from the best point'
                    endif
                    x = x_best
                    f = f_best
                    g = g_best
                    sup_norm_proj_g = sup_norm_proj_g_best
                    delta0_dir = 1.d-1*delta0_dir
                    w = -huge(1.d0)
                    m = m/5 + 1
                    k = -1
                    lj = 0
                    mj = 0
                    z_nm = min0(n_true,z)
                    w(m) = f
                    w_max = f
                    is_first_linesearch = .true.
                    restart = .true.
                    checkpoint = .false.
                    act_phase = .true.
                    
                    cycle
                    
                else
                    
                    if (verbosity.gt.1) then
                        write(*,'(a,d13.6,a)') ' function control satisfied (f = ', f, ')'
                        write(10,'(a,d13.6,a)') ' function control satisfied (f = ', f, ')'
                    endif
                    if (f.lt.f_best) then
                        if (f.le.min_f) then
                            goto 21
                        endif
                        x_best = x
                        f_best = f
                        g_best = g
                        sup_norm_proj_g_best = sup_norm_proj_g
                    endif
                    call update_w(m,f)
                    lj = k
                    
                endif
                
            endif
            
            ! check if this is the first line search
            if (is_first_linesearch) then
                f_newton_first = f
                delta_f0 = 0.d0
            endif
            
            ! line search
            if (verbosity.gt.1) then
                write(*,*) 'line search'
                write(10,*) 'line search'
            endif
            call linesearch(n_non_act,n,ind_non_act,x,l,u,f,w_max,n_f,min_stepsize,max_n_f,stepsize_exit,flag)
            
            if (flag.eq.9) then
                goto 22
            endif
            if (flag.eq.5) then
                goto 20
            endif
            
            ! check if the line search failed
            if (stepsize_exit) then
                goto 20
            endif
            
            if (verbosity.gt.1) then
                write(*,'(a,d13.6)') ' stepsize = ', stepsize
                write(10,'(a,d13.6)') ' stepsize = ', stepsize
            endif
            
            checkpoint = .true.
            
            if (is_first_linesearch) then
                delta_f0 = f_newton_first - f
                delta_dir = delta0_dir*stepsize*norm_proj_d
                is_first_linesearch = .false.   
            endif
            
            cycle
            
        else
            
            sup_norm_proj_g = huge(1.d0)
            gd = -huge(1.d0)
            norm_proj_d = huge(1.d0)
            lj = lj + 1
            cycle
            
        endif
        
22      continue ! error when computing f(x), or g(x), or H(x)*d
        x = x_best
        f = f_best
        sup_norm_proj_g = sup_norm_proj_g_best
        if (verbosity.gt.0) then
            if (flag.eq.9) then
                write(*,'(/,a)') ' error when computing the objective function, algorithm stopped'
                if (verbosity > 0) then
                    write(10,'(/,a)') ' error when computing the objective function, algorithm stopped'
                endif
            elseif (flag.eq.10) then
                write(*,'(/,a)') ' error when computing the gradient of objective function, algorithm stopped'
                if (verbosity > 0) then
                    write(10,'(/,a)') ' error when computing the gradient of objective function, algorithm stopped'
                endif
            else ! fleg.eq.11
                write(*,'(/,a)') ' error when computing the Hessian-vector product, algorithm stopped'
                if (verbosity > 0) then
                    write(10,'(/,a)') ' error when computing the Hessian-vector product, algorithm stopped'
                endif
            endif
        endif
        goto 19
        
21      continue
        it = it + 1
        flag = 8
        if (verbosity.gt.0) then
            write(*,'(/,a)') ' objective value below the minimum threshold, algorithm stopped'
            write(10,'(/,a)') ' objective value below the minimum threshold, algorithm stopped'
        endif
        goto 19
        
20      continue
        if (.not.f_computed) then
            call funct(status,n,x,f)
            if (status.ne.0) then
                flag = 9
                goto 22
            endif
            n_f = n_f + 1
            f_computed = .true.
            if (f.le.min_f) then
                goto 21
            endif
        endif
                
        if (sup_norm_proj_g.le.eps_opt) then
            if (verbosity.gt.0) then
                write(*,'(/,1x,83("="),/,a,d13.6,/,1x,83("="))') &
                        ' optimality condition satisfied: sup-norm of the projected gradient <= ', eps_opt
                write(10,'(/,1x,83("="),/,a,d13.6,/,1x,83("="))') &
                            ' optimality condition satisfied: sup-norm of the projected gradient <= ', eps_opt
            endif
            if (f.le.f_best) then
                flag = 0
            endif
        elseif (gd_exit) then
            ! active-set estimate
            act_phase = .true.
            call multipliers(n,x,l,u)
            call estimate_active_set(n,x,l,u,eps_act_set,act_phase)
            if (act_phase) then
                gd_exit = .false.
                flag = -2
            elseif (f.le.f_best) then
                flag = 1
                if (verbosity.gt.0) then
                    write(*,'(/,a)') 'directional derivative not sufficiently negative, algorithm stopped'
                    write(10,'(/,a)') 'directional derivative not sufficiently negative, algorithm stopped'
                endif
            elseif (verbosity.gt.0) then
                write(*,'(/,a)') 'directional derivative not sufficiently negative'
                write(10,'(/,a)') 'directional derivative not sufficiently negative'
            endif
        elseif (dir_exit) then
            ! active-set estimate
            act_phase = .true.
            call multipliers(n,x,l,u)
            call estimate_active_set(n,x,l,u,eps_act_set,act_phase)
            if (act_phase) then
                dir_exit = .false.
                flag = -2
            elseif (f.le.f_best) then
                flag = 2
                if (verbosity.gt.0) then
                    write(*,'(/,a)') ' norm of the projected direction too small, algorithm stopped'
                    write(10,'(/,a)') ' norm of the projected direction too small, algorithm stopped'
                endif
            elseif (verbosity.gt.0) then
                write(*,'(/,a)') ' norm of the projected direction too small'
                write(10,'(/,a)') ' norm of the projected direction too small'
            endif
        elseif (stepsize_exit) then
            ! active-set estimate
            act_phase = .true.
            call multipliers(n,x,l,u)
            call estimate_active_set(n,x,l,u,eps_act_set,act_phase)
            if (act_phase) then
                stepsize_exit = .false.
                flag = -2
            elseif (f.le.f_best) then
                flag = 3
                if (verbosity.gt.0) then
                    write(*,'(/,a)') ' stepsize too small, algorithm stopped'
                    write(10,'(/,a)') ' stepsize too small, algorithm stopped'
                endif
            elseif (verbosity.gt.0) then
                write(*,'(/,a)') ' stepsize too small'
                write(10,'(/,a)') ' stepsize too small'
            endif
        else
            if (f.gt.f_best) then
                x = x_best
                f = f_best
                sup_norm_proj_g = sup_norm_proj_g_best
            endif
            if (it.ge.max_it) then
                flag = 4
                if (verbosity.gt.0) then
                    write(*,'(/,a)') ' too many iterations, algorithm stopped'
                    write(10,'(/,a)') ' too many iterations, algorithm stopped'
                endif
            elseif (it.ge.max_n_f) then
                flag = 5
                if (verbosity.gt.0) then
                    write(*,'(/,a)') ' too function evaluations, algorithm stopped'
                    write(10,'(/,a)') ' too function evaluations, algorithm stopped'
                endif
            elseif (it.ge.max_n_g) then
                flag = 6
                if (verbosity.gt.0) then
                    write(*,'(/,a)') ' too gradient evaluations, algorithm stopped'
                    write(10,'(/,a)') ' too gradient evaluations, algorithm stopped'
                endif
            else ! n_hd >= max_n_hd
                flag = 7
                if (verbosity.gt.0) then
                    write(*,'(/,a)') ' too many Hessian-vector products, algorithm stopped'
                    write(10,'(/,a)') ' too many Hessian-vector products, algorithm stopped'
                endif
            endif
        endif
        
        if (flag.ge.0) then
            goto 19
        endif
        
        if ((it.lt.max_it).and.(n_f.lt.max_n_f).and.(n_g.lt.max_n_g).and.(n_hd.lt.max_n_hd)) then
        
            if (flag.eq.-1) then
            
                ! restart
                if (verbosity.gt.0) then
                    write(*,'(//,a)') ' restart from the best point'
                    write(10,'(//,a)') ' restart from the best point'
                endif
                x = x_best
                f = f_best
                g = g_best
                sup_norm_proj_g = sup_norm_proj_g_best
                delta0_dir = 1.d-1*delta0_dir
                w = -huge(1.d0)
                m = m/5 + 1
                k = -1
                lj = 0
                mj = 0
                z_nm = min0(n_true,z)
                w(m) = f
                w_max = f
                is_first_linesearch = .true.
                restart = .true.
                checkpoint = .false.
                act_phase = .true.
            
            endif
            
        else

            x = x_best
            f = f_best
            sup_norm_proj_g = sup_norm_proj_g_best
            
            if (it.ge.max_it) then
            
                flag = 4
                if (verbosity.gt.0) then
                    write(*,'(/,a)') 'restart not possible (too many iterations), algorithm stopped'
                    write(10,'(/,a)') 'restart not possible (too many iterations), algorithm stopped'
                endif
            
            elseif (n_f.ge.max_n_f) then
            
                flag = 5
                if (verbosity.gt.0) then
                    write(*,'(/,a)') 'restart not possible (too many function evaluations), algorithm stopped'
                    write(10,'(/,a)') 'restart not possible (too many function evaluations), algorithm stopped'
                endif
                
            elseif (n_g.ge.max_n_g) then
            
                flag = 6
                if (verbosity.gt.0) then
                    write(*,'(/,a)') 'restart not possible (too many gradient evaluations), algorithm stopped'
                    write(10,'(/,a)') 'restart not possible (too many gradient evaluations), algorithm stopped'
                endif
            
            else ! n_hd >= max_n_hd
            
                flag = 7
                if (verbosity.gt.0) then
                    write(*,'(/,a)') 'restart not possible (too many Hessian-vector products), algorithm stopped'
                    write(10,'(/,a)') 'restart not possible (too many Hessian-vector products), algorithm stopped'
                endif
                
            endif
            
        endif
        
    enddo
    
19 continue
    
    if (verbosity.gt.0) then
        write(*,'(/,a,/)') ' WARNING: using ''verbosity = 0'' may be faster'
        write(10,'(/,a)') ' WARNING: using ''verbosity = 0'' may be faster'
        close(10)
    endif
    deallocate(ind_act_l,ind_act_u,ind_non_act,ind_fix,lambda,mu,g,d,v,w) ! free allocated memory
    
! non-executable statements
31 format(/, 1x, 56('-'), //, &
          ' iteration ', i0, //, &
          ' best f = ', d13.6, /, &
          ' sup-norm of the projected gradient at the current point = ', d13.6, /, &
          ' number of function evaluations = ', i0, /, &
          ' number of gradient evaluations = ', i0, /, &
          ' number of Hessian-vector products = ', i0, /, &
          ' number of inner cg iterations = ', i0)

end subroutine
!-------------------------------------------------------------------------------------




! other subroutines
!-------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------
subroutine multipliers(n,x,l,u)
    
    use asa_bcp_active_set, only : lambda,mu,ind_fix
    use asa_bcp_gradients, only : g
    use asa_bcp_max_bounds
    
    implicit none
    
    integer, intent(in) :: n
    double precision, intent(in) :: x(n),l(n),u(n)
    
    integer :: i
    double precision :: mult1,mult2
    
    do i = 1,n
        if (.not.ind_fix(i)) then
            if (l(i).gt.l_min) then
                mult1 = (u(i)-x(i))*(u(i)-x(i))
                mult2 = (x(i)-l(i))*(x(i)-l(i))
                lambda(i) = g(i)*mult1/(mult1+mult2)
                if (u(i).lt.u_max) then
                    mu(i) = -g(i)*mult2/(mult1+mult2)
                endif
            elseif (u(i).lt.u_max) then
                mu(i) = -g(i)*(x(i)-l(i))*(x(i)-l(i))/((u(i)-x(i))*(u(i)-x(i))+(x(i)-l(i))*(x(i)-l(i)))
            endif
        endif
    enddo
    
    return

end subroutine
!-------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------
subroutine estimate_active_set(n,x,l,u,eps_act_set,act_phase)
    
    use asa_bcp_active_set
    use asa_bcp_gradients, only : g
    use asa_bcp_max_bounds
    
    implicit none
    
    integer, intent(in) :: n
    double precision, intent(in) :: x(n),l(n),u(n),eps_act_set
    logical, intent(inout) :: act_phase
    
    integer :: i,count
    logical :: est
    
    ind_act_l = .false.
    ind_act_u = .false.
    ind_non_act = .false.
    
    n_act_new = 0
    n_non_act = 0
    
    est = .true.
    count = 0
    
    do while (est.and.(count<2))
        if (act_phase) then
            do i = 1,n
                if (.not.ind_fix(i)) then
                    if ( (x(i)-l(i).le.eps_act_set*lambda(i)).and.(g(i).gt.0.d0).and. &
                         (x(i).gt.l(i)).and.(l(i).gt.l_min) ) then
                        ind_act_l(i) = .true.
                        n_act_new = n_act_new + 1
                    elseif ( (u(i)-x(i).le.eps_act_set*mu(i)).and.(g(i).lt.0.d0).and. &
                             (x(i).lt.u(i)).and.(u(i).lt.u_max) ) then
                        ind_act_u(i) = .true.
                        n_act_new = n_act_new + 1
                    endif
                endif
            enddo
            act_phase = (n_act_new>0)
            est = .not.act_phase
            count = count + 1
        else
            do i = 1,n
                if (.not.ind_fix(i)) then
                    if ( .not.((x(i).le.l(i)+eps_act_set*lambda(i)).and.(g(i).gt.0.d0).and.(l(i).gt.l_min)).and. &
                         .not.((x(i).ge.u(i)-eps_act_set*mu(i)).and.(g(i).lt.0.d0).and.(u(i).lt.u_max)) ) then
                        ind_non_act(i) = .true.
                        n_non_act = n_non_act + 1
                    endif
                endif
            enddo
            act_phase = (n_non_act==0)
            est = act_phase
            count = count + 1
        endif
    enddo
    
    return

end subroutine
!-------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------
subroutine dir(n,n_orig,ind_non_act,x,hd_exact,flag)
    
    use asa_bcp_gradients, only : g
    use asa_bcp_dir
    
    implicit none
    
    integer, intent(in) :: n,n_orig
    double precision, intent(in) :: x(n_orig)
    logical, intent(in) :: ind_non_act(n_orig),hd_exact
    integer, intent(inout) :: flag
    
    integer :: i,counter
    double precision :: d_non_act(n)
    
    call trunc_newton_dir(n,n_orig,ind_non_act,x,d_non_act,hd_exact,flag)
    
    i = 1
    counter = 1
    do while (i.le.n_orig)
        if (ind_non_act(i)) then
            d(i) = d_non_act(counter)
            counter = counter + 1
        else
            d(i) = 0.d0
        endif
        i = i + 1
    enddo
    
    return

end subroutine
!-------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------
subroutine trunc_newton_dir(n,n_orig,ind_non_act,x,d,hd_exact,flag)
    
    use asa_bcp_counters, only : k,it_cg,n_g,n_hd
    use asa_bcp_gradients, only : g
    use asa_bcp_trunc_dir
    
    implicit none
    
    integer, intent(in) :: n,n_orig
    double precision, intent(in) :: x(n_orig)
    logical, intent(in) :: ind_non_act(n_orig),hd_exact
    integer, intent(inout) :: flag
    double precision, intent(out) :: d(n)
    
    integer :: i,counter,status
    double precision :: curv,alpha,beta,norm_g_non_act,norm_d,eps_cg_tr,gp,gqp,norm_p,fquad,fquad_old,gd_old
    double precision :: g_non_act(n),hp(n),hd(n),gq(n),p(n),p_normalized(n),v1(n_orig),v2(n_orig)
    logical :: goth
    
    ! in this subroutine, 'p' are the conjugate directions and 'd' is the resulting search direction
    
    warn_grad = .false.
    warn_small = .false.
    warn_noposdef = .false.
	warn_sing = .false.
    warn_conjfail = .false.
    warn_maxcgit = .false.
    is_d_neg_curv = .false.
    
    i = 1
    counter = 1
    do while (counter.le.n)
        if (ind_non_act(i)) then
            g_non_act(counter) = g(i)
            counter = counter + 1
        endif
        i = i + 1
    enddo
    
    if (hd_exact) then
        v1 = 0.d0
        goth = .false.
    else
        v1 = x
    endif
    
    d = 0.d0
    p = -g_non_act
    gq = g_non_act
    hd = 0.d0
    
    fquad = 0.d0
    gd = 0.d0
    
    norm_g_non_act = dsqrt(dot_product(g_non_act,g_non_act))
    norm_p = norm_g_non_act
    
    it_cg = 1
    do
        
        ! normalize 'p'
        p_normalized = p/norm_p
        
        ! compute H(x)*p_normalized in the subspace of the estimated non-active variables
        !---------------------------------------------
        if (hd_exact) then
            counter = 1
            i = 1
            do while (counter.le.n)
                if (ind_non_act(i)) then
                    v1(i) = p_normalized(counter)
                    counter = counter + 1
                endif
                i = i + 1
            enddo
            call hd_prod(status,n_orig,goth,x,v1,v2)
            if (status.ne.0) then
                flag = 11
                exit
            endif
            n_hd = n_hd + 1
            goth = .true.
            counter = 1
            i = 1
            do while (counter.le.n)
                if (ind_non_act(i)) then
                    hp(counter) = v2(i)
                    counter = counter + 1
                endif
                i = i + 1
            enddo
        else
            counter = 1
            i = 1
            do while (counter.le.n)
                if (ind_non_act(i)) then
                    v1(i) = x(i) + eps_approx*p_normalized(counter)
                    counter = counter + 1
                endif
                i = i + 1
            enddo
            call grad(status,n_orig,v1,v2)
            if (status.ne.0) then
                flag = 10
                exit
            endif
            n_g = n_g + 1
            hp = 0.d0
            counter = 1
            i = 1
            do while (counter.le.n)
                if (ind_non_act(i)) then
                    hp(counter) = (v2(i)-g(i))/eps_approx
                    counter = counter + 1
                endif
                i = i + 1
            enddo
        endif
        !---------------------------------------------
        
        ! compute curvature
        curv = dot_product(p_normalized,hp)
        if (it_cg.eq.1) then
            curv_g_non_act = curv*norm_p*norm_p
        endif
        
        ! check curvature
        !---------------------------------------------
        if (curv.le.eps_cg_curv) then ! Hessian matrix not (sufficiently) positive definite
			if (curv.lt.-eps_cg_curv) then
                if (it_cg.gt.1) then
                    alpha = dot_product(-gq,p_normalized)/curv
                    d = d - alpha*p_normalized
                    gd = dot_product(g_non_act,d)
                    hd = hd - alpha*hp
                    curv = dot_product(d,hd)
                    norm_d = dot_product(d,d) ! squared norm of 'd'
                    if (curv.lt.eps_cg_curv*norm_d) then
                        is_d_neg_curv = (curv.lt.-eps_cg_curv*norm_d)
                        warn_sing = .not.is_d_neg_curv
                    endif
                else
                    d = g_non_act/curv
                    gd = norm_g_non_act*norm_g_non_act/curv
                    is_d_neg_curv = .true.
                endif
            else
                if (it_cg.eq.1) then
                    d = -g_non_act
                    gd = -norm_g_non_act*norm_g_non_act
                    warn_grad = .true.
                endif
                warn_sing = .true.
            endif
            warn_noposdef = .true.
            exit
        endif
        !---------------------------------------------

        ! check if the norm of 'p' is sufficiently large
        !---------------------------------------------
        if (norm_p.le.eps_cg_dir) then
            warn_small = .true.        
            if (it_cg.eq.1) then
                d = -g_non_act
                gd = -norm_g_non_act*norm_g_non_act
                warn_grad = .true.
            endif
            exit
        endif
        !---------------------------------------------
        
        ! update 'd'
        alpha = dot_product(-gq,p_normalized)/curv
        d = d + alpha*p_normalized
        
        hd = hd + alpha*hp
        gd_old = gd
        fquad_old = fquad
        gq = gq + alpha*hp
        gd = dot_product(g_non_act,d)
        fquad = (dot_product(gq,d)+gd)/2.d0
        norm_d = dsqrt(dot_product(d,d))
        
        ! truncated-Newton termination test
        !--------------------------------------------- 
        if (norm_d.lt.norm_g_non_act) then
            eps_cg_tr = dmin1(1.d0,1.d-1+exp(-1.d-3*dble(k)))*dmin1(1.d0,norm_d)
        else
            eps_cg_tr = dmin1(1.d0,1.d-1+exp(-1.d-3*dble(k)))
        endif
        if (dabs((-(3.d0/2.d0)*(gd-gd_old)+fquad-fquad_old)/(-(3.d0/2.d0)*gd+fquad))*dble(it_cg).le.eps_cg_tr) then
            exit
        endif
        !--------------------------------------------- 
        
        ! check if the maximum number of inner iterations has been reached
        !--------------------------------------------- 
        if (it_cg.ge.max_it_cg) then
            warn_maxcgit = .true.
            exit
        endif
        !--------------------------------------------- 
        
        ! update 'p'
        beta = dot_product(gq,hp)/curv
        p = -gq + beta*p_normalized
        
        norm_p = dsqrt(dot_product(p,p))
        gp = dot_product(g_non_act,p)
        gqp = dot_product(gq,p)
        
        ! check conjugacy
        !--------------------------------------------- 
        if ((dabs(gqp-gp).gt.(dabs(gqp)+1.d-6)).or.((gp*gqp.lt.0.d0).and.(dabs(gqp-gp).gt.1.d-6*(dabs(gqp)+1.d-9)))) then
            warn_conjfail = .true.
            exit
        endif
        !--------------------------------------------- 
        
        it_cg = it_cg + 1
        
    enddo
    
    return

end subroutine
!-------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------
subroutine linesearch(n,n_orig,ind_non_act,x,l,u,f,w_max,n_f,min_stepsize,max_n_f,ls_fail,flag)
    
    use asa_bcp_dir
    use asa_bcp_trunc_dir
    use asa_bcp_line_search
    use asa_bcp_max_bounds
    
    implicit none
    
    integer, intent(in) :: n_orig,n,max_n_f
    double precision, intent (in) :: w_max,l(n_orig),u(n_orig),min_stepsize
    logical, intent (in) :: ind_non_act(n_orig)
    integer, intent(inout) :: n_f,flag
    double precision, intent(inout) :: x(n_orig),f
    logical, intent(out) :: ls_fail
    
    double precision, parameter :: gamma = 1.d-6
    double precision, parameter :: delta = 5.d-1
    double precision, parameter :: f_dec = 1.d6
    
    integer :: i,counter,status
    double precision :: fw,fv,aa,a1
    logical :: red_par_computed
    
    ! N.B. 'x' is the current point, 'v' is the point obtained by using
    ! the search direction with unit stepsize (and projecting onto the box)
    ! and 'f' is f(x)
    
    red_par_computed = .false.
    
    if (warn_grad.and.(.not.is_d_neg_curv)) then
		stepsize = dmin1(dmax1(-gd/curv_g_non_act,1.d-3),1.d0)
        if (stepsize.ne.1.d0) then
		    i = 1
            counter = 1
            do while (counter.le.n)
                if (ind_non_act(i)) then
                    v(i) = x(i) + stepsize*d(i)
                    if (dmax1(l_min,v(i)).lt.l(i)) then
                        v(i) = l(i)
                    elseif (dmin1(u_max,v(i)).gt.u(i)) then
                        v(i) = u(i)
                    endif
                    counter = counter + 1
                endif
                i = i + 1
		    enddo
        endif
    else
		stepsize = 1.d0
    endif
	
    call funct(status,n_orig,v,fv)
    if (status.ne.0) then
        flag = 9
        return
    endif
    n_f = n_f + 1
    
    if ((fv-f_newton_first.ge.f_dec*delta_f0).or.warn_noposdef) then
        fw = f
    else
        fw = w_max
    endif
    
    ls_fail = .false.
    
    do
        
        if (fv.le.fw+gamma*stepsize*gd) then
            x = v
            f = fv
            exit
        endif
        
        if (n_f.ge.max_n_f) then
            flag = 5
            exit
        endif
        
        ! update the stepsize
        if ((fv-f)/dmin1(1.d12,dmax1(1.d-12,-gamma*stepsize*gd)).lt.1.d10) then
            stepsize = delta*stepsize
        else
            if (.not.red_par_computed) then
                a1 = 1.d0*dmax1(1.d0,dsqrt(dot_product(x,x)))/dmax1(1.d-15,norm_proj_d)
                red_par_computed = .true.
            endif
            a1 = dmin1(a1,stepsize*delta*delta)
            aa = 2.d-12*stepsize
            if (dmax1(aa,a1).gt.min_stepsize) then
                stepsize = dmax1(aa,a1)
            else
                stepsize = delta*stepsize
            endif
        endif
        
        if (stepsize.le.min_stepsize) then
            ls_fail = .true.
            exit
        endif
        
        ! compute a new point
        i = 1
        counter = 1
        do while (counter.le.n)
            if (ind_non_act(i)) then
                v(i) = x(i) + stepsize*d(i)
                if (dmax1(l_min,v(i)).lt.l(i)) then
                    v(i) = l(i)
                elseif (dmin1(u_max,v(i)).gt.u(i)) then
                    v(i) = u(i)
                endif
                counter = counter + 1
            endif
            i = i + 1
        enddo
        
        call funct(status,n_orig,v,fv)
        if (status.ne.0) then
            flag = 9
            exit
        endif
        n_f = n_f + 1
        
    enddo
    
    return

end subroutine

!-------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------
subroutine update_w(m,f)
    
    use asa_bcp_nm_ref
    
    implicit none
    
    integer, intent(in) :: m
    double precision, intent(in) :: f
    
    integer :: i
    
    i = 1
    do while (i.le.m-1)
        w(i) = w(i+1)
        i = i + 1
    enddo
    w(m) = f
    w_max = f
    mj = min0(mj+1,m-1)
    do i = 1,mj
        if (w(m-i).gt.w_max) then
            w_max = w(m-i)
        endif
    enddo
    
    return

end subroutine
!-------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------
subroutine get_sup_norm_proj_g(sup_norm_proj_grad)
    
    use asa_bcp_gradients, only : sup_norm_proj_g
    
    implicit none
    
    double precision, intent(out) :: sup_norm_proj_grad
    
    sup_norm_proj_grad = sup_norm_proj_g
    
    return

end subroutine
!-------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------
subroutine get_it(iter)
    
    use asa_bcp_counters, only : it
    
    implicit none
    
    integer, intent(out) :: iter
    
    iter = it
    
    return

end subroutine
!-------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------
subroutine get_n_f(funct_eval)
    
    use asa_bcp_counters, only : n_f
    
    implicit none
    
    integer, intent(out) :: funct_eval
    
    funct_eval = n_f
    
    return

end subroutine
!-------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------
subroutine get_n_g(grad_eval)
    
    use asa_bcp_counters, only : n_g
    
    implicit none
    
    integer, intent(out) :: grad_eval
    
    grad_eval = n_g
    
    return

end subroutine
!-------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------
subroutine get_n_hd(hd_eval)
    
    use asa_bcp_counters, only : n_hd
    
    implicit none
    
    integer, intent(out) :: hd_eval
    
    hd_eval = n_hd
    
    return

end subroutine
!-------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------
subroutine get_inner_it(inner_it_cg)
    
    use asa_bcp_counters, only : it_cg_tot
    
    implicit none
    
    integer, intent(out) :: inner_it_cg
    
    inner_it_cg = it_cg_tot
    
    return

end subroutine
!-------------------------------------------------------------------------------------