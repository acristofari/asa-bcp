! ******************************************************
!            GENERALIZED ROSENBROCK FUNCTION
! ******************************************************


! PROBLEM DIMENSION
!
! Required input/output arguments:
!   status (out) is 0 if no error occurred, >0 otherwise
!   n (out) is the problem dimension
!
!-------------------------------------------------------------------------------------
subroutine prob_dim(status,n)
   
   implicit none
   integer, intent(out) :: status,n
   
   n = 1000
   if (n > 0) then
       status = 0
   else
       status = 1
   endif
   
   return

end subroutine
!-------------------------------------------------------------------------------------


! OBJECTIVE FUNCTION
!
! Required input/output arguments:
!   status (out) is 0 if no error occurred, >0 otherwise
!   n (in) is the dimension
!   x (in) is the point where the objective function will be computed
!   f (out) is the value of objective function at x
!
!-------------------------------------------------------------------------------------
subroutine funct(status,n,x,f)
    
    implicit none
    
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)
    integer, intent(out) :: status
    double precision, intent(out) :: f
    
    integer :: i
    double precision, parameter :: c = 1.d2
    double precision :: t
    
    f = 0.d0
    do i = 1,n-1
        t = x(i)*x(i)
        f = f + (x(i)-1.d0)*(x(i)-1.d0) + c*(x(i+1)-t)*(x(i+1)-t)
    enddo
    
    status = 0
    
    return

    end subroutine
!-------------------------------------------------------------------------------------


! GRADIENT OF THE OBJECTIVE FUNCTION
!
! Required input/output arguments:
!   status (out) is 0 if no error occurred, >0 otherwise
!   n (in) is the dimension
!   x (in) is the point where the gradient of the objective function will be computed
!   g (out) is the gradient of the objective function at x
!
!-------------------------------------------------------------------------------------
subroutine grad(status,n,x,g)
    
    implicit none
    
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)
    integer, intent(out) :: status
    double precision, intent(out) :: g(n)
    
    integer :: i
    double precision, parameter :: c = 1.d2
    
    g(1) = 2.d0*(x(1)-1.d0) + 4.d0*c*x(1)*((x(1)*x(1))-x(2))
    do i = 2,n-1
        g(i) = 2.d0*(x(i)-1.d0) + 4.d0*c*x(i)*((x(i)*x(i))-x(i+1)) - 2.d0*c*(x(i-1)*x(i-1)-x(i))
    enddo
    g(n) = -2.d0*c*(x(n-1)*x(n-1)-x(n))
    
    status = 0
    
    return

    end subroutine
!-------------------------------------------------------------------------------------


! HESSIAN-VECTOR PRODUCT
!
! Required input/output arguments:
!   status (out) is 0 if no error occurred, >0 otherwise
!   n (in) is the dimension
!   goth (in) is .true. if x is the same point used in the previous call, .false. otherwise
!             (it can be used to store values of the Hessian matrix)
!   x (in) is the point where the Hessian-vector product will be computed to be multiplied with d
!   d (in) is a vector whose product with the Hessian matrix is required
!   hd (out) is the Hessian-vector product
!
! N.B. This subroutine can be ignored (in the sense that it can return any dummy value)
! if Hessian-vector products are approximated (i.e., if 'hd_exact' is equal to .false.
! in the subroutine 'asa_bcp')
!
!-------------------------------------------------------------------------------------
subroutine hd_prod(status,n,goth,x,d,hd)
    
    implicit none
    
    integer, intent(in) :: n
    double precision, intent(in) :: x(n),d(n)
    logical, intent(in) :: goth
    integer, intent(out) :: status
    double precision, intent(out) :: hd(n)
    
    integer :: i
    double precision, parameter :: c = 1.d2
    double precision :: a

    ! example of how to use 'goth'
    !
    ! if (.not.goth) then
    !     ! store values of the Hessian matrix in some module or by save variables
    ! endif
    !
    ! ! use the previously saved values to compute H(x)*d 
    
    hd = 0.d0
    a = -4.d0*c*x(1)
    hd(1) = (2.d0+12.d0*c*x(1)*x(1)-4.d0*c*x(2))*d(1) + a*d(2)
    hd(2) = a*d(1)
    do i = 2,n-1
        a = -4.d0*c*x(i)
        hd(i) = hd(i) + (2.d0*c+2.d0+12.d0*c*x(i)*x(i)-4.d0*c*x(i+1))*d(i) + a*d(i+1)
        hd(i+1) = hd(i+1) + a*d(i)
    enddo
    hd(n) = hd(n) + 2.d0*c*d(n)
    
    status = 0
    
    return

    end subroutine
!-------------------------------------------------------------------------------------


! BOUNDS
!
! Required input/output arguments:
!   status (out) is 0 if no error occurred, >0 otherwise
!   n (in) is the dimension
!   l (out) is the vector of lower bound for the variables
!   u (out) is the vector of upper bound for the variables
!
! N.B. There is no i-th lower bound if l(i) <= -1.d20 and there is no i-th upper bound
! if u(i) >= 1.d20.
!
!-------------------------------------------------------------------------------------
subroutine bounds(status,n,l,u)
    
    implicit none
    
    integer, intent(in) :: n
    integer, intent(out) :: status
    double precision, intent(out) :: l(n),u(n)
    
    integer :: i
    
    do i = 1,n,2
        l(i) = -15.d-1
        l(i+1) = 5.d-1
        u(i) = 5.d-1
        u(i+1) = 2.d0
    enddo
    
    status = 0
    
    return

end subroutine
!-------------------------------------------------------------------------------------


! STARTING POINT
!
! Required input/output arguments:
!   status (out) is 0 if no error occurred, >0 otherwise
!   n (in) is the dimension
!   x (out) is the starting point for the algorithm
!
!-------------------------------------------------------------------------------------
subroutine starting_point(status,n,x)
    
    implicit none
    
    integer, intent(in) :: n
    integer, intent(out) :: status
    double precision, intent(out) :: x(n)
    
    integer :: i
    
    do i = 1,n-1,2
       x(i) = -12.d-1
       x(i+1) = 1.d0
    enddo
    
    status = 0
    
    return

end subroutine
!-------------------------------------------------------------------------------------