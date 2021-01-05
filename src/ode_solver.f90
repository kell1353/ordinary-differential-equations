!-----------------------------------------------------------------------
!Module: ode_solver
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! The purpose of this module is solve any ordinary differential equation
!! numerically using the 4th Order Runga Kutta algorithm for given initial
!! conditions.
!!
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! solve_runge_kutta_4
!!----------------------------------------------------------------------
!
module ode_solver
use types
implicit none
private

public :: solve_runge_kutta_4

interface
    function func(r, t, work) result(f)
        use types, only : dp
        implicit none
        real(dp), intent(in) :: r(:), t, work(:)
        real(dp) :: f(1:size(r))
    end function func
end interface

contains


!-----------------------------------------------------------------------
!! Subroutine: solve_runge_kutta_4
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This subroutine takes in initial conditions, and interval steps for the ordinary 
!! differential equation and applies the 4th Order Runga Kutta algorithm to solve for the 
!! dependent variables against the independent variable (time). It also
!! takes in the work array which holds unchanging variables useful for the ODE.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! t_i 		real 		initial time
!! t_f 		real 		final time
!! q_i		real 		array of initial differential values
!! f        real        array of the function to apply the algorithm to
!! n 		real 		number of points
!!----------------------------------------------------------------------
!! Output:
!!
!! t        real        array of the time at each step 
!! q        real        matrix of the solutions for each dependent variables
!-----------------------------------------------------------------------
subroutine solve_runge_kutta_4(f, t_i, q_i, work, t_f, n_points, t, q)
    implicit none
    procedure(func) :: f 

    real(dp), intent(in) :: t_i, q_i(:), work(:), t_f
    integer, intent(in) :: n_points
    real(dp), intent(out), allocatable :: t(:), q(:, :)
    integer :: n_variables, i, n
    real(dp) :: h
    real(dp), allocatable :: k1(:), k2(:), k3(:), k4(:), q_sol(:), t_sol

    n_variables = size(q_i)

    allocate(t(n_points), q(n_variables, n_points))
    allocate(k1(n_variables), k2(n_variables), k3(n_variables), k4(n_variables), q_sol(n_variables))


    n = n_points

    h = (t_f - t_i)/n
    q_sol = q_i
    t_sol = t_i

    do i=1,n
        k1 = h*f(q_sol, t_sol, work)
        k2 = h*f(q_sol + 0.5_dp*k1, t_sol + 0.5_dp*h, work)
        k3 = h*f(q_sol + 0.5_dp*k2, t_sol + 0.5_dp*h, work)
        k4 = h*f(q_sol + k3, t_sol + h, work)

        q_sol = q_sol + (1.0_dp/6.0_dp)*(k1 + 2*k2 + 2*k3 + k4)
        t_sol = t_sol + h
        
        t(i) = t_sol
        q(:, i) = q_sol
    enddo

end subroutine solve_runge_kutta_4

end module ode_solver