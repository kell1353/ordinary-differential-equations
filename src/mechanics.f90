!-----------------------------------------------------------------------
!Module: mechanics
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This module is for setting up the differential equations for the case 
!! of orbital motion and subsequently solving for the energy and angular
!! momentum of the system based on given initial conditions.
!! 
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! calculate_energy
!! calculate_angular_momentum
!!----------------------------------------------------------------------
!! Included functions:
!!
!! planets_ode
!-----------------------------------------------------------------------
module mechanics
use types
implicit none

integer, parameter :: G = 1   !Gravitational Constant


private
public :: planets_ode, calculate_energy, calculate_angular_momentum

contains

!-----------------------------------------------------------------------
!! function: planets_ode
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This function takes the masses of each object, the intitial conditions
!! of the dependent variables and the time. It then runs the variables through
!! a set of eight equations and returns the results of each of the dependent
!! variables.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! r        real    array containing the values of the system
!! t        real    the current time value 
!! work     real    array containing mass of each object in the system
!!----------------------------------------------------------------------
!! Output:
!!
!! f        real    array containing the results from the ODE
!-----------------------------------------------------------------------
function planets_ode(r, t, work) result(f)
    implicit none
    real(dp), intent(in) :: r(:), t, work(:)
    real(dp) :: f(1:size(r))
    integer :: n

    real(dp) :: M, m_1, m_2, d_1, d_2, d_12
    ! This is the function that will be sent to 
    ! solve_runge_kutta_4 as an argument.

    ! Initial positions
    !x_1 = r(1) !y_1 = r(2)
    !x_2 = r(3) !y_2 = r(4)

    ! Initial positions
    !v_x1 = r(5) !v_y2 = r(6)
    !v_x2 = r(7) !v_y2 = r(8)

    ! Masses
    M = work(1)
    m_1 = work(2)
    m_2 = work(3)

    ! Distances from other objects
    d_1 = sqrt(r(1)**2 + r(2)**2)
    d_2 = sqrt(r(3)**2 + r(4)**2)
    d_12 = sqrt((r(1) - r(3))**2 + (r(2) - r(4))**2)

    ! System of differential equation 
    ! defined here (d/dt)
    f(1) = r(5)
    f(2) = r(6)
    f(3) = r(7)
    f(4) = r(8)

    f(5) = -r(1)*((G*M)/(d_1**3)) - (r(1) - r(3))*((G*m_2)/(d_12**3))
    f(6) = -r(2)*((G*M)/(d_1**3)) - (r(2) - r(4))*((G*m_2)/(d_12**3))
    f(7) = -r(3)*((G*M)/(d_2**3)) - (r(3) - r(1))*((G*m_1)/(d_12**3))
    f(8) = -r(4)*((G*M)/(d_2**3)) - (r(4) - r(2))*((G*m_1)/(d_12**3))

end function planets_ode


!-----------------------------------------------------------------------
!! Subroutine: calculate_energy
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This subroutine takes the masses of each object and the solutions of 
!! the dependent variables and calculates the total energy of the 
!! system at each time step.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! work          real    array containing mass of each object in the system
!! sol           real    matrix of all positions and velocities calculated at each time step
!!----------------------------------------------------------------------
!! Output:
!!
!! energy        real    array of the energies at each time step
!-----------------------------------------------------------------------
subroutine calculate_energy(work, sol, energy)
    implicit none
    real(dp), intent(in) :: work(:), sol(:,:)
    real(dp), intent(out) :: energy(:)

    real(dp) :: M, m_1, m_2!, r_1, r_2, r_12
    integer :: n_shape(1:2), n, i 

    real(dp), allocatable :: r_1(:), r_2(:), r_12(:)

    ! Initial positions
    !x_1 = sol(1, i) !y_1 = sol(2, i)
    !x_2 = sol(3, i) !y_2 = sol(4, i)

    ! Initial velocities
    !v_x1 = sol(5, i) !v_y1 = sol(6, i)
    !v_x2 = sol(7, i) !v_y2 = sol(8, i)

    n_shape = shape(sol)
    n = n_shape(2)

    allocate(r_1(n), r_2(n), r_12(n))

    ! Masses
    M = work(1)
    m_1 = work(2)
    m_2 = work(3)

    ! Distances from other objects
    r_1 = sqrt(sol(1, :)**2 + sol(2, :)**2)
    r_2 = sqrt(sol(3, :)**2 + sol(4, :)**2)
    r_12 = sqrt((sol(1, :) - sol(3, :))**2 + (sol(2, :) - sol(4, :))**2)

    ! Calculate Energy
    energy = .5*m_1*(sol(5, :)**2 + sol(6, :)**2) + .5*m_2*(sol(7, :)**2 + sol(8, :)**2) - & 
    (G*m_1*M)/r_1 - (G*m_2*M)/r_2 - (G*m_1*m_2)/r_12

    print *, 'The total energy of the system is ', energy(1)
end subroutine calculate_energy


!-----------------------------------------------------------------------
!! subroutine: calculate_angular_momentum
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This subroutine takes the masses of each object and the solutions of 
!! the dependent variables and calculates the total angular momentum of the 
!! system at each time step.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! work      real    array containing mass of each object in the system
!! sol       real    matrix of all positions and velocities calculated at each time step
!!----------------------------------------------------------------------
!! Output:
!!
!! angular_momentum        real    array of the angular momentum at each time step
!-----------------------------------------------------------------------
subroutine calculate_angular_momentum(work, sol, angular_momentum)
    implicit none
    real(dp), intent(in) :: work(:), sol(:,:)
    real(dp), intent(out) :: angular_momentum(:)

    real(dp) :: m_1, m_2
    integer :: n_shape(1:2), n, i 
    real(dp), allocatable :: p_x1(:), p_y1(:), p_x2(:), p_y2(:), l_1(:), l_2(:)

    n_shape = shape(sol)
    n = n_shape(2)

    allocate(p_x1(n), p_y1(n), p_x2(n), p_y2(n), l_1(n), l_2(n))

    ! Masses
    m_1 = work(2)
    m_2 = work(3)

    ! Linear Momentums
    p_x1 = m_1*sol(5, :)
    p_y1 = m_1*sol(6, :)
    p_x2 = m_2*sol(7, :)
    p_y2 = m_2*sol(8, :)

    ! Angular Momentums (l = r x p)
    l_1 = sol(1, :)*p_y1 - sol(2, :)*p_x1
    l_2 = sol(3, :)*p_y2 - sol(4, :)*p_x2

    angular_momentum = l_1 + l_2

end subroutine calculate_angular_momentum

end module mechanics