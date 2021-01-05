! Program: planets
! 
!-----------------------------------------------------------------------------
! By: Austin Keller
!
! This program calculates numerically the orbital motion of two masses, 
! around a central mass, using the 4th Order Runga Kutta algorithm and 
! provides output for analysis. It allows for the caluculation of multiple
! namelist files as arguments.
!-----------------------------------------------------------------------------
program planets

use types
use read_write, only : read_input, write_results
use ode_solver, only : solve_runge_kutta_4
use mechanics, only : planets_ode, calculate_energy, calculate_angular_momentum

implicit none

real(dp) :: work_array(1:3), initial_conditions(1:8)
real(dp) :: final_time, initial_time
integer :: n_steps, num_arguments, i, n
real(dp), allocatable :: time(:), solution(:,:), energy(:), angular_momentum(:)
character(len=1024) :: output_file


num_arguments = command_argument_count() 

if (num_arguments == 0) then
    n = num_arguments + 1 ! To enter the do loop when 0 args
else 
    n = num_arguments
endif


do i=1,n
    call read_input(i, num_arguments, work_array, initial_conditions, initial_time, final_time, n_steps, output_file)
    call solve_runge_kutta_4(planets_ode, initial_time, initial_conditions, work_array, final_time, n_steps, time, solution)

    ! Allocate arrays
    if(allocated(energy)) deallocate(energy)
    allocate(energy(1:n_steps))
    if(allocated(angular_momentum)) deallocate(angular_momentum)
    allocate(angular_momentum(1:n_steps))

    ! Call energy and angular momentum calculation
    call calculate_energy(work_array, solution, energy)
    call calculate_angular_momentum(work_array, solution, angular_momentum)

    ! Writes the results
    call write_results(output_file, initial_time, initial_conditions, energy, angular_momentum, time, solution)
enddo

end program planets