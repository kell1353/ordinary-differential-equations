!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! The purpose of this module to read initial conditions from one or 
!! several namelist files and write the results of the program to .dat
!! files.
!!
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! read_input
!! write_results
!!----------------------------------------------------------------------
module read_write
use types
implicit none

private
public :: read_input, write_results

contains


!-----------------------------------------------------------------------
!! Subroutine: read_input
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This subroutine takes in several inputs, sets up the structure of the 
!! namelist file to be read. Then sets default values for the initial 
!! conditions in the case no file is enetered as an argument. 
!!
!! In the other case, it then checks if the number of arguments is greater 
!! then zero, then calls for the argument number and checks if the namelist
!! file exists. 
!! If so, it reads each sections values in the namelist file and stops
!! the program if there is an error.
!!
!! After all that it outputs the values in two seperate arrays one for 
!! the initial conditions and one for the other important values for the 
!! ODE.
!!
!!----------------------------------------------------------------------
!! Input:
!! i                    integer     the argument number 
!! n_arguments          integer     the number of arguments the user inputs
!!----------------------------------------------------------------------
!! Output:
!!
!! work_array           real        array containing mass of each object in the system
!! initial_condition    real        array containing the initial conditions of the system
!! initial_time         real        initial time provided by the namelist file
!! final_time           real        final time provided by the namelist file
!! n_steps              integer     the number of time steps 
!! output_file          character   the name of the output file 
!-----------------------------------------------------------------------
subroutine read_input(i, n_arguments, work_array, initial_condition, initial_time, final_time, n_steps, output_file)
    implicit none
    integer, intent(in) :: i, n_arguments
    real(dp), intent(out) :: work_array(1:3)
    real(dp), intent(out) :: initial_condition(1:8)
    real(dp), intent(out) :: initial_time, final_time
    integer, intent(out) :: n_steps
    character(len=*), intent(out) :: output_file
    
    real(dp) :: primary_mass, planet_mass_1, planet_mass_2
    real(dp) :: initial_pos_1(1:2), initial_pos_2(1:2)
    real(dp) :: initial_vel_1(1:2), initial_vel_2(1:2)

    integer :: unit, ierror
    character(len=1024) :: namelist_file
    logical :: file_exists


    ! Must have an empty line at the end of the namelist
    namelist /masses/ primary_mass, planet_mass_1, planet_mass_2
    namelist /initial_conditions/ initial_time, initial_pos_1, initial_pos_2, initial_vel_1, initial_vel_2
    namelist /solution_parameters/ final_time, n_steps
    namelist /output/ output_file


    primary_mass  = 1._dp
    planet_mass_1 = 1._dp
    planet_mass_2 = 0._dp
    initial_pos_1 = [1._dp, 0._dp]
    initial_pos_2 = [0._dp, 1._dp]
    initial_vel_1 = [0._dp, 1._dp]
    initial_vel_2 = [0._dp, 0._dp]
    initial_time = 0._dp
    final_time = 10._dp
    n_steps = 100
    output_file = 'planetary_motion.dat'
    

    !n_arguments = command_argument_count()

    if (n_arguments > 0) then
        ! get namelist file name from command line
        call get_command_argument(i, namelist_file)
        inquire(file=trim(namelist_file), exist = file_exists)

        if (file_exists) then
            open(newunit = unit, file = namelist_file )

            ! read namelists
            read(unit, nml = masses, iostat = ierror)
            if(ierror /= 0) then
                print*, 'Error reading masses namelist'
                stop
            endif

            read(unit, nml = initial_conditions, iostat = ierror)
            if(ierror /= 0) then
                print*, 'Error reading initial_conditions namelist'
                stop
            endif

            read(unit, nml = solution_parameters, iostat = ierror)
            if(ierror /= 0) then
                print*, 'Error reading solution_parameters namelist'
                stop
            endif

            read(unit, nml = output, iostat = ierror)
            if(ierror /= 0) then
                print*, 'Error reading output namelist'
                stop
            endif
            close(unit)
        else
            print*, 'Argument, ', trim(namelist_file)
            print*, 'does not exist. Ending program'
            stop
        endif

    elseif(n_arguments /= 0) then
        print*, 'Incorrect number of arguments'
        print*, 'The program takes either 0 or 1 arguments'
        print*, 'See documentation on README.md for details'
        stop
    endif

    initial_condition = [initial_pos_1, initial_pos_2, initial_vel_1, initial_vel_2]
    work_array = [primary_mass, planet_mass_1, planet_mass_2]

end subroutine read_input


!-----------------------------------------------------------------------
!! Subroutine: read_input
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This function takes in the results from the program and writes the 
!! results into a file to be analyzed using Jupyter Notebook.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! file_name        character   the name of the file to write to
!! t_i              real        initial time provided by the namelist file
!! q_i              real        array containing the initial conditions of the system
!! energy           real        array containing the energies at each time step
!! ang_momentum     real        array containing the angular momentum at each time step
!! t                real        array containing the time at each step
!! q                real        array containing the ODE solution values at each step
!-----------------------------------------------------------------------
subroutine write_results(file_name, t_i, q_i, energy, ang_momentum, t, q)
    implicit none
    character(len=*), intent(in) :: file_name
    real(dp), intent(in) :: t_i, q_i(:), t(:), q(:,:), energy(:), ang_momentum(:)

    integer :: n, i, unit

    n = size(t)

    open(newunit = unit, file = trim(file_name))

    ! number of headers: 6, number of spaces: 18
    write(unit,*) 'The initial conditions used: '
    write(unit,*) 'Inital Time: ', t_i
    write(unit,*) 'Inital Position of Object 1: ', q_i(1), q_i(2)
    write(unit,*) 'Inital Position of Object 2: ', q_i(3), q_i(4)
    write(unit,*) 'Inital Velocity of Object 1: ', q_i(5), q_i(6)
    write(unit,*) 'Inital Velocity of Object 2: ', q_i(7), q_i(8)

    !write, q_i

    write(unit,*) ! To space out output and initial conditions

    write(unit,'(7a20)') 'time', 'x_1', 'y_1', 'x_2', 'y_2', 'energy', 'angular_momentum'
    do i=1,n
        write(unit,*) t(i), q(1, i), q(2, i), q(3, i), q(4, i), energy(i), ang_momentum(i)
    enddo

end subroutine write_results


end module read_write
