# Overview
This program calculates numerically the orbital motion of two masses around a central mass using the 4th Order Runga Kutta algorithm and 
provides output for analysis.


# Compilation Instructions
Compilation is all done using the makefile in the repository. Type `make` into your command line to compile the files.

Initially sets planetary as the executable name.
- Then creates the types object file.
- Then creates the ode_solver object file using the types object file.
- Then creates the mechanics object file using the types object file.
- Then creates the read_write object file using the types object files.
- Then creates the main object file using the types, read_write, ode_solver and mechanics object files.


# Usage Instructions 
Once you have compiled everything execute the program using ```./planetary arg1, arg2, ...```
The program can be excuted using just the executable ```./planetary``` (which will return results for default initial conditions) or it can take 
namelist files as optional arguments. Just make sure that you keep the namelist files in the directory you are working in.

The namelist files need to be structured like:

```
&masses
	primary_mass = 100.0, 
	planet_mass_1 = .1, 
	planet_mass_2 = .2
/
&initial_conditions
	initial_time = 0.0,
	initial_pos_1 = 50.0, 0.0, 
	initial_pos_2 = 200.0, 0.0,
	initial_vel_1 = 0.0, 1.414, 
	initial_vel_2 = 0.0, .747,
/
&solution_parameters
	final_time = 5000.0,
	n_steps = 10000
/
&output
	output_file = 'planetary_motion.dat'
/

```

Make sure to include an extra line at the end of the file (after the last hash).
Once created you will be able to run the Jupyter notebook the data files produced. 


# Expected Behavior
Once you have executed it will perform the calculations, print the total energies of the system for each set of initial conditions and writes the results 
of the program into .dat files.

The output files should have 7 columns. 'time', 'x_1', 'y_1', 'x_2', 'y_2', 'energy', 'angular momentum'.

- time: This column will contain the time at each step.
- x_1: This column will contain the x position of the first object at each step.
- y_1: This column will contain the y position of the first object at each step.
- x_2: This column will contain the x position of the second object at each step.
- y_2: This column will contain the y position of the second object at each step.
- energy: This column will contain the total energy of the system at each step.
- angular momentum: This column will contain the total angular momentum of the system at each step.
