! Copyright 2012-2015 Jan von Cosel
!
! This file is part of astronx.
!
! astronx if free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! astronx is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have recieved a copy of the GNU General Public License
! along with astronx. If not, see <http://www.gnu.org/licenses/>.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
program astronx
!
!  This is the main program which initializes the simulation and calls
!  the relevant routines.
!
use types
use shared_data, only: bin_trj, elapsed_time, output, steps, trajectory
use input_module, only: names, mass, total_mass, do_steps, do_texttrj, name_directory, shift_cog, shift_mom, read_input, &
                        process_cmd_arguments, verbose
use astronx_utils, only: show_input_parameters, shiftcog, shiftmom, centre_of_gravity, linear_momentum, angular_momentum
use propagation
implicit none


! "Data Dictionary": the variables used in this program:
character(len=40) :: input_file                 ! take an educated guess...
character(len=10) :: bin_traj_file="trajectory" ! the trajectory file in binary format
character(len=9) :: traj_file="text_traj"       ! file to save the trajectory
character(len=6) :: output_file="output"        ! the output file with general information about the calculation
character(len=5) :: steps_file="steps"          ! this file contains detailed info about the propagation steps
integer(st) :: i                                ! loop indices
integer(st) :: output_status                    ! describing the file opening status
real(dp),dimension(:,:),allocatable :: X        ! spatial coordinates of all objects (m)
real(dp),dimension(:,:),allocatable :: V        ! velocity components of all objects (m/s)
real(dp),dimension(3) :: ct_of_grav             ! the centre of gravity of the system (m)
real(dp),dimension(3) :: lin_mom                ! the total linear momentum of the system (kg*m/s)
real(dp),dimension(3) :: ang_mom                ! the total angular momentum of the system (kg*m^2/s)
real(dp) :: start_propcpu                       ! cpu time at the start of the propagation
real(dp) :: end_propcpu                         ! cpu time at the end of the propagation


! process the command line arguments:
call process_cmd_arguments

! read the input and move to the working directory:
call get_command_argument(1, input_file)
call read_input(input_file, X, V)
call chdir(name_directory)

! open the output file for writing:
open(unit=output,file=output_file,status="replace",action="write",iostat=output_status)
if (output_status /= 0) then
    write(*,'("ERROR: failed to open output file, status: ",i3)') output_status
    stop
endif

! open the text trajectory file if requested:
if (do_texttrj) then
    open(unit=trajectory,file=traj_file,status="replace",action="write",iostat=output_status)
    if (output_status /= 0) then
        write(*,'("ERROR: failed to open text trajectory file, status: ",i3)') output_status
        write(output,'("ERROR: failed to open text trajectory file, status: ",i3)') output_status
        stop
    endif
    write(trajectory,'("# trajectory in gnuplot-friendly text form")')
    write(trajectory,'("#")')
    write(trajectory,'("# time          ")',advance='no')
    do i = 1, size(mass)
        write(trajectory,'(a54)',advance='no') names(i)
    enddo
    write(trajectory,*)
    write(trajectory,'("#                   ")',advance='no')
    do i = 1, size(mass)
        write(trajectory,'("x [m]             y [m]             z [m]             ")',advance='no')
    enddo
    write(trajectory,*)
endif

! open the steps file if requested:
if (do_steps) then
    open(unit=steps,file=steps_file,status="replace",action="write",iostat=output_status)
    if (output_status /= 0 ) then
        write(*,'("ERROR: failed to open steps file, status: ",i3)')  output_status
        write(output,'("ERROR: failed to open steps file, status: ",i3)')  output_status
        stop
    endif
endif

! open the binary trajectory file:
open(unit=bin_trj,file=bin_traj_file,status="replace",action="write",iostat=output_status,form="unformatted")
if (output_status /= 0) then
    write(*,'("ERROR: failed to open binary trajectory file, status: ",i3)') output_status
    write(output,'("ERROR: failed to open binary trajectory file, status: ",i3)') output_status
    stop
endif


! write the header of the output file:
write(output,'("/-------------------------------------------------------------------------------\")')
write(output,'("|                                  ** AstronX **                                |")')
write(output,'("|                                                                               |")')
write(output,'("|                A program for the simulation of celestial mechanics            |")')
write(output,'("|      /\                                                              \\    // |")')
write(output,'("|     //\\                Copyright 2012-2015 Jan von Cosel             \\  //  |")')
write(output,'("|    //  \\                                                              \\//   |")')
write(output,'("|   //====\\                  Astronx is free software.                  //\\   |")')
write(output,'("|  //      \\    You can redistribute it and/or modify it under the     //  \\  |")')
write(output,'("| //        \\   terms of the GNU General Public License as published  //    \\ |")')
write(output,'("|                by the Free Software Foundation, either version 3 of           |")')
write(output,'("|                the License, or (at your option) any later version.            |")')
write(output,'("|                     Astronx comes with absolutely no warranty.                |")')
write(output,'("\-------------------------------------------------------------------------------/")')
write(output,*)
write(output,*)
write(output,*)
if (verbose) then
    write(*,'(" AstronX: A program for the simulation of celestial mechanics.")')
    write(*,'(" Copyright 2012-2015 Jan von Cosel. Astronx is free software.")')
endif

! write some information about the system to the output file:
write(output,*)
write(output,*)
write(output,'("------------------------------------")')
write(output,'("GENERAL INFORMATION ABOUT THE SYSTEM")')
write(output,'("------------------------------------")')
write(output,*)
write(output,'("Total mass:")')
write(output,'(" m = ",es15.9," kg")') total_mass
write(output,*)

call centre_of_gravity(X, ct_of_grav, mass, total_mass)
write(output,'("Location of the centre of gravity:")')
write(output,'(" x = ",es15.8," m")') ct_of_grav(1)
write(output,'(" y = ",es15.8," m")') ct_of_grav(2)
write(output,'(" z = ",es15.8," m")') ct_of_grav(3)
write(output,*)

call linear_momentum(V, lin_mom, mass)
write(output,'("Total linear momentum:")')
write(output,'(" x = ",es15.8," kg*m/s")') lin_mom(1)
write(output,'(" y = ",es15.8," kg*m/s")') lin_mom(2)
write(output,'(" z = ",es15.8," kg*m/s")') lin_mom(3)
write(output,*)

call angular_momentum(X, V, ang_mom, mass)
write(output,'("Total angular momentum:")')
write(output,'(" x = ",es15.8," kg*m^2/s")') ang_mom(1)
write(output,'(" y = ",es15.8," kg*m^2/s")') ang_mom(2)
write(output,'(" z = ",es15.8," kg*m^2/s")') ang_mom(3)
write(output,*)
write(output,*)
write(output,*)


! write the values of all parameters to the output file:
call show_input_parameters

! shift the cog if requested:
if (shift_cog) then
    call shiftcog(X)
    call centre_of_gravity(X, ct_of_grav, mass, total_mass)
    write(output,'("--------------------------------------------")')
    write(output,'("SHIFTING THE CENTRE OF GRAVITY TO THE ORIGIN")')
    write(output,'("--------------------------------------------")')
    write(output,*)
    write(output,'("Centre of gravity after shifting:")')
    write(output,'(" x = ",es15.8," m")') ct_of_grav(1)
    write(output,'(" y = ",es15.8," m")') ct_of_grav(2)
    write(output,'(" z = ",es15.8," m")') ct_of_grav(3)
    write(output,*)
    write(output,*)
    write(output,*)
endif

! shift the momentum if requested:
if (shift_mom) then
    call shiftmom(V)
    call linear_momentum(V, lin_mom, mass)
    write(output,'("--------------------------------------------")')
    write(output,'("ELIMINATING THE TOTAL MOMENTUM OF THE SYSTEM")')
    write(output,'("--------------------------------------------")')
    write(output,*)
    write(output,'("Residual linear momentum:")')
    write(output,'(" x = ",es15.8," kg*m/s")') lin_mom(1)
    write(output,'(" y = ",es15.8," kg*m/s")') lin_mom(2)
    write(output,'(" z = ",es15.8," kg*m/s")') lin_mom(3)
    write(output,*)
    write(output,*)
    write(output,*)
endif



! write the initial positions and velocities to the output file:
write(output,'("----------------------------------")')
write(output,'("INITIAL COORDINATES AND VELOCITIES")')
write(output,'("----------------------------------")')
write(output,*)
write(output,'("    name      mass (kg)       X (m)      Y (m)      Z (m)     V_x (m/s)  V_y (m/s)  V_z (m/s)")')
write(output,*)
do i = 1, size(mass)
    write(output,'(a10,"  ",es11.3,"  ",3es11.3,"  ",3es11.3)') &
        & trim(names(i)), mass(i), X(i,1), X(i,2), X(i,3), V(i,1), V(i,2), V(i,3)
enddo
write(output,*)
write(output,*)

! write the header of the binary trajectory file:
write(bin_trj) size(mass)
do i = 1, size(mass)
    write(bin_trj) mass(i)
enddo

!initialize the time before the propagation:
elapsed_time = 0.0_dp

flush(output)

! here comes the actual simulation (either rk or bs):
write(output,'("  ----------------------------------------------------------------------------------")')
write(output,'("                                  STARTING THE PROPAGATION")')
write(output,'("  ----------------------------------------------------------------------------------")')
write(output,*)

call cpu_time(start_propcpu)
call propagate(X, V)
call cpu_time(end_propcpu)

write(output,'("  ----------------------------------------------------------------------------------")')
write(output,'("                                 FINISHED THE PROPAGATION")')
write(output,'("                           total cpu time:", f12.3, " seconds")') end_propcpu - start_propcpu
write(output,'("  ----------------------------------------------------------------------------------")')

if (verbose) then
    write(*,*)
    write(*,'(" Finished the propagation using ", f0.3, " s of CPU time.")') end_propcpu - start_propcpu
endif


!#####################################################################
! TO DO:  react to the "underflow" flag from the propagation
!#####################################################################



! and now the postproduction (coordinates of last frame):
write(output,'("--------------------------------")')
write(output,'("FINAL COORDINATES AND VELOCITIES")')
write(output,'("--------------------------------")')
write(output,*)
write(output,'("    name      mass (kg)       X (m)      Y (m)      Z (m)     V_x (m/s)  V_y (m/s)  V_z (m/s)")')
write(output,*)
do i = 1, size(mass)
    write(output,'(a10,"  ",es11.3,"  ",3es11.3,"  ",3es11.3)') &
        & trim(names(i)), mass(i), X(i,1), X(i,2), X(i,3), V(i,1), V(i,2), V(i,3)
enddo
write(output,*)
write(output,*)
write(output,*)

write(output,'("----------------------------------------------------------")')
write(output,'("CENTRE OF GRAVITY AND LINEAR MOMENTUM AFTER THE SIMULATION")')
write(output,'("----------------------------------------------------------")')
write(output,*)

call centre_of_gravity(X, ct_of_grav, mass, total_mass)

write(output,'("Location of the centre of gravity:")')
write(output,'(" x = ",es15.8," m")') ct_of_grav(1)
write(output,'(" y = ",es15.8," m")') ct_of_grav(2)
write(output,'(" z = ",es15.8," m")') ct_of_grav(3)
write(output,*)

call linear_momentum(V, lin_mom, mass)

write(output,'("Total linear momentum:")')
write(output,'(" x = ",es15.8," kg*m/s")') lin_mom(1)
write(output,'(" y = ",es15.8," kg*m/s")') lin_mom(2)
write(output,'(" z = ",es15.8," kg*m/s")') lin_mom(3)
write(output,*)

call angular_momentum(X, V, ang_mom, mass)

write(output,'("Total angular momentum:")')
write(output,'(" x = ",es15.8," kg*m^2/s")') ang_mom(1)
write(output,'(" y = ",es15.8," kg*m^2/s")') ang_mom(2)
write(output,'(" z = ",es15.8," kg*m^2/s")') ang_mom(3)
write(output,*)
write(output,*)


end program astronx
