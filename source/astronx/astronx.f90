! Copyright 2012-2017 Jan von Cosel
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
use iso_fortran_env, only: int32, real64
use globalmod, only: bin_trj, elapsed_time, output, txt_trj
use counters
use input_module, only: N_obj, names, mass, total_mass, do_txttrj, shift_cog, shift_mom, &
                        read_input, process_cmd_arguments, read_input, initialize, verbose, ndigit
use astronx_utils, only: show_input_parameters, shiftcog, shiftmom, centre_of_gravity, linear_momentum, angular_momentum
use propagation
implicit none


! "Data Dictionary": the variables used in this program:
integer(int32) :: i, j                      ! loop indices
real(real64),dimension(:),allocatable :: X  ! spatial coordinates of all objects (m)
real(real64),dimension(:),allocatable :: V  ! velocity components of all objects (m/s)
real(real64),dimension(3) :: ct_of_grav     ! the centre of gravity of the system (m)
real(real64),dimension(3) :: lin_mom        ! the total linear momentum of the system (kg*m/s)
real(real64),dimension(3) :: ang_mom        ! the total angular momentum of the system (kg*m^2/s)
real(real64) :: start_propcpu               ! cpu time at the start of the propagation
real(real64) :: end_propcpu                 ! cpu time at the end of the propagation


call process_cmd_arguments
call read_input(X, V)
call initialize

! write the header to the text trajectory file if requested:
if (do_txttrj) then
    write(txt_trj,'("# trajectory in gnuplot-friendly text form")')
    write(txt_trj,'("#")')
    write(txt_trj,'("# time               ")',advance='no')
    do i = 1, N_obj - 1
        write(txt_trj,'(a30)',advance='no') names(i)
        do j = 1, (3 * ndigit) - 5
            write(txt_trj,'(" ")',advance='no')
        enddo
    enddo
    write(txt_trj,'(a30)',advance='no') names(N_obj)
    write(txt_trj,*)
    write(txt_trj,'("#                    ")',advance='no')
    do i = 1, N_obj - 1
        write(txt_trj,'("x [m]")',advance='no')
        do j = 1, ndigit + 3
            write(txt_trj,'(" ")',advance='no')
        enddo
        write(txt_trj,'("y [m]")',advance='no')
        do j = 1, ndigit + 3
            write(txt_trj,'(" ")',advance='no')
        enddo
        write(txt_trj,'("z [m]")',advance='no')
        do j = 1, ndigit + 4
            write(txt_trj,'(" ")',advance='no')
        enddo
    enddo
    write(txt_trj,'("x [m]")',advance='no')
    do j = 1, ndigit + 3
        write(txt_trj,'(" ")',advance='no')
    enddo
    write(txt_trj,'("y [m]")',advance='no')
    do j = 1, ndigit + 3
        write(txt_trj,'(" ")',advance='no')
    enddo
    write(txt_trj,'("z [m]")',advance='no')
    write(txt_trj,*)
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
    write(output_unit,'(" AstronX: A program for the simulation of celestial mechanics.")')
    write(output_unit,'(" Copyright 2012-2015 Jan von Cosel. Astronx is free software.")')
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

call centre_of_gravity(X, ct_of_grav)
write(output,'("Location of the centre of gravity:")')
write(output,'(" x = ",es15.8," m")') ct_of_grav(1)
write(output,'(" y = ",es15.8," m")') ct_of_grav(2)
write(output,'(" z = ",es15.8," m")') ct_of_grav(3)
write(output,*)

call linear_momentum(V, lin_mom)
write(output,'("Total linear momentum:")')
write(output,'(" x = ",es15.8," kg*m/s")') lin_mom(1)
write(output,'(" y = ",es15.8," kg*m/s")') lin_mom(2)
write(output,'(" z = ",es15.8," kg*m/s")') lin_mom(3)
write(output,*)

call angular_momentum(X, V, ang_mom)
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
    call centre_of_gravity(X, ct_of_grav)
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
    call linear_momentum(V, lin_mom)
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
do i = 1, N_obj
    write(output,'(a10,"  ",es11.3,"  ",3es11.3,"  ",3es11.3)') &
        & trim(names(i)), mass(i), X(3 * (i - 1) + 1), X(3 * (i - 1) + 2), X(3 * (i - 1) + 3), &
                                   V(3 * (i - 1) + 1), V(3 * (i - 1) + 2), V(3 * (i - 1) + 3)
enddo
write(output,*)
write(output,*)

! write the header of the binary trajectory file:
write(bin_trj) N_obj
do i = 1, N_obj
    write(bin_trj) mass(i)
enddo

!initialize the time before the propagation:
elapsed_time = 0.0_real64

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
do i = 1, N_obj
    write(output,'(a10,"  ",es11.3,"  ",3es11.3,"  ",3es11.3)') &
        & trim(names(i)), mass(i), X(3 * (i - 1) + 1), X(3 * (i - 1) + 2), X(3 * (i - 1) + 3), &
                                   V(3 * (i - 1) + 1), V(3 * (i - 1) + 2), V(3 * (i - 1) + 3)
enddo
write(output,*)
write(output,*)
write(output,*)

write(output,'("----------------------------------------------------------")')
write(output,'("CENTRE OF GRAVITY AND LINEAR MOMENTUM AFTER THE SIMULATION")')
write(output,'("----------------------------------------------------------")')
write(output,*)

call centre_of_gravity(X, ct_of_grav)

write(output,'("Location of the centre of gravity:")')
write(output,'(" x = ",es15.8," m")') ct_of_grav(1)
write(output,'(" y = ",es15.8," m")') ct_of_grav(2)
write(output,'(" z = ",es15.8," m")') ct_of_grav(3)
write(output,*)

call linear_momentum(V, lin_mom)

write(output,'("Total linear momentum:")')
write(output,'(" x = ",es15.8," kg*m/s")') lin_mom(1)
write(output,'(" y = ",es15.8," kg*m/s")') lin_mom(2)
write(output,'(" z = ",es15.8," kg*m/s")') lin_mom(3)
write(output,*)

call angular_momentum(X, V, ang_mom)

write(output,'("Total angular momentum:")')
write(output,'(" x = ",es15.8," kg*m^2/s")') ang_mom(1)
write(output,'(" y = ",es15.8," kg*m^2/s")') ang_mom(2)
write(output,'(" z = ",es15.8," kg*m^2/s")') ang_mom(3)
write(output,*)
write(output,*)

write(output,'("   routine         number of calls")')
write(output,'(" ---------------------------------")')
write(output,'(" acceleration    ",i12)') N_acceleration
write(output,'(" bs_largestep    ",i12)') N_bs_largestep
write(output,'(" bs_onestep      ",i12)') N_bs_onestep
write(output,'(" bs_substeps     ",i12)') N_bs_substeps
write(output,'(" extrapolate     ",i12)') N_extrapolate

end program astronx
