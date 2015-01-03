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
module propagate
!
! this module contains the basic propagating routines
!
contains

!##################################################################################################
!##################################################################################################

subroutine propagate_bs(X, V)
!
! this routine will propagate the system using the bulirsch-stoer algorithm
!
use types
use shared_data, only: elapsed_time, restart, restart_file, output, steps, underflow
use input_module, only: N_obj, mass, do_restart, tfinal, do_steps, init_step, verbose
use astronx_utils, only: write_to_trj
use bulirsch_stoer, only: bs_largestep
implicit none


! arguments to the subroutine:
real(dp),dimension(N_obj,3),intent(inout) :: X  ! spatial coordinates of all objects (m)
real(dp),dimension(N_obj,3),intent(inout) :: V  ! velocity components of all objects (m/s)

! internal variables:
integer(st) :: restart_status   ! IO status of the restart file
integer(st) :: i                ! loop index
integer(st) :: N_ok             ! number of successful BS steps per large step
integer(st) :: N_ok_tot         ! total number of successful BS steps
integer(st) :: N_fail           ! number of BS steps per large step, where the stepsize had to be reduced
integer(st) :: N_fail_tot       ! total number of failed BS steps
integer(st) :: N_bssteps        ! total number of BS steps per large step
integer(st) :: N_bstotal        ! overall number of BS steps
integer(st) :: N_smallsteps     ! number of BS substeps per large step
integer(st) :: N_smalltotal     ! total number of BS substeps
real(dp) :: timestep            ! initial stepsize for the next large step (s)
real(dp) :: next_step           ! suggested next stepsize returned by bs_largestep (s)
real(dp) :: start_stepcpu       ! cpu time before one large step
real(dp) :: end_stepcpu         ! cpu time after one large step
real(dp) :: start_totcpu        ! cpu time at the beginning of the calculation
real(dp) :: end_totcpu          ! cpu time at the end of the calculation


timestep = init_step
N_ok_tot = 0
N_fail_tot = 0
N_bstotal = 0
N_smalltotal = 0

write(output,*) "                  ***************************************"
write(output,*) "                  * USING THE BULIRSCH-STOER INTEGRATOR *"
write(output,*) "                  ***************************************"
write(output,*) ""
write(output,*) ""
write(output,*) "--------------------------"
write(output,*) "DETAILS OF THE PROPAGATION"
write(output,*) "--------------------------"
write(output,*) ""
write(output,*) "elapsed time    large steps    BS steps    small steps    cpu time [ms]"
write(output,*) "                good    bad"

flush(output)

! do the propagation:

call cpu_time(start_totcpu)

do
    call write_to_trj(elapsed_time, X, V)

    if (do_restart) then
        open(unit=restart,file=restart_file,status="replace",action="write",iostat=restart_status)
        write(restart,*) elapsed_time
        do i = 1, N_obj
            write(restart,*) mass(i), X(i,1), X(i,2), X(i,3), V(i,1), V(i,2), V(i,3)
        enddo
        close(unit=restart)
    endif

    if ((elapsed_time >= tfinal) .or. (tfinal - elapsed_time < 1.0_dp)) exit

    if (do_steps) then
        write(steps,100) elapsed_time
        100 format ("# elapsed time:    ", es16.9)
        write(steps,101) timestep
        101 format ("# trying propagation with step:   ", es16.9)
    endif

    call cpu_time(start_stepcpu)
    call bs_largestep(timestep, next_step, X, V, N_ok, N_fail, N_bssteps, N_smallsteps)
    call cpu_time(end_stepcpu)

    if (underflow) exit

    timestep = next_step
    N_ok_tot = N_ok_tot + N_ok
    N_fail_tot = N_fail_tot + N_fail
    N_bstotal = N_bstotal + N_bssteps
    N_smalltotal = N_smalltotal + N_smallsteps

    if (do_steps) then
        write(steps,'("# successful / failed steps: ", i4, ",", i4)') N_ok, N_fail
        write(steps,'("#")')
        write(steps,'("######################################################################")')
        write(steps,'("#")')
    endif

    if (verbose) then
    endif

    write(output,102) elapsed_time, N_ok, N_fail, N_bssteps, N_smallsteps, (end_stepcpu-start_stepcpu)*1000.0_dp
    102 format  (' ', es11.5, "   ", i5, "  ", i5, "     ", i5, "      ", i7, "       ", f8.1)
    flush(output)
enddo

call cpu_time(end_totcpu)

write(output,*) ""
write(output,*) ""
write(output,*) ""
write(output,*) "--------------------------"
write(output,*) "SUMMARY OF THE CALCULATION"
write(output,*) "--------------------------"
write(output,*) ""
write(output,*) "total time      large steps    BS steps    small steps    cpu time used [s]"
write(output,*) "                good    bad"

write(output,131) elapsed_time, N_ok_tot, N_fail_tot, N_bstotal, N_smalltotal, (end_totcpu-start_totcpu)
131 format (' ', es10.3, "   ", i7, i7, "  ", i8, "     ", i8, "     ", f12.3)

write(output,*) ""
write(output,*) ""
write(output,*) ""

if (do_steps) then
    write(steps,*) "# Simulation done!"
    write(steps,*) "# Total number of successful / failed steps: ", N_ok_tot, N_fail_tot
endif

end subroutine propagate_bs

!##################################################################################################
!##################################################################################################

!subroutine propagate_rk(X, V)
!
! this routine will propagate the system using the runge-kutta algorithm of fourth order
!
!use types
!use shared_data, only: N_obj, elapsed_time, restart, restart_file, mass, output, steps
!use input_parameters, only: do_restart, tfinal, do_steps, init_step
!use astronx_utils, only: write_to_trj
!use bulirsch_stoer, only: bs_largestep
!implicit none


! arguments to the subroutine:
!real(dp),dimension(N_obj,3),intent(inout) :: X  ! spatial coordinates of all objects (m)
!real(dp),dimension(N_obj,3),intent(inout) :: V  ! velocity components of all objects (m/s)

! internal variables:
!integer(st) :: restart_status
!integer(st) :: i
!integer(st) :: N_ok
!integer(st) :: N_ok_tot
!integer(st) :: N_fail
!integer(st) :: N_fail_tot
!integer(st) :: N_bssteps
!integer(st) :: N_bstotal
!integer(st) :: N_smallsteps
!integer(st) :: N_smalltotal
!real(dp) :: timestep
!real(dp) :: next_step
!real(dp) :: start_stepcpu
!real(dp) :: end_stepcpu
!real(dp) :: start_totcpu
!real(dp) :: end_totcpu
!
!
!timestep = init_step
!N_ok_tot = 0
!N_fail_tot = 0
!N_bstotal = 0
!N_smalltotal = 0
!
!write(output,*) " Using the Runge-Kutta integrator."
!write(output,*) "-----------------------------------------------------------------------------"
!write(output,*) ""
!write(output,*) "elapsed time    large steps    BS steps    small steps    cpu time [ms]"
!write(output,*) "                good    bad"


! doing the propagation:
!call cpu_time(start_totcpu)
!do
!    call write_to_trj(elapsed_time, X, V)
!    if (do_restart) then
!        open(unit=restart,file=restart_file,status="replace",action="write",iostat=restart_status)
!        write(restart,*) elapsed_time
!        do i = 1, N_obj
!            write(restart,*) mass(i), X(i,1), X(i,2), X(i,3), V(i,1), V(i,2), V(i,3)
!        enddo
!        close(unit=restart)
!    endif
!    if (do_steps) then
!        write(steps,100) elapsed_time
!        100 format (" elapsed time:    ", es16.9)
!        write(steps,101) timestep
!        101 format (" trying propagation with step:   ", es16.9)
!    endif
!    if ((elapsed_time >= tfinal) .or. (tfinal - elapsed_time < 1.0_dp)) exit
!    call cpu_time(start_stepcpu)
!    call bs_largestep(timestep, next_step, X, V, N_ok, N_fail, N_bssteps, N_smallsteps)
!    call cpu_time(end_stepcpu)
!    timestep = next_step
!    N_ok_tot = N_ok_tot + N_ok
!    N_fail_tot = N_fail_tot + N_fail
!    N_bstotal = N_bstotal + N_bssteps
!    N_smalltotal = N_smalltotal + N_smallsteps
!    if (do_steps) then
!        write(steps,*) "main: successful / failed steps: ", N_ok, N_fail
!        write(steps,*) ""
!        write(steps,*) "######################################################################"
!        write(steps,*) ""
!    endif
!    write(output,102) elapsed_time, N_ok, N_fail, N_bssteps, N_smallsteps, (end_stepcpu-start_stepcpu)*1000.0_dp
!    102 format  (' ', es12.5, "   ", i5, "  ", i5, "     ", i5, "      ", i7, "       ", f8.1)
!enddo
!
!call cpu_time(end_totcpu)
!
!write(output,*) "----------------------------------------------------------------------------"
!write(output,*) ""
!write(output,*) "Summary of the calculation:"
!write(output,*) "total time      large steps    BS steps    small steps    cpu time used [ms]"
!write(output,*) "                good    bad"
!write(output,131) elapsed_time, N_ok_tot, N_fail_tot, N_bstotal, N_smalltotal, (end_totcpu-start_totcpu)*1000.0_dp
!131 format (' ', es10.3, "   ", i7, i7, "  ", i8, "     ", i8, "     ", f10.1)
!
!
!if (do_steps) then
!    write(steps,*) "Simulation done!"
!    write(steps,*) "Total number of successful / failed steps: ", N_ok_tot, N_fail_tot
!endif
!end subroutine propagate_rk

end module propagate
