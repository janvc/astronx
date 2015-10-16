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
module propagation
!
! this module contains the basic propagating routines
!
contains

!##################################################################################################
!##################################################################################################

subroutine propagate(X, V)
!
! this routine will propagate the system with the chosen integration algoithm
!
use types
use shared_data, only: elapsed_time, restart, restart_file, output, steps, underflow
use input_module, only: N_obj, mass, do_restart, tfinal, do_steps, init_step, do_bs, do_rk4mid, verbose
use astronx_utils, only: write_to_trj
use bulirsch_stoer, only: bs_largestep
use runge_kutta_4_nr, only: rk4nr_largestep
implicit none

! arguments to the subroutine:
real(dp),dimension(3,N_obj),intent(inout) :: X  ! spatial coordinates of all objects (m)
real(dp),dimension(3,N_obj),intent(inout) :: V  ! velocity components of all objects (m/s)

! internal variables:
integer(st) :: restart_status   ! IO status of the restart file
integer(st) :: i                ! loop index
integer(st) :: N_ok             ! number of successful BS steps per large step
integer(st) :: N_ok_tot         ! total number of successful BS steps
integer(st) :: N_fail           ! number of BS steps per large step, where the stepsize had to be reduced
integer(st) :: N_fail_tot       ! total number of failed BS steps
integer(st) :: N_rksteps        ! total number of RK steps per large step
integer(st) :: N_bssteps        ! total number of BS steps per large step
integer(st) :: N_bstotal        ! overall number of BS steps
integer(st) :: N_rktotal        ! overall number of RK steps
integer(st) :: N_smallsteps     ! number of BS substeps per large step
integer(st) :: N_smalltotal     ! total number of BS substeps
real(dp) :: timestep            ! initial stepsize for the next large step (s)
real(dp) :: next_step           ! suggested next stepsize returned by bs_largestep (s)
real(dp) :: start_stepcpu       ! cpu time before one large step
real(dp) :: end_stepcpu         ! cpu time after one large step


timestep = init_step
N_ok_tot = 0
N_fail_tot = 0
N_bstotal = 0
N_rktotal = 0
N_smalltotal = 0

if (do_bs) then
    write(output,'("                          ***************************************")')
    write(output,'("                          * USING THE BULIRSCH-STOER INTEGRATOR *")')
    write(output,'("                          ***************************************")')
    write(output,*)
    write(output,'("       elapsed time      large steps      BS steps      small steps      cpu time [ms]")')
    write(output,'("                         good    bad")')
    if (verbose) then
        write(*,'(" Starting propagation with the Bulirsch-Stoer integrator")')
    endif
else
    write(output,'("                     *************************************************")')
    write(output,'("                     * USING THE RUNGE-KUTTA INTEGRATOR OF 4TH ORDER *")')
    if (do_rk4mid) then
        write(output,'("                     *               WITH MIDPOINT STEP              *")')
    else
        write(output,'("                     *                WITH THIRDS STEP               *")')
    endif
    write(output,'("                     *************************************************")')
    write(output,*)
    write(output,'("       elapsed time      large steps              RK steps               cpu time [ms]")')
    write(output,'("                         good    bad")')
    if (verbose) then
        write(*,'(" Starting propagation with the Runge-Kutta integrator of fourth order")')
    endif
endif
flush(output)

! here comes the real propagation loop:
do
    call write_to_trj(elapsed_time, X, V)

    if (verbose) then
        call write_status(real(elapsed_time, dp), real(tfinal, dp))
    endif

    if (do_restart) then
        open(unit=restart,file=restart_file,status="replace",action="write",iostat=restart_status)
        write(restart,*) elapsed_time
        do i = 1, N_obj
            write(restart,*) mass(i), X(1,i), X(2,i), X(3,i), V(1,i), V(2,i), V(3,i)
        enddo
        close(unit=restart)
    endif

    if ((elapsed_time >= tfinal) .or. (tfinal - elapsed_time < 1.0_dp)) exit

    if (do_steps) then
        write(steps,'("# elapsed time:    ", es16.9)') elapsed_time
        write(steps,'("# trying propagation with step:   ", es16.9)') timestep
        if (.not. do_bs) then
            write(steps,'("#        Time             Stepsize           Error")')
        endif
    endif

    call cpu_time(start_stepcpu)
    if (do_bs) then
        call bs_largestep(timestep, next_step, X, V, N_ok, N_fail, N_bssteps, N_smallsteps)
    else
        call rk4nr_largestep(timestep, next_step, X, V, N_ok, N_fail, N_rksteps)
    endif
    call cpu_time(end_stepcpu)

    if (underflow) exit

    timestep = next_step
    N_ok_tot = N_ok_tot + N_ok
    N_fail_tot = N_fail_tot + N_fail
    N_bstotal = N_bstotal + N_bssteps
    N_smalltotal = N_smalltotal + N_smallsteps
    N_rktotal = N_rktotal + N_rksteps

    if (do_steps) then
        write(steps,'("# successful / failed steps: ", i4, ",", i4)') N_ok, N_fail
        write(steps,'("#")')
        write(steps,'("######################################################################")')
        write(steps,'("#")')
    endif

    if (do_bs) then
        write(output,'("       ", es11.5, "      ", i5, "  ", i5, "         ", i5, "         ", i7, "       ", f8.1)') &
            elapsed_time, N_ok, N_fail, N_bssteps, N_smallsteps, (end_stepcpu-start_stepcpu)*1000.0_dp
    else
        write(output,'("       ", es11.5, "      ", i5, "  ", i5, "                 ", i5, "                ", f8.1)') &
            elapsed_time, N_ok, N_fail, N_rksteps, (end_stepcpu-start_stepcpu)*1000.0_dp
    endif
    flush(output)
enddo

write(output,*)
write(output,'("       **************************************************************")')
write(output,'("       *                 SUMMARY OF THE CALCULATION                 *")')
write(output,'("       *                                                            *")')
if (do_bs) then
    write(output,'("       * total time      large steps      BS steps      small steps *")')
    write(output,'("       *                 good    bad                                *")')

    write(output,'("       *", es10.3, "    ", i7, i7, "      ", i8, "        ", i8, "  *")') &
        elapsed_time, N_ok_tot, N_fail_tot, N_bstotal, N_smalltotal
else
    write(output,'("       * total time      large steps              RK steps          *")')
    write(output,'("       *                 good    bad                                *")')

    write(output,'("       *", es10.3, "    ", i7, i7, "              ", i8, "          *")') &
        elapsed_time, N_ok_tot, N_fail_tot, N_rktotal
endif

write(output,'("       **************************************************************")')
write(output,*)

if (do_steps) then
    write(steps,'("# Simulation done!")')
    write(steps,'("# Total number of successful / failed steps: ",2i7)') N_ok_tot, N_fail_tot
endif

end subroutine propagate

end module propagation
