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
module astronx_utils
!
! This module contains utility subroutines of astronx:
!
!      routine                  purpose
!     ---------                ---------
!  scale_error           scale the absolute error of the extrapolation with the radius of gyration
!  write_to_trj          write the current state of the system to the binary trajectory
!  show_input_parameters write all input parameters and their values to the output file
!  shiftcog              shift the centre of gravity into the origin
!  shiftmom              shift the total linear momentum to zero
!
contains

!##################################################################################################
!##################################################################################################

subroutine scale_error(X_new, V_new, dX, dV, delta)
!
! This subroutine will calculate the normalized error estimate of the current Bulirsch-Stoer
! Iteration based on the absolute error estimate returned from the extrapolation routine. The
! error will be scaled based on the radius of gyration. This acts as an approximate measure
! for the size of the system.
!
use types
use shared_data, only: small
use common_utils, only: centre_of_gravity, radius_of_gyration
implicit none


! arguments to the routine:
real(ep),intent(in),dimension(:,:) :: X_new
real(ep),intent(in),dimension(:,:) :: V_new
real(ep),intent(in),dimension(:,:) :: dX
real(ep),intent(in),dimension(:,:) :: dV
real(ep),intent(out) :: delta

! internal variables:
real(ep),dimension(size(X_new,1),3) :: dX_scal
real(ep),dimension(size(X_new,1),3) :: dV_scal
real(ep) :: gyrate
real(ep) :: V_avg


call radius_of_gyration(X_new, gyrate)

V_avg = sum(abs(V_new)) / real(3*size(X_new,1),ep)

dX_scal = abs(dX / gyrate)
dV_scal = abs(dV / V_avg)

delta = (sum(dX_scal) + sum(dV_scal)) / real(6*size(X_new,1),ep)


end subroutine scale_error

!##################################################################################################
!##################################################################################################

subroutine write_to_trj(time, X, V)
!
! This subroutine will write one frame to the trajectory consisting of the elapsed time and
! the contents of the arrays X and V.
!
use types
use shared_data, only: bin_trj, trajectory
use input_module, only: do_texttrj
implicit none


! arguments to the routine:
real(dp),intent(in) :: time                 ! the current time
real(dp),intent(in),dimension(:,:) :: X     ! the positions
real(dp),intent(in),dimension(:,:) :: V     ! the velocities

! internal variables:
integer(st) :: i                            ! loop index

! write to the text trajectory if requested:
if (do_texttrj) then
    write(trajectory,100) time, (X(i,1), X(i,2), X(i,3), i=1,size(X,1))
    100 format (' ', 300es22.14)
endif


write(bin_trj) time
do i = 1, size(X,1)
    write(bin_trj) X(i,1), X(i,2), X(i,3)
enddo
do i = 1, size(X,1)
    write(bin_trj) V(i,1), V(i,2), V(i,3)
enddo


end subroutine write_to_trj

!##################################################################################################
!##################################################################################################

subroutine show_input_parameters
!
! this routine will write the values of all input parameters to the output file
!
use input_module
use shared_data, only: output
implicit none


write(output,*) "---------------------"
write(output,*) "SIMULATION PARAMETERS"
write(output,*) "---------------------"
write(output,*) ""
write(output,111) eps
111 format (' ', "eps         ", es11.3)
write(output,112) tfinal
112 format (' ', "tfinal      ", es11.3)
write(output,113) write_step
113 format (' ', "tout        ", es11.3)
write(output,114) init_step
114 format (' ', "init_step   ", es11.3)
write(output,115) maxsubstep
115 format (' ', "maxsubstep  ", i3)
write(output,116) thres
116 format (' ', "inc_thres   ", i3)
write(output,117) min_step
117 format (' ', "min_step    ", es11.3)
write(output,118) maxinc
118 format (' ', "maxinc      ", es11.3)
write(output,119) redmin
119 format (' ', "redmin      ", es11.3)
write(output,120) redmax
120 format (' ', "redmax      ", es11.3)
if (shift_cog) then
    write(output,*) "shift_cog     yes"
else
    write(output,*) "shift_cog     no"
endif
if (shift_mom) then
    write(output,*) "shift_mom     yes"
else
    write(output,*) "shift_mom     no"
endif
if (do_restart) then
    write(output,*) "restart       yes"
else
    write(output,*) "restart       no"
endif
if (do_steps) then
    write(output,*) "steps         yes"
else
    write(output,*) "steps         no"
endif
if (do_texttrj) then
    write(output,*) "text_trj      yes"
else
    write(output,*) "text_trj      no"
endif
if (do_unrestrictedprop) then
    write(output,*) "prop_type     unrestricted"
else
    write(output,*) "prop_type     normal"
endif

write(output,*) ""
write(output,*) ""
write(output,*) ""

end subroutine show_input_parameters

!##################################################################################################
!##################################################################################################

subroutine shiftcog(X)
!
! this routine will the positions of all objects by an equal amount, so that the centre of
! gravity coincides with the origin of the coordinate system
!
use types
use input_module, only: mass, total_mass
use common_utils, only: centre_of_gravity
implicit none


! arguments to the routine:
real(dp),dimension(:,:),intent(inout) :: X  ! the positions that will be shifted

! internal variables:
real(dp),dimension(3) :: cog    ! the position of the centre of gravity
integer(st) :: i, k             ! loop indices


call centre_of_gravity(X, cog, mass, total_mass)

forall ( i=1:size(X,1), k=1:3 )
    X(i,k) = X(i,k) - cog(k)
end forall


end subroutine shiftcog

!##################################################################################################
!##################################################################################################

subroutine shiftmom(V)
!
! this routine will add an additional contribution to all velocity components
! to make the total linear momentum of the system vanish
!
use types
use input_module, only: mass, total_mass
use common_utils, only: linear_momentum
implicit none


! arguments to the routine:
real(dp),dimension(:,:),intent(inout) :: V      ! the velocities to be shifted

! internal variables:
real(dp),dimension(3) :: mom    ! the linear momentum vector
integer(st) :: i, k             ! loop indices


call linear_momentum(V, mom, mass)

forall ( i=1:size(mass), k=1:3 )
    V(i,k) = V(i,k) - (mom(k) / total_mass)
end forall


end subroutine shiftmom

end module astronx_utils
