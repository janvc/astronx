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
module enprop_module
!
! This module contains routines that will calculate energetic properties of a
! system from a trajectory file
!
contains

!##################################################################################################
!##################################################################################################

subroutine enprop(N_frames)
!
! This is the main routine which organizes the calculation
!
use types
use extract_share
use extractor_utils, only: number_of_frames
implicit none


! arguments to the routine:
integer(st),intent(in) :: N_frames          ! number of (intact) frames in the trajectory

! internal variables:
integer(st) :: N_obj                        ! number of objects in the system
integer(st) :: i, j                         ! loop indices
integer(st) :: ioerror                      ! status variable for the file opening operation
integer(st),parameter :: ener = 3           ! IO unit of the energy file
character(len=20) :: energy_file            ! output file with the energy data
real(dp) :: E_kin                           ! total kinetic energy of the system
real(dp) :: E_pot                           ! total potential energy of the system
real(dp) :: elapsed_time                    ! the time
real(dp),dimension(:),allocatable :: mass   ! array with masses of the objects
real(dp),dimension(:,:),allocatable :: X    ! the positions of the objects
real(dp),dimension(:,:),allocatable :: V    ! the velocities of the objects


! start reading the trajectory file. First some information about the system:
rewind(trj)
read(trj) N_obj         ! Number of objects
allocate(mass(N_obj))
do i = 1, N_obj         ! read the masses. Now the record counter should stand at the first
    read(trj) mass(i)   ! frame of the actual trajectory.
enddo

allocate(X(N_obj,3))
allocate(V(N_obj,3))

! open the output file for writing:
write(*,*) ""
write(*,*) "    Enter the name of the file to write the energies to (max 20 characters)"
do
    write(*,'(">>> ")',advance="no")
    read(*,*) energy_file
    open(unit=ener,file=energy_file,status="new",action="write",iostat=ioerror)
    if (ioerror == 0) exit
    write(*,*) ""
    write(*,*) "    ERROR opening the the file '", trim(energy_file), "'. Maybe it already exists."
    write(*,*) "    Please select another file."
enddo


! this is the main processing loop:
do i = 1, N_frames
    read(trj,iostat=ioerror) elapsed_time
!    if (ioerror == -1) exit
    do j = 1, N_obj
        read(trj) X(j,1), X(j,2), X(j,3)
    enddo
    do j = 1, N_obj
        read(trj) V(j,1), V(j,2), V(j,3)
    enddo
    call pot_ener(N_obj, mass, X, E_pot)
    call kin_ener(N_obj, mass, V, E_kin)
    write(ener,*) elapsed_time, E_pot, E_kin, E_pot + E_kin
enddo

close(unit=ener)

deallocate(mass)
deallocate(X)
deallocate(V)
write(*,*) ""
write(*,*) "    The energy has been calculated."

end subroutine enprop

!##################################################################################################
!##################################################################################################

subroutine pot_ener(N_obj, mass, X, E_pot)
!
! This routine calculates the potential energy of a configuration from the positions
!
use types
use shared_data, only: G
implicit none


! arguments to the routine:
integer(st),intent(in) :: N_obj                 ! the number of objects
real(dp),dimension(N_obj),intent(in) :: mass    ! the masses of the objects
real(dp),dimension(N_obj,3),intent(in) :: X     ! the positions of the objects
real(dp),intent(out) :: E_pot                   ! the total potential energy

! internal variables:
integer(st) :: i, j                     ! loop indices
real(dp) :: E_pot_local                 ! contribution of a single object pair
real(dp),dimension(N_obj,N_obj) :: R    ! distances between all objects


E_pot = 0.0_dp

! calculate the distances:
forall(i=1:N_obj, j=1:N_obj, j>=i+1)
    R(i,j) = sqrt(((X(j,1)-X(i,1))*(X(j,1)-X(i,1))) &
                + ((X(j,2)-X(i,2))*(X(j,2)-X(i,2))) &
                + ((X(j,3)-X(i,3))*(X(j,3)-X(i,3))))
    R(j,i) = R(i,j)
end forall

! calculate the potential energy:
do i = 1, N_obj
    do j = i+1, N_obj
        E_pot_local = -G * (mass(i) * mass(j)) / R(i,j)
        E_pot = E_pot + E_pot_local
    enddo
enddo

end subroutine pot_ener

!##################################################################################################
!##################################################################################################

subroutine kin_ener(N_obj, mass, V, E_kin)
!
! This routine calculates the kinetic energy of a configuration from the velocities
!
use types
implicit none

! arguments to the routine:
integer(st),intent(in) :: N_obj                 ! the number of objects
real(dp),dimension(N_obj),intent(in) :: mass    ! the masses of the objects
real(dp),dimension(N_obj,3),intent(in) :: V     ! the velocities of the objects
real(dp),intent(out) :: E_kin                   ! the total kinetic energy

! internal variables:
integer(st) :: i, j         ! loop indices
real(dp) :: E_kin_local     ! contribution of one object


E_kin = 0.0_dp

do i = 1, N_obj
    do j = 1, 3
        E_kin_local = 0.5_dp * mass(i) * (V(i,j)*V(i,j))
        E_kin = E_kin + E_kin_local
    enddo
enddo

end subroutine kin_ener

end module enprop_module
