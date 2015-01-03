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
module mkdist_module
!
! This module contains the routines to generate object distributions
!
contains

!##################################################################################################
!##################################################################################################

subroutine mkdist(dist_shape, pos_type, vel_type, mass_type, coord_sys, N_obj, max_mass, centre, X, V, mass)
!
! This routine organises the calculation and calls the respective subroutines.
!
use types
use common_utils, only: random_num
implicit none


! arguments to the routine:
integer(st),intent(in) :: dist_shape            ! the shape of the distribution
integer(st),intent(in) :: pos_type              ! the type of the position distribution
integer(st),intent(in) :: vel_type              ! the type of the velocity distribution
integer(st),intent(in) :: mass_type             ! the type of the mass distribution
integer(st),intent(in) :: coord_sys             ! the type of coordinate system to be used
integer(st),intent(in) :: N_obj                 ! the number of objects
real(dp),intent(in) :: max_mass                 ! the upper bound for the object mass
real(dp),intent(in),dimension(3) :: centre      ! the position of the centre of the distribution
real(dp),intent(in),dimension(N_obj,3) :: X     ! the position array
real(dp),intent(in),dimension(N_obj,3) :: V     ! the velocity array
real(dp),intent(in),dimension(N_obj) :: mass    ! the mass array


! internal variables:
integer(st) :: i                    ! loop index
real(dp),dimension(3) :: dist_size  ! the size of a cuboidal configuration


select case (dist_shape)
    case (1) ! cuboidal distribution
        ! get the size of the configuration:
        write(*,*) ""
        write(*,*) "    Enter extension of the configuration in the x-, y-, and z-direction:"
        write(*,*) ""
        write(*,'(">>> ")',advance="no")
        read(*,*) dist_size(1), dist_size(2), dist_size(3)

        ! create the mass distribution:
        do i = 1, N_obj
            mass(i) = random_num(mass_type) * max_mass
        enddo

        ! create the position distribution:
        do i = 1, N_obj
            do j = 1, 3
                X(i,j) = centre(j) + dist_size(j) * (random_num(pos_type) - 0.5_dp)
            enddo
        enddo

        ! create the velocity distribution:
        do i = 1, N_obj
            do j = 1, 3
                V(i,j) = centre(j) + dist_size(j) * (random_num(vel_type) - 0.5_dp)
            enddo
        enddo
    case (2) ! spherical distribution
        call mk_spher_dist
    case (3) ! rectangular distribution
        call mk_rec_dist
    case (4) ! circular distribution
        call mk_circ_dist(coordinates, N_obj, centre)
end select

end subroutine mkdist



end module mkdist_module
