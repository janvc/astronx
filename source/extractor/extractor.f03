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
program extractor
!
! This program will read the trajectory file written by AstronX
! and extract various physical properties of the system from it.
! This file is basically just a collection of calls to the respective subroutines.
!
use types
use extract_share
use extractor_utils, only: number_of_frames
use bin2text_module
use enprop_module, only: enprop
implicit none


! "Data Dictionary":
integer(st) :: action           ! to define the action to be taken with the trajectory
integer(st) :: trj_stat         ! IO status of the trajectory
integer(st) :: N_obj            ! number of objects
integer(st) :: N_frames         ! number of frames
character(len=40) :: trj_file   ! name of the binary trajectory file


! greet the user:
!write(*,*) ""
!write(*,*) "--------------------------------------------------------------------------------"
write(*,*) ""
write(*,*) "    Welcome to the Extractor!"
write(*,*) ""
write(*,*) "    The extractor is part of the Astronx package."
write(*,*) "    Copyright 2012-2013 by Jan von Cosel"
write(*,*) "    Astronx is free software."
write(*,*) "    You can redistribute it and/or modify it under the terms of the GNU General"
write(*,*) "    Public License as published by the Free Software Foundation, either version 3"
write(*,*) "    of the License, or (at your option) any later version."
write(*,*) ""
write(*,*) "    Astronx comes with absolutely no warranty."

! open the trajectory:
write(*,*) "    --------------------------------------------------------------------------------"
write(*,*) ""
write(*,*) "    Enter the name of the trajectory file to analyze"
write(*,*) ""
write(*,'(">>> ")',advance="no")
!101 format (">>> ")

do
    read(*,*) trj_file
    open(unit=trj,file=trj_file,status="old",action="read",form="unformatted",iostat=trj_stat)
    if (trj_stat == 0) exit
    write(*,*) "    ERROR: failed to open trajectory file '", trim(trj_file), "'. Enter the correct filename."
    write(*,*) ""
    write(*,'(">>> ")',advance="no")
enddo


read(trj) N_obj
call number_of_frames(N_frames)

write(*,102) N_frames, N_obj
102 format ("     The trajectory contains ", i5, " frames and ", i3, " objects.")
write(*,*) ""

! here comes the main interactive loop:
do
    write(*,*) "    ---------------------------------------------------------------------------------"
    write(*,*) ""
    write(*,*) "    What would you like to do?"
    write(*,*) ""
    write(*,*) "    (0)  Exit"
    write(*,*) "    (1)  Translate a binary trajectory into a textfile"
    write(*,*) "    (2)  Compare two trajectories"
    write(*,*) "    (3)  Calculate energetic properties from a trajectory"
    write(*,*) "    (4)  Calculate configurational properties from a trajectory"
    write(*,*) "    (5)  Translate the trajectory into a MOLDEN/VMD compatible format"
    write(*,*) ""
    write(*,'(">>> ")',advance="no")

    read(*,*) action

    if (action == 0) exit

    select case (action)
        case default
            write(*,*) "    ERROR: illegal input!"
        case (1)
            call bin2text  ! this routine will write (selected) objects to a textfile.
        case (2)
            write(*,*) "Sorry, this feature is not implemented yet"
        case (3)
            call enprop(N_frames)
        case (4)
            write(*,*) "Sorry, this feature is not implemented yet"
!            call confprop
        case (5)
            write(*,*) "Sorry, this feature is not implemented yet"
!            call trans_vmdmol
    end select
enddo

end program extractor
