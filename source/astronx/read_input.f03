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
module input_module
use types
implicit none


! first the data to be declared:
integer(st),protected :: N_obj                          ! number of objects to be simulated
integer(st),protected :: maxsubstep                     ! maximum number of substeps in one BS-step
integer(st),protected :: thres                          ! number of substeps below which the stepsize will be increased
real(dp),protected :: tfinal                            ! the total length of the simulation (s)
real(dp),protected :: write_step                        ! intervall of successive writes to the trajectory
real(dp),protected :: eps                               ! the error tolerance for the propagation
real(dp),protected :: min_step                          ! minimum timestep
real(dp),protected :: maxinc                            ! max. factor by which to increase stepsize
real(dp),protected :: redmin                            ! minimum factor for stepsize reduction
real(dp),protected :: redmax                            ! maximum reduction factor
real(dp),protected :: init_step                         ! the initial value for the timestep
real(dp),protected :: total_mass                        ! the total mass of the system (kg)
real(dp),dimension(:),allocatable,protected :: mass     ! masses of the objects (kg)
real(dp),dimension(:,:),allocatable,protected :: mass_2 ! mass-products of all object-pairs (kg^2)
character(len=100),protected :: name_directory          ! directory that will be created for the output
logical,protected :: do_restart                         ! create a restart file at every write step
logical,protected :: do_steps                           ! create a file containing information about the steps
logical,protected :: do_texttrj                         ! write a trajectory file in text form
logical,protected :: do_unrestrictedprop                ! propagate without caring about the write step
logical,protected :: do_bs                              ! use BS instead of RK as integrator
logical,protected :: do_overwrite                       ! overwrite the contents of the name-directory
logical,protected :: shift_cog                          ! wether or not to shift the cog at the beginning
logical,protected :: shift_mom                          ! wether or not to compensate for cog motion
logical,protected :: verbose                            ! wether of not to write info about the propagation to terminal

contains

subroutine process_cmd_arguments
!
! This routine reads the command line arguments of the program and
! sets the flags accordingly
!
use types
!use input_module, only: do_overwrite
implicit none

! internal variables:
integer(st) :: i                ! loop index
character(len=128) :: argument    ! command line argument template


! see if there are any arguments:
if (command_argument_count() == 0) then
    write(*,*) "ERROR: No arguments given. We need an input file."
    stop
endif

! loop over the arguments
do i = 1, command_argument_count()
    call get_command_argument(i, argument)
    if (trim(argument) == "-w") then
        do_overwrite = .true.
    else
        ! check, if the argument is the input file...
    endif
enddo

end subroutine process_cmd_arguments

!##################################################################################################
!##################################################################################################

subroutine read_input(input_file, X, V)
!
! the purpose of this subroutine is to read the simulation parameters from the free-format input
! file and setting the corresponding variables in the module shared_data. All input parameters,
! which are not stored in shared_data are made available to the main program via the calling arguments.
!
use types
use shared_data, only: input, output
implicit none


! arguments to the routine:
character(len=*),intent(in) :: input_file             ! name of the input file to be read
real(dp),dimension(:,:),allocatable,intent(out) :: X  ! spatial coordinates of all objects (m)
real(dp),dimension(:,:),allocatable,intent(out) :: V  ! velocity components of all objects (m/s)

! internal variables:
character(len=20) :: keyword            ! dummy variable for the first column of the input file
character(len=100) :: input_buffer      ! one line of the input file, that is checked for keywords
character(len=100) :: raw_buffer        ! input buffer without comments removed
integer(st) :: input_status             ! IO status of the input file ( if /= 0 -> trouble!!)
integer(st) :: number_of_lines          ! number of lines in the input file
integer(st) :: readstatus               ! IO status for the read operation
integer(st) :: i, k                     ! loop counting indices
integer(st) :: coord_start, coord_end   ! first and last line of the coordinates-block in the input file
logical :: test_name = .false.          !\
logical :: test_N_obj = .false.         ! | logical variables to test for the existence of the corresponding
logical :: test_tfinal = .false.        ! | parameters in the input file
logical :: test_tout = .false.          !/
logical :: test_tinit = .false.         !


! open the input file:
open(unit=input,file=input_file,status="old",action="read",iostat=input_status)
if (input_status /= 0) then
    write(*,*) "ERROR: failed to open input file: '", trim(input_file), "', status: ", input_status
    stop
endif

! determine the number of lines in the input file:
number_of_lines = 0
do
    read(input,100,iostat=readstatus) input_buffer
    100 format (a100)
    if(readstatus /= 0) exit
    number_of_lines = number_of_lines + 1
enddo
rewind(input)

! set the default values for all input parameters:
maxsubstep = 12
thres = 8
eps = 1.0e-6_dp
min_step = 100.0_dp
maxinc = 10.0_dp
redmin = 0.9_dp
redmax = 0.01_dp
shift_cog = .false.
shift_mom = .false.
do_restart = .false.
do_steps = .false.
do_texttrj = .false.
do_unrestrictedprop = .false.
do_bs = .true.

! find and read the parameters from the input file and set the test-variables:
do i=1 , number_of_lines

    ! first, read one line of the input file:
    read(input,100) raw_buffer

    ! check for a comment sign and trim the string:
    if (index(raw_buffer,"#") /= 0) then
        input_buffer = raw_buffer(1:(index(raw_buffer,"#")-1))
    else
        input_buffer = raw_buffer
    endif

    ! get the parameters:
    if (index(input_buffer,"name_dir") /= 0) then
        test_name = .true.
        read(input_buffer,*) keyword, name_directory
    else if (index(input_buffer,"n_obj") /= 0) then
        test_N_obj = .true.
        read(input_buffer,*) keyword, N_obj
    else if (index(input_buffer,"eps") /= 0 .and. index(input_buffer,"steps") == 0) then
        read(input_buffer,*) keyword, eps
    else if (index(input_buffer,"tfinal") /= 0) then
        test_tfinal = .true.
        read(input_buffer,*) keyword, tfinal
    else if (index(input_buffer,"tout") /= 0) then
        test_tout = .true.
        read(input_buffer,*) keyword, write_step
    else if (index(input_buffer,"init_step") /= 0) then
        test_tinit = .true.
        read(input_buffer,*) keyword, init_step
    else if (index(input_buffer,"maxsubstep") /= 0) then
        read(input_buffer,*) keyword, maxsubstep
    else if (index(input_buffer,"inc_thres") /= 0) then
        read(input_buffer,*) keyword, thres
    else if (index(input_buffer,"min_step") /= 0) then
        read(input_buffer,*) keyword, min_step
    else if (index(input_buffer,"maxinc") /= 0) then
        read(input_buffer,*) keyword, maxinc
    else if (index(input_buffer,"redmin") /= 0) then
        read(input_buffer,*) keyword, redmin
    else if (index(input_buffer,"redmax") /= 0) then
        read(input_buffer,*) keyword, redmax
    else if (index(input_buffer,"shift_cog") /= 0 .and. index(input_buffer,"yes") /= 0) then
        shift_cog = .true.
    else if (index(input_buffer,"shift_mom") /= 0 .and. index(input_buffer,"yes") /= 0) then
        shift_mom = .true.
    else if (index(input_buffer,"restart") /= 0 .and. index(input_buffer,"yes") /= 0) then
        do_restart = .true.
    else if (index(input_buffer,"steps") /= 0 .and. index(input_buffer,"yes") /= 0) then
        do_steps = .true.
    else if (index(input_buffer,"text_trj") /= 0 .and. index(input_buffer,"yes") /= 0) then
        do_texttrj = .true.
    else if (index(input_buffer,"prop_type") /= 0) then
        if (index(input_buffer,"unrestricted") /= 0) then
            do_unrestrictedprop = .true.
        else if (index(input_buffer,"normal") /= 0) then
            do_unrestrictedprop = .false.
        else
            write(output,*) "Invalid value for 'prop_type', using default: normal"
            do_unrestrictedprop = .false.
        endif
    else if (index(input_buffer,"integrator") /= 0) then
        if (index(input_buffer,"bs") /= 0) then
            do_bs = .true.
        else if (index(input_buffer,"rk") /= 0) then
            do_bs = .false.
        else
            write(output,*) "Invalid value for 'integrator', using default: bs"
            do_bs = .true.
        endif
!    else   !  need to improve this!
!        write(*,*) "Warning: unrecognized keyword in line ", i, " : ", input_buffer
    endif
enddo
rewind(input)

! check if all mandatory parameters were there:
if (.not. test_name) then
    write(*,*) "INPUT ERROR: No name directory specified. Exiting"
    stop
endif
if (.not. test_N_obj) then
    write(*,*) "INPUT ERROR: Number of objects not specified. Exiting"
    stop
endif
if (.not. test_tfinal) then
    write(*,*) "INPUT ERROR: Propagation time not specified. Exiting"
    stop
endif
if (.not. test_tout) then
    write(*,*) "INPUT ERROR: Printout step not specified. Exiting"
    stop
endif

! if the initial step has not been specified, set to write_step as the default:
if (.not. test_tinit) then
    init_step = write_step
endif

! find the first and last line of the coordinate block:
coord_start = 0
coord_end = 0
do
    read(input,100) input_buffer
    coord_start = coord_start + 1
    if (index(input_buffer,"begin_coords") /= 0) exit
enddo
rewind(input)
do
    read(input,100) input_buffer
    coord_end = coord_end + 1
    if (index(input_buffer,"end_coords") /= 0) exit
enddo
rewind(input)

! check for consistency:
if ((coord_end-coord_start-N_obj-1) /= 0) then
    write(*,*) "INPUT ERROR: Number of objects (", N_obj, ") does not match number of lines in coordinate block. Exiting."
    stop
endif


! if everything is fine, read the initial coordinates:

! allocate the memory:
allocate(mass(N_obj))
allocate(mass_2(N_obj,N_obj))
allocate(X(N_obj,3))
allocate(V(N_obj,3))

! move to start of coordinate block:
do i=1 , coord_start
    read(input,*)
enddo

! read coordinates:
do i=1 , N_obj
    read(input,*,iostat=readstatus) mass(i), X(i,1), X(i,2), X(i,3), V(i,1), V(i,2), V(i,3)
    if (readstatus /= 0) then
        write(*,*) "Error reading coordinate no. ", i, " Exiting."
        stop
    endif
enddo

close(input)

! set up the mass-squared-array to save a multiplication in each force calculation:
forall ( i=1:N_obj, k=1:N_obj, i/=k )
    mass_2(i,k) = mass(i) * mass(k)
end forall

! calculate the total system mass:
total_mass = 0.0_dp
do i=1, N_obj
    total_mass = total_mass + mass(i)
enddo

! create the name directory:
call system(("mkdir "//name_directory), input_status)
if (input_status /= 0) then
    write(*,*) "ERROR: Could not create directory '", trim(name_directory), "'. Exiting."
    stop
endif

end subroutine read_input

end module input_module