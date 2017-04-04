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
module input_module
use iso_fortran_env, only: int32, real64, error_unit
use globalmod, only: input
implicit none


! first the data to be declared:
integer(int32),protected :: N_obj                               ! number of objects to be simulated
integer(int32),protected :: maxsubstep                          ! maximum number of substeps in one BS-step
integer(int32),protected :: inc_thres                           ! number of substeps below which the stepsize will be increased
integer(int32),protected :: ndigit                              ! number of significant digits in text trajectory
integer(int32),protected :: nstep                               ! number of steps in rk4fix
integer(int32),protected :: int_type                            ! bs (1), rkqs (2), rk4fixm (3), rk4fixt (4)
real(real64),protected :: tfinal                                ! the total length of the simulation (s)
real(real64),protected :: tout                                  ! intervall of successive writes to the trajectory
real(real64),protected :: eps                                   ! the error tolerance for the propagation
real(real64),protected :: eps_thres                             ! error scaling threshold for increasing the stepsize in rk-runs
real(real64),protected :: min_step                              ! minimum timestep
real(real64),protected :: maxinc                                ! max. factor by which to increase stepsize
real(real64),protected :: redmin                                ! minimum factor for stepsize reduction
real(real64),protected :: redmax                                ! maximum reduction factor
real(real64),protected :: init_step                             ! the initial value for the timestep
real(real64),protected :: total_mass                            ! the total mass of the system (kg)
real(real64),dimension(:),allocatable,protected :: mass         ! masses of the objects (kg)
character(len=30),dimension(:),allocatable,protected :: names   ! the names of the objects
character(len=:),allocatable,protected :: name_directory        ! directory that will be created for the output
character(len=:),allocatable,protected :: input_file            ! name of the input file
logical,protected :: do_restart                                 ! create a restart file at every write step
logical,protected :: do_steps                                   ! create a file containing information about the steps
logical,protected :: do_txttrj                                  ! write a trajectory file in text form
logical,protected :: do_unResProp                               ! propagate without caring about the write step
logical,protected :: do_overwrite                               ! overwrite the contents of the name-directory
logical,protected :: shift_cog                                  ! wether or not to shift the cog at the beginning
logical,protected :: shift_mom                                  ! wether or not to compensate for cog motion
logical,protected :: verbose                                    ! wether of not to write info about the propagation to terminal
logical,protected :: acc_test                                   ! wether we are actually using acc_test
character(len=14),protected :: format_string                    ! format string to use for the text trajectory

contains

!##################################################################################################

subroutine process_cmd_arguments
!
! This routine reads the command line arguments of the program and
! sets the flags accordingly
!
implicit none

! internal variables:
integer(int32) :: i             ! loop index
integer(int32) :: name_length   ! length of the input file name
character(len=256) :: argument  ! command line argument template
logical :: file_exists          ! to check if the file exists


! see if there are any arguments:
if (command_argument_count() == 0) then
    write(error_unit,*) "ERROR: No arguments given. We need an input file."
    stop
endif

! loop over the arguments
do i = 1, command_argument_count()
    call get_command_argument(i, argument)
    if (trim(argument) == "-w") then
        do_overwrite = .true.
    elseif (trim(argument) == "-v") then
        verbose = .true.
    elseif (trim(argument) == "-t") then
        acc_test = .true.
    else
        if (.not. allocated(input_file)) then
            inquire(file=trim(argument), exist=file_exists)
            if (.not. file_exists) then
                write(error_unit,*) "ERROR: Input file ", trim(argument), " not found."
                stop
            else
                name_length = len(trim(argument))
                allocate(character(len=name_length) :: input_file)
                input_file = trim(argument)
            endif
        endif
    endif
enddo

end subroutine process_cmd_arguments

!##################################################################################################
!##################################################################################################

subroutine read_input(X, V)
!
! the purpose of this subroutine is to read the simulation parameters from the free-format input
! file and setting the corresponding variables in the module shared_data. All input parameters,
! which are not stored in shared_data are made available to the main program via the calling arguments.
!
use globalmod, only: input
implicit none


! arguments to the routine:
real(real64),dimension(:),allocatable,intent(out) :: X  ! spatial coordinates of all objects (m)
real(real64),dimension(:),allocatable,intent(out) :: V  ! velocity components of all objects (m/s)

! internal variables:
integer(int32) :: fw                        ! field width for writing to the text trajectory
character(len=40) :: keyword                ! dummy variable for the first column of the input file
character(len=40) :: paraValue              ! the value of the parameter to be read
character(len=200) :: input_buffer          ! one line of the input file, that is checked for keywords
character(len=200) :: raw_buffer            ! input buffer without comments removed
integer(int32) :: input_status              ! IO status of the input file ( if /= 0 -> trouble!!)
integer(int32) :: Nlines                    ! number of lines in the input file
integer(int32) :: readstatus                ! IO status for the read operation
integer(int32) :: i                         ! loop counting index
integer(int32) :: coord_start, coord_end    ! first and last line of the coordinates-block in the input file
integer(int32) :: namedir_length            ! length of the string describing the name directory
logical :: test_name = .false.              !
logical :: test_tfinal = .false.            ! parameters in the input file
logical :: test_tout = .false.              !
logical :: test_tinit = .false.             !


! open the input file:
open(newunit=input,file=input_file,status="old",action="read",iostat=input_status)
if (input_status /= 0) then
    write(error_unit,*) "ERROR: failed to open input file: '", trim(input_file), "', status: ", input_status
    stop
endif

! determine the number of lines in the input file:
Nlines = 0
do
    read(input,'(a256)',iostat=readstatus) input_buffer
    if(readstatus /= 0) exit
    Nlines = Nlines + 1
enddo
rewind(input)

! set the default values for all input parameters:
maxsubstep = 12
inc_thres = 8
ndigit = 10
nstep = 5
int_type = 1
eps = 1.0e-6_real64
eps_thres = 0.9_real64
min_step = 100.0_real64
maxinc = 10.0_real64
redmin = 0.9_real64
redmax = 0.01_real64
shift_cog = .false.
shift_mom = .false.
do_restart = .false.
do_steps = .false.
do_txttrj = .false.
do_unResProp = .false.
format_string = '(1x, 3es18.10)'

! find and read the parameters from the input file and set the test-variables:
do i=1 , Nlines

    ! first, read one line of the input file:
    read(input,'(a256)') raw_buffer

    ! check for a comment sign and trim the string:
    if (index(raw_buffer,"#") /= 0) then
        input_buffer = raw_buffer(1:(index(raw_buffer,"#")-1))
    else
        input_buffer = raw_buffer
    endif

    ! check whether we have a boolean keyword or a parameter specification:
    if (index(input_buffer, "=") /= 0) then     ! parameter
        keyword = trim(input_buffer(1:index(input_buffer, "=") - 1))
        paraValue = trim(input_buffer(index(input_buffer, "=") + 1:))

        ! get the parameters:
        if (index(keyword, "name_dir") /= 0) then
            test_name = .true.
            namedir_length = len(trim(paraValue))
            allocate(character(len=namedir_length) :: name_directory)
            read(paraValue,*) name_directory
        else if (index(keyword, "tfinal") /= 0) then
            test_tfinal = .true.
            read(paraValue,*) tfinal
        else if (index(keyword, "tout") /= 0) then
            test_tout = .true.
            read(paraValue,*) tout
        else if (index(keyword, "eps") /= 0) then
            read(paraValue,*) eps
        else if (index(keyword, "eps_thres") /= 0) then
            read(paraValue,*) eps_thres
        else if (index(keyword, "init_step") /= 0) then
            test_tinit = .true.
            read(paraValue,*) init_step
        else if (index(keyword, "maxsubstep") /= 0) then
            read(paraValue,*) maxsubstep
        else if (index(keyword, "inc_thres") /= 0) then
            read(paraValue,*) inc_thres
        else if (index(keyword, "min_step") /= 0) then
            read(paraValue,*) min_step
        else if (index(keyword, "maxinc") /= 0) then
            read(paraValue,*) maxinc
        else if (index(keyword, "redmin") /= 0) then
            read(paraValue,*) redmin
        else if (index(keyword, "redmax") /= 0) then
            read(paraValue,*) redmax
        else if (index(keyword, "nstep") /= 0) then
            read(paraValue,*) nstep
        else if (index(keyword, "prop_type") /= 0) then
            if (index(paraValue, "unrestricted") /= 0) then
                do_unResProp = .true.
            else if (index(paraValue, "normal") /= 0) then
                do_unResProp = .false.
            else
                write(error_unit,*) "Invalid value for 'prop_type', using default: normal"
                do_unResProp = .false.
            endif
        else if (index(keyword, "int_type") /= 0) then
            if (index(paraValue, "bs") /= 0) then
                int_type = 1
            else if (index(paraValue, "rkqs") /= 0) then
                int_type = 2
            else if (index(paraValue, "rk4fixm") /= 0) then
                int_type = 3
            else if (index(paraValue, "rk4fixt") /= 0) then
                int_type = 4
            else
                write(error_unit,*) "Invalid value for 'integrator', using default: bs"
                int_type = 1
            endif
        else if (index(keyword, "ndigit") /= 0) then
            read(paravalue,*) ndigit
            fw = ndigit + 8
            if (fw < 10) then
                write(format_string,'("(1x,  3es",i1,".",i1,")")') fw, ndigit
            else if (ndigit > 10) then
                write(format_string,'("(1x,3es",i2,".",i2,")")') fw, ndigit
            else
                write(format_string,'("(1x, 3es",i2,".",i1,")")') fw, ndigit
            endif
        endif
    else    ! boolean
        keyword = trim(input_buffer)
        if (index(keyword, "shift_cog") /= 0) then
            shift_cog = .true.
        else if (index(keyword, "shift_mom") /= 0) then
            shift_mom = .true.
        else if (index(keyword, "restart") /= 0) then
            do_restart = .true.
        else if (index(keyword, "steps") /= 0) then
            do_steps = .true.
        else if (index(keyword, "text_trj") /= 0) then
            do_txttrj = .true.
        endif
    endif

enddo
rewind(input)

! check if all mandatory parameters were there:
if (.not. test_name) then
    write(error_unit,*) "INPUT ERROR: No name directory specified. Exiting"
    stop
endif
if (.not. test_tfinal) then
    write(error_unit,*) "INPUT ERROR: Propagation time not specified. Exiting"
    stop
endif
if (.not. test_tout) then
    write(error_unit,*) "INPUT ERROR: Printout step not specified. Exiting"
    stop
endif

if ((int_type == 3 .or. int_type == 4) .and. do_steps) then
    write(error_unit,*) "Steps file makes no sense with a fixed-step integrator..."
    do_steps = .false.
endif

! if the initial step has not been specified, set to tout as the default:
if (.not. test_tinit) then
    init_step = tout
endif

! find the first and last line of the coordinate block:
coord_start = 0
coord_end = 0
do
    read(input,'(a256)') input_buffer
    coord_start = coord_start + 1
    if (index(input_buffer,"begin_coords") /= 0) exit
enddo
rewind(input)
do
    read(input,'(a256)') input_buffer
    coord_end = coord_end + 1
    if (index(input_buffer,"end_coords") /= 0) exit
enddo
rewind(input)
N_obj = coord_end - coord_start - 1

! if everything is fine, read the initial coordinates:

! allocate the memory:
allocate(names(N_obj))
allocate(mass(N_obj))
allocate(X(3 * N_obj))
allocate(V(3 * N_obj))

! move to start of coordinate block:
do i=1 , coord_start
    read(input,*)
enddo

! read coordinates:
do i=1, N_obj
    read(input,*,iostat=readstatus) names(i), mass(i), &
                                    X(3 * (i-1) + 1), X(3 * (i-1) + 2), X(3 * (i-1) + 3), &
                                    V(3 * (i-1) + 1), V(3 * (i-1) + 2), V(3 * (i-1) + 3)
    if (readstatus /= 0) then
        write(error_unit,*) "Error reading coordinate of object ", i, " Exiting."
        stop
    endif
enddo

close(input)


end subroutine read_input

!##################################################################################################

subroutine initialize
!
! this routine creates the name directory, opens the requested output files
! and initializes the counters
!
use globalmod, only: bin_trj, bin_trj_file, txt_trj, txt_trj_file, output, output_file, &
                     restart, restart_file, steps, steps_file
use counters
implicit none


integer(int32) :: i     ! loop index
integer(int32) :: stat  ! status of system operation
logical :: test_name    ! to test for the name directory


! calculate the total system mass:
total_mass = 0.0_real64
do i=1, N_obj
    total_mass = total_mass + mass(i)
enddo

! create the name directory and move there:
inquire(file=name_directory, exist=test_name)
if (test_name) then
    if (.not. do_overwrite) then
        write(error_unit,*) "ERROR: The name directory '", trim(name_directory), "' already exists."
        write(error_unit,*) "Use the option '-w' to overwrite previous results."
        stop
    endif
else
    call system(("mkdir "//name_directory), stat)
    if (stat /= 0) then
        write(error_unit,*) "ERROR: Could not create directory '", trim(name_directory), "'. Exiting."
        stop
    endif
endif
call chdir(name_directory)

! open the output files:
open(newunit=output,file=output_file,status="replace",action="write",iostat=stat)
if (stat /= 0) then
    write(error_unit,'("ERROR: failed to open output file, status: ",i3)') stat
    stop
endif

open(unit=bin_trj,file=bin_trj_file,status="replace",action="write",iostat=stat,form="unformatted")
if (stat /= 0) then
    write(error_unit,'("ERROR: failed to open binary trajectory file, status: ",i3)') stat
    write(output,'("ERROR: failed to open binary trajectory file, status: ",i3)') stat
    stop
endif

if (do_txttrj) then
    open(newunit=txt_trj,file=txt_trj_file,status="replace",action="write",iostat=stat)
    if (stat /= 0) then
        write(error_unit,'("ERROR: failed to open text trajectory file, status: ",i3)') stat
        write(output,'("ERROR: failed to open text trajectory file, status: ",i3)') stat
        stop
    endif
endif

if (do_restart) then
    open(newunit=restart,file=restart_file,status="replace",action="write",iostat=stat)
    if (stat /= 0) then
        write(error_unit,'("ERROR: failed to open restart file, status: ",i3)') stat
        write(output,'("ERROR: failed to open restart file, status: ",i3)') stat
        stop
    endif
endif

if (do_steps) then
    open(newunit=steps,file=steps_file,status="replace",action="write",iostat=stat)
    if (stat /= 0) then
        write(error_unit,'("ERROR: failed to open steps file, status: ",i3)') stat
        write(output,'("ERROR: failed to open steps file, status: ",i3)') stat
        stop
    endif
endif

! initialize routine counters:
N_acceleration = 0
N_bs_largestep = 0
N_bs_onestep = 0
N_bs_substeps = 0
N_extrapolate = 0

end subroutine initialize

end module input_module
