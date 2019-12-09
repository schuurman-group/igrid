!######################################################################
! igrid: A program to sample initial positions and momenta using the
!        Wigner distribution function and a wavefunction corresponding
!        to a Hartree product of the eigenfunctions of 1D Hamiltonians.
!
!        Works in terms of mass- and frequency-scaled normal modes
!        and assumes a separable total Hamiltonian.
!
!        1D eigenfunctions are calculated using the Fourier grid
!        Hamiltonian (FGH) method of Balint-Kurti [J. Chem. Phys., 91,
!        3571 (1989)].
!######################################################################

program igrid

  use ioqc
  use igridglobal
  
  implicit none

!----------------------------------------------------------------------
! Open the input and log files
!----------------------------------------------------------------------
  call open_files_igrid

!----------------------------------------------------------------------
! Read the input file
!----------------------------------------------------------------------
  call read_input_igrid

!----------------------------------------------------------------------
! Determine the normal mode file type
!----------------------------------------------------------------------
  call freqtype

!----------------------------------------------------------------------
! Determine the no. atoms and allocate xcoo0 and related arrays
!----------------------------------------------------------------------
  call getdim

!----------------------------------------------------------------------
! Read the Cartesian coordinates
!----------------------------------------------------------------------
  call getxcoo0

!----------------------------------------------------------------------
! Determine the number of normal modes from the moment of intertia
! tensor and allocate associated arrays
!----------------------------------------------------------------------
  call getnmodes

!----------------------------------------------------------------------
! Read the normal modes, frequencies and symmetry labels
!----------------------------------------------------------------------
  call getmodes

!----------------------------------------------------------------------
! Create the transformation matrices
!----------------------------------------------------------------------
  call nm2xmat

!----------------------------------------------------------------------
! Read the names of the quantum chemistry output files
!----------------------------------------------------------------------
  call rdqcfilenames
  
contains

!######################################################################

  subroutine open_files_igrid
    
    use constants
    use channels
    use iomod

    implicit none

    integer :: ilbl
    logical :: found

!----------------------------------------------------------------------
! Exit if no input file has been given
!----------------------------------------------------------------------
    if (iargc().eq.0) then
       write(6,'(/,2x,a,/)') 'No input file has been given'
       stop
    endif
    
!----------------------------------------------------------------------
! Read the name of the input file and set the names of the log file
! operator file, and binary file
!----------------------------------------------------------------------
    call getarg(1,ain)

    ilbl=index(ain,'.inp')

    if (ilbl.ne.0) then
       alog=ain(1:ilbl-1)//'.log'
    else
       alog=trim(ain)//'.log'
       ain=trim(ain)//'.inp'
    endif

!----------------------------------------------------------------------
! Exit if the input file does not exist
!----------------------------------------------------------------------
    inquire(file=trim(ain),exist=found)
    
    if (.not.found) then
       write(6,'(/,2x,a,/)') 'The file '//trim(ain)//' does not exist'
       stop
    endif

!----------------------------------------------------------------------
! Open the input, log, operator and binary files
!----------------------------------------------------------------------
    iin=1
    open(iin,file=ain,form='formatted',status='old')
    
    ilog=2
    open(ilog,file=alog,form='formatted',status='unknown')

    return
    
  end subroutine open_files_igrid

!######################################################################

  subroutine read_input_igrid

    use constants
    use channels
    use iomod
    use parsemod
    use sysinfo
    use igridglobal
    
    implicit none

    integer :: i,k,k1,k2
    
!----------------------------------------------------------------------
! Set defaults
!----------------------------------------------------------------------
    ! Frequency file
    freqfile=''

    ! Set file
    setfile=''
    lsetfile=.false.

!----------------------------------------------------------------------
! Read the input file
!----------------------------------------------------------------------
    rewind(iin)

15  continue
    call rdinp(iin)

    i=0
    if (.not.lend) then
    
20     continue
       i=i+1

       if (keyword(i).eq.'$freqfile') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             freqfile=keyword(i)
          else
             goto 100
          endif

       else if (keyword(i).eq.'$qcfiles') then
          if (keyword(i+1).eq.'=') then
             ! Filenames are to be read from a set file
             lsetfile=.true.
             i=i+2
             setfile=keyword(i)
          else
             ! Filenames are given directly in the input file
             ! skip past for now and read these later
25           call rdinp(iin)
             if (keyword(1).ne.'$end') goto 25
          endif

       else
          ! Exit if the keyword is not recognised
          errmsg='Unknown keyword: '//trim(keyword(i))
          call error_control
       endif

       ! If there are more keywords to be read on the current line,
       ! then read them, else read the next line
       if (i.lt.inkw) then
          goto 20
       else
          goto 15
       endif
       
       ! Exit if a required argument has not been given with a keyword
100    continue
       errmsg='No argument given with the keyword '//trim(keyword(i))
       call error_control
       
    endif

!----------------------------------------------------------------------
! Make sure that all the required information has been given
!----------------------------------------------------------------------
    if (freqfile.eq.'') then
       errmsg='The name of the frequency calculation file has not &
            been given'
       call error_control
    endif
    
    return
    
  end subroutine read_input_igrid

!######################################################################

  subroutine rdqcfilenames

    use igridglobal
    
    implicit none

    if (lsetfile) then
       ! Read the QC ouput filenames from a set files
       call rdqcfilenames_setfile
    else
       ! Read the QC output filenames from the input file
       call rdqcfilenames_inpfile
    endif
    
    return
    
  end subroutine rdqcfilenames

!######################################################################

  subroutine rdqcfilenames_setfile

    use constants
    use channels
    use parsemod
    use iomod
    use igridglobal
    
    implicit none

    integer :: unit,ierr,i

!----------------------------------------------------------------------
! Open the set file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=setfile,form='formatted',status='old',iostat=ierr)

    if (ierr.ne.0) then
       errmsg='Error opening the set file: '//trim(setfile)
       call error_control
    endif

!----------------------------------------------------------------------
! First pass: determine the no. files and allocate arrays
!----------------------------------------------------------------------
    ! Determine the no. files
    nfiles=0
5   call rdinp(unit)
    if (.not.lend) then
       nfiles=nfiles+1
       goto 5
    endif

    ! Allocate arrays
    allocate(qcfiles(nfiles))
    qcfiles=''

!----------------------------------------------------------------------
! Second pass: read in the filenames
!----------------------------------------------------------------------
    rewind(unit)

    do i=1,nfiles
       call rdinp(unit)
       qcfiles(i)=keyword(1)
    enddo

!----------------------------------------------------------------------
! Close the set file
!----------------------------------------------------------------------
    close(unit)
    
    return
    
  end subroutine rdqcfilenames_setfile

!######################################################################

  subroutine rdqcfilenames_inpfile

    use constants
    use channels
    use parsemod
    use iomod
    use igridglobal
    
    implicit none

    integer :: i

!----------------------------------------------------------------------
! Read to the qcfiles section
!----------------------------------------------------------------------
    rewind(iin)
5   call rdinp(iin)
    if (lend) goto 100
    if (keyword(1).ne.'$qcfiles') goto 5

!----------------------------------------------------------------------
! First pass: determine the no. files and allocate arrays
!----------------------------------------------------------------------
    ! Determine the no. files
    nfiles=0
10  call rdinp(iin)
    if (keyword(1).ne.'$end') then
       nfiles=nfiles+1
       goto 10
    endif

    ! Allocate arrays
    allocate(qcfiles(nfiles))
    qcfiles=''

!----------------------------------------------------------------------
! Second pass: read in the filenames
!----------------------------------------------------------------------
    do i=1,nfiles+1
       backspace(iin)
    enddo

    do i=1,nfiles
       call rdinp(iin)
       qcfiles(i)=keyword(1)
    enddo

    return

100 continue
    errmsg='The $qcfiles section could not be found'
    call error_control
    
  end subroutine rdqcfilenames_inpfile
    
!######################################################################
  
end program igrid
