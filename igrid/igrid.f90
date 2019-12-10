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
  use eigenmod
  use interpolation
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

!----------------------------------------------------------------------
! Determine the quantum chemistry calculation type.
! Note that we are assuming here that the same level of theory was
! used for all points, but that we are not checking this...
!----------------------------------------------------------------------
  call entype(qctyp,qcfiles(1))

!----------------------------------------------------------------------
! Parse the quantum chemistry output files
!----------------------------------------------------------------------
  call parse_qcfiles

!----------------------------------------------------------------------
! Parse the eigenfunction indices in the input file
!----------------------------------------------------------------------
  call parse_initwf
  
!----------------------------------------------------------------------
! If necessary, interpolate the 1D potentials
!----------------------------------------------------------------------
  if (interpolate) call interpolate_potentials
  
!----------------------------------------------------------------------
! Calculate the eigenfunctions of the 1D Hamiltonians
!----------------------------------------------------------------------
  call eigen1d
  
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

    ! Eigenfunction numbers to be read from the input file
    eiginp=.false.
    
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

       else if (keyword(i).eq.'$initwf') then
          eiginp=.true.
          ! Read to the end of the initwf section: this will
          ! be parsed later once the number of normal modes
          ! is known
30        call rdinp(iin)
          if (keyword(1).ne.'$end') goto 30
          
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

  subroutine parse_initwf

    use constants
    use channels
    use iomod
    use parsemod
    use sysinfo
    use igridglobal
    
    implicit none

    integer          :: n
    character(len=4) :: an
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(eigindx(nmodes))
    allocate(ngrid(nmodes))
    
!----------------------------------------------------------------------
! Default values: 1 <-> ground states of the 1D Hamiltonians
!----------------------------------------------------------------------
    eigindx=1

!----------------------------------------------------------------------
! If the initwf section was present in the input file, then parse it
! here
!----------------------------------------------------------------------
    if (eiginp) then
       
       rewind(iin)

       ! Read to the initwf section
5      call rdinp(iin)
       if (keyword(1).ne.'$initwf') goto 5

       ! Read the user-specified 1-mode eigenstate indices
10     call rdinp(iin)
       if (keyword(1).ne.'$end') then
          read(keyword(1),*) n
          read(keyword(2),*) ngrid(n)
          read(keyword(3),*) eigindx(n)
          goto 10
       endif
       
    endif

!----------------------------------------------------------------------
! Check that the numbers of grid points are all odd
!----------------------------------------------------------------------
    do n=1,nmodes
       if (mod(ngrid(n),2).eq.0) then
          write(an,'(i4)') n
          errmsg='Error: even number of grid points for mode '&
               //trim(adjustl(an))
          call error_control
       endif
    enddo

!----------------------------------------------------------------------
! Set the interpolate flag to true if the grids do not correspond to
! those read from the QC output
!----------------------------------------------------------------------
    interpolate=.false.
    do n=1,nmodes
       if (ngrid(n).ne.npnts(n)) then
          interpolate=.true.
          exit
       endif
    enddo
    
    return
    
  end subroutine parse_initwf
  
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

  subroutine parse_qcfiles

    use constants
    use channels
    use iomod
    use sysinfo
    use ioqc
    use igridglobal
    
    implicit none

    integer          :: n
    character(len=3) :: an
    
!----------------------------------------------------------------------
! First pass: determine the number of points along a normal
! mode (npnts) and allocate arrays
!----------------------------------------------------------------------
    ! Get the numbers of grid points for each mode
    call get_npnts

    ! Maximum number of grid points
    maxpnts=maxval(npnts)

    ! Allocate arrays
    allocate(pot(maxpnts,nmodes))
    allocate(qgrid(maxpnts,nmodes))
        
!----------------------------------------------------------------------
! Check to make sure that all grids have an odd number of points
!----------------------------------------------------------------------
    do n=1,nmodes
       if (mod(npnts(n),2).eq.0) then
          write(an,'(i3)') n
          errmsg='Error: even number of grid points for mode '&
               //trim(adjustl(an))
          call error_control
       endif
    enddo
    
!----------------------------------------------------------------------
! Second pass: fill in the potential and normal mode displacement
! arrays
!----------------------------------------------------------------------
    call get_grids
    
    return
    
  end subroutine parse_qcfiles

!######################################################################

  subroutine get_npnts

    use constants
    use channels
    use iomod
    use sysinfo
    use ioqc
    use igridglobal
    
    implicit none

    integer             :: i,n,ndisp,indx
    real(dp)            :: q(nmodes)
    real(dp)            :: x(ncoo)
    real(dp), parameter :: thrsh=1e-5_dp
    
!----------------------------------------------------------------------
! Allocate the npnts array
!----------------------------------------------------------------------
    allocate(npnts(nmodes))
    npnts=0
    
!----------------------------------------------------------------------
! Determine the number of points for each normal mode
!----------------------------------------------------------------------
    ! Loop over files
    do i=1,nfiles

       ! Read the Cartesian coordinates (in Bohr)
       call getxcoo(qcfiles(i),x)

       ! Calculate the normal mode coordinates
       q=matmul(coonm,(x-xcoo0)/ang2bohr)

       ! Determine the index of the displaced mode
       ndisp=0
       indx=0
       do n=1,nmodes
          if (abs(q(n)).gt.thrsh) then
             ndisp=ndisp+1
             indx=n
          endif
       enddo
       
       ! Exit if more than a single mode has been displaced
       if (ndisp.gt.1) then
          errmsg='More than 1 normal mode is displaced in ' &
               //trim(qcfiles(i))
          call error_control
       endif

       ! Update npnts for the displaced mode
       npnts(indx)=npnts(indx)+1
       
    enddo

    ! Include the Q0 reference point
    npnts=npnts+1
    
!----------------------------------------------------------------------
! Ouput some information to the log file
!----------------------------------------------------------------------
    write(ilog,'(/,20a)') ('-',i=1,20)
    write(ilog,'(a)') ' Mode | Grid points'
    write(ilog,'(20a)') ('-',i=1,20)

    do n=1,nmodes
       write(ilog,'(x,i3,2x,a,x,i4)') n,'|',npnts(n)
    enddo
    
    write(ilog,'(20a)') ('-',i=1,20)
    
    return
    
  end subroutine get_npnts

!######################################################################

  subroutine get_grids

    use constants
    use channels
    use iomod
    use sysinfo
    use ioqc
    use utils
    use igridglobal
    
    implicit none

    integer             :: i,n,indx
    integer             :: icount(nmodes)
    integer             :: gindx(maxpnts)
    real(dp), parameter :: thrsh=1e-5_dp 
    real(dp)            :: q(nmodes)
    real(dp)            :: x(ncoo)
    real(dp)            :: v(1),v0 
    real(dp)            :: swap(maxpnts)
    
!----------------------------------------------------------------------
! Read the normal mode coordinate and potential values on the 1D grids
!----------------------------------------------------------------------
    ! Initialise the mode counter
    icount=0
    
    ! Loop over files
    do i=1,nfiles

       ! Read the Cartesian coordinates (in Bohr)
       call getxcoo(qcfiles(i),x)

       ! Calculate the normal mode coordinates
       q=matmul(coonm,(x-xcoo0)/ang2bohr)

       ! Read the potential value
       call geten(v,1,qcfiles(i))
       
       ! Determine the index of the displaced mode
       indx=0
       do n=1,nmodes
          if (abs(q(n)).gt.thrsh) indx=n
       enddo

       ! Normal mode and potential value
       if (indx.eq.0) then
          ! Reference point (Q0): contributes to all grids
          do n=1,nmodes
             icount(n)=icount(n)+1
             qgrid(icount(n),n)=0.0d0
             pot(icount(n),n)=v(1)
          enddo
          v0=v(1)
       else
          ! Single displaced mode
          icount(indx)=icount(indx)+1
          qgrid(icount(indx),indx)=q(indx)
          pot(icount(indx),indx)=v(1)
       endif

    enddo

!----------------------------------------------------------------------
! Sort the grids
!----------------------------------------------------------------------
    ! Loop over modes
    do n=1,nmodes

       ! Sort the grid points in order of increasing coordinate
       ! value
       call dsortindxa1('A',npnts(n),qgrid(1:npnts(n),n),gindx)
       
       ! Sort the grid points
       do i=1,npnts(n)
          swap(i)=qgrid(gindx(i),n)
       enddo
       qgrid(1:npnts(n),n)=swap(1:npnts(n))

       ! Sort the potential values
       do i=1,npnts(n)
          swap(i)=pot(gindx(i),n)
       enddo
       pot(1:npnts(n),n)=swap(1:npnts(n))
       
    enddo
    
!----------------------------------------------------------------------
! Shift the potentials by V0
!----------------------------------------------------------------------
    do n=1,nmodes
       pot(1:npnts(n),n)=pot(1:npnts(n),n)-v0
    enddo
    
    return
    
  end subroutine get_grids
    
!######################################################################
  
end program igrid
