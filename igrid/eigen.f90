module eigenmod

contains
  
!######################################################################

  subroutine eigen1d

    use constants
    use channels
    use sysinfo
    use igridglobal
    
    implicit none

    integer :: n
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(eigval1d(maxpnts,nmodes))
    allocate(eigvec1d(maxpnts,maxpnts,nmodes))
    eigval1d=0.0d0
    eigvec1d=0.0d0

!----------------------------------------------------------------------
! Calculate the eigenpairs of the one-mode Hamiltonians using the
!----------------------------------------------------------------------
    do n=1,nmodes
       call fgh1d(n)
    enddo

!----------------------------------------------------------------------
! Ouput the eigenvalues of the eigenstates of interest
!----------------------------------------------------------------------
    call wreigenvalues

!----------------------------------------------------------------------
! Ouput the eigenvectors for inspection
!----------------------------------------------------------------------
    call wrwfs
    
    return

  end subroutine eigen1d

!######################################################################

  subroutine fgh1d(m)

    use constants
    use channels
    use sysinfo
    use igridglobal
    
    implicit none

    integer, intent(in) :: m
    integer             :: dim
    real(dp)            :: hmat(maxpnts,maxpnts)

!----------------------------------------------------------------------
! Number of grid points for the current mode
!----------------------------------------------------------------------
    dim=npnts(m)
    
!----------------------------------------------------------------------
! Compute the FGH Hamiltonian matrix
!----------------------------------------------------------------------
    call calc_hamiltonian(m,hmat(1:dim,1:dim),dim)

!----------------------------------------------------------------------
! Diagonalise the FGH Hamiltonian matrix
!----------------------------------------------------------------------
    call diag_hamiltonian(m,hmat(1:dim,1:dim),dim)
    
    return
    
  end subroutine fgh1d

!######################################################################

  subroutine calc_hamiltonian(m,hmat,dim)

    use constants
    use channels
    use sysinfo
    use igridglobal
    
    implicit none

    integer, intent(in)     :: m,dim
    integer                 :: i,j,l,n
    real(dp), intent(inout) :: hmat(dim,dim)
    real(dp)                :: Tl((dim-1)/2)
    real(dp)                :: dQ,omega,ftmp
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    hmat=0.0d0

!----------------------------------------------------------------------
! Potential contribution
!----------------------------------------------------------------------
    do i=1,dim
       hmat(i,i)=pot(i,m)
    enddo

!----------------------------------------------------------------------
! Kinetic energy operator contribution
!----------------------------------------------------------------------
    n=(dim-1)/2

    ! Coordinate grid spacing
    dQ=qgrid(2,m)-qgrid(1,m)
    
    ! Mode frequency
    omega=freq(m)/eh2ev

    ! Precalculation of the T_l terms
    do l=1,n
       Tl(l)=2.0d0*omega*(pi*l/(dim*dQ))**2
    enddo

    ! Kinetic energy contribution
    do i=1,dim
       do j=i,dim
          ftmp=0.0d0
          do l=1,n
             ftmp=ftmp+cos(l*2.0d0*pi*(i-j)/dim)*Tl(l)
          enddo
          hmat(i,j)=hmat(i,j)+ftmp*2.0d0/dim
          hmat(j,i)=hmat(i,j)
       enddo
    enddo
    
    return

  end subroutine calc_hamiltonian

!######################################################################

  subroutine diag_hamiltonian(m,hmat,dim)

    use constants
    use channels
    use iomod
    use sysinfo
    use igridglobal
    
    implicit none

    integer, intent(in)     :: m,dim
    integer                 :: workdim,error,i
    real(dp), intent(inout) :: hmat(dim,dim)
    real(dp)                :: work(3*dim)

!----------------------------------------------------------------------
! Diagonalise the Hamiltonian matrix
!----------------------------------------------------------------------
    eigvec1d(1:dim,1:dim,m)=hmat
    
    call dsyev('V','U',dim,eigvec1d(1:dim,1:dim,m),dim,&
         eigval1d(1:dim,m),work,3*dim,error)

    if (error.ne.0) then
       errmsg='Error in the diagonalisation of the Hamiltonian matrix'
       call error_control
    endif

    return
    
  end subroutine diag_hamiltonian

!######################################################################

  subroutine wreigenvalues

    use constants
    use channels
    use iomod
    use sysinfo
    use igridglobal
    
    implicit none

    integer  :: n,k,indx
    real(dp) :: eharm,e

!----------------------------------------------------------------------
! Table header
!----------------------------------------------------------------------
    write(ilog,'(/,51a)') ('-',k=1,51)
    write(ilog,'(a)') ' Mode | Eigenstate | Eigenvalue | &
         Harmonic Approx.'
    write(ilog,'(51a)') ('-',k=1,51)

!----------------------------------------------------------------------
! Eigenvalues
!----------------------------------------------------------------------
    ! Loop over modes
    do n=1,nmodes

       ! Index of the eigenstate for the current mode
       indx=eigindx(n)

       ! Harmonic approximation value
       eharm=(indx-1+0.5d0)*freq(n)
       
       ! Actual value
       e=eigval1d(indx,n)*eh2ev
       
       write(ilog,'(x,i3,2x,a,x,i2,9x,a,x,F6.4,x,a,2x,a,x,F6.4,x,a)') &
            n,'|',indx,'|',e,'eV','|',eharm,'eV'
            
    enddo

    write(ilog,'(51a)') ('-',k=1,51)
    
    return
    
  end subroutine wreigenvalues

!######################################################################

  subroutine wrwfs

    use constants
    use channels
    use iomod
    use sysinfo
    use igridglobal
    
    implicit none

    integer           :: i,n,unit
    character(len=60) :: filename
    logical           :: found
    
!----------------------------------------------------------------------
! Create or clean up the eigen directory
!----------------------------------------------------------------------
    ! Works with ifort
    !inquire(directory='eigen',exist=found)

    ! Works with gfortran
    inquire(file='eigen/.',exist=found)

    if (found) then
       call system('rm -rf eigen/*')
    else
       call system('mkdir eigen')
    endif

!----------------------------------------------------------------------
! Write the eigenfunction (WFs at the grid points) to file
!----------------------------------------------------------------------
    call freeunit(unit)
    
    ! Loop over modes
    do n=1,nmodes

       ! Open the wf file
       write(filename,'(a,i0,a,i0,a)') 'eigen/q',n,'_eig',eigindx(n),'.dat'
       open(unit,file=filename,form='formatted',status='unknown')

       ! Write the wf file
       do i=1,npnts(n)
          write(unit,*) qgrid(i,n),eigvec1d(i,eigindx(n),n)
       enddo
       
       ! Close the wf file
       close(unit)
       
    enddo
    
    return
    
  end subroutine wrwfs
  
!######################################################################
  
end module eigenmod
