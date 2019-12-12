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
! Get the grids in the momentum representation
!----------------------------------------------------------------------
    call get_momentum_grids

!----------------------------------------------------------------------
! Compute the 1D eigenstates in the momentum representation
!----------------------------------------------------------------------
    call eigen_momrep
    
!----------------------------------------------------------------------
! Ouput the eigenvalues of the eigenstates of interest
!----------------------------------------------------------------------
    call wreigenvalues

!----------------------------------------------------------------------
! Ouput the eigenvectors for inspection
!----------------------------------------------------------------------
    call wrwfs

!----------------------------------------------------------------------
! Ouput the first and last grid populations
!----------------------------------------------------------------------
    call wrgridpop
    
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

    subroutine get_momentum_grids

    use constants
    use channels
    use iomod
    use sysinfo
    use igridglobal
    
    implicit none

    integer  :: n,k,nk,ik
    real(dp) :: dq
    real(dp) :: dk(nmodes)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(pgrid(maxpnts,nmodes))
    pgrid=0.0d0
    
!----------------------------------------------------------------------
! Loop over modes and construct the momentum representation grids for
! each
!----------------------------------------------------------------------
    ! Loop over modes
    do n=1,nmodes

       ! Position grid spacing
       dq=qgrid(2,n)-qgrid(1,n)

       ! Momentum grid spacing
       dk(n)=2.0d0*pi/(npnts(n)*dq)

       ! Momentum grids
       nk=(npnts(n)-1)/2
       ik=0
       do k=-nk,nk
          ik=ik+1
          pgrid(ik,n)=dk(n)*k
       enddo
       
    enddo
    
    return
    
  end subroutine get_momentum_grids

!######################################################################

  subroutine eigen_momrep

    use constants
    use channels
    use iomod
    use sysinfo
    use igridglobal

    implicit none

    integer              :: k,l,n,i,nk,il,np
    real(dp)             :: dk,dq,qovrlp,povrlp
    complex, allocatable :: F(:,:)
    
!----------------------------------------------------------------------
! Allocate and initialise arrays
!----------------------------------------------------------------------
    allocate(peigvec1d(maxpnts,maxpnts,nmodes))
    allocate(peigval1d(maxpnts,nmodes))
    peigvec1d=czero
    peigval1d=czero
    
!----------------------------------------------------------------------
! Discrete Fourier transform of the position representation 1D
! eigenfunctions to get the momentum representation eigenfunctions.
! This, of course, would be faster to compute using an FFT, but our
! grids are small enough to not care.
!----------------------------------------------------------------------
    ! Loop over modes
    do n=1,nmodes

       ! No. grid points for the current mode
       np=npnts(n)

       ! Momentum grid bounds
       nk=(np-1)/2

       ! Position grid spacing
       dq=qgrid(2,n)-qgrid(1,n)
       
       ! Compute the transformation matrix F
       allocate(F(np,np))
       F=czero
       do i=1,np
          il=0
          do l=-nk,nk
             il=il+1
             F(il,i)=exp(-ci*2.0d0*pi*i*l/(np-1))
          enddo
       enddo

       ! Transform the position representation eigenvectors to the
       ! momentum representation
       do k=1,np
          peigvec1d(:,k,n)=matmul(F,eigvec1d(:,k,n))
       enddo

       ! Include the dk factor in the momentum representation
       ! eigenstates
       do k=1,np
          peigvec1d(:,k,n)=peigvec1d(:,k,n)&
               /sqrt(dot_product(peigvec1d(:,k,n),peigvec1d(:,k,n)))
       enddo
       
       ! Deallocate the transformation matrix F
       deallocate(F)
       
    enddo

    return
    
  end subroutine eigen_momrep
  
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

  subroutine wrgridpop

    use constants
    use channels
    use sysinfo
    use igridglobal
    
    implicit none

    integer :: k,n

!----------------------------------------------------------------------
! Write the populations of the first and last grid points to the log
! file
!----------------------------------------------------------------------
    write(ilog,'(/,34a)') ('-',k=1,34)
    write(ilog,'(8x,a)') 'Grid Populations'
    write(ilog,'(34a)') ('-',k=1,34)
    write(ilog,'(a)') ' Mode |    First    |    Last'
    write(ilog,'(34a)') ('-',k=1,34)

    do n=1,nmodes
       write(ilog,'(x,i3,2x,a,x,ES11.4,x,a,x,ES11.4)') &
            n,'|',eigvec1d(1,eigindx(n),n),'|',&
            eigvec1d(npnts(n),eigindx(n),n)
    enddo
    
    write(ilog,'(34a)') ('-',k=1,34)
    
    return
    
  end subroutine wrgridpop
    
!######################################################################
  
end module eigenmod
