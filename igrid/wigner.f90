module wigner

contains

!######################################################################

  subroutine sample_wigner

    use constants
    use igridglobal

    implicit none

!!----------------------------------------------------------------------
!! Get the quadrature points
!! *** We cannot do this yet: the quadrature points will vary with the
!! sampled positions ***
!!----------------------------------------------------------------------
!    call init_quadrature

!----------------------------------------------------------------------
! Sample the Wigner distribution
!----------------------------------------------------------------------
    call wigner_distribution
    
    return
    
  end subroutine sample_wigner

!######################################################################

  subroutine init_quadrature

    use constants
    use channels
    use iomod
    use sysinfo
    use splinemod
    use igridglobal
    
    implicit none

    integer               :: n,i,indx,np
    real(dp)              :: ds,q,dq
    real(dp), allocatable :: ypp(:,:)
    real(dp)              :: ypval,yppval
    
!----------------------------------------------------------------------
! For now, we will hard-wire the no. quadrature points
!----------------------------------------------------------------------
    nquad=501

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(qquad(nquad,nmodes))
    allocate(fquad(nquad,nmodes))
    qquad=0.0d0
    fquad=0.0d0

    allocate(ypp(maxpnts,nmodes))
    ypp=0.0d0
    
!----------------------------------------------------------------------
! Get the quadrature points
!----------------------------------------------------------------------
    ! Loop over modes
    do n=1,nmodes

       ! Quadrature grid spacing
       ds=(qgrid(npnts(n),n)-qgrid(1,n))/(nquad-1)

       ! Loop over quadrature points
       do i=1,nquad
          qquad(i,n)=qgrid(1,n)+(i-1)*ds
       enddo
       
    enddo

!----------------------------------------------------------------------
! Get the value of the 1D eigenfunctions at quadrature points using
! cubic spline interpolation.
! Note that we have to account for the sqrt(dQ) normalisation factors
! that were folded into the eigenfunction values.
!----------------------------------------------------------------------
    ! Loop over modes
    do n=1,nmodes

       np=npnts(n)
       indx=eigindx(n)
       dq=qgrid(2,n)-qgrid(1,n)
       
       ! Calculate the 2nd derivatives at the knot points
       call spline_cubic_set(np,qgrid(1:np,n),&
            eigvec1d(1:np,indx,n)/sqrt(dq),0,0.0d0,0,0.0d0,ypp(1:np,n))

       ! Calculate the interpolated eigenfunction values at the
       ! quadrature points
       do i=1,nquad
          call spline_cubic_val(np,qgrid(1:np,n),&
               eigvec1d(1:np,indx,n)/sqrt(dq),ypp(1:np,n),&
               qquad(i,n),fquad(i,n),ypval,yppval)
       enddo
       
    enddo

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ypp)
    
    return
    
  end subroutine init_quadrature

!######################################################################

  subroutine wigner_distribution

    use constants
    use channels
    use iomod
    use sysinfo
    use igridglobal
    
    implicit none

    integer            :: j,m,n
    integer, parameter :: maxtry=100
    real(dp)           :: q(nmodes),p(nmodes)
    real(dp)           :: fw
    logical            :: ok
    
!----------------------------------------------------------------------
! Sample the Wigner distribution
!----------------------------------------------------------------------
    ! Loop over samples
    do j=1,nsample

       ok=.false.
       do m=1,maxtry

          ! Random position and momentum vectors
          call qp_rand(q,p)

          ! Calculate the Wigner distribution value at the
          ! current point in phase space
          call calc_wigner(q,p,fw)
          
          ! Accept or reject the current phase space vector
          
       enddo
          
    enddo

    return
    
  end subroutine wigner_distribution

!######################################################################

  subroutine qp_rand(q,p)

    use constants
    use sysinfo
    use igridglobal
    
    implicit none

    integer  :: n,i1,i2
    real(dp) :: q(nmodes),p(nmodes)
    real(dp) :: rand
    logical  :: ok

    ! Loop over modes
    do n=1,nmodes

       ! Position
       i1=qbounds(1,n)
       i2=qbounds(2,n)
       call random_number(rand)
       q(n)=qgrid(i1,n)+rand*(qgrid(i2,n)-qgrid(i1,n))

       ! Momentum
       i1=pbounds(1,n)
       i2=pbounds(2,n)
       call random_number(rand)
       p(n)=pgrid(i1,n)+rand*(pgrid(i2,n)-pgrid(i1,n))
       
    enddo
    
    return
    
  end subroutine qp_rand

!######################################################################

  subroutine calc_wigner(q,p,fw)

    use constants
    use sysinfo
    use igridglobal
    
    implicit none

    integer  :: n
    real(dp) :: q(nmodes),p(nmodes)
    real(dp) :: fw,fw1mode

    ! Initialisation
    fw=1.0d0
    
    ! Loop over modes
    do n=1,nmodes

       ! Calculate the contribution from the current mode
       call calc_wigner_1mode(n,q(n),p(n),fw1mode)
       fw=fw*fw1mode
       
    enddo

    ! Normalisation factor
    ! ??? We need to figure out what this is ???
    
    return
    
  end subroutine calc_wigner

!######################################################################

  subroutine calc_wigner_1mode(n,q,p,fw1mode)

    use constants
    use sysinfo
    use igridglobal
    
    implicit none

    integer,  intent(in)    :: n
    real(dp), intent(in)    :: q,p
    real(dp), intent(inout) :: fw1mode
    real(dp)                :: dq
    real(dp)                :: qwf(maxpnts)
    
!----------------------------------------------------------------------
! 'Un-normalise' the qgrid array to get the position representaion
! wavefunction values at the grid points
!----------------------------------------------------------------------
    dq=qgrid(2,n)-qgrid(1,n)
    qwf=eigvec1d(:,eigindx(n),n)/sqrt(dq)
    
    ! Do the same for the peigvec1d array to get the values of the
    ! momentum representation wavefunction at the grid points
    
!    print*,"FINISH WRITING CALC_WIGNER_1MODE!"
!    stop
    
    return
    
  end subroutine calc_wigner_1mode
    
!######################################################################
  
end module wigner
