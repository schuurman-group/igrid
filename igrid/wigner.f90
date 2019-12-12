module wigner

contains

!######################################################################

  subroutine sample_wigner

    use constants
    use igridglobal

    implicit none

!----------------------------------------------------------------------
! Get the quadrature points
!----------------------------------------------------------------------
    call init_quadrature

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

    integer  :: n
    real(dp) :: q(nmodes),p(nmodes)
    real(dp) :: dq,dp
    logical  :: ok

    ! Loop over modes
    do n=1,nmodes
       
       
       
    enddo
    
    return
    
  end subroutine qp_rand
  
!######################################################################
  
end module wigner
