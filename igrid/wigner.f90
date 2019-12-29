module wigner

contains

!######################################################################

  subroutine sample_wigner

    use constants
    use igridglobal

    implicit none

!----------------------------------------------------------------------
! Precalculate the 2nd derivatives of the 1D position representation
! eigenfunctions at the grid points (i.e., at the knot points in the
! cubic spline interpolation scheme)
!----------------------------------------------------------------------
    call precalc_wfpp

!----------------------------------------------------------------------
! Sample the Wigner distribution
!----------------------------------------------------------------------
    call wigner_distribution
    
    return
    
  end subroutine sample_wigner

!######################################################################

  subroutine precalc_wfpp

    use constants
    use channels
    use sysinfo
    use splinemod
    use igridglobal
    
    implicit none

    integer :: n,np,indx
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(wfpp(maxpnts,nmodes))
    wfpp=0.0d0

!----------------------------------------------------------------------
! Calculate the 2nd derivatives of the wavefunctions at the knot points
!----------------------------------------------------------------------
    ! Loop over modes
    do n=1,nmodes

       ! No. knot, i.e., grid, points for the current mode
       np=npnts(n)

       ! Eigenfunction index for the current mode
       indx=eigindx(n)

       ! Calculate the knot point second derivatives
       call spline_cubic_set(np,qgrid(1:np,n),&
            eigvec1d(1:np,indx,n),0,0.0d0,0,0.0d0,&
            wfpp(1:np,n))
       
    enddo
       
    return
    
  end subroutine precalc_wfpp
    
!######################################################################
!
!  subroutine init_quadrature
!
!    use constants
!    use channels
!    use iomod
!    use sysinfo
!    use splinemod
!    use igridglobal
!    
!    implicit none
!
!    integer               :: n,i,indx,np
!    real(dp)              :: ds,q
!    real(dp), allocatable :: ypp(:,:)
!    real(dp)              :: ypval,yppval
!    
!!----------------------------------------------------------------------
!! For now, we will hard-wire the no. quadrature points
!!----------------------------------------------------------------------
!    nquad=501
!
!!----------------------------------------------------------------------
!! Allocate arrays
!!----------------------------------------------------------------------
!    allocate(qquad(nquad,nmodes))
!    allocate(fquad(nquad,nmodes))
!    qquad=0.0d0
!    fquad=0.0d0
!
!    allocate(ypp(maxpnts,nmodes))
!    ypp=0.0d0
!    
!!----------------------------------------------------------------------
!! Get the quadrature points
!!----------------------------------------------------------------------
!    ! Loop over modes
!    do n=1,nmodes
!
!       ! Quadrature grid spacing
!       ds=(qgrid(npnts(n),n)-qgrid(1,n))/(nquad-1)
!
!       ! Loop over quadrature points
!       do i=1,nquad
!          qquad(i,n)=qgrid(1,n)+(i-1)*ds
!       enddo
!       
!    enddo
!
!!----------------------------------------------------------------------
!! Get the value of the 1D eigenfunctions at quadrature points using
!! cubic spline interpolation.
!! Note that we have to account for the sqrt(dQ) normalisation factors
!! that were folded into the eigenfunction values.
!!----------------------------------------------------------------------
!    ! Loop over modes
!    do n=1,nmodes
!
!       np=npnts(n)
!       indx=eigindx(n)
!              
!       ! Calculate the 2nd derivatives at the knot points
!       call spline_cubic_set(np,qgrid(1:np,n),&
!            eigvec1d(1:np,indx,n),0,0.0d0,0,0.0d0,&
!            ypp(1:np,n))
!
!       ! Calculate the interpolated eigenfunction values at the
!       ! quadrature points
!       do i=1,nquad
!          call spline_cubic_val(np,qgrid(1:np,n),&
!               eigvec1d(1:np,indx,n),ypp(1:np,n),&
!               qquad(i,n),fquad(i,n),ypval,yppval)
!       enddo
!       
!    enddo
!
!!----------------------------------------------------------------------
!! Deallocate arrays
!!----------------------------------------------------------------------
!    deallocate(ypp)
!    
!    return
!    
!  end subroutine init_quadrature
!
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
          !**************************************************
          ! Note that we need to first determine the maximum
          ! value of W(q,p) so that we can normalise it...
          !**************************************************
          
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
       !call calc_wigner_1mode(n,q(n),p(n),fw1mode)

       call calc_wigner_1mode(n,1.0d0,1.0d0,fw1mode)
       
       ! Accumulate the result
       fw=fw*fw1mode
       
    enddo

    return
    
  end subroutine calc_wigner

!######################################################################

  subroutine calc_wigner_1mode(n,q,p,wqp)

    use constants
    use sysinfo
    use splinemod
    use igridglobal
    
    implicit none

    integer,  intent(in)    :: n
    integer                 :: ia,ib,i,ii,ic,ip,im
    integer                 :: np,indx
    real(dp), intent(in)    :: q,p
    real(dp), intent(inout) :: wqp
    real(dp)                :: qa,qb,deltaq,ds,s
    real(dp)                :: dista,distb
    real(dp)                :: psi1,psi2,tmpp,tmppp
    
!**********************************************************************
! Compute the integral
!
! W(q,p) = 1/(2pi) int_{-infty}^{+infty} exp(ips) psi(q+s/2) psi(q-s/2) ds
!
! for the given values of q and p using the interpolated grid
! representation of the position wavefunction psi(q).
!
! For now, we will just use trapazoidal rule quadrature with a
! suitably large number of quadrature points. We can figure out a
! better quadrature scheme later if this proves to be too slow.
!
! Note that the position representation wavefunction is real,
! allowing to write W(q,p) as follows:
!
! W(q,p) = 1/pi int_{0}^{+infty} cos(ps) psi(q+s/2) psi(q-s/2) ds,
!
! which will be our working equation.
!**********************************************************************

!----------------------------------------------------------------------
! Integration bounds.
!----------------------------------------------------------------------
! We only need to consider quadrature points in
! the interval [q-deltaq,q+deltaq], where deltaq is the shortest
! distance from q to either the effective start or end of the grid.
!----------------------------------------------------------------------
    ia=qbounds(1,n)
    ib=qbounds(2,n)

    qa=qgrid(ia,n)
    qb=qgrid(ib,n)

    dista=abs(q-qa)
    distb=abs(q-qb)

    if (dista.lt.distb) then
       deltaq=dista
    else
       deltaq=distb
    endif
    
!----------------------------------------------------------------------
! Quadrature point spacing
!----------------------------------------------------------------------
    ds=deltaq/((nquad-1)/2-1)
    
!----------------------------------------------------------------------
! Quadrature points
!----------------------------------------------------------------------
    ! LHS points
    do i=1,(nquad-1)/2
       qquad(i)=q-deltaq+(i-1)*ds
    enddo

    ! Centre point
    qquad((nquad-1)/2+1)=q
    
    ! RHS points
    ii=0
    do i=(nquad-1)/2+2,nquad
       ii=ii+1
       qquad(i)=q+ii*ds
    enddo

!----------------------------------------------------------------------
! Evalute the integral W(q,p)
!----------------------------------------------------------------------
    np=npnts(n)
    indx=eigindx(n)

    ! Initialisation
    wqp=0.0d0

    ! Loop over quadrature points, start from the middle of the qquad
    ! array and working outwards in each direction simulataneously
    ic=(nquad-1)/2+1
    do i=0,(nquad-1)/2

       ! psi(q-s/2)
       im=ic-i
       call spline_cubic_val(np,qgrid(1:np,n),&
            eigvec1d(1:np,indx,n),wfpp(1:np,n),&
            qquad(im),psi1,tmpp,tmppp)
       
       ! psi(q+s/2)
       ip=ic+i
       call spline_cubic_val(np,qgrid(1:np,n),&
            eigvec1d(1:np,indx,n),wfpp(1:np,n),&
            qquad(ip),psi2,tmpp,tmppp)

       ! Current value of s
       s=2.0d0*(qquad(ip)-q)
       
       ! Contribution to W(q,p)
       wqp=wqp+cos(s*p)*psi1*psi2

    enddo
    
!----------------------------------------------------------------------
! Prefactors
!----------------------------------------------------------------------
    wqp=wqp/(2.0d0*pi)
    wqp=wqp*ds

    print*,wqp
    
    return
    
  end subroutine calc_wigner_1mode
    
!######################################################################
  
end module wigner
