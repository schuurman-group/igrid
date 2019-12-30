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
! Get the effective grids for use in the evaluation of the 1D Wigner
! distribution functions
!----------------------------------------------------------------------
  call effective_grids
    
!----------------------------------------------------------------------
! Determine the maximum of the Wigner distribution
!----------------------------------------------------------------------
  call get_maxfw
    
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
  
  subroutine get_maxfw

    use constants
    use channels
    use iomod
    use sysinfo
    use igridglobal
    
    implicit none

    integer           :: iq,ip,n,unit
    real(dp)          :: int,fw
    character(len=60) :: acheck
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(maxfw1m(nmodes))
    maxfw1m=0.0d0
    
!----------------------------------------------------------------------
! Compute the 1D Wigner distributions at the grid points
!----------------------------------------------------------------------
    do n=1,nmodes

       int=0

       do iq=qbounds(1,n),qbounds(2,n)
          do ip=pbounds(1,n),pbounds(2,n)
          
             call calc_wigner_1mode(n,qgrid(iq,n),pgrid(ip,n),fw)

             if (abs(fw).gt.abs(maxfw1m(n))) maxfw1m(n)=fw
             
             int=int+fw*dq(n)*dk(n)
          
          enddo
       enddo

    enddo

!----------------------------------------------------------------------
! Set the maximum value of W(q,p).
! Here we are:
! (i)  making use of the fact that we are dealing with a
!      single Hartree product wavefunction, and;
! (ii) assuming that the maximum values of the 1D |W(q_n,p_n)|'s will
!      correspond to positive values of W(q_n,p_n).
!----------------------------------------------------------------------
    maxfw=1.0d0
    do n=1,nmodes
       maxfw=maxfw*maxfw1m(n)
    enddo

    return
    
  end subroutine get_maxfw

!######################################################################
  
  subroutine effective_grids

    use constants
    use channels
    use iomod
    use sysinfo
    use igridglobal
    
    implicit none

    integer             :: n,i,i1,np
    real(dp)            :: pop
    real(dp), parameter :: thrsh=1e-8_dp

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(qbounds(2,nmodes))
    allocate(pbounds(2,nmodes))
    qbounds=0
    pbounds=0

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    do n=1,nmodes
       np=npnts(n)
       pbounds(1,n)=1
       pbounds(2,n)=np
       qbounds(1,n)=1
       qbounds(2,n)=np
    enddo
    
!----------------------------------------------------------------------
! Determine the bounds of the effective grids that will be used in
! the position and momentum sampling. Here, we will discard the grid
! points at the edges for which the grid populations are vanishingly
! small
!----------------------------------------------------------------------
    ! Loop over modes
    do n=1,nmodes

       np=npnts(n)
       
       ! Position grid: LHS grid points
       do i=1,(np-1)/2
          pop=abs(eigvec1d(i,eigindx(n),n))
          if (pop.ge.thrsh) then
             qbounds(1,n)=i
             exit
          endif
       enddo

       ! Position grid: RHS grid points
       do i1=1,(np-1)/2
          i=np-i1+1
          pop=abs(eigvec1d(i,eigindx(n),n))
          if (pop.ge.thrsh) then
             qbounds(2,n)=i
             exit
          endif
       enddo

       ! Momentum grid: LHS grid points
       do i=1,(np-1)/2
          pop=abs(peigvec1d(i,eigindx(n),n))
          if (pop.ge.thrsh) then
             pbounds(1,n)=i
             exit
          endif
       enddo

       ! Momentum grid: RHS grid points
       do i1=1,(np-1)/2
          i=np-i1+1
          pop=abs(peigvec1d(i,eigindx(n),n))
          if (pop.ge.thrsh) then
             pbounds(2,n)=i
             exit
          endif
       enddo

    enddo
    
    return
    
  end subroutine effective_grids
  
!######################################################################

  subroutine wigner_distribution

    use constants
    use channels
    use iomod
    use sysinfo
    use igridglobal
    
    implicit none

    integer            :: j,m,n
    integer, parameter :: maxtry=50000
    real(dp)           :: q(nmodes),p(nmodes)
    real(dp)           :: fw,fwnorm,rand
    logical            :: ok

    character*1 cr
    cr = char(13)
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(qsample(nmodes,nsample))
    qsample=0.0d0

    allocate(psample(nmodes,nsample))
    psample=0.0d0
    
!----------------------------------------------------------------------
! Sample the Wigner distribution
!----------------------------------------------------------------------
    ! Loop over samples
    do j=1,nsample

       ! Output our progress
       write(*,100,advance='no') int(real(j)/nsample*100),cr
100    format ('Sampling progress ','[ ', I0, '% ]',A)
       
       ok=.false.
       do m=1,maxtry

          ! Random position and momentum vectors
          call qp_rand(q,p)

          ! Calculate the Wigner distribution value at the
          ! current point in phase space
          call calc_wigner(q,p,fw)
          
          ! Accept or reject the current phase space vector
          fwnorm=fw/maxfw
          call random_number(rand)
          if (rand.lt.fwnorm) then
             ok=.true.
             exit
          endif
          
       enddo

       ! Exit here if the sampling was not successful
       if (.not.ok) then
          errmsg='Unsuccessful sampling in subroutine &
               wigner_distribution. Quitting.'
          call error_control
       endif

       ! Save the sampled position and momentum vector
       qsample(:,j)=q
       psample(:,j)=p
       
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
    real(dp) :: b1,b2
    
    do n=1,nmodes

       ! Position
       b1=qexp(n)-sqrt(2.0d0)*qvar(n)
       b2=qexp(n)+sqrt(2.0d0)*qvar(n)
       call random_number(rand)
       q(n)=b1+rand*(b2-b1)

       ! Momentum
       b1=pexp(n)-sqrt(2.0d0)*pvar(n)
       b2=pexp(n)+sqrt(2.0d0)*pvar(n)
       call random_number(rand)
       p(n)=b1+rand*(b2-b1)

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
! W(q,p) = 1/(pi) int_{-infty}^{+infty} exp(2ips) psi(q+s) psi(q-s) ds
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
! W(q,p) = 2/pi int_{0}^{+infty} cos(2ps) psi(q+s) psi(q-s) ds,
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
    ds=deltaq/((nquad-1)/2)
    
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

       ! psi(q-s)
       im=ic-i
       call spline_cubic_val(np,qgrid(1:np,n),&
            eigvec1d(1:np,indx,n),wfpp(1:np,n),&
            qquad(im),psi1,tmpp,tmppp)
       
       ! psi(q+s)
       ip=ic+i
       call spline_cubic_val(np,qgrid(1:np,n),&
            eigvec1d(1:np,indx,n),wfpp(1:np,n),&
            qquad(ip),psi2,tmpp,tmppp)

       ! Current value of s
       s=qquad(ip)-q
       
       ! Contribution to W(q,p)
       wqp=wqp+cos(2.0d0*s*p)*psi1*psi2

    enddo
    
!----------------------------------------------------------------------
! Prefactors
!----------------------------------------------------------------------
    ! 1/pi prefactor in W(q,p)
    wqp=wqp/pi

    ! Multiplication by the quadrature point spacing
    wqp=wqp*ds*2.0d0

    return
    
  end subroutine calc_wigner_1mode
    
!######################################################################
  
end module wigner
