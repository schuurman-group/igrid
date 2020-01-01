module wigner_harmonic

contains

!######################################################################

  subroutine sample_wigner_harmonic

    use constants
    use channels
    use iomod
    use sysinfo
    use igridglobal
    
    implicit none

    integer            :: j,m,n
    integer, parameter :: maxtry=75000
    real(dp)           :: q(nmodes),p(nmodes)
    real(dp)           :: fwnorm,rand
    character(len=1)   :: cr
    logical            :: ok

    cr = char(13)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(qsample(nmodes,nsample))
    qsample=0.0d0

    allocate(psample(nmodes,nsample))
    psample=0.0d0

!----------------------------------------------------------------------
! Sample the Wigner distribution within the harmonic approximation
!----------------------------------------------------------------------
    ! Loop over samples
    do j=1,nsample

       ! Output our progress
       write(*,100,advance='no') int(real(j-1)/nsample*100),cr
100    format ('Sampling progress ','[ ', i0, '% ]',a)

       ok=.false.
       do m=1,maxtry

          ! Random position and momentum vectors
          call qp_rand_harmonic(q,p)

          ! Calculate the Wigner distribution value within the
          ! harmonic approximation at the current point in phase
          ! space
          ! Note that here we are calculating the *normalised*
          ! Wigner distribution function value
          call calc_wigner_harmonic(q,p,fwnorm)
          
          ! Accept or reject the current phase space vector
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
    
  end subroutine sample_wigner_harmonic

!######################################################################

  subroutine qp_rand_harmonic(q,p)

    use constants
    use sysinfo
    use igridglobal
    
    implicit none

    integer             :: n
    real(dp)            :: q(nmodes),p(nmodes)
    real(dp)            :: rand
    real(dp)            :: b1,b2
    real(dp), parameter :: sigma=1.0d0/sqrt(2.0d0)

    ! Loop over modes
    do n=1,nmodes

       ! Position vector
       b1=-sigma
       b2=+sigma
       call random_number(rand)
       q(n)=b1+rand*(b2-b1)

       ! Momentum vector
       b1=-sigma
       b2=+sigma
       call random_number(rand)
       p(n)=b1+rand*(b2-b1)

    enddo
    
    return
    
  end subroutine qp_rand_harmonic
    
!######################################################################

  subroutine calc_wigner_harmonic(q,p,fwnorm)

    use constants
    use sysinfo
    use igridglobal
    
    implicit none

    integer  :: n
    real(dp) :: q(nmodes),p(nmodes)
    real(dp) :: fwnorm,fw1m

    !***************************************************************
    !**** Note that the normalised Wigner distribution function ****
    !**** is being calculated here. Hence no division by 1/pi^N ****
    !***************************************************************

    ! Initialisation
    fwnorm=1.0d0
    
    ! Loop over modes
    do n=1,nmodes

       ! Calculate the contribution from the current mode
       fw1m=exp(-q(n)**2)*exp(-p(n)**2)

       ! Accumulate the result
       fwnorm=fwnorm*fw1m

    enddo

    return
    
  end subroutine calc_wigner_harmonic
    
!######################################################################
    
end module wigner_harmonic
