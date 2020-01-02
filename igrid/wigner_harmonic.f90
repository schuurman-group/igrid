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

    integer             :: j,k,m,n
    integer, parameter  :: maxtry=1000
    real(dp)            :: q(nmodes),p(nmodes)
    real(dp)            :: fwscaled,rand,b1,b2
    real(dp), parameter :: sigma=1.0d0/sqrt(2.0d0)
    character(len=1)    :: cr
    logical             :: ok

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

       ! Loop over modes
       do n=1,nmodes

          ok=.false.
          do k=1,maxtry

             ! Generate a random position and momentum for the
             ! current mode
             b1=-3.0d0*sigma
             b2=+3.0d0*sigma
             call random_number(rand)
             q(n)=b1+rand*(b2-b1)
             call random_number(rand)
             p(n)=b1+rand*(b2-b1)

             ! Compute the value of the (scaled) 1-mode Wigner
             ! distribution function
             fwscaled=exp(-q(n)**2)*exp(-p(n)**2)

             ! Accept or reject the current 1-mode phase space point
             call random_number(rand)
             if (rand.lt.fwscaled) then
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
          
       enddo

       ! Save the sampled position and momentum vector
       qsample(:,j)=q
       psample(:,j)=p

    enddo

    return
    
  end subroutine sample_wigner_harmonic

!######################################################################
    
end module wigner_harmonic
