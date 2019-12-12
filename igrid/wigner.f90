module wigner

contains

!######################################################################

  subroutine sample_wigner

    use constants
    use igridglobal

    implicit none

!----------------------------------------------------------------------
! Get the grids in the momentum representation
!----------------------------------------------------------------------
    call get_momentum_grids
    
!----------------------------------------------------------------------
! Get the quadrature points and weights
!----------------------------------------------------------------------
    call init_quadrature
    
    return
    
  end subroutine sample_wigner

!######################################################################

  subroutine get_momentum_grids

    use constants
    use channels
    use iomod
    use sysinfo
    use igridglobal
    
    implicit none

    integer   :: n,k,nk,ik
    real(dp)  :: dq
    real(dp)  :: dk(nmodes)

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

  subroutine init_quadrature

    use constants
    use channels
    use iomod
    use sysinfo
    use igridglobal
    
    implicit none

    integer :: n
    
!----------------------------------------------------------------------
! Determine 'effective' position and momenta grids where the 1D
! wavefunctions have non-negligible values
!----------------------------------------------------------------------
    
    return
    
  end subroutine init_quadrature
    
!######################################################################
  
end module wigner
