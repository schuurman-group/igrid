module wigner

contains

!######################################################################

  subroutine sample_wigner

    use constants
    use igridglobal

    implicit none

!----------------------------------------------------------------------
! Get the quadrature points and weights
!----------------------------------------------------------------------
    call init_quadrature
    
    return
    
  end subroutine sample_wigner

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
