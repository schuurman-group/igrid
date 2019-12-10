module interpolation

contains

!######################################################################

  subroutine interpolate_potentials

    use constants
    use channels
    use sysinfo
    use igridglobal
    
    implicit none

    integer :: n
    
!----------------------------------------------------------------------
! Re-allocate the grid and potential arrays
!----------------------------------------------------------------------
    call realloc_grids

!----------------------------------------------------------------------
! Modify the grids and calculate interpolated potential values
!----------------------------------------------------------------------
    ! Loop over modes
    do n=1,nmodes

       ! Modify the grid for the current mode if the no. grid points
       ! is different
       if (npnts(n).ne.ngrid(n)) call modify_grid(n)

       ! Reset npnts
       npnts(n)=ngrid(n)
       
    enddo

    return
    
  end subroutine interpolate_potentials

!######################################################################

  subroutine realloc_grids

    use constants
    use channels
    use sysinfo
    use igridglobal
    
    implicit none

    integer :: n,maxpnts1
    real(dp), allocatable :: qgrid1(:,:),pot1(:,:)

!----------------------------------------------------------------------
! Allocate the temporary swap arrays
!----------------------------------------------------------------------
    maxpnts1=maxpnts
    allocate(qgrid1(maxpnts1,nmodes))
    allocate(pot1(maxpnts1,nmodes))
    qgrid1=0.0d0
    pot1=0.0d0
    
!----------------------------------------------------------------------
! Copy the old grid and potential arrays
!----------------------------------------------------------------------
    qgrid1=qgrid
    pot1=pot

!----------------------------------------------------------------------
! Re-dimension the grid and potential arrays
!----------------------------------------------------------------------
    deallocate(qgrid)
    deallocate(pot)
    maxpnts=maxval(ngrid)
    allocate(qgrid(maxpnts,nmodes))
    allocate(pot(maxpnts,nmodes))

!----------------------------------------------------------------------
! Fill in the grid and potential arrays with their original values
! (note that not all grids may change)
!----------------------------------------------------------------------
    qgrid=0.0d0
    pot=0.0d0
    qgrid(1:maxpnts1,:)=qgrid1(1:maxpnts1,:)
    pot(1:maxpnts1,:)=pot1(1:maxpnts1,:)  
    
!----------------------------------------------------------------------
! Deallocate the temporary swap arrays
!----------------------------------------------------------------------
    deallocate(qgrid1)
    deallocate(pot1)
    
    return
    
  end subroutine realloc_grids

!######################################################################

  subroutine modify_grid(n)

    use constants
    use channels
    use sysinfo
    use splinemod
    use igridglobal
    
    implicit none

    integer, intent(in)   :: n
    integer               :: i,nnew,n1,n2
    real(dp)              :: range,q1
    real(dp), allocatable :: qgrid1(:),pot1(:),ypp(:)
    real(dp)              :: ypval,yppval
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(qgrid1(npnts(n)))
    allocate(pot1(npnts(n)))
    allocate(ypp(npnts(n)))
    
!----------------------------------------------------------------------
! Copies of the orginal grid and potential values for use in the
! interpolation procedure
!----------------------------------------------------------------------
    qgrid1(1:npnts(n))=qgrid(1:npnts(n),n)
    pot1(1:npnts(n))=pot(1:npnts(n),n)
    
!----------------------------------------------------------------------
! Modify the grid points
!----------------------------------------------------------------------
    nnew=ngrid(n)
    range=qgrid(npnts(n),n)-qgrid(1,n)
    q1=qgrid(1,n)

    do i=1,ngrid(n)
       qgrid(i,n)=q1+(i-1)*range/(nnew-1)
    enddo

!----------------------------------------------------------------------
! Compute the interpolated potential values at the new grid points
!----------------------------------------------------------------------
    n1=npnts(n)
    n2=ngrid(n)

    call spline_cubic_set(n1,qgrid1,pot1,0,0.0d0,0,0.0d0,ypp)

    do i=1,n2
       call spline_cubic_val(n1,qgrid1,pot1,ypp,qgrid(i,n),pot(i,n),&
            ypval,yppval)
    enddo
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(qgrid1)
    deallocate(pot1)
    deallocate(ypp)
    
    return
    
  end subroutine modify_grid
    
!######################################################################
  
end module interpolation
