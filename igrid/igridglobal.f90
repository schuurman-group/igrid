module igridglobal

  use constants
  
  implicit none

  ! Number of QC output files
  integer                         :: nfiles 

  ! Filenames
  character(len=120), allocatable :: qcfiles(:)
  character(len=120)              :: setfile
  logical                         :: lsetfile

  ! QC type
  integer                         :: qctyp
  
  ! Dimensions
  integer                         :: maxpnts
  integer, allocatable            :: npnts(:)
  
  ! Potential values
  real(dp), allocatable           :: pot(:,:)
  
  ! Normal mode coordinates
  real(dp), allocatable           :: qgrid(:,:)
  real(dp), allocatable           :: dq(:)
  
end module igridglobal
