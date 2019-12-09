module igridglobal

  use constants
  
  implicit none

  ! Number of QC output files
  integer                         :: nfiles 

  ! Filenames
  character(len=120), allocatable :: qcfiles(:)
  character(len=120)              :: setfile
  logical                         :: lsetfile

  ! Potential values
  real(dp), allocatable           :: pot(:,:)
  
  ! Normal mode coordinates
  real(dp), allocatable           :: qvec(:,:)
  
end module igridglobal
