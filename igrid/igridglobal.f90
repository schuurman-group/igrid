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
  integer, allocatable            :: npnts(:),ngrid(:)
  
  ! Potential values
  real(dp), allocatable           :: pot(:,:)
  
  ! Normal mode coordinates
  real(dp), allocatable           :: qgrid(:,:)
  real(dp), allocatable           :: dq(:)
  logical                         :: interpolate

  ! Momentum grids
  real(dp), allocatable           :: pgrid(:,:)
  real(dp), allocatable           :: dk(:)
  
  ! Eigenpairs of the 1-mode Hamiltonians
  integer, allocatable            :: eigindx(:)
  real(dp), allocatable           :: eigvec1d(:,:,:)
  real(dp), allocatable           :: eigval1d(:,:)
  real(dp), allocatable           :: qexp(:)
  real(dp), allocatable           :: qvar(:)
  real(dp), allocatable           :: pexp(:)
  real(dp), allocatable           :: pvar(:)
  complex(dp), allocatable        :: peigvec1d(:,:,:)
  logical                         :: eiginp

  ! Wigner distribution sampling
  integer                         :: nsample
  integer, allocatable            :: qbounds(:,:)
  integer, allocatable            :: pbounds(:,:)
  integer, parameter              :: nquad=501 ! (This has to be odd)
  real(dp), dimension(nquad)      :: qquad
  real(dp), allocatable           :: wfpp(:,:)
  real(dp), allocatable           :: maxfw1m(:)
  real(dp)                        :: maxfw

  ! Sampled positions and momenta
  real(dp), allocatable           :: qsample(:,:)
  real(dp), allocatable           :: psample(:,:)
  
end module igridglobal
