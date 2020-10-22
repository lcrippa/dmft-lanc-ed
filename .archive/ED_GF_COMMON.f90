MODULE ED_GF_COMMON
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,free_unit,reg,free_units,txtfy
  USE SF_LINALG,  only: inv,eigh,eye
  USE SF_SP_LINALG, only: sp_lanc_tridiag
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_IO                     !< this contains the routine to print GF,Sigma and G0
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_BATH_FUNCTIONS
  USE ED_SETUP
  USE ED_HAMILTONIAN
  USE ED_AUX_FUNX



  real(8),dimension(:),pointer                :: state_cvec
  real(8)                                     :: state_e
  complex(8),allocatable,dimension(:,:)       :: auxGmats,auxGreal



END MODULE ED_GF_COMMON




  ! USE SF_CONSTANTS, only:one,xi,zero,pi
  ! USE SF_TIMER  
  ! USE SF_IOTOOLS, only: str,free_unit,reg,free_units,txtfy
  ! USE SF_LINALG,  only: inv,eigh,eye
  ! USE SF_SP_LINALG, only: sp_lanc_tridiag
  ! USE ED_INPUT_VARS
  ! USE ED_VARS_GLOBAL
  ! USE ED_IO                     !< this contains the routine to print GF,Sigma and G0
  ! USE ED_EIGENSPACE
  ! USE ED_BATH
  ! USE ED_BATH_FUNCTIONS
  ! USE ED_SETUP
  ! USE ED_HAMILTONIAN
  ! USE ED_AUX_FUNX
  ! !
  ! implicit none

  ! !Lanczos shared variables
  ! !=========================================================
  ! real(8),dimension(:),pointer                :: state_cvec
  ! real(8)                                     :: state_e


  ! !AUX GF
  ! !=========================================================
  ! complex(8),allocatable,dimension(:,:)       :: auxGmats,auxGreal
