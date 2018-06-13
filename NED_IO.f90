MODULE ED_IO
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE SF_LINALG
  USE SF_ARRAYS, only: linspace,arange
  USE SF_IOTOOLS, only: str,reg,free_unit,splot,sread
  implicit none
  private

  !Retrieve self-energy through routines:
  interface ed_get_Smats
     module procedure ed_get_sigma_matsubara_1
     module procedure ed_get_sigma_matsubara_2
     module procedure ed_get_sigma_matsubara_lattice_1
     module procedure ed_get_sigma_matsubara_lattice_2
     !
     module procedure ed_get_all_sigma_matsubara_1
     module procedure ed_get_all_sigma_matsubara_2
     module procedure ed_get_all_sigma_matsubara_lattice_1
     module procedure ed_get_all_sigma_matsubara_lattice_2
  end interface ed_get_Smats

  interface ed_get_Sreal
     module procedure ed_get_sigma_real_1
     module procedure ed_get_sigma_real_2
     module procedure ed_get_sigma_real_lattice_1
     module procedure ed_get_sigma_real_lattice_2
     !
     module procedure ed_get_all_sigma_real_1
     module procedure ed_get_all_sigma_real_2
     module procedure ed_get_all_sigma_real_lattice_1
     module procedure ed_get_all_sigma_real_lattice_2
  end interface ed_get_Sreal

  interface ed_get_SAmats
     module procedure ed_get_self_matsubara_1
     module procedure ed_get_self_matsubara_2
     module procedure ed_get_self_matsubara_lattice_1
     module procedure ed_get_self_matsubara_lattice_2
  end interface ed_get_SAmats

  interface ed_get_SAreal
     module procedure ed_get_self_real_1
     module procedure ed_get_self_real_2
     module procedure ed_get_self_real_lattice_1
     module procedure ed_get_self_real_lattice_2
  end interface ed_get_SAreal





  interface ed_get_sigma_matsubara
     module procedure ed_get_sigma_matsubara_1
     module procedure ed_get_sigma_matsubara_2
     module procedure ed_get_sigma_matsubara_lattice_1
     module procedure ed_get_sigma_matsubara_lattice_2
     !
     module procedure ed_get_all_sigma_matsubara_1
     module procedure ed_get_all_sigma_matsubara_2
     module procedure ed_get_all_sigma_matsubara_lattice_1
     module procedure ed_get_all_sigma_matsubara_lattice_2
  end interface ed_get_sigma_matsubara

  interface ed_get_self_matsubara
     module procedure ed_get_self_matsubara_1
     module procedure ed_get_self_matsubara_2
     module procedure ed_get_self_matsubara_lattice_1
     module procedure ed_get_self_matsubara_lattice_2
  end interface ed_get_self_matsubara




  interface ed_get_sigma_real
     module procedure ed_get_sigma_real_1
     module procedure ed_get_sigma_real_2
     module procedure ed_get_sigma_real_lattice_1
     module procedure ed_get_sigma_real_lattice_2
     !
     module procedure ed_get_all_sigma_real_1
     module procedure ed_get_all_sigma_real_2
     module procedure ed_get_all_sigma_real_lattice_1
     module procedure ed_get_all_sigma_real_lattice_2
  end interface ed_get_sigma_real

  interface ed_get_self_real
     module procedure ed_get_self_real_1
     module procedure ed_get_self_real_2
     module procedure ed_get_self_real_lattice_1
     module procedure ed_get_self_real_lattice_2
  end interface ed_get_self_real




  !Retrieve imp GF through routines.
  interface ed_get_Gmats
     module procedure ed_get_gimp_matsubara_1
     module procedure ed_get_gimp_matsubara_2
     module procedure ed_get_gimp_matsubara_lattice_1
     module procedure ed_get_gimp_matsubara_lattice_2
     !
     module procedure ed_get_all_gimp_matsubara_1
     module procedure ed_get_all_gimp_matsubara_2
     module procedure ed_get_all_gimp_matsubara_lattice_1
     module procedure ed_get_all_gimp_matsubara_lattice_2
  end interface ed_get_Gmats

  interface ed_get_Fmats
     module procedure ed_get_fimp_matsubara_1
     module procedure ed_get_fimp_matsubara_2
     module procedure ed_get_fimp_matsubara_lattice_1
     module procedure ed_get_fimp_matsubara_lattice_2
  end interface ed_get_Fmats



  interface ed_get_Greal
     module procedure ed_get_gimp_real_1
     module procedure ed_get_gimp_real_2
     module procedure ed_get_gimp_real_lattice_1
     module procedure ed_get_gimp_real_lattice_2
     !
     module procedure ed_get_all_gimp_real_1
     module procedure ed_get_all_gimp_real_2
     module procedure ed_get_all_gimp_real_lattice_1
     module procedure ed_get_all_gimp_real_lattice_2
  end interface ed_get_Greal

  interface ed_get_Freal
     module procedure ed_get_fimp_real_1
     module procedure ed_get_fimp_real_2
     module procedure ed_get_fimp_real_lattice_1
     module procedure ed_get_fimp_real_lattice_2
  end interface ed_get_Freal




  interface ed_get_gimp_matsubara
     module procedure ed_get_gimp_matsubara_1
     module procedure ed_get_gimp_matsubara_2
     module procedure ed_get_gimp_matsubara_lattice_1
     module procedure ed_get_gimp_matsubara_lattice_2
     !
     module procedure ed_get_all_gimp_matsubara_1
     module procedure ed_get_all_gimp_matsubara_2
     module procedure ed_get_all_gimp_matsubara_lattice_1
     module procedure ed_get_all_gimp_matsubara_lattice_2
  end interface ed_get_gimp_matsubara

  interface ed_get_fimp_matsubara
     module procedure ed_get_fimp_matsubara_1
     module procedure ed_get_fimp_matsubara_2
     module procedure ed_get_fimp_matsubara_lattice_1
     module procedure ed_get_fimp_matsubara_lattice_2
  end interface ed_get_fimp_matsubara




  interface ed_get_gimp_real
     module procedure ed_get_gimp_real_1
     module procedure ed_get_gimp_real_2
     module procedure ed_get_gimp_real_lattice_1
     module procedure ed_get_gimp_real_lattice_2
     !
     module procedure ed_get_all_gimp_real_1
     module procedure ed_get_all_gimp_real_2
     module procedure ed_get_all_gimp_real_lattice_1
     module procedure ed_get_all_gimp_real_lattice_2
  end interface ed_get_gimp_real

  interface ed_get_fimp_real
     module procedure ed_get_fimp_real_1
     module procedure ed_get_fimp_real_2
     module procedure ed_get_fimp_real_lattice_1
     module procedure ed_get_fimp_real_lattice_2
  end interface ed_get_fimp_real


  !Retrieve static common observables  
  interface ed_get_dens
     module procedure ed_get_dens_1
     module procedure ed_get_dens_2
     module procedure ed_get_dens_lattice_1
     module procedure ed_get_dens_lattice_2
  end interface ed_get_dens

  interface ed_get_mag
     module procedure ed_get_mag_1
     module procedure ed_get_mag_2
     module procedure ed_get_mag_lattice_1
     module procedure ed_get_mag_lattice_2
  end interface ed_get_mag

  interface ed_get_docc
     module procedure ed_get_docc_1
     module procedure ed_get_docc_2
     module procedure ed_get_docc_lattice_1
     module procedure ed_get_docc_lattice_2
  end interface ed_get_docc

  interface ed_get_phisc
     module procedure ed_get_phisc_1
     module procedure ed_get_phisc_2
     module procedure ed_get_phisc_lattice_1
     module procedure ed_get_phisc_lattice_2
  end interface ed_get_phisc

  interface ed_get_eimp
     module procedure :: ed_get_eimp_
     module procedure :: ed_get_eimp_lattice
  end interface ed_get_eimp

  interface ed_get_epot
     module procedure :: ed_get_epot_
     module procedure :: ed_get_epot_lattice
  end interface ed_get_epot

  interface ed_get_eint
     module procedure :: ed_get_eint_
     module procedure :: ed_get_eint_lattice
  end interface ed_get_eint

  interface ed_get_ehartree
     module procedure :: ed_get_ehartree_
     module procedure :: ed_get_ehartree_lattice
  end interface ed_get_ehartree

  interface ed_get_eknot
     module procedure :: ed_get_eknot_
     module procedure :: ed_get_eknot_lattice
  end interface ed_get_eknot

  interface ed_get_doubles
     module procedure :: ed_get_doubles_
     module procedure :: ed_get_doubles_lattice
  end interface ed_get_doubles

  interface ed_get_dust
     module procedure :: ed_get_dust_
     module procedure :: ed_get_dust_lattice
  end interface ed_get_dust

  interface ed_get_dund
     module procedure :: ed_get_dund_
     module procedure :: ed_get_dund_lattice
  end interface ed_get_dund

  interface ed_get_dse
     module procedure :: ed_get_dse_
     module procedure :: ed_get_dse_lattice
  end interface ed_get_dse

  interface ed_get_dph
     module procedure :: ed_get_dph_
     module procedure :: ed_get_dph_lattice
  end interface ed_get_dph

  interface ed_get_density_matrix
     module procedure :: ed_get_density_matrix_single
     module procedure :: ed_get_density_matrix_lattice
  end interface ed_get_density_matrix

  public :: ed_get_Smats
  public :: ed_get_SAmats
  public :: ed_get_Sreal
  public :: ed_get_SAreal

  public :: ed_get_sigma_matsubara
  public :: ed_get_self_matsubara
  public :: ed_get_sigma_real
  public :: ed_get_self_real


  public :: ed_get_Gmats
  public :: ed_get_Fmats
  public :: ed_get_Greal
  public :: ed_get_Freal

  public :: ed_get_gimp_matsubara
  public :: ed_get_fimp_matsubara
  public :: ed_get_gimp_real
  public :: ed_get_fimp_real

  public :: ed_get_dens
  public :: ed_get_mag
  public :: ed_get_docc
  public :: ed_get_phisc

  public :: ed_get_eimp
  public :: ed_get_epot
  public :: ed_get_eint 
  public :: ed_get_ehartree
  public :: ed_get_eknot

  public :: ed_get_doubles
  public :: ed_get_dust
  public :: ed_get_dund
  public :: ed_get_dse
  public :: ed_get_dph

  public :: ed_get_density_matrix
  public :: ed_get_quantum_SOC_operators_single
  public :: ed_get_quantum_SOC_operators_lattice

  public :: ed_read_impSigma_single
  public :: ed_read_impSigma_lattice

  !****************************************************************************************!
  !****************************************************************************************!

  public :: ed_print_impSigma
  public :: ed_print_impG
  public :: ed_print_impG0
  public :: ed_print_impChi


  !****************************************************************************************!
  !****************************************************************************************!



  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable :: wm,tau,wr,vm
  character(len=64)                :: suffix





contains




  !+------------------------------------------------------------------+
  !PURPOSE  : Print impurity Functions case:
  ! - impSigma
  ! - impG
  ! - impG0
  !+------------------------------------------------------------------+
  include "ED_IO/print_impSigma.f90"
  subroutine ed_print_impSigma
    call print_impSigma_normal
  end subroutine ed_print_impSigma


  include "ED_IO/print_impG.f90"
  subroutine ed_print_impG
    call print_impG_normal
  end subroutine ed_print_impG


  include "ED_IO/print_impG0.f90"
  subroutine ed_print_impG0
    call print_impG0_normal
  end subroutine ed_print_impG0


  include "ED_IO/print_impChi.f90"
  subroutine ed_print_impChi
    call print_chi_spin
    call print_chi_dens
    call print_chi_dens_mix
    call print_chi_dens_tot
    call print_chi_pair
  end subroutine ed_print_impChi




  !+-----------------------------------------------------------------------------+!



  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the impurity self-energy 
  !+-----------------------------------------------------------------------------+!
  include "ED_IO/get_all_sigma_matsubara.f90"
  include "ED_IO/get_all_sigma_realaxis.f90"
  include "ED_IO/get_sigma_matsubara.f90"
  include "ED_IO/get_self_matsubara.f90"
  include "ED_IO/get_sigma_realaxis.f90"
  include "ED_IO/get_self_realaxis.f90"


  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the impurity green's functions 
  !+-----------------------------------------------------------------------------+!
  include "ED_IO/get_all_gimp_matsubara.f90"
  include "ED_IO/get_all_gimp_realaxis.f90"
  include "ED_IO/get_gimp_matsubara.f90"
  include "ED_IO/get_fimp_matsubara.f90"
  include "ED_IO/get_gimp_realaxis.f90"
  include "ED_IO/get_fimp_realaxis.f90"


  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the local observables
  !+-----------------------------------------------------------------------------+!
  include "ED_IO/get_dens.f90"
  include "ED_IO/get_mag.f90"
  include "ED_IO/get_docc.f90"
  include "ED_IO/get_phisc.f90"
  include "ED_IO/get_eimp.f90"
  include "ED_IO/get_doubles.f90"
  !
  include "ED_IO/get_imp_dm.f90"
  include "ED_IO/get_imp_SOC_op.f90"






  ! PURPOSE: Read self-energy function(s) - also for inequivalent sites.
  !+-----------------------------------------------------------------------------+!
  include "ED_IO/read_impSigma.f90"
  subroutine ed_read_impSigma_single
    !
    if(allocated(impSmats))deallocate(impSmats)
    if(allocated(impSreal))deallocate(impSreal)
    if(allocated(impSAmats))deallocate(impSAmats)
    if(allocated(impSAreal))deallocate(impSAreal)
    allocate(impSmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impSreal(Nspin,Nspin,Norb,Norb,Lreal))
    allocate(impSAmats(Nspin,Nspin,Norb,Norb,Lmats)) !THIS SHOULD NOT DEPEND ON SPIN: NSPIN=>1
    allocate(impSAreal(Nspin,Nspin,Norb,Norb,Lreal)) !THIS SHOULD NOT DEPEND ON SPIN: NSPIN=>1
    impSmats=zero
    impSreal=zero
    impSAmats=zero
    impSAreal=zero
    !
    call read_impSigma_normal
  end subroutine ed_read_impSigma_single

  subroutine ed_read_impSigma_lattice(Nineq)
    integer :: Nineq
    integer :: ilat
    !
    if(allocated(Smatsii))deallocate(Smatsii)
    if(allocated(Srealii))deallocate(Srealii)
    if(allocated(SAmatsii))deallocate(SAmatsii)
    if(allocated(SArealii))deallocate(SArealii)
    allocate(Smatsii(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Srealii(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(SAmatsii(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(SArealii(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    Smatsii  = zero 
    Srealii  = zero 
    SAmatsii = zero 
    SArealii = zero
    !
    do ilat=1,Nineq
       ed_file_suffix=reg(ineq_site_suffix)//str(ilat,site_indx_padding)
       call ed_read_impSigma_single
       Smatsii(ilat,:,:,:,:,:)  = impSmats
       Srealii(ilat,:,:,:,:,:)  = impSreal
       SAmatsii(ilat,:,:,:,:,:) = impSAmats
       SArealii(ilat,:,:,:,:,:) = impSAreal
    enddo
    ed_file_suffix=""
  end subroutine ed_read_impSigma_lattice


END MODULE ED_IO







