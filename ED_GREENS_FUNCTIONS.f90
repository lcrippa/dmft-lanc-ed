MODULE ED_GREENS_FUNCTIONS
  USE ED_GF_SHARED
  USE ED_GF_NORMAL
  USE ED_GF_PHONON
  USE ED_GF_CHISPIN
  USE ED_GF_CHIDENS
  ! USE ED_CHI_PAIR
  !
  implicit none
  private 

  public :: buildGf_impurity
  public :: buildChi_impurity

contains



  !+------------------------------------------------------------------+
  ! GF CALCULATIONS
  !+------------------------------------------------------------------+
  subroutine buildGF_impurity()
    !
    call allocate_grids
    !
    impGmats=zero
    impGreal=zero
    !
    impSmats = zero
    impSreal = zero
    !
    impG0mats=zero
    impG0real=zero
    !
    impDmats_ph=zero
    impDreal_ph=zero
    !
    write(LOGfile,"(A)")"Get impurity Greens functions:"
    call build_gf_normal()
    if(DimPh>1) call build_gf_phonon()
    call build_sigma_normal()
    !
    if(MPIMASTER)then
       if(ed_print_Sigma)call ed_print_impSigma()
       if(ed_print_G) then
          call ed_print_impG()
          if(DimPh>1)call ed_print_impD()
       endif
       if(ed_print_G0)call ed_print_impG0()
    endif
    !
    call deallocate_grids
    !
  end subroutine buildgf_impurity








  !+------------------------------------------------------------------+
  ! SUSCEPTIBILITY CALCULATIONS
  !+------------------------------------------------------------------+
  subroutine buildChi_impurity()
    !
    call allocate_grids
    !
    !
    !BUILD SPIN SUSCEPTIBILITY
    spinChi_tau=zero
    spinChi_w=zero
    spinChi_iv=zero
    if(chispin_flag)call build_chi_spin()
    !
    !
    ! !BUILD CHARGE SUSCEPTIBILITY
    densChi_tau=zero
    densChi_w=zero
    densChi_iv=zero
    if(chidens_flag)call build_chi_dens()
    !
    !
    ! !BUILD PAIR SUSCEPTIBILITY
    ! pairChi_tau=zero
    ! pairChi_w=zero
    ! pairChi_iv=zero
    !if(chipair_flag)call build_chi_pair()
    !
    !
    !PRINTING:
    if(MPIMASTER.AND.(any([chispin_flag,chidens_flag,chipair_flag])))call ed_print_impChi()
    !
    call deallocate_grids
    !
  end subroutine buildChi_impurity




end MODULE ED_GREENS_FUNCTIONS
