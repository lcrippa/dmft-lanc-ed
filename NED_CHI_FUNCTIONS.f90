MODULE ED_CHI_FUNCTIONS
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,free_unit,reg,free_units,txtfy,splot
  USE SF_ARRAYS,  only: arange,linspace
  USE SF_LINALG,  only: inv,inv_sym,inv_her,eye
  USE SF_SP_LINALG, only: sp_lanc_tridiag
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_IO                     !< this contains the routine to print Chi
  USE ED_EIGENSPACE
  USE ED_BATH_FUNCTIONS
  USE ED_SETUP
  USE ED_HAMILTONIAN_MATVEC
  USE ED_AUX_FUNX
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  !
  implicit none
  private 



  public :: buildChi_impurity

  public :: ed_chi_functions_set_MPI

  public :: ed_chi_functions_del_MPI



  !Lanczos shared variables
  !=========================================================
  real(8),dimension(:),pointer                :: state_vec
  complex(8),dimension(:),pointer             :: state_cvec
  real(8)                                     :: state_e

  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable            :: wm,tau,wr,vm


#ifdef _MPI
  integer                                     :: MpiComm=MPI_UNDEFINED
#else
  integer                                     :: MpiComm=0
#endif
  logical                                     :: MpiStatus=.false.  
  integer                                     :: MPI_RANK=0
  integer                                     :: MPI_SIZE=1
  logical                                     :: MPI_MASTER=.true.
  integer                                     :: mpi_ierr



contains


  subroutine ed_chi_functions_set_MPI(comm)
#ifdef _MPI
    integer :: comm
    MpiComm  = comm
    MpiStatus = .true.
    MPI_RANK = get_Rank_MPI(MpiComm)
    MPI_SIZE = get_Size_MPI(MpiComm)
    MPI_MASTER=get_Master_MPI(MpiComm)
#else
    integer,optional :: comm
#endif
  end subroutine ed_chi_functions_set_MPI


  subroutine ed_chi_functions_del_MPI()
#ifdef _MPI
    MpiComm  = MPI_UNDEFINED
    MpiStatus = .false.
    MPI_RANK=0
    MPI_SIZE=1
    MPI_MASTER=.true.
#endif
  end subroutine ed_chi_functions_del_MPI





  !+------------------------------------------------------------------+
  ! SUSCEPTIBILITY CALCULATIONS
  !+------------------------------------------------------------------+
  subroutine buildChi_impurity()
    integer :: i
    !
    call allocate_grids
    !
    !
    !BUILD SPIN SUSCEPTIBILITY
    if(.not.allocated(spinChi_tau)) stop "buildChi_impurity: spinChi_tau not allocated"
    if(.not.allocated(spinChi_w))  stop "buildChi_impurity: spinChi_w not allocated"
    if(.not.allocated(spinChi_iv)) stop "buildChi_impurity: spinChi_iv not allocated"
    spinChi_tau=zero
    spinChi_w=zero
    spinChi_iv=zero
    call build_chi_spin()
    !
    !BUILD CHARGE SUSCEPTIBILITY
    if(.not.allocated(densChi_tau)) stop "buildChi_impurity: densChi_tau not allocated"
    if(.not.allocated(densChi_w))  stop "buildChi_impurity: densChi_w not allocated"
    if(.not.allocated(densChi_iv)) stop "buildChi_impurity: densChi_iv not allocated"
    if(.not.allocated(densChi_mix_tau))stop "buildChi_impurity: densChi_mix_tau not allocated"
    if(.not.allocated(densChi_mix_w))  stop "buildChi_impurity: densChi_mix_w not allocated"
    if(.not.allocated(densChi_mix_iv)) stop "buildChi_impurity: densChi_mix_iv not allocated"
    if(.not.allocated(densChi_tot_tau))stop "buildChi_impurity: densChi_tot_tau not allocated"
    if(.not.allocated(densChi_tot_w))  stop "buildChi_impurity: densChi_tot_w not allocated"
    if(.not.allocated(densChi_tot_iv)) stop "buildChi_impurity: densChi_tot_iv not allocated"
    densChi_tau=zero
    densChi_w=zero
    densChi_iv=zero
    densChi_mix_tau=zero
    densChi_mix_w=zero
    densChi_mix_iv=zero
    densChi_tot_tau=zero
    densChi_tot_w=zero
    densChi_tot_iv=zero
    call build_chi_dens()
    !
    !BUILD PAIR SUSCEPTIBILITY
    if(.not.allocated(pairChi_tau))stop "buildChi_impurity: pairChi_tau not allocated"
    if(.not.allocated(pairChi_w))  stop "buildChi_impurity: pairChi_w not allocated"
    if(.not.allocated(pairChi_iv)) stop "buildChi_impurity: pairChi_iv not allocated"
    pairChi_tau=zero
    pairChi_w=zero
    pairChi_iv=zero
    call build_chi_pair()
    !
    !PRINTING:
    if(MPI_MASTER)call ed_print_impChi()
    !
    !
    call deallocate_grids
  end subroutine buildChi_impurity

  !+------------------------------------------------------------------+
  !                    SUSCEPTIBILITIES (SPIN, CHARGE, PAIR)
  !+------------------------------------------------------------------+
  ! include 'ED_GREENS_FUNCTIONS/build_chi_spin.f90'
  ! include 'ED_GREENS_FUNCTIONS/build_chi_dens.f90'
  ! include 'ED_GREENS_FUNCTIONS/build_chi_pair.f90'




  !+------------------------------------------------------------------+
  !                            SPIN
  !+------------------------------------------------------------------+
  !PURPOSE  : Evaluate the Spin susceptibility \Chi_spin for a 
  ! single orbital: \chi = <S_a(\tau)S_a(0)>
  ! note: as S_a is hermitian particle and holes contributions (isign=1,-1)
  ! are identical so work out only one lanczos tridiag. work out the 
  ! reduction for both values of isign in the same call.
  !+------------------------------------------------------------------+
  subroutine build_chi_spin()
    integer :: iorb
    write(LOGfile,"(A)")"Get impurity spin Chi:"
    do iorb=1,Norb
       write(LOGfile,"(A)")"Get Chi_spin_l"//reg(txtfy(iorb))
       if(MPI_MASTER)call start_timer()
       call lanc_ed_build_spinChi_c(iorb)
       if(MPI_MASTER)call stop_timer(LOGfile)
    enddo
    if(Norb>1)then
       write(LOGfile,"(A)")"Get Chi_spin_tot"
       if(MPI_MASTER)call start_timer()
       call lanc_ed_build_spinChi_tot_c()
       if(MPI_MASTER)call stop_timer(LOGfile)
    endif
    spinChi_tau = SpinChi_tau/zeta_function
    spinChi_w   = spinChi_w/zeta_function
    spinChi_iv  = spinChi_iv/zeta_function
  end subroutine build_chi_spin

  subroutine lanc_ed_build_spinChi_c(iorb)
    integer                          :: iorb,isite,isector,izero
    integer                          :: numstates
    integer                          :: nlanc,idim
    integer                          :: iup0,idw0,isign
    integer                          :: ib(Nlevels)
    integer                          :: m,i,j,r
    real(8)                          :: norm0,sgn
    real(8),allocatable              :: alfa_(:),beta_(:)
    complex(8),allocatable           :: vvinit(:)
    integer                          :: Nitermax
    type(sector_map) :: HI    !map of the Sector S to Hilbert space H
    !
    !
    
    !
    do izero=1,state_list%size
       isector     =  es_return_sector(state_list,izero)
       idim      =  getdim(isector)
       state_e    =  es_return_energy(state_list,izero)
       state_cvec => es_return_cvector(state_list,izero)
       norm0=sqrt(dot_product(state_cvec,state_cvec))
       if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
       allocate(vvinit(idim))
       if(ed_verbose==3)write(LOGfile,"(A,2I3)")'Apply Sz:',getnup(isector),getndw(isector)
       call build_sector(isector,HI)
       vvinit=0.d0
       do m=1,idim                     !loop over |gs> components m
          i=HI%map(m)
          ib = bdecomp(i,2*Ns)
          sgn = dble(ib(iorb))-dble(ib(iorb+Ns))
          vvinit(m) = 0.5d0*sgn*state_cvec(m)   !build the cdg_up|gs> state
       enddo
       deallocate(HI%map)
       norm0=sqrt(dot_product(vvinit,vvinit))
       vvinit=vvinit/norm0
       !
       call setup_Hv_sector(isector)
       if(ed_sparse_H)call ed_buildH_c()
       !
       nlanc=min(idim,lanc_nGFiter)
       allocate(alfa_(nlanc),beta_(nlanc))
#ifdef _MPI
       if(MpiStatus)then
          call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
       else
          call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
       endif
#else
       call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
#endif
       !particles
       isign=1
       call add_to_lanczos_spinChi(norm0,state_e,alfa_,beta_,isign,iorb)
       !holes
       isign=-1
       call add_to_lanczos_spinChi(norm0,state_e,alfa_,beta_,isign,iorb)
       !
       call delete_Hv_sector()
       !
       deallocate(vvinit,alfa_,beta_)
       if(spH0%status)call sp_delete_matrix(spH0)
       nullify(state_cvec)
    enddo
    
  end subroutine lanc_ed_build_spinChi_c

  subroutine lanc_ed_build_spinChi_tot_c()
    integer                          :: iorb,isite,isector,izero
    integer                          :: numstates
    integer                          :: nlanc,idim
    integer                          :: iup0,idw0,isign
    integer                          :: ib(Nlevels)
    integer                          :: m,i,j,r
    real(8)                          :: norm0,sgn
    real(8),allocatable              :: alfa_(:),beta_(:)
    complex(8),allocatable           :: vvinit(:)
    integer                          :: Nitermax
    type(sector_map) :: HI    !map of the Sector S to Hilbert space H
    !
    !
    
    !
    do izero=1,state_list%size
       isector     =  es_return_sector(state_list,izero)
       idim       =  getdim(isector)
       state_e    =  es_return_energy(state_list,izero)
       state_cvec => es_return_cvector(state_list,izero)
       norm0=sqrt(dot_product(state_cvec,state_cvec))
       if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
       allocate(vvinit(idim))
       if(ed_verbose==3)write(LOGfile,"(A,2I3)")'Apply Sz:',getnup(isector),getndw(isector)
       call build_sector(isector,HI)
       vvinit=0.d0
       do m=1,idim  
          i=HI%map(m)
          ib = bdecomp(i,2*Ns)
          sgn = sum(dble(ib(1:Norb)))-sum(dble(ib(Ns+1:Ns+Norb)))
          vvinit(m) = 0.5d0*sgn*state_cvec(m) 
       enddo
       deallocate(HI%map)
       norm0=sqrt(dot_product(vvinit,vvinit))
       vvinit=vvinit/norm0
       !
       call setup_Hv_sector(isector)
       if(ed_sparse_H)call ed_buildH_c()
       !
       nlanc=min(idim,lanc_nGFiter)
       allocate(alfa_(nlanc),beta_(nlanc))
#ifdef _MP
       if(MpiStatus)then
          call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
       else
          call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
       endif
#else
       call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
#endif
       !particles
       isign=1
       call add_to_lanczos_spinChi(norm0,state_e,alfa_,beta_,isign,Norb+1)
       !holes
       isign=-1
       call add_to_lanczos_spinChi(norm0,state_e,alfa_,beta_,isign,Norb+1)
       !
       call delete_Hv_sector()
       !
       deallocate(vvinit,alfa_,beta_)
       if(spH0%status)call sp_delete_matrix(spH0)
       nullify(state_cvec)
    enddo
    
  end subroutine lanc_ed_build_spinChi_tot_c

  subroutine add_to_lanczos_spinChi(vnorm,Ei,alanc,blanc,isign,iorb)
    real(8)                                    :: vnorm,Ei,Ej,Egs,pesoF,pesoAB,pesoBZ,de,peso
    integer                                    :: nlanc
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: isign,iorb
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw,chisp
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    pesoF  = vnorm**2/zeta_function 
    pesoBZ = 1d0
    if(finiteT)pesoBZ = exp(-beta*(Ei-Egs))
    !
    diag             = 0.d0
    subdiag          = 0.d0
    Z                = eye(Nlanc)
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call tql2(Nlanc,diag,subdiag,Z,ierr)
    !
    select case(isign)
    case (1)
       do j=1,nlanc
          Ej     = diag(j)
          dE     = Ej-Ei
          pesoAB = Z(1,j)*Z(1,j)
          peso   = pesoF*pesoAB*pesoBZ
          if(beta*dE < 1d-1)then     !abs(X - (1-exp(-X)) is about 5*10^-3 for X<10^-1 this is a satisfactory bound
             spinChi_iv(iorb,0)=spinChi_iv(iorb,0) + peso*beta
          else
             spinChi_iv(iorb,0)=spinChi_iv(iorb,0) + peso*(1d0-exp(-beta*dE))/dE 
          endif
          do i=1,Lmats
             spinChi_iv(iorb,i)=spinChi_iv(iorb,i) + peso*(exp(-beta*dE)-1d0)/(dcmplx(0d0,vm(i)) - dE)
          enddo
          do i=0,Ltau
             spinChi_tau(iorb,i)=spinChi_tau(iorb,i) + peso*exp(-tau(i)*de)
          enddo
          do i=1,Lreal
             spinChi_w(iorb,i)=spinChi_w(iorb,i) + peso*(exp(-beta*dE)-1.d0)/(dcmplx(wr(i),eps) - dE)
          enddo
       enddo
    case (-1)
       do j=1,nlanc
          Ej     = diag(j)
          dE     = Ej-Ei
          pesoAB = Z(1,j)*Z(1,j)
          peso   = pesoF*pesoAB*pesoBZ
          if(beta*dE < 1d-1)then     !abs(X - (1-exp(-X)) is about 5*10^-3 for X<10^-1 this is a satisfactory bound
             spinChi_iv(iorb,0)=spinChi_iv(iorb,0) + peso*beta
          else
             spinChi_iv(iorb,0)=spinChi_iv(iorb,0) + peso*(1d0-exp(-beta*dE))/dE 
          endif
          do i=1,Lmats
             spinChi_iv(iorb,i)=spinChi_iv(iorb,i) + peso*(1d0-exp(-beta*dE))/(dcmplx(0d0,vm(i)) + dE)
          enddo
          do i=0,Ltau
             spinChi_tau(iorb,i)=spinChi_tau(iorb,i) + peso*exp(-(beta-tau(i))*dE)
          enddo
          do i=1,Lreal
             spinChi_w(iorb,i)=spinChi_w(iorb,i) + peso*(1d0-exp(-beta*dE))/(dcmplx(wr(i),eps) + dE)
          enddo
       enddo
    case default
       stop "add_to_lanczos_spinChi: isign not in {-1,1}"
    end select
  end subroutine add_to_lanczos_spinChi














  !+------------------------------------------------------------------+
  !                            CHARGE
  !+------------------------------------------------------------------+
  !PURPOSE  : Evaluate Charge-Charge Susceptibility <n_a(tau)n_b(0)>
  !+------------------------------------------------------------------+
  subroutine build_chi_dens()
    integer :: iorb,jorb
    write(LOGfile,"(A)")"Get impurity dens Chi:"
    do iorb=1,Norb
       write(LOGfile,"(A)")"Get Chi_dens_diag_l"//reg(txtfy(iorb))
       if(MPI_MASTER)call start_timer()
       call lanc_ed_build_densChi_diag_c(iorb)
       if(MPI_MASTER)call stop_timer(LOGfile)
    enddo
    !
    !
    if(Norb>1)then
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             write(LOGfile,"(A)")"Get Chi_dens_offdiag_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
             if(MPI_MASTER)call start_timer()
             call lanc_ed_build_densChi_offdiag_c(iorb,jorb)
             if(MPI_MASTER)call stop_timer(LOGfile)
          end do
       end do
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             denschi_w(iorb,jorb,:) = 0.5d0*( denschi_w(iorb,jorb,:) - (one+xi)*denschi_w(iorb,iorb,:) - (one+xi)*denschi_w(jorb,jorb,:))
          enddo
       enddo
       !
       do iorb=1,Norb
          do jorb=1,Norb
             write(LOGfile,"(A)")"Get Chi_dens_offdiag_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
             if(MPI_MASTER)call start_timer()
             call lanc_ed_build_densChi_mix_c(iorb,jorb)
             if(MPI_MASTER)call stop_timer(LOGfile)
          end do
       end do
       !
       write(LOGfile,"(A)")"Get Chi_dens_tot"
       if(MPI_MASTER)call start_timer()
       call lanc_ed_build_densChi_tot_c()
       if(MPI_MASTER)call stop_timer(LOGfile)
    endif
    !
    denschi_tau = Denschi_tau/zeta_function
    denschi_w   = denschi_w/zeta_function
    denschi_iv  = denschi_iv/zeta_function
    !
  end subroutine build_chi_dens

  !+------------------------------------------------------------------+
  !PURPOSE  : Evaluate the Charge-Charge susceptibility \Chi_dens for  
  ! the orbital diagonal case: \chi_dens_aa = <N_a(\tau)N_a(0)>
  !+------------------------------------------------------------------+
  subroutine lanc_ed_build_densChi_diag_c(iorb)
    integer                :: iorb,isite,isector,izero
    integer                :: numstates
    integer                :: nlanc,idim
    integer                :: iup0,idw0,isign
    integer                :: ib(Nlevels)
    integer                :: m,i,j,r
    real(8)                :: norm0,sgn
    complex(8)             :: cnorm2
    real(8),allocatable    :: alfa_(:),beta_(:)
    complex(8),allocatable :: vvinit(:)
    integer                :: Nitermax
    type(sector_map)       :: HI    !map of the Sector S to Hilbert space H
    !
    !
    
    !
    do izero=1,state_list%size
       isector    =  es_return_sector(state_list,izero)
       idim       =  getdim(isector)
       state_e    =  es_return_energy(state_list,izero)
       state_cvec => es_return_cvector(state_list,izero)
       norm0=sqrt(dot_product(state_cvec,state_cvec))
       if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
       allocate(vvinit(idim))
       if(ed_verbose==3)write(LOGfile,"(A,2I3)")'Apply N:',getnup(isector),getndw(isector)
       call build_sector(isector,HI)
       vvinit=zero
       do m=1,idim                     !loop over |gs> components m
          i=HI%map(m)
          ib = bdecomp(i,2*Ns)
          sgn = dble(ib(iorb))+dble(ib(iorb+Ns))
          vvinit(m) = sgn*state_cvec(m)   !build the cdg_up|gs> state
       enddo
       deallocate(HI%map)
       norm0=dot_product(vvinit,vvinit)
       vvinit=vvinit/sqrt(norm0)
       !
       call setup_Hv_sector(isector)
       if(ed_sparse_H)call ed_buildH_c()
       !
       nlanc=min(idim,lanc_nGFiter)
       allocate(alfa_(nlanc),beta_(nlanc))
#ifdef _MPI
       if(MpiStatus)then
          call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
       else
          call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
       endif
#else
       call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
#endif
       cnorm2=one*norm0
       isign=1
       call add_to_lanczos_densChi(cnorm2,state_e,alfa_,beta_,isign,iorb,iorb)
       isign=-1
       call add_to_lanczos_densChi(cnorm2,state_e,alfa_,beta_,isign,iorb,iorb)
       !
       call delete_Hv_sector()
       !
       deallocate(vvinit,alfa_,beta_)
       if(spH0%status)call sp_delete_matrix(spH0)
       nullify(state_cvec)
    enddo
    
  end subroutine lanc_ed_build_densChi_diag_c

  !+------------------------------------------------------------------+
  !PURPOSE  : Evaluate the TOTAL Charge-Charge susceptibility \Chi_dens  
  ! \chi_dens_tot = <N(\tau)N(0)>, N=sum_a N_a
  !+------------------------------------------------------------------+
  subroutine lanc_ed_build_densChi_tot_c()
    integer                :: iorb,isite,isector,izero
    integer                :: numstates
    integer                :: nlanc,idim
    integer                :: iup0,idw0,isign
    integer                :: ib(Nlevels)
    integer                :: m,i,j,r
    complex(8)             :: cnorm2
    real(8)                :: norm0,sgn
    real(8),allocatable    :: alfa_(:),beta_(:)
    complex(8),allocatable :: vvinit(:)
    integer                :: Nitermax
    type(sector_map)       :: HI    !map of the Sector S to Hilbert space H
    !
    
    !
    do izero=1,state_list%size
       isector    =  es_return_sector(state_list,izero)
       idim       =  getdim(isector)
       state_e    =  es_return_energy(state_list,izero)
       state_cvec => es_return_cvector(state_list,izero)
       norm0=sqrt(dot_product(state_cvec,state_cvec))
       if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
       allocate(vvinit(idim))
       if(ed_verbose==3)write(LOGfile,"(A,2I3)")'Apply N:',getnup(isector),getndw(isector)
       call build_sector(isector,HI)
       vvinit=zero
       do m=1,idim
          i=HI%map(m)
          ib = bdecomp(i,2*Ns)
          sgn = sum(dble(ib(1:Norb)))+sum(dble(ib(Ns+1:Ns+Norb)))
          vvinit(m) = sgn*state_cvec(m) 
       enddo
       deallocate(HI%map)
       norm0=dot_product(vvinit,vvinit)
       vvinit=vvinit/sqrt(norm0)
       !
       call setup_Hv_sector(isector)
       if(ed_sparse_H)call ed_buildH_c()
       !
       nlanc=min(idim,lanc_nGFiter)
       allocate(alfa_(nlanc),beta_(nlanc))
#ifdef _MPI
       if(MpiStatus)then
          call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
       else
          call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
       endif
#else
       call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
#endif
       cnorm2=one*norm0
       isign=1
       call add_to_lanczos_densChi_tot(cnorm2,state_e,alfa_,beta_,isign)
       isign=-1
       call add_to_lanczos_densChi_tot(cnorm2,state_e,alfa_,beta_,isign)
       !
       call delete_Hv_sector()
       !
       deallocate(vvinit,alfa_,beta_)
       if(spH0%status)call sp_delete_matrix(spH0)
       nullify(state_cvec)
    enddo
    
  end subroutine lanc_ed_build_densChi_tot_c

  !+------------------------------------------------------------------+
  !PURPOSE  : Evaluate the Charge-Charge susceptibility \Chi_dens for
  ! the orbital off-diagonal case: \chi_dens_ab = <N_a(\tau)N_b(0)>
  !+------------------------------------------------------------------+
  subroutine lanc_ed_build_densChi_offdiag_c(iorb,jorb)
    integer                :: iorb,jorb,isite,isector,izero,isign
    integer                :: numstates
    integer                :: nlanc,idim
    integer                :: iup0,idw0
    integer                :: ib(Nlevels)
    integer                :: m,i,j,r
    complex(8)             :: cnorm2
    real(8)                :: norm0,sgn
    real(8),allocatable    :: alfa_(:),beta_(:)
    complex(8),allocatable :: vvinit(:)
    complex(8),allocatable :: cvinit(:)
    integer                :: Nitermax
    type(sector_map)       :: HI    !map of the Sector S to Hilbert space H
    !
    
    !
    do izero=1,state_list%size
       ! properties of the ground states
       isector     =  es_return_sector(state_list,izero)
       idim        = getdim(isector)
       state_e     =  es_return_energy(state_list,izero)
       state_cvec  => es_return_cvector(state_list,izero)
       norm0=sqrt(dot_product(state_cvec,state_cvec))
       if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
       allocate(vvinit(idim),cvinit(idim))
       call build_sector(isector,HI)
       !
       !build the (N_iorb+N_jorb)|gs> state
       if(ed_verbose==3)write(LOGfile,"(A)")'Apply N_iorb + N_jorb:'
       vvinit=zero
       do m=1,idim                     !loop over |gs> components m
          i=HI%map(m)
          ib = bdecomp(i,2*Ns)
          sgn = dble(ib(iorb))+dble(ib(iorb+Ns))
          vvinit(m) = sgn*state_cvec(m)   
          !
          sgn = dble(ib(jorb))+dble(ib(jorb+Ns))
          vvinit(m) = vvinit(m) + sgn*state_cvec(m)   
          !
       enddo
       norm0=dot_product(vvinit,vvinit)
       vvinit=vvinit/sqrt(norm0)
       !
       call setup_Hv_sector(isector)
       if(ed_sparse_H)call ed_buildH_c()
       !
       nlanc=min(idim,lanc_nGFiter)
       allocate(alfa_(nlanc),beta_(nlanc))
#ifdef _MPI
       if(MpiStatus)then
          call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
       else
          call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
       endif
#else
       call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
#endif
       cnorm2=one*norm0
       !particle and holes excitations all at once
       isign=1                    !<---
       call add_to_lanczos_densChi(cnorm2,state_e,alfa_,beta_,isign,iorb,jorb)
       isign=-1                   !<---
       call add_to_lanczos_densChi(cnorm2,state_e,alfa_,beta_,isign,iorb,jorb)
       !
       call delete_Hv_sector()
       !
       !
       !build the (N_iorb - xi*N_jorb)|gs> state
       if(ed_verbose==3)write(LOGfile,"(A)")'Apply N_iorb + xi*N_jorb:'
       cvinit=zero
       do m=1,idim
          i=HI%map(m)
          ib = bdecomp(i,2*Ns)
          sgn = dble(ib(iorb))+dble(ib(iorb+Ns))
          cvinit(m) = sgn*state_cvec(m)   
          !
          sgn = dble(ib(jorb))+dble(ib(jorb+Ns))
          cvinit(m) = cvinit(m) - xi*sgn*state_cvec(m)   
          !
       enddo
       norm0=dot_product(cvinit,cvinit)
       cvinit=cvinit/sqrt(norm0)
       !
       call setup_Hv_sector(isector)
       if(ed_sparse_H)call ed_buildH_c()
       !
       nlanc=min(idim,lanc_nGFiter)
       allocate(alfa_(nlanc),beta_(nlanc))
#ifdef _MPI
       if(MpiStatus)then
          call sp_lanc_tridiag(MpiComm,spHtimesV_cc,cvinit,alfa_,beta_)
       else
          call sp_lanc_tridiag(spHtimesV_cc,cvinit,alfa_,beta_)
       endif
#else
       call sp_lanc_tridiag(spHtimesV_cc,cvinit,alfa_,beta_)
#endif     
       cnorm2=xi*norm0
       isign=1
       call add_to_lanczos_densChi(cnorm2,state_e,alfa_,beta_,isign,iorb,jorb)
       !
       call delete_Hv_sector()
       !
       !
       !build the (N_iorb + xi*N_jorb)|gs> state
       if(ed_verbose==3)write(LOGfile,"(A)")'Apply N_iorb + xi*N_jorb:'
       cvinit=zero
       do m=1,idim
          i=HI%map(m)
          ib = bdecomp(i,2*Ns)
          sgn = dble(ib(iorb))+dble(ib(iorb+Ns))
          cvinit(m) = sgn*state_cvec(m)   
          !
          sgn = dble(ib(jorb))+dble(ib(jorb+Ns))
          cvinit(m) = cvinit(m) + xi*sgn*state_cvec(m)   
          !
       enddo
       norm0=dot_product(cvinit,cvinit)
       cvinit=cvinit/sqrt(norm0)
       !
       call setup_Hv_sector(isector)
       if(ed_sparse_H)call ed_buildH_c()
       !
       nlanc=min(idim,lanc_nGFiter)
       allocate(alfa_(nlanc),beta_(nlanc))
#ifdef _MPI
       if(MpiStatus)then
          call sp_lanc_tridiag(MpiComm,spHtimesV_cc,cvinit,alfa_,beta_)
       else
          call sp_lanc_tridiag(spHtimesV_cc,cvinit,alfa_,beta_)
       endif
#else
       call sp_lanc_tridiag(spHtimesV_cc,cvinit,alfa_,beta_)
#endif
       cnorm2=xi*norm0
       isign=-1
       call add_to_lanczos_densChi(cnorm2,state_e,alfa_,beta_,isign,iorb,jorb)
       !
       call delete_Hv_sector()
       !
       deallocate(cvinit,vvinit,alfa_,beta_)
       deallocate(HI%map)
       if(spH0%status)call sp_delete_matrix(spH0)
       nullify(state_cvec)
    enddo
    
  end subroutine lanc_ed_build_densChi_offdiag_c

  !+------------------------------------------------------------------+
  !PURPOSE  : Evaluate the inter-orbital charge susceptibility \Chi_mix 
  ! \chi_mix = <C^+_a(\tau)N_a(0)>
  !+------------------------------------------------------------------+
  subroutine lanc_ed_build_densChi_mix_c(iorb,jorb)
    integer             :: iorb,jorb,ispin
    complex(8),allocatable :: vvinit(:),vvinit_tmp(:)
    real(8),allocatable :: alfa_(:),beta_(:)
    integer             :: isite,jsite,istate
    integer             :: isector,jsector,ksector
    integer             :: idim,jdim,kdim
    type(sector_map)    :: HI,HJ,HK
    integer             :: ib(Nlevels)
    integer             :: m,i,j,r,numstates
    real(8)             :: sgn,norm2,norm0
    complex(8)          :: cnorm2
    integer             :: Nitermax,Nlanc
    !
    !   
    
    !
    do istate=1,state_list%size
       isector     =  es_return_sector(state_list,istate)
       idim        = getdim(isector)
       state_e     =  es_return_energy(state_list,istate)
       state_cvec  => es_return_cvector(state_list,istate)
       norm0=sqrt(dot_product(state_cvec,state_cvec))
       if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
       idim  = getdim(isector)
       call build_sector(isector,HI)
       !
       !+- Apply Sum_ispin c^dg_{jorb,ispin} c_{iorb,ispin} -+!
       do ispin=1,Nspin
          isite=impIndex(iorb,ispin)
          jsector = getCsector(ispin,isector)
          if(jsector/=0)then
             jdim  = getdim(jsector)
             allocate(vvinit_tmp(jdim))
             call build_sector(jsector,HJ)
             vvinit_tmp=zero
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(isite)==1)then
                   call c(isite,i,r,sgn)
                   j=binary_search(HJ%map,r)
                   vvinit_tmp(j) = sgn*state_cvec(m)
                end if
             enddo
          endif
          jsite = impIndex(jorb,ispin)
          ksector = getCDGsector(ispin,jsector)
          if(ksector/=0) then       
             kdim  = getdim(ksector)
             allocate(vvinit(kdim)) !<==== ACTHUNG! 
             call build_sector(ksector,HK)
             vvinit=zero              !<==== ACTHUNG! 
             do m=1,jdim
                i=HJ%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(jsite)==0)then
                   call cdg(jsite,i,r,sgn)
                   j=binary_search(HK%map,r)
                   vvinit(j) = sgn*vvinit_tmp(m)
                endif
             enddo
          end if
          deallocate(HJ%map,HK%map,vvinit_tmp)
          !
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)
          !
          call setup_Hv_sector(ksector)
          if(ed_sparse_H)call ed_buildH_c()
          !
          nlanc=min(kdim,lanc_nGFiter)
          allocate(alfa_(nlanc),beta_(nlanc))
#ifdef _MPI
          if(MpiStatus)then
             call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
          else
             call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
          endif
#else
          call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
#endif
          cnorm2=one*norm2
          call add_to_lanczos_densChi_mix(cnorm2,state_e,alfa_,beta_,1,iorb,jorb)
          !
          call delete_Hv_sector()
          !
          deallocate(vvinit,alfa_,beta_)
       enddo
       !
       !
       !+- Apply Sum_ispin c^dg_{iorb,ispin} c_{jorb,ispin} -+!
       do ispin=1,Nspin
          jsite=impIndex(jorb,ispin)
          jsector = getCsector(ispin,isector)
          if(jsector/=0)then
             jdim  = getdim(jsector)
             allocate(vvinit_tmp(jdim))
             call build_sector(jsector,HJ)
             vvinit_tmp=zero
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(jsite)==1)then
                   call c(jsite,i,r,sgn)
                   j=binary_search(HJ%map,r)
                   vvinit_tmp(j) = sgn*state_cvec(m)
                endif
             enddo
          endif
          isite = impIndex(iorb,ispin)
          ksector = getCDGsector(ispin,jsector)
          if(ksector/=0) then       
             kdim  = getdim(ksector)
             allocate(vvinit(kdim))
             call build_sector(ksector,HK)
             vvinit=zero
             do m=1,jdim
                i=HJ%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(isite)==0)then
                   call cdg(isite,i,r,sgn)
                   j=binary_search(HK%map,r)
                   vvinit(j) = sgn*vvinit_tmp(m)
                endif
             enddo
          end if
          deallocate(HJ%map,HK%map,vvinit_tmp)
          !
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)
          !
          call setup_Hv_sector(ksector)
          if(ed_sparse_H)call ed_buildH_c()
          !
          nlanc=min(kdim,lanc_nGFiter)
          allocate(alfa_(nlanc),beta_(nlanc))
#ifdef _MPI
          if(MpiStatus)then
             call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
          else
             call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
          endif
#else
          call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
#endif
          cnorm2=one*norm2
          call add_to_lanczos_densChi_mix(cnorm2,state_e,alfa_,beta_,-1,iorb,jorb)
          !
          call delete_Hv_sector()
          !
          deallocate(vvinit,alfa_,beta_)
       enddo
       !
       nullify(state_cvec)
       deallocate(HI%map)
       !
    enddo
    
  end subroutine lanc_ed_build_densChi_mix_c

  subroutine add_to_lanczos_densChi(vnorm2,Ei,alanc,blanc,isign,iorb,jorb)
    integer                                    :: iorb,jorb,isign
    complex(8)                                 :: pesoF,pesoAB,pesoBZ,peso,vnorm2  
    real(8)                                    :: Ei,Ej,Egs,de
    integer                                    :: nlanc
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw,chisp
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    pesoF  = vnorm2/zeta_function 
    pesoBZ = 1d0
    if(finiteT)pesoBZ = exp(-beta*(Ei-Egs))
    !
    diag             = 0.d0
    subdiag          = 0.d0
    Z                = eye(Nlanc)
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call tql2(Nlanc,diag,subdiag,Z,ierr)
    !
    select case(isign)
    case (1)
       do j=1,nlanc
          Ej     = diag(j)
          dE     = Ej-Ei
          pesoAB = Z(1,j)*Z(1,j)
          peso   = pesoF*pesoAB*pesoBZ
          if(beta*dE < 1d-1)then     !abs(X - (1-exp(-X)) is about 5*10^-3 for X<10^-1 this is a satisfactory bound
             densChi_iv(iorb,jorb,0)=densChi_iv(iorb,jorb,0) - peso*beta
          else
             densChi_iv(iorb,jorb,0)=densChi_iv(iorb,jorb,0) + peso*(exp(-beta*dE)-1d0)/dE 
          endif
          do i=1,Lmats
             densChi_iv(iorb,jorb,i)=densChi_iv(iorb,jorb,i) + peso*(exp(-beta*dE)-1d0)/(dcmplx(0d0,vm(i)) - dE)
          enddo
          do i=0,Ltau
             densChi_tau(iorb,jorb,i)=densChi_tau(iorb,jorb,i) + peso*exp(-tau(i)*de)
          enddo
          do i=1,Lreal
             densChi_w(iorb,jorb,i)=densChi_w(iorb,jorb,i) + peso*(exp(-beta*dE)-1.d0)/(dcmplx(wr(i),eps) - dE)
          enddo
       enddo
    case (-1)
       do j=1,nlanc
          Ej     = diag(j)
          dE     = Ej-Ei
          pesoAB = Z(1,j)*Z(1,j)
          peso   = pesoF*pesoAB*pesoBZ
          if(beta*dE < 1d-1)then     !abs(X - (1-exp(-X)) is about 5*10^-3 for X<10^-1 this is a satisfactory bound
             densChi_iv(iorb,jorb,0)=densChi_iv(iorb,jorb,0) + peso*beta
          else
             densChi_iv(iorb,jorb,0)=densChi_iv(iorb,jorb,0) + peso*(1d0-exp(-beta*dE))/dE 
          endif
          do i=1,Lmats
             densChi_iv(iorb,jorb,i)=densChi_iv(iorb,jorb,i) + peso*(1d0-exp(-beta*dE))/(dcmplx(0d0,vm(i)) + dE)
          enddo
          do i=0,Ltau
             densChi_tau(iorb,jorb,i)=densChi_tau(iorb,jorb,i) + peso*exp(-(beta-tau(i))*dE)
          enddo
          do i=1,Lreal
             densChi_w(iorb,jorb,i)=densChi_w(iorb,jorb,i) + peso*(1d0-exp(-beta*dE))/(dcmplx(wr(i),eps) + dE)
          enddo
       enddo
    case default
       stop "add_to_lanczos_densChi: isign not in {-1,1}"
    end select
  end subroutine add_to_lanczos_densChi

  subroutine add_to_lanczos_densChi_mix(vnorm2,Ei,alanc,blanc,isign,iorb,jorb)
    integer                                    :: iorb,jorb,isign
    complex(8)                                 :: pesoF,pesoAB,pesoBZ,peso,vnorm2  
    real(8)                                    :: Ei,Ej,Egs,de
    integer                                    :: nlanc
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw,chisp
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    pesoF  = vnorm2/zeta_function 
    pesoBZ = 1d0
    if(finiteT)pesoBZ = exp(-beta*(Ei-Egs))
    !
    diag             = 0.d0
    subdiag          = 0.d0
    Z                = eye(Nlanc)
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call tql2(Nlanc,diag,subdiag,Z,ierr)
    !
    select case(isign)
    case (1)
       do j=1,nlanc
          Ej     = diag(j)
          dE     = Ej-Ei
          pesoAB = Z(1,j)*Z(1,j)
          peso   = pesoF*pesoAB*pesoBZ
          if(beta*dE < 1d-1)then     !abs(X - (1-exp(-X)) is about 5*10^-3 for X<10^-1 this is a satisfactory bound
             densChi_mix_iv(iorb,jorb,0)=densChi_mix_iv(iorb,jorb,0) - peso*beta
          else
             densChi_mix_iv(iorb,jorb,0)=densChi_mix_iv(iorb,jorb,0) + peso*(exp(-beta*dE)-1d0)/dE 
          endif
          do i=1,Lmats
             densChi_mix_iv(iorb,jorb,i)=densChi_mix_iv(iorb,jorb,i) + peso*(exp(-beta*dE)-1d0)/(dcmplx(0d0,vm(i)) - dE)
          enddo
          do i=0,Ltau
             densChi_mix_tau(iorb,jorb,i)=densChi_mix_tau(iorb,jorb,i) + peso*exp(-tau(i)*de)
          enddo
          do i=1,Lreal
             densChi_mix_w(iorb,jorb,i)=densChi_mix_w(iorb,jorb,i) + peso*(exp(-beta*dE)-1.d0)/(dcmplx(wr(i),eps) - dE)
          enddo
       enddo
    case (-1)
       do j=1,nlanc
          Ej     = diag(j)
          dE     = Ej-Ei
          pesoAB = Z(1,j)*Z(1,j)
          peso   = pesoF*pesoAB*pesoBZ
          if(beta*dE < 1d-1)then     !abs(X - (1-exp(-X)) is about 5*10^-3 for X<10^-1 this is a satisfactory bound
             densChi_mix_iv(iorb,jorb,0)=densChi_mix_iv(iorb,jorb,0) + peso*beta
          else
             densChi_mix_iv(iorb,jorb,0)=densChi_mix_iv(iorb,jorb,0) + peso*(1d0-exp(-beta*dE))/dE 
          endif
          do i=1,Lmats
             densChi_mix_iv(iorb,jorb,i)=densChi_mix_iv(iorb,jorb,i) + peso*(1d0-exp(-beta*dE))/(dcmplx(0d0,vm(i)) + dE)
          enddo
          do i=0,Ltau
             densChi_mix_tau(iorb,jorb,i)=densChi_mix_tau(iorb,jorb,i) + peso*exp(-(beta-tau(i))*dE)
          enddo
          do i=1,Lreal
             densChi_mix_w(iorb,jorb,i)=densChi_mix_w(iorb,jorb,i) + peso*(1d0-exp(-beta*dE))/(dcmplx(wr(i),eps) + dE)
          enddo
       enddo
    case default
       stop "add_to_lanczos_densChi_mix: isign not in {-1,1}"
    end select
  end subroutine add_to_lanczos_densChi_mix

  subroutine add_to_lanczos_densChi_tot(vnorm2,Ei,alanc,blanc,isign)
    complex(8)                                 :: pesoF,pesoAB,pesoBZ,peso,vnorm2  
    real(8)                                    :: Ei,Ej,Egs,de
    integer                                    :: nlanc,isign
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw,chisp
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    pesoF  = vnorm2/zeta_function 
    pesoBZ = 1d0
    if(finiteT)pesoBZ = exp(-beta*(Ei-Egs))
    !
    diag             = 0.d0
    subdiag          = 0.d0
    Z                = eye(Nlanc)
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call tql2(Nlanc,diag,subdiag,Z,ierr)
    !
    select case(isign)
    case (1)
       do j=1,nlanc
          Ej     = diag(j)
          dE     = Ej-Ei
          pesoAB = Z(1,j)*Z(1,j)
          peso   = pesoF*pesoAB*pesoBZ
          if(beta*dE < 1d-1)then     !abs(X - (1-exp(-X)) is about 5*10^-3 for X<10^-1 this is a satisfactory bound
             densChi_tot_iv(0)=densChi_tot_iv(0) - peso*beta
          else
             densChi_tot_iv(0)=densChi_tot_iv(0) + peso*(exp(-beta*dE)-1d0)/dE 
          endif
          do i=1,Lmats
             densChi_tot_iv(i)=densChi_tot_iv(i) + peso*(exp(-beta*dE)-1d0)/(dcmplx(0d0,vm(i)) - dE)
          enddo
          do i=0,Ltau
             densChi_tot_tau(i)=densChi_tot_tau(i) + peso*exp(-tau(i)*de)
          enddo
          do i=1,Lreal
             densChi_tot_w(i)=densChi_tot_w(i) + peso*(exp(-beta*dE)-1.d0)/(dcmplx(wr(i),eps) - dE)
          enddo
       enddo
    case (-1)
       do j=1,nlanc
          Ej     = diag(j)
          dE     = Ej-Ei
          pesoAB = Z(1,j)*Z(1,j)
          peso   = pesoF*pesoAB*pesoBZ
          if(beta*dE < 1d-1)then     !abs(X - (1-exp(-X)) is about 5*10^-3 for X<10^-1 this is a satisfactory bound
             densChi_tot_iv(0)=densChi_tot_iv(0) + peso*beta
          else
             densChi_tot_iv(0)=densChi_tot_iv(0) + peso*(1d0-exp(-beta*dE))/dE 
          endif
          do i=1,Lmats
             densChi_tot_iv(i)=densChi_tot_iv(i) + peso*(1d0-exp(-beta*dE))/(dcmplx(0d0,vm(i)) + dE)
          enddo
          do i=0,Ltau
             densChi_tot_tau(i)=densChi_tot_tau(i) + peso*exp(-(beta-tau(i))*dE)
          enddo
          do i=1,Lreal
             densChi_tot_w(i)=densChi_tot_w(i) + peso*(1d0-exp(-beta*dE))/(dcmplx(wr(i),eps) + dE)
          enddo
       enddo
    case default
       stop "add_to_lanczos_densChi_tot: isign not in {-1,1}"
    end select
  end subroutine add_to_lanczos_densChi_tot



















  !+------------------------------------------------------------------+
  !                            PAIR
  !+------------------------------------------------------------------+
  !PURPOSE  : Evaluate Pair Susceptibility using Lanczos algorithm
  !+------------------------------------------------------------------+
  subroutine build_chi_pair()
    integer :: iorb
    write(LOGfile,"(A)")"Get impurity pair Chi:"
    do iorb=1,Norb
       write(LOGfile,"(A)")"Get Chi_pair_l"//reg(txtfy(iorb))
       if(MPI_MASTER)call start_timer()
       call lanc_ed_build_pairChi_c(iorb)
       if(MPI_MASTER)call stop_timer(LOGfile)
    enddo
    pairChi_tau = PairChi_tau/zeta_function
    pairChi_w   = pairChi_w/zeta_function
    pairChi_iv  = pairChi_iv/zeta_function
  end subroutine build_chi_pair

  !+------------------------------------------------------------------+
  !PURPOSE  : Evaluate the Pair susceptibility \Chi_pair for a 
  ! single orbital: \chi = <Phi_a(\tau)Phi_a(0)>
  !+------------------------------------------------------------------+
  subroutine lanc_ed_build_pairChi_c(iorb)
    integer                          :: iorb,isite,isector,izero
    integer                          :: numstates
    integer                          :: nlanc,idim
    integer                          :: iup0,idw0,isign
    integer                          :: ib(Nlevels)
    integer                          :: m,i,i1,i2,j,r
    real(8)                          :: norm0,sgn,sgn1,sgn2
    real(8),allocatable              :: alfa_(:),beta_(:)
    complex(8),allocatable           :: vvinit(:)
    integer                          :: Nitermax
    type(sector_map) :: HI    !map of the Sector S to Hilbert space H
    !
    
    !
    do izero=1,state_list%size
       isector    =  es_return_sector(state_list,izero)
       idim       =  getdim(isector)
       state_e    =  es_return_energy(state_list,izero)
       state_cvec => es_return_cvector(state_list,izero)
       norm0=sqrt(dot_product(state_cvec,state_cvec))
       if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
       call build_sector(isector,HI)
       !
       !Build the C_{iorb,up}C_{iorb,dw}|eigvec> 
       if(ed_verbose==3)write(LOGfile,"(A,2I3)")'Apply C_{iorb,up}C_{iorb,dw}:',getsz(isector)
       allocate(vvinit(idim))
       vvinit=0.d0
       do m=1,idim
          i=HI%map(m)
          ib = bdecomp(i,2*Ns)
          if(ib(iorb+Ns)==0.OR.ib(iorb)==0)cycle
          call c(iorb+Ns,i,i1,sgn1)
          call c(iorb,i1,i2,sgn2)
          j = binary_search(HI%map,i2)
          vvinit(j) = sgn1*sgn2*state_cvec(m)
       enddo
       norm0=sqrt(dot_product(vvinit,vvinit))
       vvinit=vvinit/norm0
       !
       call setup_Hv_sector(isector)
       if(ed_sparse_H)call ed_buildH_c()
       !
       nlanc=min(idim,lanc_nGFiter)
       allocate(alfa_(nlanc),beta_(nlanc))
#ifdef _MPI
       if(MpiStatus)then
          call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
       else
          call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
       endif
#else
       call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
#endif
       isign=-1 !<== ACTHUNG!!!! check this is the correct value of isign
       call add_to_lanczos_pairChi(norm0,state_e,alfa_,beta_,isign,iorb)
       !
       call delete_Hv_sector()
       !
       if(spH0%status)call sp_delete_matrix(spH0)
       deallocate(vvinit,alfa_,beta_)
       !
       !Build the CDG_{iorb,dw}CDG_{iorb,up}|eigvec> 
       if(ed_verbose==3)write(LOGfile,"(A,2I3)")'Apply CDG_{iorb,dw}CDG_{iorb,up}:',getsz(isector)
       allocate(vvinit(idim))
       vvinit=0.d0
       do m=1,idim
          i=HI%map(m)
          ib = bdecomp(i,2*Ns)
          if(ib(iorb+Ns)==1.OR.ib(iorb)==1)cycle
          call cdg(iorb,i,i1,sgn1)
          call cdg(iorb+Ns,i1,i2,sgn2)
          j = binary_search(HI%map,i2)
          vvinit(j) = sgn1*sgn2*state_cvec(m)
       enddo
       norm0=sqrt(dot_product(vvinit,vvinit))
       vvinit=vvinit/norm0
       !
       call setup_Hv_sector(isector)
       if(ed_sparse_H)call ed_buildH_c()
       !
       nlanc=min(idim,lanc_nGFiter)
       allocate(alfa_(nlanc),beta_(nlanc))
#ifdef _MPI
       if(MpiStatus)then
          call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
       else
          call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
       endif
#else
       call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
#endif
       isign=1 !<== ACTHUNG!!!! check this is the correct value of isign
       call add_to_lanczos_pairChi(norm0,state_e,alfa_,beta_,isign,iorb)
       !
       call delete_Hv_sector()
       !
       if(spH0%status)call sp_delete_matrix(spH0)
       deallocate(vvinit,alfa_,beta_)
       deallocate(HI%map)
       nullify(state_cvec)
    enddo
    
  end subroutine lanc_ed_build_pairChi_c

  subroutine add_to_lanczos_pairChi(vnorm,Ei,alanc,blanc,isign,iorb)
    real(8)                                    :: vnorm,Ei,Ej,Egs,pesoF,pesoAB,pesoBZ,de,peso
    integer                                    :: nlanc
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: isign,iorb
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw,chisp
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    pesoF  = vnorm**2/zeta_function 
    pesoBZ = 1d0
    if(finiteT)pesoBZ = exp(-beta*(Ei-Egs))
    !
    diag             = 0.d0
    subdiag          = 0.d0
    Z                = eye(Nlanc)
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call tql2(Nlanc,diag,subdiag,Z,ierr)
    !
    select case(isign)
    case (1)
       do j=1,nlanc
          Ej     = diag(j)
          dE     = Ej-Ei
          pesoAB = Z(1,j)*Z(1,j)
          peso   = pesoF*pesoAB*pesoBZ
          if(beta*dE < 1d-1)then     !abs(X - (1-exp(-X)) is about 5*10^-3 for X<10^-1 this is a satisfactory bound
             pairChi_iv(iorb,0)=pairChi_iv(iorb,0) + peso*beta
          else
             pairChi_iv(iorb,0)=pairChi_iv(iorb,0) + peso*(1d0-exp(-beta*dE))/dE 
          endif
          do i=1,Lmats
             pairChi_iv(iorb,i)=pairChi_iv(iorb,i) + peso*(exp(-beta*dE)-1d0)/(dcmplx(0d0,vm(i)) - dE)
          enddo
          do i=0,Ltau
             pairChi_tau(iorb,i)=pairChi_tau(iorb,i) + peso*exp(-tau(i)*de)
          enddo
          do i=1,Lreal
             pairChi_w(iorb,i)=pairChi_w(iorb,i) + peso*(exp(-beta*dE)-1.d0)/(dcmplx(wr(i),eps) - dE)
          enddo
       enddo
    case (-1)
       do j=1,nlanc
          Ej     = diag(j)
          dE     = Ej-Ei
          pesoAB = Z(1,j)*Z(1,j)
          peso   = pesoF*pesoAB*pesoBZ
          if(beta*dE < 1d-1)then     !abs(X - (1-exp(-X)) is about 5*10^-3 for X<10^-1 this is a satisfactory bound
             pairChi_iv(iorb,0)=pairChi_iv(iorb,0) + peso*beta
          else
             pairChi_iv(iorb,0)=pairChi_iv(iorb,0) + peso*(1d0-exp(-beta*dE))/dE 
          endif
          do i=1,Lmats
             pairChi_iv(iorb,i)=pairChi_iv(iorb,i) + peso*(1d0-exp(-beta*dE))/(dcmplx(0d0,vm(i)) + dE)
          enddo
          do i=0,Ltau
             pairChi_tau(iorb,i)=pairChi_tau(iorb,i) + peso*exp(-(beta-tau(i))*dE)
          enddo
          do i=1,Lreal
             pairChi_w(iorb,i)=pairChi_w(iorb,i) + peso*(1d0-exp(-beta*dE))/(dcmplx(wr(i),eps) + dE)
          enddo
       enddo
    case default
       stop "add_to_lanczos_pairChi: isign not in {-1,1}"
    end select
  end subroutine add_to_lanczos_pairChi
























  !+------------------------------------------------------------------+
  !PURPOSE  : Allocate arrays and setup frequencies and times
  !+------------------------------------------------------------------+
  subroutine allocate_grids
    integer :: i
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(vm))allocate(vm(0:Lmats))          !bosonic frequencies
    if(.not.allocated(wr))allocate(wr(Lreal))
    if(.not.allocated(tau))allocate(tau(0:Ltau))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    do i=0,Lmats
       vm(i) = pi/beta*2.d0*dble(i)
    enddo
    wr     = linspace(wini,wfin,Lreal)
    tau(0:)= linspace(0.d0,beta,Ltau+1)
  end subroutine allocate_grids


  subroutine deallocate_grids
    if(allocated(wm))deallocate(wm)
    if(allocated(vm))deallocate(vm)
    if(allocated(tau))deallocate(tau)
    if(allocated(wr))deallocate(wr)
  end subroutine deallocate_grids








  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !++++++++++++++++++COMPUTATIONAL ROUTINE: TQL2++++++++++++++++++++++++ 
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !---------------------------------------------------------------------
  ! PURPOSE computes all eigenvalues/vectors, real symmetric tridiagonal matrix.
  !    This subroutine finds the eigenvalues and eigenvectors of a symmetric
  !    tridiagonal matrix by the QL method.  The eigenvectors of a full
  !    symmetric matrix can also be found if TRED2 has been used to reduce this
  !    full matrix to tridiagonal form.
  !  Parameters:
  !    Input, integer ( kind = 4 ) N, the order of the matrix.
  !
  !    Input/output, real ( kind = 8 ) D(N).  On input, the diagonal elements of
  !    the matrix.  On output, the eigenvalues in ascending order.  If an error
  !    exit is made, the eigenvalues are correct but unordered for indices
  !    1,2,...,IERR-1.
  !
  !    Input/output, real ( kind = 8 ) E(N).  On input, E(2:N) contains the
  !    subdiagonal elements of the input matrix, and E(1) is arbitrary.
  !    On output, E has been destroyed.
  !
  !    Input, real ( kind = 8 ) Z(N,N).  On input, the transformation matrix
  !    produced in the reduction by TRED2, if performed.  If the eigenvectors of
  !    the tridiagonal matrix are desired, Z must contain the identity matrix.
  !    On output, Z contains the orthonormal eigenvectors of the symmetric
  !    tridiagonal (or full) matrix.  If an error exit is made, Z contains
  !    the eigenvectors associated with the stored eigenvalues.
  !
  !    Output, integer ( kind = 4 ) IERR, error flag.
  !    0, normal return,
  !    J, if the J-th eigenvalue has not been determined after
  !    30 iterations.
  !
  !---------------------------------------------------------------------
  subroutine tql2 ( n, d, e, z, ierr )
    integer :: n
    real(8) :: c
    real(8) :: c2
    real(8) :: c3
    real(8) :: d(n)
    real(8) :: dl1
    real(8) :: e(n)
    real(8) :: el1
    real(8) :: f
    real(8) :: g
    real(8) :: h
    integer ( kind = 4 ) i
    integer ( kind = 4 ) ierr
    integer ( kind = 4 ) ii
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    integer ( kind = 4 ) l
    integer ( kind = 4 ) l1
    integer ( kind = 4 ) l2
    integer ( kind = 4 ) m
    integer ( kind = 4 ) mml
    real(8) :: p
    real(8) :: r
    real(8) :: s
    real(8) :: s2
    real(8) :: tst1
    real(8) :: tst2
    real(8) :: z(n,n)
    ierr = 0
    if ( n == 1 ) then
       return
    end if
    do i = 2, n
       e(i-1) = e(i)
    end do
    f = 0.0D+00
    tst1 = 0.0D+00
    e(n) = 0.0D+00
    do l = 1, n
       j = 0
       h = abs ( d(l) ) + abs ( e(l) )
       tst1 = max ( tst1, h )
       !
       !  Look for a small sub-diagonal element.
       !
       do m = l, n
          tst2 = tst1 + abs ( e(m) )
          if ( tst2 == tst1 ) then
             exit
          end if
       end do
       if ( m == l ) then
          go to 220
       end if
130    continue
       if ( 30 <= j ) then
          ierr = l
          return
       end if
       j = j + 1
       !
       !  Form shift.
       !
       l1 = l + 1
       l2 = l1 + 1
       g = d(l)
       p = ( d(l1) - g ) / ( 2.0D+00 * e(l) )
       r = pythag ( p, 1.0D+00 )
       d(l) = e(l) / ( p + sign ( r, p ) )
       d(l1) = e(l) * ( p + sign ( r, p ) )
       dl1 = d(l1)
       h = g - d(l)
       d(l2:n) = d(l2:n) - h
       f = f + h
       !
       !  QL transformation.
       !
       p = d(m)
       c = 1.0D+00
       c2 = c
       el1 = e(l1)
       s = 0.0D+00
       mml = m - l
       do ii = 1, mml
          c3 = c2
          c2 = c
          s2 = s
          i = m - ii
          g = c * e(i)
          h = c * p
          r = pythag ( p, e(i) )
          e(i+1) = s * r
          s = e(i) / r
          c = p / r
          p = c * d(i) - s * g
          d(i+1) = h + s * ( c * g + s * d(i) )
          !
          !  Form vector.
          !
          do k = 1, n
             h = z(k,i+1)
             z(k,i+1) = s * z(k,i) + c * h
             z(k,i) = c * z(k,i) - s * h
          end do
       end do
       p = - s * s2 * c3 * el1 * e(l) / dl1
       e(l) = s * p
       d(l) = c * p
       tst2 = tst1 + abs ( e(l) )
       if ( tst2 > tst1 ) then
          go to 130
       end if
220    continue
       d(l) = d(l) + f
    end do
    !
    !  Order eigenvalues and eigenvectors.
    !
    do ii = 2, n
       i = ii - 1
       k = i
       p = d(i)
       do j = ii, n
          if ( d(j) < p ) then
             k = j
             p = d(j)
          end if
       end do
       if ( k /= i ) then
          d(k) = d(i)
          d(i) = p
          do j = 1, n
             call r8_swap ( z(j,i), z(j,k) )
          end do
       end if
    end do
    return
  end subroutine tql2


  !---------------------------------------------------------------------
  ! PURPOSE: computes SQRT ( A * A + B * B ) carefully.
  !    The formula
  !    PYTHAG = sqrt ( A * A + B * B )
  !    is reasonably accurate, but can fail if, for example, A**2 is larger
  !    than the machine overflow.  The formula can lose most of its accuracy
  !    if the sum of the squares is very large or very small.
  !  Parameters:
  !    Input, real(8) :: A, B, the two legs of a right triangle.
  !    Output, real(8) :: PYTHAG, the length of the hypotenuse.
  !---------------------------------------------------------------------
  function pythag ( a, b )
    implicit none
    real(8) :: a
    real(8) :: b
    real(8) :: p
    real(8) :: pythag
    real(8) :: r
    real(8) :: s
    real(8) :: t
    real(8) :: u
    p = max ( abs ( a ), abs ( b ) )
    if ( p /= 0.0D+00 ) then
       r = ( min ( abs ( a ), abs ( b ) ) / p )**2
       do
          t = 4.0D+00 + r
          if ( t == 4.0D+00 ) then
             exit
          end if
          s = r / t
          u = 1.0D+00 + 2.0D+00 * s
          p = u * p
          r = ( s / u )**2 * r
       end do
    end if
    pythag = p
    return
  end function pythag

  !---------------------------------------------------------------------
  ! PURPOSE: swaps two R8's.
  !  Parameters:
  !    Input/output, real(8) :: X, Y.  On output, the values of X and
  !    Y have been interchanged.
  !---------------------------------------------------------------------
  subroutine r8_swap ( x, y )
    real(8) :: x
    real(8) :: y
    real(8) :: z
    z = x
    x = y
    y = z
    return
  end subroutine r8_swap



end MODULE ED_CHI_FUNCTIONS
