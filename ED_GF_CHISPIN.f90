MODULE ED_GF_CHISPIN
  USE ED_GF_SHARED
  implicit none
  private


  public :: build_chi_spin

contains


  !+------------------------------------------------------------------+
  !                            SPIN
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






  !################################################################
  !################################################################
  !################################################################
  !################################################################






  subroutine lanc_ed_build_spinChi_c(iorb)
    integer                :: iorb,isite,isector,istate
    integer                :: numstates
    integer                :: nlanc
    integer                :: isign
    integer                :: idim,idimUP,idimDW
    integer                :: jdim,jdimUP,jdimDW
    integer                :: nup(Ns),ndw(Ns)
    integer                :: m,i,j,r
    integer                :: iup,idw,jup,jdw,mup,mdw
    real(8)                :: norm2,sgn
    real(8),allocatable    :: alfa_(:),beta_(:)
    complex(8),allocatable :: vvinit(:)
    integer                :: Nitermax
    type(sector_map)       :: HI(2)    !map of the Sector S to Hilbert space H
    !
    !
    !
    do istate=1,state_list%size
       isector     =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
       ! state_cvec => es_return_cvector(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          state_cvec => es_return_cvector(MpiComm,state_list,istate)
       else
          state_cvec => es_return_cvector(state_list,istate)
       endif
#else
       state_cvec => es_return_cvector(state_list,istate)
#endif
       !
       !
       idim    = getdim(isector)
       idimUp  = getDimUp(isector)
       idimDw  = getDimDw(isector)
       call build_sector(isector,HI)
       !
       if(ed_verbose==3)write(LOGfile,"(A,2I3)")'Apply Sz:',getnup(isector),getndw(isector)
       !
       allocate(vvinit(idim));vvinit=0.d0
       do iup=1,idimUP
          mup = HI(1)%map(iup)
          nup = bdecomp(mup,Ns)
          !
          do idw=1,idimDW
             mdw = HI(2)%map(idw)
             ndw = bdecomp(mdw,Ns)
             !
             i = iup + (idw-1)*idimUP
             !
             sgn = nup(iorb)-ndw(iorb)
             !
             vvinit(i) = 0.5d0*sgn*state_cvec(i)   !build the cdg_up|gs> state
          enddo
       enddo
       call delete_sector(isector,HI)
       !
       norm2=dot_product(vvinit,vvinit)
       vvinit=vvinit/sqrt(norm2)
       !
       nlanc=min(jdim,lanc_nGFiter)
       allocate(alfa_(nlanc),beta_(nlanc))
       !
       call build_Hv_sector(isector)
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
       call add_to_lanczos_spinChi(norm2,state_e,alfa_,beta_,1,iorb)
       !holes
       call add_to_lanczos_spinChi(norm2,state_e,alfa_,beta_,-1,iorb)
       !
       call delete_Hv_sector()
       !
       deallocate(vvinit,alfa_,beta_)
       nullify(state_cvec)
    enddo
    return
  end subroutine lanc_ed_build_spinChi_c



  !################################################################
  !################################################################
  !################################################################
  !################################################################




  subroutine lanc_ed_build_spinChi_tot_c()
    integer                :: iorb,isite,isector,istate
    integer                :: numstates
    integer                :: nlanc
    integer                :: isign
    integer                :: idim,idimUP,idimDW
    integer                :: jdim,jdimUP,jdimDW
    integer                :: nup(Ns),ndw(Ns)
    integer                :: m,i,j,r
    integer                :: iup,idw,jup,jdw,mup,mdw
    real(8)                :: norm2,sgn,Sup,Sdw
    real(8),allocatable    :: alfa_(:),beta_(:)
    complex(8),allocatable :: vvinit(:)
    integer                :: Nitermax
    type(sector_map)       :: HI(2)    !map of the Sector S to Hilbert space H
    !
    !
    do istate=1,state_list%size
       isector     =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
       ! state_cvec => es_return_cvector(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          state_cvec => es_return_cvector(MpiComm,state_list,istate)
       else
          state_cvec => es_return_cvector(state_list,istate)
       endif
#else
       state_cvec => es_return_cvector(state_list,istate)
#endif
       !
       !
       idim  = getdim(isector)
       idimUp  = getDimUp(isector)
       idimDw  = getDimDw(isector)
       call build_sector(isector,HI)
       !
       if(ed_verbose==3)write(LOGfile,"(A,2I3)")'Apply Sz:',getnup(isector),getndw(isector)
       !
       allocate(vvinit(idim));vvinit=0.d0
       do iup=1,idimUP
          mup = HI(1)%map(iup)
          nup = bdecomp(mup,Ns)
          !
          do idw=1,idimDW
             mdw = HI(2)%map(idw)
             ndw = bdecomp(mdw,Ns)
             !
             i = iup + (idw-1)*idimUP
             !
             Sup = sum(nup(1:Norb))
             Sdw = sum(ndw(1:Norb))
             sgn = Sup - Sdw
             vvinit(m) = 0.5d0*sgn*state_cvec(m) 
          enddo
       enddo
       call delete_sector(isector,HI)
       !
       norm2=dot_product(vvinit,vvinit)
       vvinit=vvinit/sqrt(norm2)
       !
       nlanc=min(jdim,lanc_nGFiter)
       allocate(alfa_(nlanc),beta_(nlanc))
       !
       call build_Hv_sector(isector)
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
       call add_to_lanczos_spinChi(norm2,state_e,alfa_,beta_,1,Norb+1)
       !holes
       call add_to_lanczos_spinChi(norm2,state_e,alfa_,beta_,-1,Norb+1)
       !
       call delete_Hv_sector()
       !
       deallocate(vvinit,alfa_,beta_)
       nullify(state_cvec)
    enddo
    return
  end subroutine lanc_ed_build_spinChi_tot_c





  !################################################################
  !################################################################
  !################################################################
  !################################################################




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





END MODULE ED_GF_CHISPIN
