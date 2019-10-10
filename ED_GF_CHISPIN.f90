MODULE ED_GF_CHISPIN
  USE ED_GF_SHARED
  implicit none
  private


  public :: build_chi_spin

  integer             :: istate,iorb,jorb,ispin,jspin
  integer             :: isector,jsector
  integer             :: idim,idimUP,idimDW
  integer             :: jdim,jdimUP,jdimDW
  real(8),allocatable :: vvinit(:),vvloc(:)
  real(8),allocatable :: alfa_(:),beta_(:)
  integer             :: ialfa,ibeta
  integer             :: iorb1,jorb1
  integer             :: r
  integer             :: i,iup,idw
  integer             :: j,jup,jdw  
  integer             :: m,mup,mdw
  real(8)             :: sgn,norm2,norm0
  integer             :: Nitermax,Nlanc,vecDim


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
       if(MPIMASTER)call start_timer()
       select case(ed_diag_type)
       case default
          call lanc_ed_build_spinChi_main(iorb)
       case ("full")
          call full_ed_build_spinChi_main(iorb)
       end select
       if(MPIMASTER)call stop_timer(LOGfile)
    enddo
    if(Norb>1)then
       write(LOGfile,"(A)")"Get Chi_spin_tot"
       if(MPIMASTER)call start_timer()
       select case(ed_diag_type)
       case default
          call lanc_ed_build_spinChi_tot_main()
       case ("full")
       end select
       if(MPIMASTER)call stop_timer(LOGfile)
    endif
    spinChi_tau = SpinChi_tau/zeta_function
    spinChi_w   = spinChi_w/zeta_function
    spinChi_iv  = spinChi_iv/zeta_function
  end subroutine build_chi_spin






  !################################################################
  !################################################################
  !################################################################
  !################################################################






  subroutine lanc_ed_build_spinChi_main(iorb)
    integer                     :: iorb
    integer,dimension(2*Ns_Ud)  :: Indices
    integer,dimension(2*Ns_Ud)  :: Jndices
    integer,dimension(Ns_Ud)    :: iDimUps,iDimDws
    integer,dimension(Ns_Ud)    :: jDimUps,jDimDws
    integer,dimension(2,Ns_Orb) :: Nud
    integer                     :: Iud(2)
    type(sector_map)            :: HI(2*Ns_Ud),HJ(2*Ns_Ud)
    !
    if(ed_total_ud)then
       ialfa = 1
       iorb1 = iorb
    else
       ialfa = iorb
       iorb1 = 1
    endif
    ibeta  = ialfa + (ispin-1)*Ns_Ud
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
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
       idim = getdim(isector)
       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       !
       if(MpiMaster)then
          !
          if(ed_verbose==3)write(LOGfile,"(A,I12)")'Apply Sz:',isector
          !
          allocate(vvinit(idim));vvinit=0.d0
          !
          call build_sector(isector,HI)
          do i=1,iDim
             call state2indices(i,[iDimUps,iDimDws],Indices)
             iud(1)   = HI(ialfa)%map(Indices(ialfa))
             iud(2)   = HI(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
             nud(1,:) = Bdecomp(iud(1),Ns_Orb)
             nud(2,:) = Bdecomp(iud(2),Ns_Orb)
             !
             sgn = dble(nud(1,iorb1))-dble(nud(2,iorb1))
             !
             vvinit(i) = 0.5d0*sgn*state_cvec(i)
          enddo
          call delete_sector(isector,HI)
          !
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)
       endif
       !
       nlanc=min(idim,lanc_nGFiter)
       allocate(alfa_(nlanc),beta_(nlanc))
       !
       call build_Hv_sector(isector)
#ifdef _MPI
       if(MpiStatus)then
          call Bcast_MPI(MpiComm,norm2)
          vecDim = vecDim_Hv_sector(isector)
          allocate(vvloc(vecDim))
          call scatter_vector_MPI(MpiComm,vvinit,vvloc)
          call sp_lanc_tridiag(MpiComm,spHtimesV_p,vvloc,alfa_,beta_)
       else
          call sp_lanc_tridiag(spHtimesV_p,vvinit,alfa_,beta_)
       endif
#else
       call sp_lanc_tridiag(spHtimesV_p,vvinit,alfa_,beta_)
#endif
       !particles
       call add_to_lanczos_spinChi(norm2,state_e,alfa_,beta_,1,iorb)
       !holes
       call add_to_lanczos_spinChi(norm2,state_e,alfa_,beta_,-1,iorb)
       !
       call delete_Hv_sector()
       !
       deallocate(alfa_,beta_)
       if(allocated(vvinit))deallocate(vvinit)
       if(allocated(vvloc))deallocate(vvloc)
       nullify(state_cvec)
    enddo
    return
  end subroutine lanc_ed_build_spinChi_main




  !################################################################




  subroutine lanc_ed_build_spinChi_tot_main()
    integer,dimension(2*Ns_Ud)  :: Indices
    integer,dimension(2*Ns_Ud)  :: Jndices
    integer,dimension(Ns_Ud)    :: iDimUps,iDimDws
    integer,dimension(Ns_Ud)    :: jDimUps,jDimDws
    integer,dimension(2,Ns_Orb) :: Nud
    integer                     :: Iud(2)
    real(8)                     :: Sup,Sdw
    type(sector_map)            :: HI(2*Ns_Ud)    !map of the Sector S to Hilbert space H
    !
    if(ed_total_ud)then
       ialfa = 1
       iorb1 = iorb
    else
       ialfa = iorb
       iorb1 = 1
    endif
    ibeta  = ialfa + (ispin-1)*Ns_Ud
    !
    do istate=1,state_list%size
       isector     =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
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
       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       !
       ! if(MpiMaster)then
       call build_sector(isector,HI)
       !
       if(MpiMaster)then
          if(ed_verbose==3)write(LOGfile,"(A,I15)")'Apply Sz:',isector
          !
          allocate(vvinit(idim));vvinit=0.d0
          !
          call build_sector(isector,HI)
          do i=1,iDim
             call state2indices(i,[iDimUps,iDimDws],Indices)
             iud(1)   = HI(ialfa)%map(Indices(ialfa))
             iud(2)   = HI(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
             nud(1,:) = Bdecomp(iud(1),Ns_Orb)
             nud(2,:) = Bdecomp(iud(2),Ns_Orb)
             !
             Sup = sum(nud(1,1:Norb))
             Sdw = sum(nud(2,1:Norb))
             sgn = Sup - Sdw
             !
             vvinit(i) = sgn*state_cvec(i)
          enddo
          call delete_sector(isector,HI)
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)
       endif
       !
       nlanc=min(idim,lanc_nGFiter)
       allocate(alfa_(nlanc),beta_(nlanc))
       !
       call build_Hv_sector(isector)
#ifdef _MP
       if(MpiStatus)then
          call Bcast_MPI(MpiComm,norm2)
          vecDim = vecDim_Hv_sector(isector)
          allocate(vvloc(vecDim))
          call scatter_vector_MPI(MpiComm,vvinit,vvloc)
          call sp_lanc_tridiag(MpiComm,spHtimesV_p,vvloc,alfa_,beta_)
       else
          call sp_lanc_tridiag(spHtimesV_p,vvinit,alfa_,beta_)
       endif
#else
       call sp_lanc_tridiag(spHtimesV_p,vvinit,alfa_,beta_)
#endif
       !particles
       call add_to_lanczos_spinChi(norm2,state_e,alfa_,beta_,1,Norb+1)
       !holes
       call add_to_lanczos_spinChi(norm2,state_e,alfa_,beta_,-1,Norb+1)
       !
       call delete_Hv_sector()
       !
       deallocate(alfa_,beta_)
       if(allocated(vvinit))deallocate(vvinit)
       if(allocated(vvloc))deallocate(vvloc)
       nullify(state_cvec)
    enddo
    return
  end subroutine lanc_ed_build_spinChi_tot_main






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
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,alanc)
       call Bcast_MPI(MpiComm,blanc)
    endif
#endif
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
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




  !################################################################
  !################################################################
  !################################################################
  !################################################################



  subroutine full_ed_build_spinChi_main(iorb)
    integer                     :: iorb
    integer,dimension(2*Ns_Ud)  :: Indices
    integer,dimension(2*Ns_Ud)  :: Jndices
    integer,dimension(Ns_Ud)    :: iDimUps,iDimDws
    integer,dimension(Ns_Ud)    :: jDimUps,jDimDws
    integer,dimension(2,Ns_Orb) :: Nud
    integer                     :: Iud(2)
    type(sector_map)            :: HI(2*Ns_Ud)
    real(8)                     :: chij,spin
    integer                     :: i,j,ll,isector
    integer                     :: idim,ia
    real(8)                     :: Ei,Ej,cc,peso,pesotot
    real(8)                     :: expterm,de,w0,it
    complex(8)                  :: iw 
    !
    !
    !Spin susceptibility \X(tau). |<i|S_z|j>|^2
    !
    if(ed_total_ud)then
       ialfa = 1
       iorb1 = iorb
    else
       ialfa = iorb
       iorb1 = 1
    endif
    !
    do isector=1,Nsectors !loop over <i| total particle number
       call eta(isector,Nsectors,LOGfile)
       idim = getdim(isector)
       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       call build_sector(isector,HI)
       !
       do i=1,idim 
          do j=1,idim
             chij=0.d0
             expterm=exp(-beta*espace(isector)%e(j))
             if(expterm<cutoff)cycle
             do ll=1,idim
                call state2indices(ll,[iDimUps,iDimDws],Indices)
                iud(1)  = HI(ialfa)%map(Indices(ialfa))
                iud(2)  = HI(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
                nud(1,:)= Bdecomp(iud(1),Ns_Orb)
                nud(2,:)= Bdecomp(iud(2),Ns_Orb)
                !
                spin    = nud(1,iorb1) - nud(2,iorb1)
                chij    = chij + espace(isector)%M(ll,i)*spin*espace(isector)%M(ll,j)
             enddo
             Ei=espace(isector)%e(i)
             Ej=espace(isector)%e(j)
             de=Ei-Ej
             peso=chij/zeta_function
             !
             !Matsubara (bosonic) frequency
             !abs(X - (1-exp(-X)) is about 5*10^-3 for X<10^-1 this is a satisfactory bound
             if(beta*dE < 1d-1)then 
                spinChi_iv(iorb,0)=spinChi_iv(iorb,0) + peso*beta
             else
                spinChi_iv(iorb,0)=spinChi_iv(iorb,0) + peso*(1d0-exp(-beta*dE))/dE 
             endif
             do m=1,Lmats
                iw=xi*vm(m)
                spinChi_iv(iorb,m)=spinChi_iv(iorb,m) + peso*(1d0-exp(-beta*dE))/(iw - de)
             enddo
             ! if(de>cutoff)spinChi_iv(iorb,0)=spinChi_iv(iorb,0)-peso*exp(-beta*Ej)*(exp(-beta*de)-1.d0)/de
             ! do m=1,Lmats
             !    iw=xi*vm(m)
             !    spinChi_iv(iorb,m)=spinChi_iv(iorb,m)+peso*exp(-beta*Ej)*(exp(-beta*de)-1.d0)/(iw-de)
             ! enddo
             !
             !Real-frequency: Retarded = Commutator = response function
             do m=1,Lreal
                w0=wr(m);iw=dcmplx(w0,eps)
                spinChi_w(iorb,m)=spinChi_w(iorb,m)+peso*exp(-beta*Ej)*(exp(-beta*de)-1.d0)/(iw-de)
             enddo
             !
             !Imaginary time:
             do m=0,Ltau 
                it=tau(m)
                spinChi_tau(iorb,m)=spinChi_tau(iorb,m) + exp(-it*Ei)*exp(-(beta-it)*Ej)*peso
             enddo
             !
          enddo
       enddo
    enddo
  end subroutine full_ed_build_spinChi_main



END MODULE ED_GF_CHISPIN
























