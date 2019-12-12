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
  integer             :: ialfa
  integer             :: jalfa
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
  ! note: as S_a is hermitian particle and holes contributions
  ! are identical so work out only one lanczos tridiag. work out the 
  ! reduction for both values of isign in the same call.
  !+------------------------------------------------------------------+
  subroutine build_chi_spin()
    write(LOGfile,"(A)")"Get impurity spin Chi:"
    do iorb=1,Norb
       write(LOGfile,"(A)")"Get Chi_spin_l"//reg(txtfy(iorb))
       if(MPIMASTER)call start_timer()
       select case(ed_diag_type)
       case default
          call lanc_ed_build_spinChi_main(iorb)
       case ("full")
          call full_ed_build_spinChi_main(iorb,iorb)
       end select
       if(MPIMASTER)call stop_timer(LOGfile)
    enddo
    !
    if(Norb>1)then
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             write(LOGfile,"(A)")"Get Chi_spin_mix_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
             if(MPIMASTER)call start_timer()
             select case(ed_diag_type)
             case default
                call lanc_ed_build_spinChi_mix_main(iorb,jorb)
             case ("full")
                call full_ed_build_spinChi_main(iorb,jorb)
             end select
             if(MPIMASTER)call stop_timer(LOGfile)
          end do
       end do
       !
       write(LOGfile,"(A)")"Get Chi_spin_tot"
       if(MPIMASTER)call start_timer()
       select case(ed_diag_type)
       case default
          call lanc_ed_build_spinChi_tot_main()
       case ("full")
       ! Chi_spin_tot not implemented yet in the FULL ED case
       end select
       if(MPIMASTER)call stop_timer(LOGfile)
       !
       !
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             select case(ed_diag_type)
             case default
                spinChi_w(iorb,jorb,:)   = 0.5d0*(spinChi_w(iorb,jorb,:) - spinChi_w(iorb,iorb,:) - spinChi_w(jorb,jorb,:))
                spinChi_tau(iorb,jorb,:) = 0.5d0*(spinChi_tau(iorb,jorb,:) - spinChi_tau(iorb,iorb,:) - spinChi_tau(jorb,jorb,:))
                spinChi_iv(iorb,jorb,:)  = 0.5d0*(spinChi_iv(iorb,jorb,:) - spinChi_iv(iorb,iorb,:) - spinChi_iv(jorb,jorb,:))
                !
             case ("full")
             ! The previous calculation is not needed in the FULL ED case
             end select
             !
             spinChi_w(jorb,iorb,:)   = spinChi_w(iorb,jorb,:)
             spinChi_tau(jorb,iorb,:) = spinChi_tau(iorb,jorb,:)
             spinChi_iv(jorb,iorb,:)  = spinChi_iv(iorb,jorb,:)
          enddo
       enddo
    endif
    !
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
       else
          allocate(vvinit(0))
       endif
       !
       nlanc=min(idim,lanc_nGFiter)
       allocate(alfa_(nlanc),beta_(nlanc));alfa_=0d0;beta_=0d0
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
       call add_to_lanczos_spinChi(norm2,state_e,alfa_,beta_,iorb,iorb)
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
             sgn = 0.5d0*(Sup - Sdw)
             !
             vvinit(i) = sgn*state_cvec(i)
          enddo
          call delete_sector(isector,HI)
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)
       else
          allocate(vvinit(0))
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
       call add_to_lanczos_spinChi(norm2,state_e,alfa_,beta_,Norb+1,Norb+1)
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


  subroutine lanc_ed_build_spinChi_mix_main(iorb,jorb)
    integer                     :: iorb,jorb
    integer,dimension(2*Ns_Ud)  :: Indices
    integer,dimension(2*Ns_Ud)  :: Jndices
    integer,dimension(Ns_Ud)    :: iDimUps,iDimDws
    integer,dimension(Ns_Ud)    :: jDimUps,jDimDws
    integer,dimension(2,Ns_Orb) :: Nud
    integer                     :: Iud(2)
    real(8)                     :: Siorb,Sjorb
    type(sector_map)            :: HI(2*Ns_Ud)    !map of the Sector S to Hilbert space H
    !
    if(ed_total_ud)then
       ialfa = 1
       jalfa = 1
       iorb1 = iorb
       jorb1 = jorb
    else
       ialfa = iorb
       jalfa = jorb
       iorb1 = 1
       jorb1 = 1
    endif
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
       idim  = getdim(isector)
       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       !
       !EVALUATE (Sz_jorb + Sz_iorb)|gs> = Sz_jorb|gs> + Sz_iorb|gs>
       if(MpiMaster)then
          if(ed_verbose==3)write(LOGfile,"(A,I15)")'Apply Na*Nb:',isector
          !
          allocate(vvinit(idim));vvinit=0.d0
          !
          call build_sector(isector,HI)
          do i=1,iDim
             call state2indices(i,[iDimUps,iDimDws],Indices)
             !
             iud(1)   = HI(ialfa)%map(Indices(ialfa))
             iud(2)   = HI(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
             nud(1,:) = Bdecomp(iud(1),Ns_Orb)
             nud(2,:) = Bdecomp(iud(2),Ns_Orb)
             Siorb    = nud(1,iorb1) - nud(2,iorb1)
             !
             iud(1)   = HI(jalfa)%map(Indices(jalfa))
             iud(2)   = HI(jalfa+Ns_Ud)%map(Indices(jalfa+Ns_Ud))
             nud(1,:) = Bdecomp(iud(1),Ns_Orb)
             nud(2,:) = Bdecomp(iud(2),Ns_Orb)
             Sjorb    = nud(1,jorb1) - nud(2,jorb1)
             !
             sgn       = 0.5d0*Siorb + 0.5d0*Sjorb
             vvinit(i) = sgn*state_cvec(i)
          enddo
          call delete_sector(isector,HI)
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)
       else
          allocate(vvinit(0))
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
       call add_to_lanczos_spinChi(norm2,state_e,alfa_,beta_,iorb,jorb)
       !
       call delete_Hv_sector()
       !
       deallocate(alfa_,beta_)
       if(allocated(vvinit))deallocate(vvinit)
       if(allocated(vvloc))deallocate(vvloc)
       nullify(state_cvec)
    enddo
    return
  end subroutine lanc_ed_build_spinChi_mix_main




  !################################################################


  subroutine add_to_lanczos_spinChi(vnorm2,Ei,alanc,blanc,iorb,jorb)
    real(8)                                    :: vnorm2,Ei,Ej,Egs,pesoF,pesoAB,pesoBZ,de,peso
    integer                                    :: nlanc
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: iorb,jorb
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
    do j=1,nlanc
       Ej     = diag(j)
       dE     = Ej-Ei
       pesoAB = Z(1,j)*Z(1,j)
       peso   = pesoF*pesoAB*pesoBZ
       ! the correct behavior for beta*dE << 1 is recovered only by assuming that v_n is still finite
       ! beta*dE << v_n for v_n--> 0 slower. First limit beta*dE--> 0 and only then v_n -->0.
       ! This ensures that the correct null contribution is obtained.
       ! So we impose that: if (beta*dE is larger than a small qty) we sum up the contribution, else
       ! we do not include the contribution (because we are in the situation described above).
       ! For the real-axis case this problem is circumvented by the usual i*0+ = xi*eps
       if(beta*dE > 1d-3)spinChi_iv(iorb,jorb,0)=spinChi_iv(iorb,jorb,0) + peso*2*(1d0-exp(-beta*dE))/dE 
       do i=1,Lmats
          spinChi_iv(iorb,jorb,i)=spinChi_iv(iorb,jorb,i) + peso*(1d0-exp(-beta*dE))*2d0*dE/(vm(i)**2+dE**2)
       enddo
       do i=0,Ltau
          spinChi_tau(iorb,jorb,i)=spinChi_tau(iorb,jorb,i) + exp(-tau(i)*dE)*peso
       enddo
       do i=1,Lreal
          spinChi_w(iorb,jorb,i)=spinChi_w(iorb,jorb,i) - peso*(1d0-exp(-beta*dE))*(1d0/(dcmplx(vr(i),eps) - dE) - 1d0/(dcmplx(vr(i),eps) + dE))
       enddo
    enddo
  end subroutine add_to_lanczos_spinChi




  !################################################################
  !################################################################
  !################################################################
  !################################################################



  subroutine full_ed_build_spinChi_main(iorb,jorb)
    integer                     :: iorb,jorb
    integer,dimension(2*Ns_Ud)  :: Indices
    integer,dimension(2*Ns_Ud)  :: Jndices
    integer,dimension(Ns_Ud)    :: iDimUps,iDimDws
    integer,dimension(Ns_Ud)    :: jDimUps,jDimDws
    integer,dimension(2,Ns_Orb) :: Nud
    integer                     :: Iud(2)
    type(sector_map)            :: HI(2*Ns_Ud)
    real(8)                     :: Chiorb,Chjorb,Siorb,Sjorb
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
       jalfa = 1
       iorb1 = iorb
       jorb1 = jorb
    else
       ialfa = iorb
       jalfa = jorb
       iorb1 = 1
       jorb1 = 1
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
             Chiorb=0d0
             Chjorb=0d0
             expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(isector)%e(j))
             if(expterm<cutoff)cycle
             do ll=1,idim
                call state2indices(ll,[iDimUps,iDimDws],Indices)
                iud(1)   = HI(ialfa)%map(Indices(ialfa))
                iud(2)   = HI(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
                nud(1,:) = Bdecomp(iud(1),Ns_Orb)
                nud(2,:) = Bdecomp(iud(2),Ns_Orb)
                Siorb    = nud(1,iorb1) - nud(2,iorb1)
                Chiorb   = Chiorb + espace(isector)%M(ll,i)*0.5d0*Siorb*espace(isector)%M(ll,j)
                !
                iud(1)   = HI(jalfa)%map(Indices(jalfa))
                iud(2)   = HI(jalfa+Ns_Ud)%map(Indices(jalfa+Ns_Ud))
                nud(1,:) = Bdecomp(iud(1),Ns_Orb)
                nud(2,:) = Bdecomp(iud(2),Ns_Orb)
                Sjorb    = nud(1,jorb1) - nud(2,jorb1)
                Chjorb   = Chjorb + espace(isector)%M(ll,i)*0.5d0*Sjorb*espace(isector)%M(ll,j)
             enddo
             Ei=espace(isector)%e(i)
             Ej=espace(isector)%e(j)
             de=Ei-Ej
             peso = Chiorb*Chjorb/zeta_function
             !
             !Matsubara (bosonic) frequency
             if(beta*dE > 1d-3)spinChi_iv(iorb,jorb,0)=spinChi_iv(iorb,jorb,0) + peso*2*exp(-beta*Ej)*(1d0-exp(-beta*dE))/dE
             do m=1,Lmats
                spinChi_iv(iorb,jorb,m)=spinChi_iv(iorb,jorb,m)+ peso*exp(-beta*Ej)*2*dE/(vm(m)**2 + de**2)
             enddo
             !
             !Imaginary time: V
             do m=0,Ltau 
                it=tau(m)
                spinChi_tau(iorb,jorb,m)=spinChi_tau(iorb,jorb,m) + exp(-it*Ei)*exp(-(beta-it)*Ej)*Chiorb*Chjorb
             enddo
             !
             !Real-frequency: Retarded = Commutator = response function
             do m=1,Lreal
                iw=dcmplx(vr(m),eps)
                spinChi_w(iorb,jorb,m)=spinChi_w(iorb,jorb,m)-peso*(exp(-beta*Ei) - exp(-beta*Ej))/(iw+de)
             enddo
             !
          enddo
       enddo
    enddo
  end subroutine full_ed_build_spinChi_main



END MODULE ED_GF_CHISPIN
























