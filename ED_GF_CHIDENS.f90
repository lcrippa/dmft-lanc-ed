MODULE ED_GF_CHIDENS
  USE ED_GF_SHARED
  implicit none
  private


  public :: build_chi_dens

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
  !                            DENS
  !PURPOSE  : Evaluate the Dens susceptibility \Chi_dens for a 
  ! \chi_ab = <n_a(\tau)n_b(0)>
  !+------------------------------------------------------------------+
  subroutine build_chi_dens()
    integer :: iorb
    if(ed_diag_type=='full')return
    write(LOGfile,"(A)")"Get impurity dens Chi:"
    do iorb=1,Norb
       write(LOGfile,"(A)")"Get Chi_dens_l"//reg(txtfy(iorb))
       if(MPIMASTER)call start_timer()
       call lanc_ed_build_densChi_main(iorb)
       if(MPIMASTER)call stop_timer(LOGfile)
    enddo
    !
    if(Norb>1)then
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             write(LOGfile,"(A)")"Get Chi_dens_mix_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
             if(MPIMASTER)call start_timer()
             call lanc_ed_build_densChi_mix_main(iorb,jorb)
             if(MPIMASTER)call stop_timer(LOGfile)
          end do
       end do
       !
       write(LOGfile,"(A)")"Get Chi_dens_tot"
       if(MPIMASTER)call start_timer()
       call lanc_ed_build_densChi_tot_main()
       if(MPIMASTER)call stop_timer(LOGfile)
       !
       !
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             densChi_w(iorb,jorb,:)   = 0.5d0*(densChi_w(iorb,jorb,:) - densChi_w(iorb,iorb,:) - densChi_w(jorb,jorb,:))
             densChi_tau(iorb,jorb,:) = 0.5d0*(densChi_tau(iorb,jorb,:) - densChi_tau(iorb,iorb,:) - densChi_tau(jorb,jorb,:))
             densChi_iv(iorb,jorb,:)  = 0.5d0*(densChi_iv(iorb,jorb,:) - densChi_iv(iorb,iorb,:) - densChi_iv(jorb,jorb,:))
             !
             densChi_w(jorb,iorb,:)   = densChi_w(iorb,jorb,:)
             densChi_tau(jorb,iorb,:) = densChi_tau(iorb,jorb,:)
             densChi_iv(jorb,iorb,:)  = densChi_iv(iorb,jorb,:)
          enddo
       enddo
    endif
    !
    !> Guess it is wrong!
    densChi_tau = densChi_tau/zeta_function
    densChi_w   = densChi_w/zeta_function
    densChi_iv  = densChi_iv/zeta_function
    !
  end subroutine build_chi_dens






  !################################################################
  !################################################################
  !################################################################
  !################################################################






  subroutine lanc_ed_build_densChi_main(iorb)
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
          if(ed_verbose==3)write(LOGfile,"(A,I12)")'Apply N:',isector
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
             sgn = nud(1,iorb1)+nud(2,iorb1)
             !
             vvinit(i) = sgn*state_cvec(i)
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
       call add_to_lanczos_densChi(norm2,state_e,alfa_,beta_,1,iorb,iorb)
       !holes
       call add_to_lanczos_densChi(norm2,state_e,alfa_,beta_,-1,iorb,iorb)
       !
       call delete_Hv_sector()
       !
       deallocate(alfa_,beta_)
       if(allocated(vvinit))deallocate(vvinit)
       if(allocated(vvloc))deallocate(vvloc)
       nullify(state_cvec)
    enddo
    return
  end subroutine lanc_ed_build_densChi_main




  !################################################################




  subroutine lanc_ed_build_densChi_tot_main()
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
             sgn = sum(nud(1,1:Norb)) + sum(nud(2,1:Norb))
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
       call add_to_lanczos_densChi(norm2,state_e,alfa_,beta_,1,Norb+1,Norb+1)
       !holes
       call add_to_lanczos_densChi(norm2,state_e,alfa_,beta_,-1,Norb+1,Norb+1)
       !
       call delete_Hv_sector()
       !
       deallocate(alfa_,beta_)
       if(allocated(vvinit))deallocate(vvinit)
       if(allocated(vvloc))deallocate(vvloc)
       nullify(state_cvec)
    enddo
    return
  end subroutine lanc_ed_build_densChi_tot_main



  !################################################################


  subroutine lanc_ed_build_densChi_mix_main(iorb,jorb)
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
       !EVALUATE (N_jorb + N_iorb)|gs> = N_jorb|gs> + N_iorb|gs>
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
             Siorb    = nud(1,iorb1) + nud(2,iorb1)
             !
             iud(1)   = HI(jalfa)%map(Indices(jalfa))
             iud(2)   = HI(jalfa+Ns_Ud)%map(Indices(jalfa+Ns_Ud))
             nud(1,:) = Bdecomp(iud(1),Ns_Orb)
             nud(2,:) = Bdecomp(iud(2),Ns_Orb)
             Sjorb    = nud(1,jorb1) + nud(2,jorb1)
             !
             sgn       = Siorb + Sjorb
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
       call add_to_lanczos_densChi(norm2,state_e,alfa_,beta_,1,iorb,jorb)
       !holes
       call add_to_lanczos_densChi(norm2,state_e,alfa_,beta_,-1,iorb,jorb)
       !
       call delete_Hv_sector()
       !
       deallocate(alfa_,beta_)
       if(allocated(vvinit))deallocate(vvinit)
       if(allocated(vvloc))deallocate(vvloc)
       nullify(state_cvec)
    enddo
    return
  end subroutine lanc_ed_build_densChi_mix_main






  !################################################################


  subroutine add_to_lanczos_densChi(vnorm2,Ei,alanc,blanc,isign,iorb,jorb)
    integer                                    :: iorb,jorb,isign
    real(8)                                    :: pesoF,pesoAB,pesoBZ,peso,vnorm2  
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






END MODULE ED_GF_CHIDENS
























