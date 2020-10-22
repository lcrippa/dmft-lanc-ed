MODULE ED_CHI_PAIR
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

  implicit none
  private


  public :: build_chi_pair

  integer             :: istate,iorb,jorb,ispin,jspin
  integer             :: isector,jsector,ksector
  integer             :: idim,idimUP,idimDW
  integer             :: jdim,jdimUP,jdimDW
  integer             :: kdim,kdimUP,kdimDW
  real(8),allocatable :: vvinit(:),vvloc(:),vvinit_tmp(:)
  real(8),allocatable :: alfa_(:),beta_(:)
  integer             :: ialfa,ibeta
  integer             :: jalfa,jbeta
  integer             :: iorb1,jorb1
  integer             :: r
  integer             :: i,iup,idw
  integer             :: j,jup,jdw  
  integer             :: m,mup,mdw
  integer             :: k,kup,kdw
  integer             :: iph,i_el
  integer             :: jph,j_el
  real(8)             :: sgn,norm2,norm0
  integer             :: Nitermax,Nlanc,vecDim

  real(8),dimension(:),pointer                :: state_cvec
  real(8)                                     :: state_e

contains


  !+------------------------------------------------------------------+
  !                            PAIR
  !PURPOSE  : Evaluate the pair susceptibility \Chi_pair for a 
  ! \chi_ab = <Delta*_a(\tau)Delta_b(0)>
  !+------------------------------------------------------------------+
  subroutine build_chi_pair()
    write(LOGfile,"(A)")"Get impurity pair Chi:"
    do iorb=1,Norb
       write(LOGfile,"(A)")"Get Chi_pair_l"//reg(txtfy(iorb))
       if(MPIMASTER)call start_timer()
       select case(ed_diag_type)
       case default
          call lanc_ed_build_pairChi_diag(iorb)
       case ("full")
          write(LOGfile,"(A)")"Chi_pair not available in Full ED"
       end select
       if(MPIMASTER)call stop_timer(unit=LOGfile)
    enddo

    if(Norb>1)then
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             write(LOGfile,"(A)")"Get Chi_pair_mix_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
             if(MPIMASTER)call start_timer()
             select case(ed_diag_type)
             case default
                call lanc_ed_build_pairChi_mix(iorb,jorb)
             case ("full")
                write(LOGfile,"(A)")"Chi_pair not available in Full ED"
             end select
             if(MPIMASTER)call stop_timer(unit=LOGfile)
          end do
       end do
       !
       !
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             select case(ed_diag_type)
             case default
                pairChi_w(iorb,jorb,:)   = 0.5d0*(pairChi_w(iorb,jorb,:) - pairChi_w(iorb,iorb,:) - pairChi_w(jorb,jorb,:))
                pairChi_tau(iorb,jorb,:) = 0.5d0*(pairChi_tau(iorb,jorb,:) - pairChi_tau(iorb,iorb,:) - pairChi_tau(jorb,jorb,:))
                pairChi_iv(iorb,jorb,:)  = 0.5d0*(pairChi_iv(iorb,jorb,:) - pairChi_iv(iorb,iorb,:) - pairChi_iv(jorb,jorb,:))
                !
             case ("full")
                !
             end select
             !
             pairChi_w(jorb,iorb,:)   = pairChi_w(iorb,jorb,:)
             pairChi_tau(jorb,iorb,:) = pairChi_tau(iorb,jorb,:)
             pairChi_iv(jorb,iorb,:)  = pairChi_iv(iorb,jorb,:)
          enddo
       enddo
    endif
    !

  end subroutine build_chi_pair




  ! \chi_aa = <Delta*_a(\tau)Delta_a(0)>
  !         = <[C^+_a(\tau)C^+_a(\tau)][C_a(0)C_a(0)]>
  subroutine lanc_ed_build_pairChi_diag(iorb)
    integer                     :: iorb
    integer,dimension(2*Ns_Ud)  :: Indices
    integer,dimension(2*Ns_Ud)  :: Jndices
    integer,dimension(2*Ns_Ud)  :: Kndices
    integer,dimension(Ns_Ud)    :: iDimUps,iDimDws
    integer,dimension(Ns_Ud)    :: jDimUps,jDimDws
    integer,dimension(Ns_Ud)    :: kDimUps,kDimDws
    integer,dimension(2,Ns_Orb) :: Nud
    integer                     :: Iud(2)
    type(sector_map)            :: HI(2*Ns_Ud),HJ(2*Ns_Ud),HK(2*Ns_Ud)
    integer                     :: Nups(Ns_Ud)
    integer                     :: Ndws(Ns_Ud)
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
       call get_Nup(isector,Nups)
       call get_Ndw(isector,Ndws)
       if(MpiMaster.AND.ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")'From sector:',isector,Nups,Ndws
       !
       !
       idim  = getdim(isector)
       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       !
       jsector = getCsector(ialfa,ispin,isector)
       jdim    = getdim(jsector)
       call get_DimUp(jsector,jDimUps)
       call get_DImDw(jsector,jDimDws)
       jDimUp = product(jDimUps)
       jDimDw = product(jDimDws)
       !
       ksector = getCsector(ialfa,ispin,jsector)
       kdim    = getdim(ksector)
       call get_DimUp(ksector,kDimUps)
       call get_DImDw(ksector,kDimDws)
       kDimUp = product(kDimUps)
       kDimDw = product(kDimDws)
       !
       if(jsector/=0.AND.ksector/=0)then
          !
          if(MpiMaster)then
             call build_sector(isector,HI)
             call build_sector(jsector,HJ)
             call build_sector(ksector,HK)
             !
             if(ed_verbose==3)write(LOGfile,"(A,I12)")'Apply Cup*Cdw:',isector
             !DW
             ispin  = 2
             ibeta  = ialfa + (ispin-1)*Ns_Ud
             allocate(vvinit_tmp(jdim)) ;  vvinit_tmp=0d0
             do i=1,iDim
                iph = (i-1)/(iDimUp*iDimDw) + 1
                i_el = mod(i-1,iDimUp*iDimDw) + 1
                !
                call state2indices(i_el,[iDimUps,iDimDws],Indices)
                iud(1)   = HI(ialfa)%map(Indices(ialfa))
                iud(2)   = HI(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
                nud(1,:) = Bdecomp(iud(1),Ns_Orb)
                nud(2,:) = Bdecomp(iud(2),Ns_Orb)
                if(nud(ispin,iorb1)/=1)cycle
                call c(iorb1,iud(ispin),r,sgn)
                Jndices        = Indices
                Jndices(ibeta) = binary_search(HJ(ibeta)%map,r)
                call indices2state(Jndices,[jDimUps,jDimDws],j)
                j = j + (iph-1)*jDimUp*jDimDw
                !
                vvinit_tmp(j) = sgn*state_cvec(i)
             enddo
             !UP
             ispin  = 1
             ibeta  = ialfa + (ispin-1)*Ns_Ud
             allocate(vvinit(kdim)) ;  vvinit=0d0
             do j=1,jDim
                jph = (j-1)/(jDimUp*jDimDw) + 1
                j_el = mod(j-1,jDimUp*jDimDw) + 1
                !
                call state2indices(j_el,[jDimUps,jDimDws],Jndices)
                iud(1)   = HJ(ialfa)%map(Jndices(ialfa))
                iud(2)   = HJ(ialfa+Ns_Ud)%map(Jndices(ialfa+Ns_Ud))
                nud(1,:) = Bdecomp(iud(1),Ns_Orb)
                nud(2,:) = Bdecomp(iud(2),Ns_Orb)
                if(nud(ispin,iorb1)/=1)cycle
                call c(iorb1,iud(ispin),r,sgn)
                Kndices        = Jndices
                Kndices(ibeta) = binary_search(HK(ibeta)%map,r)
                call indices2state(Kndices,[kDimUps,kDimDws],k)
                k = k + (jph-1)*kDimUp*kDimDw
                !
                vvinit(k) = sgn*vvinit_tmp(j)
             enddo
             deallocate(vvinit_tmp)
             call delete_sector(isector,HI)
             call delete_sector(jsector,HJ)
             call delete_sector(ksector,HK)             
             !
             norm2=dot_product(vvinit,vvinit)
             vvinit=vvinit/sqrt(norm2)
          else
             allocate(vvinit(1));vvinit=0.d0
          endif
          !
          nlanc=min(idim,lanc_nGFiter)
          allocate(alfa_(nlanc),beta_(nlanc))
          alfa_=0.d0
          beta_=0.d0
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
          call delete_Hv_sector()
          !
          call add_to_lanczos_pairChi(norm2,state_e,alfa_,beta_,iorb,iorb)
          !
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
          if(allocated(vvloc))deallocate(vvloc)
       endif
       !
#ifdef _MPI
       if(MpiStatus)then
          if(associated(state_cvec))deallocate(state_cvec)
       else
          if(associated(state_cvec))nullify(state_cvec)
       endif
#else
       if(associated(state_cvec))nullify(state_cvec)
#endif
       !
    enddo
    return
  end subroutine lanc_ed_build_pairChi_diag






  ! \chi_ab = <Delta*_a(\tau)Delta_b(0)>
  !         = <[C^+_a(\tau)C^+_a(\tau)][C_b(0)C_b(0)]>
  !from aux: <[C^+_a C^+_a + C^+_b C^+_b][C_a C_a + C_b C_b]>  
  subroutine lanc_ed_build_pairChi_mix(iorb,jorb)
    integer                     :: iorb,jorb
    integer,dimension(2*Ns_Ud)  :: Indices
    integer,dimension(2*Ns_Ud)  :: Jndices
    integer,dimension(2*Ns_Ud)  :: Kndices
    integer,dimension(Ns_Ud)    :: iDimUps,iDimDws
    integer,dimension(Ns_Ud)    :: jDimUps,jDimDws
    integer,dimension(Ns_Ud)    :: kDimUps,kDimDws
    integer,dimension(2,Ns_Orb) :: Nud
    integer                     :: Iud(2)
    type(sector_map)            :: HI(2*Ns_Ud),HJ(2*Ns_Ud),HK(2*Ns_Ud)
    integer                     :: Nups(Ns_Ud)
    integer                     :: Ndws(Ns_Ud)
    !
    if(ed_total_ud)then
       ialfa = 1
       jalfa = 1
       iorb1 = iorb
       jorb1 = jorb
    else
       ! ialfa = iorb; jalfa = jorb; iorb1 = 1; jorb1 = 1
       write(LOGfile,"(A)")"ED_CHI_PAIR warning: can not evaluate \Chi_pair_ab with ed_total_ud=F"
       return
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
       ! --> Apply [C_b C_b + C_a C_a]|state>
       ! ialfa=jalfa => jsector,ksector are the same for iorb and jorb
       ! i.e. the output vector dimension vvinit is the same.
       jsector = getCsector(ialfa,ispin,isector)
       jdim    = getdim(jsector)
       call get_DimUp(jsector,jDimUps)
       call get_DImDw(jsector,jDimDws)
       jDimUp = product(jDimUps)
       jDimDw = product(jDimDws)
       !
       ksector = getCsector(ialfa,ispin,jsector)
       kdim    = getdim(ksector)
       call get_DimUp(ksector,kDimUps)
       call get_DImDw(ksector,kDimDws)
       kDimUp = product(kDimUps)
       kDimDw = product(kDimDws)
       !
       if(jsector/=0.AND.ksector/=0)then
          !
          if(MpiMaster)then
             call build_sector(isector,HI)
             call build_sector(jsector,HJ)
             call build_sector(ksector,HK)
             !
             if(ed_verbose==3)write(LOGfile,"(A,I12)")'Apply C_a,up*C_a,dw:',isector
             !DW
             ispin  = 2
             ibeta  = ialfa + (ispin-1)*Ns_Ud
             allocate(vvinit_tmp(jdim)) ;  vvinit_tmp=0d0
             do i=1,iDim
                iph = (i-1)/(iDimUp*iDimDw) + 1
                i_el = mod(i-1,iDimUp*iDimDw) + 1
                !
                call state2indices(i_el,[iDimUps,iDimDws],Indices)
                iud(1)   = HI(ialfa)%map(Indices(ialfa))
                iud(2)   = HI(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
                nud(1,:) = Bdecomp(iud(1),Ns_Orb)
                nud(2,:) = Bdecomp(iud(2),Ns_Orb)
                if(nud(ispin,iorb1)/=1)cycle
                call c(iorb1,iud(ispin),r,sgn)
                Jndices        = Indices
                Jndices(ibeta) = binary_search(HJ(ibeta)%map,r)
                call indices2state(Jndices,[jDimUps,jDimDws],j)
                j = j + (iph-1)*jDimUp*jDimDw
                !
                vvinit_tmp(j) = sgn*state_cvec(i)
             enddo
             !
             !UP
             ispin  = 1
             ibeta  = ialfa + (ispin-1)*Ns_Ud
             allocate(vvinit(kdim)) ;  vvinit=0d0
             do j=1,jDim
                jph = (j-1)/(jDimUp*jDimDw) + 1
                j_el = mod(j-1,jDimUp*jDimDw) + 1
                !
                call state2indices(j_el,[jDimUps,jDimDws],Jndices)
                iud(1)   = HJ(ialfa)%map(Jndices(ialfa))
                iud(2)   = HJ(ialfa+Ns_Ud)%map(Jndices(ialfa+Ns_Ud))
                nud(1,:) = Bdecomp(iud(1),Ns_Orb)
                nud(2,:) = Bdecomp(iud(2),Ns_Orb)
                if(nud(ispin,iorb1)/=1)cycle
                call c(iorb1,iud(ispin),r,sgn)
                Kndices        = Jndices
                Kndices(ibeta) = binary_search(HK(ibeta)%map,r)
                call indices2state(Kndices,[kDimUps,kDimDws],k)
                k = k + (jph-1)*kDimUp*kDimDw
                !
                vvinit(k) = sgn*vvinit_tmp(j)
             enddo
             deallocate(vvinit_tmp)
             !
             !
             if(ed_verbose==3)write(LOGfile,"(A,I12)")'Apply + C_b,up*C_b,dw:',isector
             !DW
             ispin  = 2
             ibeta  = jalfa + (ispin-1)*Ns_Ud
             allocate(vvinit_tmp(jdim)) ;  vvinit_tmp=0d0
             do i=1,iDim
                iph = (i-1)/(iDimUp*iDimDw) + 1
                i_el = mod(i-1,iDimUp*iDimDw) + 1
                !
                call state2indices(i_el,[iDimUps,iDimDws],Indices)
                iud(1)   = HI(jalfa)%map(Indices(jalfa))
                iud(2)   = HI(jalfa+Ns_Ud)%map(Indices(jalfa+Ns_Ud))
                nud(1,:) = Bdecomp(iud(1),Ns_Orb)
                nud(2,:) = Bdecomp(iud(2),Ns_Orb)
                if(nud(ispin,jorb1)/=1)cycle
                call c(jorb1,iud(ispin),r,sgn)
                Jndices        = Indices
                Jndices(ibeta) = binary_search(HJ(ibeta)%map,r)
                call indices2state(Jndices,[jDimUps,jDimDws],j)
                j = j + (iph-1)*jDimUp*jDimDw
                !
                vvinit_tmp(j) = sgn*state_cvec(i)
             enddo
             !
             !UP
             ispin  = 1
             ibeta  = jalfa + (ispin-1)*Ns_Ud
             do j=1,jDim
                jph = (j-1)/(jDimUp*jDimDw) + 1
                j_el = mod(j-1,jDimUp*jDimDw) + 1
                !
                call state2indices(j_el,[jDimUps,jDimDws],Jndices)
                iud(1)   = HJ(jalfa)%map(Jndices(jalfa))
                iud(2)   = HJ(jalfa+Ns_Ud)%map(Jndices(jalfa+Ns_Ud))
                nud(1,:) = Bdecomp(iud(1),Ns_Orb)
                nud(2,:) = Bdecomp(iud(2),Ns_Orb)
                if(nud(ispin,jorb1)/=1)cycle
                call c(jorb1,iud(ispin),r,sgn)
                Kndices        = Jndices
                Kndices(ibeta) = binary_search(HK(ibeta)%map,r)
                call indices2state(Kndices,[kDimUps,kDimDws],k)
                k = k + (jph-1)*kDimUp*kDimDw
                !
                vvinit(k) = vvinit(k) + sgn*vvinit_tmp(j)
             enddo
             deallocate(vvinit_tmp)
             call delete_sector(isector,HI)
             call delete_sector(jsector,HJ)
             call delete_sector(ksector,HK)             
             !
             norm2=dot_product(vvinit,vvinit)
             vvinit=vvinit/sqrt(norm2)
          else
             allocate(vvinit(1));vvinit=0.d0
          endif
          !
          nlanc=min(idim,lanc_nGFiter)
          allocate(alfa_(nlanc),beta_(nlanc))
          alfa_=0.d0
          beta_=0.d0
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
          call delete_Hv_sector()
          !
          call add_to_lanczos_pairChi(norm2,state_e,alfa_,beta_,iorb,jorb)
          !
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
          if(allocated(vvloc))deallocate(vvloc)
       endif
       !
#ifdef _MPI
       if(MpiStatus)then
          if(associated(state_cvec))deallocate(state_cvec)
       else
          if(associated(state_cvec))nullify(state_cvec)
       endif
#else
       if(associated(state_cvec))nullify(state_cvec)
#endif
       !
    enddo
    return
  end subroutine lanc_ed_build_pairChi_mix

  !################################################################
  !################################################################
  !################################################################
  !################################################################






  subroutine add_to_lanczos_pairChi(vnorm2,Ei,alanc,blanc,iorb,jorb)
    integer                                    :: iorb,jorb
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
       if(beta*dE > 1d-3)pairChi_iv(iorb,jorb,0)=pairChi_iv(iorb,jorb,0) + peso*2*(1d0-exp(-beta*dE))/dE 
       do i=1,Lmats
          pairChi_iv(iorb,jorb,i)=pairChi_iv(iorb,jorb,i) + peso*(1d0-exp(-beta*dE))*2d0*dE/(vm(i)**2+dE**2)
       enddo
       do i=0,Ltau
          pairChi_tau(iorb,jorb,i)=pairChi_tau(iorb,jorb,i) + exp(-tau(i)*dE)*peso
       enddo
       do i=1,Lreal
          pairChi_w(iorb,jorb,i)=pairChi_w(iorb,jorb,i) - &
               peso*(1d0-exp(-beta*dE))*(1d0/(dcmplx(vr(i),eps) - dE) - 1d0/(dcmplx(vr(i),eps) + dE))
       enddo
    enddo
    !
  end subroutine add_to_lanczos_pairChi



END MODULE ED_CHI_PAIR
























