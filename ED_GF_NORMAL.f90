MODULE ED_GF_NORMAL
  USE ED_GF_SHARED
  implicit none
  private


  public :: build_gf_normal
  public :: build_sigma_normal


  integer             :: istate
  integer             :: isector,jsector
  integer             :: idim,idimUP,idimDW
  integer             :: jdim,jdimUP,jdimDW
  real(8),allocatable :: vvinit(:),vvloc(:)
  real(8),allocatable :: alfa_(:),beta_(:)
  integer             :: ialfa,ibeta
  integer             :: jalfa,jbeta
  integer             :: iorb1,jorb1
  integer             :: r
  integer             :: i,iup,idw
  integer             :: j,jup,jdw  
  integer             :: m,mup,mdw
  integer             :: iph,i_el
  real(8)             :: sgn,norm2,norm0
  integer             :: Nitermax,Nlanc,vecDim


contains



  !+------------------------------------------------------------------+
  !                        NORMAL
  !+------------------------------------------------------------------+
  subroutine build_gf_normal()
    integer                                     :: iorb,jorb,ispin,i
    logical                                     :: MaskBool
    logical(8),dimension(Nspin,Nspin,Norb,Norb) :: Hmask
    !

    do ispin=1,Nspin
       do iorb=1,Norb
          write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_s"//str(ispin)
          if(MPIMASTER)call start_timer
          select case(ed_diag_type)
          case default
             call lanc_build_gf_normal_main(iorb,ispin)
          case ("full")
             call full_build_gf_normal_main(iorb,ispin)
          end select
          if(MPIMASTER)call stop_timer(unit=LOGfile)
       enddo
    enddo
    !
    if(offdiag_gf_flag)then
       Hmask=mask_hloc(impHloc,wdiag=.true.,uplo=.true.)
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=iorb+1,Norb
                MaskBool=.true.   
                if(bath_type=="replica")MaskBool=Hmask(ispin,ispin,iorb,jorb)
                if(.not.MaskBool)cycle
                !
                write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_m"//str(jorb)//"_s"//str(ispin)
                if(MPIMASTER)call start_timer
                select case(ed_diag_type)
                case default
                   call lanc_build_gf_normal_mix_main(iorb,jorb,ispin)
                case ("full")
                   call full_build_gf_normal_mix_main(iorb,jorb,ispin)
                end select
                if(MPIMASTER)call stop_timer(unit=LOGfile)
             enddo
          enddo
       enddo
       !
       !
       !Put here off-diagonal manipulation by symmetry:
       select case(ed_diag_type)
       case default
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=iorb+1,Norb
                   !if(hybrid)always T; if(replica)T iff following condition is T
                   MaskBool=.true.   
                   if(bath_type=="replica")MaskBool=Hmask(ispin,ispin,iorb,jorb)
                   !
                   if(.not.MaskBool)cycle
                   impGmats(ispin,ispin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,ispin,iorb,jorb,:) &
                        - impGmats(ispin,ispin,iorb,iorb,:) - impGmats(ispin,ispin,jorb,jorb,:))
                   impGreal(ispin,ispin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,ispin,iorb,jorb,:) &
                        - impGreal(ispin,ispin,iorb,iorb,:) - impGreal(ispin,ispin,jorb,jorb,:))
                   impGmats(ispin,ispin,jorb,iorb,:) = impGmats(ispin,ispin,iorb,jorb,:)
                   impGreal(ispin,ispin,jorb,iorb,:) = impGreal(ispin,ispin,iorb,jorb,:)
                enddo
             enddo
          enddo
       case ("full")
          !>>ACTHUNG: this relation might not be true, it depends on the value of the impHloc_ij
          ! if impHloc_ij is REAL then it is true. if CMPLX hermiticity must be ensured
          impGmats(ispin,ispin,jorb,iorb,:) = impGmats(ispin,ispin,iorb,jorb,:)
          impGreal(ispin,ispin,jorb,iorb,:) = impGreal(ispin,ispin,iorb,jorb,:)
       end select
    end if
    !
  end subroutine build_gf_normal






  !################################################################
  !################################################################
  !################################################################
  !################################################################






  subroutine lanc_build_gf_normal_main(iorb,ispin)
    integer,intent(in)          :: iorb,ispin
    integer,dimension(2*Ns_Ud)  :: Indices
    integer,dimension(2*Ns_Ud)  :: Jndices
    integer,dimension(Ns_Ud)    :: iDimUps,iDimDws
    integer,dimension(Ns_Ud)    :: jDimUps,jDimDws
    integer,dimension(2,Ns_Orb) :: Nud
    integer,dimension(2)        :: Iud,Jud
    type(sector_map)            :: HI(2*Ns_Ud),HJ(2*Ns_Ud)
    !
    integer                     :: Nups(Ns_Ud)
    integer                     :: Ndws(Ns_Ud)
    integer                     :: Mups(Ns_Ud)
    integer                     :: Mdws(Ns_Ud)
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
       idim  = getdim(isector)
       call get_Nup(isector,Nups)
       call get_Ndw(isector,Ndws)
       if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")'From sector:',isector,Nups,Ndws
       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       call build_sector(isector,HI)
       !
       !ADD ONE PARTICLE:
       jsector = getCDGsector(ialfa,ispin,isector)
       if(jsector/=0)then 
          !
          jdim   = getdim(jsector)
          call get_Nup(jsector,nups)
          call get_Ndw(jsector,ndws)
          call get_DimUp(jsector,jDimUps)
          call get_DImDw(jsector,jDimDws)
          jDimUp = product(jDimUps)
          jDimDw = product(jDimDws)
          !The Op|gs> is worked out by the master only:
          if(MpiMaster)then
             if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")' add particle:',jsector,Nups,Ndws
             !
             allocate(vvinit(jdim)) ; vvinit=zero
             !
             call build_sector(jsector,HJ)
             do i=1,iDim
                iph = (i-1)/(iDimUp*iDimDw) + 1
                i_el = mod(i-1,iDimUp*iDimDw) + 1
                !
                call state2indices(i_el,[iDimUps,iDimDws],Indices)
                iud(1)   = HI(ialfa)%map(Indices(ialfa))
                iud(2)   = HI(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
                nud(1,:) = Bdecomp(iud(1),Ns_Orb)
                nud(2,:) = Bdecomp(iud(2),Ns_Orb)
                !
                if(Nud(ispin,iorb1)/=0)cycle
                call cdg(iorb1,iud(ispin),r,sgn)
                Jndices        = Indices
                Jndices(ibeta) = binary_search(HJ(ibeta)%map,r)
                call indices2state(Jndices,[jDimUps,jDimDws],j)
                !
                j = j + (iph-1)*jDimUp*jDimDw
                !
                vvinit(j) = sgn*state_cvec(i)
             enddo
             call delete_sector(jsector,HJ)
             !
             norm2=dot_product(vvinit,vvinit)
             vvinit=vvinit/sqrt(norm2)
          else
             allocate(vvinit(1));vvinit=0.d0
          endif

          !
          nlanc=min(jdim,lanc_nGFiter)
          allocate(alfa_(nlanc),beta_(nlanc))
          alfa_=0.d0
          beta_=0.d0
          !
          call build_Hv_sector(jsector)
#ifdef _MPI
          if(MpiStatus)then
             call bcast_MPI(MpiComm,norm2)

             vecDim = vecDim_Hv_sector(jsector)
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
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,1,iorb,iorb,ispin)
          !
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
          if(allocated(vvloc))deallocate(vvloc)
       endif
       !
       !REMOVE ONE PARTICLE:
       jsector = getCsector(ialfa,ispin,isector)
       if(jsector/=0)then
          !            
          jdim   = getdim(jsector)
          call get_Nup(jsector,nups)
          call get_Ndw(jsector,ndws)
          call get_DimUp(jsector,jDimUps)
          call get_DImDw(jsector,jDimDws)
          jDimUp = product(jDimUps)
          jDimDw = product(jDimDws)
          !
          if(MpiMaster)then
             if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")' del particle:',jsector,Nups,Ndws
             allocate(vvinit(jdim)) ; vvinit=zero
             !
             call build_sector(jsector,HJ)
             do i=1,iDim
                iph = (i-1)/(iDimUp*iDimDw) + 1
                i_el = mod(i-1,iDimUp*iDimDw) + 1
                !
                call state2indices(i_el,[iDimUps,iDimDws],Indices)
                iud(1)   = HI(ialfa)%map(Indices(ialfa))
                iud(2)   = HI(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
                nud(1,:) = Bdecomp(iud(1),Ns_Orb)
                nud(2,:) = Bdecomp(iud(2),Ns_Orb)
                !
                if(nud(ispin,iorb1)/=1)cycle
                call c(iorb1,iud(ispin),r,sgn)
                Jndices        = Indices
                Jndices(ibeta) = binary_search(HJ(ibeta)%map,r)
                call indices2state(Jndices,[jDimUps,jDimDws],j)
                !
                j = j + (iph-1)*jDimUp*jDimDw
                !
                vvinit(j) = sgn*state_cvec(i)
             enddo
             call delete_sector(jsector,HJ)
             !
             norm2=dot_product(vvinit,vvinit)
             vvinit=vvinit/sqrt(norm2)
          else
             allocate(vvinit(1));vvinit=0.d0
          endif
          !
          nlanc=min(jdim,lanc_nGFiter)
          allocate(alfa_(nlanc),beta_(nlanc))
          alfa_=0.d0
          beta_=0.d0
          !
          call build_Hv_sector(jsector)
#ifdef _MPI
          if(MpiStatus)then
             call bcast_MPI(MpiComm,norm2)
             vecDim = vecDim_Hv_sector(jsector)
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
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,-1,iorb,iorb,ispin)
          !
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
          if(allocated(vvloc))deallocate(vvloc)
       endif
       !
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
       call delete_sector(isector,HI)
       !
    enddo
    return
  end subroutine lanc_build_gf_normal_main





  !################################################################






  subroutine lanc_build_gf_normal_mix_main(iorb,jorb,ispin)
    integer                     :: iorb,jorb,ispin,istate
    integer,dimension(2*Ns_Ud)  :: Indices
    integer,dimension(2*Ns_Ud)  :: Jndices
    integer,dimension(Ns_Ud)    :: iDimUps,iDimDws
    integer,dimension(Ns_Ud)    :: jDimUps,jDimDws
    integer,dimension(2,Ns_Orb) :: Nud
    integer,dimension(2)        :: iud,jud
    type(sector_map)            :: HI(2*Ns_Ud),HJ(2*Ns_Ud)
    !
    !
    if(ed_total_ud)then
       ialfa = 1
       jalfa = ialfa               !this is the condition to evaluate G_ab: ialfa=jalfa
       iorb1 = iorb
       jorb1 = jorb
    else
       ! ialfa = iorb; jalfa = jorb; iorb1 = 1; jorb1 = 1
       write(LOGfile,"(A)")"ED_GF_NORMAL warning: can not evaluate GF_ab with ed_total_ud=F"
       return
    endif
    ibeta  = ialfa + (ispin-1)*Ns_Ud
    jbeta  = jalfa + (ispin-1)*Ns_Ud
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
       call build_sector(isector,HI)
       !
       !
       !EVALUATE (c^+_iorb + c^+_jorb)|gs>
       jsector = getCDGsector(ialfa,ispin,isector)
       if(jsector/=0)then
          !
          jdim   = getdim(jsector)
          call get_DimUp(jsector,jDimUps)
          call get_DImDw(jsector,jDimDws)
          jDimUp = product(jDimUps)
          jDimDw = product(jDimDws)
          !
          if(MpiMaster)then
             if(ed_verbose>=3)write(LOGfile,"(A,I15)")' add particle:',jsector
             allocate(vvinit(jdim)) ; vvinit=zero
             !
             call build_sector(jsector,HJ)
             do i=1,iDim
                iph = (i-1)/(iDimUp*iDimDw) + 1
                i_el = mod(i-1,iDimUp*iDimDw) + 1
                !
                call state2indices(i_el,[iDimUps,iDimDws],Indices)
                iud(1)   = HI(ialfa)%map(Indices(ialfa))
                iud(2)   = HI(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
                nud(1,:) = Bdecomp(iud(1),Ns_Orb)
                nud(2,:) = Bdecomp(iud(2),Ns_Orb)
                !
                if(nud(ispin,iorb1)/=0)cycle
                call cdg(iorb1,iud(ispin),r,sgn)
                !
                Jndices        = Indices
                Jndices(ibeta) = binary_search(HJ(ibeta)%map,r)
                call indices2state(Jndices,[jDimUps,jDimDws],j)
                !
                j = j + (iph-1)*jDimUp*jDimDw
                !
                vvinit(j) = sgn*state_cvec(i)
             enddo
             do i=1,iDim
                iph = (i-1)/(iDimUp*iDimDw) + 1
                i_el = mod(i-1,iDimUp*iDimDw) + 1
                !
                call state2indices(i_el,[iDimUps,iDimDws],Indices)
                iud(1)   = HI(jalfa)%map(Indices(jalfa))
                iud(2)   = HI(jalfa+Ns_Ud)%map(Indices(jalfa+Ns_Ud))
                nud(1,:) = Bdecomp(iud(1),Ns_Orb)
                nud(2,:) = Bdecomp(iud(2),Ns_Orb)
                !
                if(nud(ispin,jorb1)/=0)cycle
                call cdg(jorb1,iud(ispin),r,sgn)
                !
                Jndices        = Indices
                Jndices(jbeta) = binary_search(HJ(jbeta)%map,r)
                call indices2state(Jndices,[jDimUps,jDimDws],j)
                !
                j = j + (iph-1)*jDimUp*jDimDw
                !
                vvinit(j) = vvinit(j) + sgn*state_cvec(i)
             enddo
             call delete_sector(jsector,HJ)
             !
             norm2=dot_product(vvinit,vvinit)
             vvinit=vvinit/sqrt(norm2)
          else
             allocate(vvinit(1));vvinit=0.d0
          endif
          !
          nlanc=min(jdim,lanc_nGFiter)
          allocate(alfa_(nlanc),beta_(nlanc))
          alfa_=0.d0
          beta_=0.d0
          !           
          call build_Hv_sector(jsector)
#ifdef _MPI
          if(MpiStatus)then
             call bcast_MPI(MpiComm,norm2)
             vecDim = vecDim_Hv_sector(jsector)
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
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,1,iorb,jorb,ispin)
          !
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
          if(allocated(vvloc))deallocate(vvloc)
       endif
       !
       !EVALUATE (c_iorb + c_jorb)|gs>
       jsector = getCsector(ialfa,ispin,isector)
       if(jsector/=0)then
          !
          jdim   = getdim(jsector)
          call get_DimUp(jsector,jDimUps)
          call get_DImDw(jsector,jDimDws)
          jDimUp = product(jDimUps)
          jDimDw = product(jDimDws)
          !
          if(MpiMaster)then
             if(ed_verbose>=3)write(LOGfile,"(A,I15)")' del particle:',jsector
             allocate(vvinit(jdim)) ; vvinit=zero
             !
             call build_sector(jsector,HJ)
             do i=1,iDim
                iph = (i-1)/(iDimUp*iDimDw) + 1
                i_el = mod(i-1,iDimUp*iDimDw) + 1
                !
                call state2indices(i_el,[iDimUps,iDimDws],Indices)
                iud(1)   = HI(ialfa)%map(Indices(ialfa))
                iud(2)   = HI(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
                nud(1,:) = Bdecomp(iud(1),Ns_Orb)
                nud(2,:) = Bdecomp(iud(2),Ns_Orb)
                !
                if(nud(ispin,iorb1)/=1)cycle
                call c(iorb1,iud(ispin),r,sgn)
                !
                Jndices        = Indices
                Jndices(ibeta) = binary_search(HJ(ibeta)%map,r)
                call indices2state(Jndices,[jDimUps,jDimDws],j)
                !
                j = j + (iph-1)*jDimUp*jDimDw
                !
                vvinit(j) = sgn*state_cvec(i)
             enddo
             do i=1,iDim
                iph = (i-1)/(iDimUp*iDimDw) + 1
                i_el = mod(i-1,iDimUp*iDimDw) + 1
                !
                call state2indices(i,[iDimUps,iDimDws],Indices)
                iud(1)   = HI(jalfa)%map(Indices(jalfa))
                iud(2)   = HI(jalfa+Ns_Ud)%map(Indices(jalfa+Ns_Ud))
                nud(1,:) = Bdecomp(iud(1),Ns_Orb)
                nud(2,:) = Bdecomp(iud(2),Ns_Orb)
                !
                if(nud(ispin,jorb1)/=1)cycle
                call c(jorb1,iud(ispin),r,sgn)
                !
                Jndices        = Indices
                Jndices(jbeta) = binary_search(HJ(jbeta)%map,r)
                call indices2state(Jndices,[jDimUps,jDimDws],j)
                !
                j = j + (iph-1)*jDimUp*jDimDw
                !
                vvinit(j) = vvinit(j) + sgn*state_cvec(i)
             enddo
             call delete_sector(jsector,HJ)
             !
             norm2=dot_product(vvinit,vvinit)
             vvinit=vvinit/sqrt(norm2)
          else
             allocate(vvinit(1));vvinit=0.d0
          endif
          !
          nlanc=min(jdim,lanc_nGFiter)
          allocate(alfa_(nlanc),beta_(nlanc))
          alfa_=0.d0
          beta_=0.d0
          !
          call build_Hv_sector(jsector)
#ifdef _MPI
          if(MpiStatus)then
             call bcast_MPI(MpiComm,norm2)
             vecDim = vecDim_Hv_sector(jsector)
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
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin)
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
       call delete_sector(isector,HI)
       !
    enddo
    return
  end subroutine lanc_build_gf_normal_mix_main





  !################################################################




  subroutine add_to_lanczos_gf_normal(vnorm2,Ei,alanc,blanc,isign,iorb,jorb,ispin)
    complex(8)                                 :: vnorm2,pesoBZ,peso
    real(8)                                    :: Ei,Egs,de
    integer                                    :: nlanc,itype
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: isign,iorb,jorb,ispin
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    if((finiteT).and.(beta*(Ei-Egs).lt.200))then
       pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
    elseif(.not.finiteT)then
       pesoBZ = vnorm2/zeta_function
    else
       pesoBZ=0.d0
    endif
    !
    !pesoBZ = vnorm2/zeta_function
    !if(finiteT)pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
    !
    !Only the nodes in Mpi_Comm_Group did get the alanc,blanc.
    !However after delete_sector_Hv MpiComm returns to be the global one
    !so we can safely Bcast the alanc,blanc (known only to the operative group)
    !to every nodes. The master is in charge of this (as a
    !participant of the operative group)
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,alanc)
       call Bcast_MPI(MpiComm,blanc)
    endif
#endif
    !
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    do j=1,nlanc
       de = diag(j)-Ei
       peso = pesoBZ*Z(1,j)*Z(1,j)
       do i=1,Lmats
          iw=xi*wm(i)
          impGmats(ispin,ispin,iorb,jorb,i)=impGmats(ispin,ispin,iorb,jorb,i) + peso/(iw-isign*de)
       enddo
       do i=1,Lreal
          iw=dcmplx(wr(i),eps)
          impGreal(ispin,ispin,iorb,jorb,i)=impGreal(ispin,ispin,iorb,jorb,i) + peso/(iw-isign*de)
       enddo
    enddo
  end subroutine add_to_lanczos_gf_normal






  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################





  subroutine full_build_gf_normal_main(iorb,ispin)
    integer                     :: iorb,ispin
    integer,dimension(2*Ns_Ud)  :: Indices
    integer,dimension(2*Ns_Ud)  :: Jndices
    integer,dimension(Ns_Ud)    :: iDimUps,iDimDws
    integer,dimension(Ns_Ud)    :: jDimUps,jDimDws
    integer,dimension(2,Ns_Orb) :: Nud
    integer,dimension(2)        :: Iud,Jud
    type(sector_map)            :: HI(2*Ns_Ud),HJ(2*Ns_Ud)
    real(8)                     :: op_mat(2)
    real(8)                     :: spectral_weight
    real(8)                     :: sgn_cdg,sgn_c
    integer                     :: m,i,j,r,k,p,li,rj
    real(8)                     :: Ei,Ej
    real(8)                     :: expterm,peso,de,w0
    complex(8)                  :: iw
    !    
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
    do isector=1,Nsectors
       jsector=getCDGsector(ialfa,ispin,isector)
       if(jsector==0)cycle
       !
       idim=getdim(isector)     !i-th sector dimension

       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       call build_sector(isector,HI)
       !
       jdim=getdim(jsector)     !j-th sector dimension
       call get_DimUp(jsector,jDimUps)
       call get_DImDw(jsector,jDimDws)
       jDimUp = product(jDimUps)
       jDimDw = product(jDimDws)
       call build_sector(jsector,HJ)
       !
       do i=1,iDim          !loop over the states in the i-th sect.
          do j=1,jDim       !loop over the states in the j-th sect.
             !
             expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(jsector)%e(j))
             if(expterm < cutoff)cycle
             !
             op_mat=0d0
             !
             do li=1,idim              !loop over the component of |I> (IN state!)
                iph = (li-1)/(iDimUp*iDimDw) + 1
                i_el = mod(li-1,iDimUp*iDimDw) + 1
                !
                call state2indices(i_el,[iDimUps,iDimDws],Indices)
                iud(1)   = HI(ialfa)%map(Indices(ialfa))
                iud(2)   = HI(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
                nud(1,:) = Bdecomp(iud(1),Ns_Orb)
                nud(2,:) = Bdecomp(iud(2),Ns_Orb)
                !
                if(Nud(ispin,iorb1)/=0)cycle
                call cdg(iorb1,iud(ispin),k,sgn_cdg)
                Jndices        = Indices
                Jndices(ibeta) = binary_search(HJ(ibeta)%map,k)
                call indices2state(Jndices,[jDimUps,jDimDws],rj)
                !
                rj = rj + (iph-1)*jDimUp*jDimDw
                !
                ! op_mat(1)=op_mat(1) + conjg(espace(jsector)%M(rj,j))*sgn_cdg*espace(isector)%M(li,i)
                op_mat(1)=op_mat(1) + (espace(jsector)%M(rj,j))*sgn_cdg*espace(isector)%M(li,i)
             enddo
             !
             do rj=1,jdim
                iph = (rj-1)/(jDimUp*jDimDw) + 1
                i_el = mod(rj-1,jDimUp*jDimDw) + 1
                !
                call state2indices(i_el,[jDimUps,jDimDws],Jndices)
                iud(1)   = HJ(ialfa)%map(Jndices(ialfa))
                iud(2)   = HJ(ialfa+Ns_Ud)%map(Jndices(ialfa+Ns_Ud))
                nud(1,:) = Bdecomp(iud(1),Ns_Orb)
                nud(2,:) = Bdecomp(iud(2),Ns_Orb)
                if(nud(ispin,iorb1)/=1)cycle
                call c(iorb1,iud(ispin),k,sgn_c)
                Indices        = Jndices
                Indices(ibeta) = binary_search(HI(ibeta)%map,k)
                call indices2state(Indices,[iDimUps,iDimDws],li)
                !
                li = li + (iph-1)*iDimUp*iDimDw
                !
                ! op_mat(2)=op_mat(2) + conjg(espace(isector)%M(li,i))*sgn_c*espace(jsector)%M(rj,j)
                op_mat(2)=op_mat(2) + (espace(isector)%M(li,i))*sgn_c*espace(jsector)%M(rj,j)
             enddo
             !
             Ei=espace(isector)%e(i)
             Ej=espace(jsector)%e(j)
             de=Ej-Ei
             peso=expterm/zeta_function
             spectral_weight=peso*product(op_mat)
             !
             do m=1,Lmats
                iw=xi*wm(m)
                impGmats(ispin,ispin,iorb,iorb,m)=impGmats(ispin,ispin,iorb,iorb,m)+spectral_weight/(iw-de)
             enddo
             !
             do m=1,Lreal
                w0=wr(m);iw=cmplx(w0,eps)
                impGreal(ispin,ispin,iorb,iorb,m)=impGreal(ispin,ispin,iorb,iorb,m)+spectral_weight/(iw-de)
             enddo
             !
          enddo
       enddo
       call delete_sector(isector,HI)
       call delete_sector(jsector,HJ)
    enddo
  end subroutine full_build_gf_normal_main





  subroutine full_build_gf_normal_mix_main(iorb,jorb,ispin)
    integer                     :: iorb,jorb,ispin
    integer,dimension(2*Ns_Ud)  :: Indices
    integer,dimension(2*Ns_Ud)  :: Jndices
    integer,dimension(Ns_Ud)    :: iDimUps,iDimDws
    integer,dimension(Ns_Ud)    :: jDimUps,jDimDws
    integer,dimension(2,Ns_Orb) :: Nud
    integer,dimension(2)        :: Iud,Jud
    type(sector_map)            :: HI(2*Ns_Ud),HJ(2*Ns_Ud)
    complex(8)                  :: op_mat(2)
    complex(8)                  :: spectral_weight
    real(8)                     :: sgn_cdg,sgn_c
    integer                     :: m,i,j,r,k,p,li,rj
    integer :: i_iDim,j_jDim
    real(8)                     :: Ei,Ej
    real(8)                     :: expterm,peso,de,w0
    complex(8)                  :: iw
    !
    if(ed_total_ud)then
       ialfa = 1
       jalfa = ialfa               !this is the condition to evaluate G_ab: ialfa=jalfa
       iorb1 = iorb
       jorb1 = jorb
    else
       ! ialfa = iorb; jalfa = jorb; iorb1 = 1; jorb1 = 1
       write(LOGfile,"(A)")"ED_GF_NORMAL warning: can not evaluate GF_ab with ed_total_ud=F"
       return
    endif
    ibeta  = ialfa + (ispin-1)*Ns_Ud
    jbeta  = jalfa + (ispin-1)*Ns_Ud
    !
    do isector=1,Nsectors
       jsector=getCDGsector(ialfa,ispin,isector)
       if(jsector==0)cycle
       !
       idim=getdim(isector)     !i-th sector dimension
       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       call build_sector(isector,HI)
       !
       jdim=getdim(jsector)     !j-th sector dimension
       call get_DimUp(jsector,jDimUps)
       call get_DImDw(jsector,jDimDws)
       jDimUp = product(jDimUps)
       jDimDw = product(jDimDws)
       call build_sector(jsector,HJ)
       !
       !
       do i=1,iDim          !loop over the states in the i-th sect.
          do j=1,jDim       !loop over the states in the j-th sect.
             !
             expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(jsector)%e(j))
             if(expterm < cutoff)cycle
             !
             op_mat=0d0
             !
             do li=1,idim              !loop over the component of |I> (IN state!)
                iph = (li-1)/(iDimUp*iDimDw) + 1
                i_el = mod(li-1,iDimUp*iDimDw) + 1
                !
                call state2indices(i_el,[iDimUps,iDimDws],Indices)
                iud(1)   = HI(ialfa)%map(Indices(ialfa))
                iud(2)   = HI(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
                nud(1,:) = Bdecomp(iud(1),Ns_Orb)
                nud(2,:) = Bdecomp(iud(2),Ns_Orb)
                !
                if(Nud(ispin,iorb1)/=0)cycle
                call cdg(iorb1,iud(ispin),k,sgn_cdg)
                Jndices        = Indices
                Jndices(ibeta) = binary_search(HJ(ibeta)%map,k)
                call indices2state(Jndices,[jDimUps,jDimDws],rj)
                !
                rj = rj + (iph-1)*jDimUp*jDimDw
                !
                ! op_mat(1)=op_mat(1) + conjg(espace(jsector)%M(rj,j))*sgn_cdg*espace(isector)%M(li,i)
                op_mat(1)=op_mat(1) + (espace(jsector)%M(rj,j))*sgn_cdg*espace(isector)%M(li,i)
             enddo
             !
             do rj=1,jdim
                iph = (rj-1)/(jDimUp*jDimDw) + 1
                i_el = mod(rj-1,jDimUp*jDimDw) + 1
                !
                call state2indices(i_el,[jDimUps,jDimDws],Jndices)
                iud(1)   = HJ(jalfa)%map(Jndices(jalfa))
                iud(2)   = HJ(jalfa+Ns_Ud)%map(Jndices(jalfa+Ns_Ud))
                nud(1,:) = Bdecomp(iud(1),Ns_Orb)
                nud(2,:) = Bdecomp(iud(2),Ns_Orb)
                if(nud(ispin,jorb1)/=1)cycle
                call c(jorb1,iud(ispin),k,sgn_c)
                Indices        = Jndices
                Indices(jbeta) = binary_search(HI(jbeta)%map,k)
                call indices2state(Indices,[iDimUps,iDimDws],li)
                !
                li = li + (iph-1)*iDimUp*iDimDw
                !
                ! op_mat(2)=op_mat(2) + conjg(espace(isector)%M(li,i))*sgn_c*espace(jsector)%M(rj,j)
                op_mat(2)=op_mat(2) + (espace(isector)%M(li,i))*sgn_c*espace(jsector)%M(rj,j)
             enddo
             !
             Ei=espace(isector)%e(i)
             Ej=espace(jsector)%e(j)
             de=Ej-Ei
             peso=expterm/zeta_function
             spectral_weight=peso*product(op_mat)
             !
             do m=1,Lmats
                iw=xi*wm(m)
                impGmats(ispin,ispin,iorb,jorb,m)=impGmats(ispin,ispin,iorb,jorb,m)+spectral_weight/(iw-de)
             enddo
             !
             do m=1,Lreal
                w0=wr(m);iw=cmplx(w0,eps)
                impGreal(ispin,ispin,iorb,jorb,m)=impGreal(ispin,ispin,iorb,jorb,m)+spectral_weight/(iw-de)
             enddo
             !
          enddo
       enddo
       call delete_sector(isector,HI)
       call delete_sector(jsector,HJ)
    enddo
  end subroutine full_build_gf_normal_mix_main




  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################






  subroutine build_sigma_normal
    integer                                           :: i,ispin,iorb
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: invG0mats,invGmats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: invG0real,invGreal
    complex(8),dimension(Norb,Norb)                   :: invGimp
    !
    invG0mats = zero
    invGmats  = zero
    invG0real = zero
    invGreal  = zero
    !
    !Get G0^-1
    invG0mats(:,:,:,:,:) = invg0_bath_function(dcmplx(0d0,wm(:)),dmft_bath)
    invG0real(:,:,:,:,:) = invg0_bath_function(dcmplx(wr(:),eps),dmft_bath)
    !
    select case(bath_type)
    case default                !Diagonal in both spin and orbital
       !
       !Get Gimp^-1
       do ispin=1,Nspin
          do iorb=1,Norb
             invGmats(ispin,ispin,iorb,iorb,:) = one/impGmats(ispin,ispin,iorb,iorb,:)
             invGreal(ispin,ispin,iorb,iorb,:) = one/impGreal(ispin,ispin,iorb,iorb,:)
          enddo
       enddo
       !Get Sigma functions: Sigma= G0^-1 - G^-1
       impSmats=zero
       impSreal=zero
       do ispin=1,Nspin
          do iorb=1,Norb
             impSmats(ispin,ispin,iorb,iorb,:) = invG0mats(ispin,ispin,iorb,iorb,:) - invGmats(ispin,ispin,iorb,iorb,:)
             impSreal(ispin,ispin,iorb,iorb,:) = invG0real(ispin,ispin,iorb,iorb,:) - invGreal(ispin,ispin,iorb,iorb,:)
          enddo
       enddo
       !
    case ("hybrid","replica")   !Diagonal in spin only. Full Orbital structure
       !
       !Get Gimp^-1
       do ispin=1,Nspin
          do i=1,Lmats
             invGimp = impGmats(ispin,ispin,:,:,i)
             call inv(invGimp)
             invGmats(ispin,ispin,:,:,i)=invGimp
          enddo
          !
          do i=1,Lreal
             invGimp = impGreal(ispin,ispin,:,:,i)
             call inv(invGimp)
             invGreal(ispin,ispin,:,:,i)=invGimp
          enddo
       enddo
       !Get Sigma functions: Sigma= G0^-1 - G^-1
       impSmats=zero
       impSreal=zero
       do ispin=1,Nspin
          impSmats(ispin,ispin,:,:,:) = invG0mats(ispin,ispin,:,:,:) - invGmats(ispin,ispin,:,:,:)
          !
          impSreal(ispin,ispin,:,:,:) = invG0real(ispin,ispin,:,:,:) - invGreal(ispin,ispin,:,:,:)
       enddo
       !
    end select
    !
    !Get G0and:
    impG0mats(:,:,:,:,:) = g0and_bath_function(dcmplx(0d0,wm(:)),dmft_bath)
    impG0real(:,:,:,:,:) = g0and_bath_function(dcmplx(wr(:),eps),dmft_bath)
    !!
    !
  end subroutine build_sigma_normal


END MODULE ED_GF_NORMAL











