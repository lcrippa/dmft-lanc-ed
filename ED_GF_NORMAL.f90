MODULE ED_GF_NORMAL
  USE ED_GF_SHARED
  implicit none
  private


  public :: build_gf_normal
  public :: build_sigma_normal


contains



  !+------------------------------------------------------------------+
  !                        NORMAL
  !+------------------------------------------------------------------+
  subroutine build_gf_normal()
    integer :: iorb,jorb,ispin,i
    logical :: MaskBool
    !
    !NORMAL: (default)
    do ispin=1,Nspin
       do iorb=1,Norb
          if(MPI_MASTER)write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_s"//str(ispin)
          if(MPI_MASTER)call start_timer
          call lanc_build_gf_normal_c(iorb,ispin)
          if(MPI_MASTER)call stop_timer(LOGfile)
       enddo
    enddo
    !
    !HYBRID:
    if(bath_type/="normal")then
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=iorb+1,Norb
                if(MPI_MASTER)write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_m"//str(jorb)//"_s"//str(ispin)
                MaskBool=.true.   
                if(bath_type=="replica")MaskBool=&
                     (dmft_bath%mask(ispin,ispin,iorb,jorb,1)).OR.(dmft_bath%mask(ispin,ispin,iorb,jorb,2))
                if(.not.MaskBool)cycle
                if(MPI_MASTER)call start_timer
                call lanc_build_gf_normal_mix_c(iorb,jorb,ispin)
                if(MPI_MASTER)call stop_timer(LOGfile)
             enddo
          enddo
       enddo
       !Put here off-diagonal manipulation by symmetry:
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=iorb+1,Norb
                !
                !if(hybrid)always T; if(replica)T iff following condition is T
                MaskBool=.true.   
                if(bath_type=="replica")MaskBool=&
                     (dmft_bath%mask(ispin,ispin,iorb,jorb,1)).OR.(dmft_bath%mask(ispin,ispin,iorb,jorb,2))
                !
                if(.not.MaskBool)cycle
                impGmats(ispin,ispin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,ispin,iorb,jorb,:) &
                     - (one-xi)*impGmats(ispin,ispin,iorb,iorb,:) - (one-xi)*impGmats(ispin,ispin,jorb,jorb,:))
                impGreal(ispin,ispin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,ispin,iorb,jorb,:) &
                     - (one-xi)*impGreal(ispin,ispin,iorb,iorb,:) - (one-xi)*impGreal(ispin,ispin,jorb,jorb,:))
                !>>ACTHUNG: this relation might not be true, it depends on the value of the impHloc_ij
                ! if impHloc_ij is REAL then it is true. if CMPLX hermiticity must be ensured
                impGmats(ispin,ispin,jorb,iorb,:) = impGmats(ispin,ispin,iorb,jorb,:)
                impGreal(ispin,ispin,jorb,iorb,:) = impGreal(ispin,ispin,iorb,jorb,:)
             enddo
          enddo
       enddo
    endif
    !
  end subroutine build_gf_normal






  !################################################################
  !################################################################
  !################################################################
  !################################################################






  subroutine lanc_build_gf_normal_c(iorb,ispin)
    complex(8),allocatable           :: vvinit(:)
    real(8),allocatable              :: alfa_(:),beta_(:)
    integer                          :: iorb,ispin,istate
    integer                          :: isector,jsector
    integer                          :: idim,idimUP,idimDW
    integer                          :: jdim,jdimUP,jdimDW
    integer                          :: nud(2,Ns),iud(2),jud(2)
    integer                          :: m,i,j,r
    integer                          :: iup,idw,jup,jdw,mup,mdw
    real(8)                          :: sgn,norm2,norm0
    complex(8)                       :: cnorm2
    integer                          :: Nitermax,Nlanc
    type(sector_map)                 :: HI(2),HJ(2)
    !
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
       idimUp  = getDimUp(isector)
       idimDw  = getDimDw(isector)
       call build_sector(isector,HI)
       !
       !
       !
       !ADD ONE PARTICLE:
       jsector = getCDGsector(ispin,isector)
       if(jsector/=0)then 
          if(ed_verbose==3)write(LOGfile,"(A,2I3)")' add particle:',getnup(jsector),getndw(jsector)
          !
          jdim   = getdim(jsector)
          jdimUp = getDimUp(jsector)
          jdimDw = getDimDw(jsector)
          allocate(vvinit(jdim)) ; vvinit=zero
          !
          call build_sector(jsector,HJ)
          do iup=1,idimUP
             iud(1)   = HI(1)%map(iup)
             nud(1,:) = bdecomp(iud(1),Ns)
             !
             do idw=1,idimDW
                iud(2)   = HI(2)%map(idw)
                nud(2,:) = bdecomp(iud(2),Ns)
                !
                i = iup + (idw-1)*idimUP
                !
                if(nud(ispin,iorb)/=0)cycle
                call cdg(iorb,iud(ispin),r,sgn)
                !
                jud = [iup,idw]
                jud(ispin) = binary_search(HJ(ispin)%map,r)
                !
                j = jud(1) + (jud(2)-1)*jdimUP
                vvinit(j) = sgn*state_cvec(i)
                !
             enddo
          enddo
          call delete_sector(jsector,HJ)
          !
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)
          !
          nlanc=min(jdim,lanc_nGFiter)
          allocate(alfa_(nlanc),beta_(nlanc))
          !
          call setup_Hv_sector(jsector)
          if(ed_sparse_H)call ed_buildH_c()
#ifdef _MPI
          if(MpiStatus)then
             call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
          else
             call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
          endif
#else
          call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
#endif
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,1,iorb,iorb,ispin)
          call delete_Hv_sector()
          !
          deallocate(vvinit,alfa_,beta_)
          if(spH0%status)call sp_delete_matrix(spH0)
       endif
       !
       !REMOVE ONE PARTICLE:
       jsector = getCsector(ispin,isector)
       if(jsector/=0)then
          if(ed_verbose==3)write(LOGfile,"(A,2I3)")' del particle:',getnup(jsector),getndw(jsector)
          !            
          jdim   = getdim(jsector)
          jdimUp = getDimUp(jsector)
          jdimDw = getDimDw(jsector)
          allocate(vvinit(jdim)) ; vvinit=zero
          !
          call build_sector(jsector,HJ)
          do iup=1,idimUP
             iud(1)   = HI(1)%map(iup)
             nud(1,:) = bdecomp(iud(1),Ns)
             !
             do idw=1,idimDW
                iud(2)   = HI(2)%map(idw)
                nud(2,:) = bdecomp(iud(2),Ns)
                !
                i = iup + (idw-1)*idimUP
                !
                if(nud(ispin,iorb)/=1)cycle
                call c(iorb,iud(ispin),r,sgn)
                !
                jud = [iup,idw]
                jud(ispin) = binary_search(HJ(ispin)%map,r)
                !
                j = jud(1) + (jud(2)-1)*jdimUP
                vvinit(j) = sgn*state_cvec(i)
                !
             enddo
          enddo
          call delete_sector(jsector,HJ)
          !
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)
          !
          nlanc=min(jdim,lanc_nGFiter)
          allocate(alfa_(nlanc),beta_(nlanc))
          !
          call setup_Hv_sector(jsector)
          if(ed_sparse_H)call ed_buildH_c()
#ifdef _MPI        
          if(MpiStatus)then
             call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
          else
             call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
          endif
#else
          call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
#endif
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,-1,iorb,iorb,ispin)
          !
          call delete_Hv_sector()
          !
          deallocate(vvinit,alfa_,beta_)
          if(spH0%status)call sp_delete_matrix(spH0)
       endif
       !
       !
       nullify(state_cvec)
       call delete_sector(isector,HI)
       !
    enddo
    return
  end subroutine lanc_build_gf_normal_c






  !################################################################
  !################################################################
  !################################################################
  !################################################################






  subroutine lanc_build_gf_normal_mix_c(iorb,jorb,ispin)
    integer                          :: iorb,jorb,ispin,istate
    integer                          :: isector,jsector
    integer                          :: idim,idimUP,idimDW
    integer                          :: jdim,jdimUP,jdimDW
    integer                          :: nud(2,Ns),iud(2),jud(2)
    integer                          :: iup,idw,jup,jdw,mup,mdw
    integer                          :: m,i,j,r,numstates
    real(8)                          :: sgn,norm2,norm0
    complex(8)                       :: cnorm2
    complex(8),allocatable           :: vvinit(:)
    real(8),allocatable              :: alfa_(:),beta_(:)
    integer                          :: Nitermax,Nlanc
    type(sector_map)                 :: HI(2),HJ(2)
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
       !
       idim  = getdim(isector)
       idimUp  = getDimUp(isector)
       idimDw  = getDimDw(isector)
       call build_sector(isector,HI)
       !
       !
       !
       !EVALUATE (c^+_iorb + c^+_jorb)|gs>
       jsector = getCDGsector(ispin,isector)
       if(jsector/=0)then
          if(ed_verbose==3)write(LOGfile,"(A,2I3,I15)")' add particle:',getnup(jsector),getndw(jsector),jdim
          !
          jdim   = getdim(jsector)
          jdimUp = getDimUp(jsector)
          jdimDw = getDimDw(jsector)
          allocate(vvinit(jdim)) ; vvinit=zero
          !
          call build_sector(jsector,HJ)
          do iup=1,idimUP
             iud(1)   = HI(1)%map(iup)
             nud(1,:) = bdecomp(iud(1),Ns)
             !
             do idw=1,idimDW
                iud(2)   = HI(2)%map(idw)
                nud(2,:) = bdecomp(iud(2),Ns)
                !
                i = iup + (idw-1)*idimUP
                !
                !
                !
                if(nud(ispin,iorb)/=0)cycle
                call cdg(iorb,iud(ispin),r,sgn)
                !
                jud = [iup,idw]
                jud(ispin) = binary_search(HJ(ispin)%map,r)
                !
                j = jud(1) + (jud(2)-1)*jdimUP
                vvinit(j) = sgn*state_cvec(i)
             enddo
          enddo
          do iup=1,idimUP
             iud(1)   = HI(1)%map(iup)
             nud(1,:) = bdecomp(iud(1),Ns)
             !
             do idw=1,idimDW
                iud(2)   = HI(2)%map(idw)
                nud(2,:) = bdecomp(iud(2),Ns)
                !
                i = iup + (idw-1)*idimUP
                !
                !
                if(nud(ispin,jorb)/=0)cycle
                call cdg(jorb,iud(ispin),r,sgn)
                !
                jud = [iup,idw]
                jud(ispin) = binary_search(HJ(ispin)%map,r)
                !
                j = jud(1) + (jud(2)-1)*jdimUP
                vvinit(j) = vvinit(j) + sgn*state_cvec(i)
             enddo
          enddo
          call delete_sector(jsector,HJ)
          !
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)
          !
          nlanc=min(jdim,lanc_nGFiter)
          allocate(alfa_(nlanc),beta_(nlanc))
          !           
          call setup_Hv_sector(jsector)
          if(ed_sparse_H)call ed_buildH_c()
#ifdef _MPI
          if(MpiStatus)then
             call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
          else
             call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
          endif
#else
          call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
#endif
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,1,iorb,jorb,ispin)
          !
          call delete_Hv_sector()
          !
          deallocate(vvinit,alfa_,beta_)
          if(spH0%status)call sp_delete_matrix(spH0)
       endif
       !
       !EVALUATE (c_iorb + c_jorb)|gs>
       jsector = getCsector(ispin,isector)
       if(jsector/=0)then
          if(ed_verbose==3)write(LOGfile,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
          !
          jdim   = getdim(jsector)
          jdimUp = getDimUp(jsector)
          jdimDw = getDimDw(jsector)
          allocate(vvinit(jdim)) ; vvinit=zero
          !
          call build_sector(jsector,HJ)
          do iup=1,idimUP
             iud(1)   = HI(1)%map(iup)
             nud(1,:) = bdecomp(iud(1),Ns)
             !
             do idw=1,idimDW
                iud(2)   = HI(2)%map(idw)
                nud(2,:) = bdecomp(iud(2),Ns)
                !
                i = iup + (idw-1)*idimUP
                !
                !
                if(nud(ispin,iorb)/=1)cycle
                call c(iorb,iud(ispin),r,sgn)
                !
                jud = [iup,idw]
                jud(ispin) = binary_search(HJ(ispin)%map,r)
                !
                j = jud(1) + (jud(2)-1)*jdimUP
                vvinit(j) = sgn*state_cvec(i)
                !
             enddo
          enddo
          do iup=1,idimUP
             iud(1)   = HI(1)%map(iup)
             nud(1,:) = bdecomp(iud(1),Ns)
             !
             do idw=1,idimDW
                iud(2)   = HI(2)%map(idw)
                nud(2,:) = bdecomp(iud(2),Ns)
                !
                i = iup + (idw-1)*idimUP
                !
                !
                if(nud(ispin,jorb)/=1)cycle
                call c(jorb,iud(ispin),r,sgn)
                !
                jud = [iup,idw]
                jud(ispin) = binary_search(HJ(ispin)%map,r)
                !
                j = jud(1) + (jud(2)-1)*jdimUP
                vvinit(j) = vvinit(j) + sgn*state_cvec(i)
                !
             enddo
          enddo
          call delete_sector(jsector,HJ)
          !
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)
          !
          nlanc=min(jdim,lanc_nGFiter)
          allocate(alfa_(nlanc),beta_(nlanc))
          !
          call setup_Hv_sector(jsector)
          if(ed_sparse_H)call ed_buildH_c()
#ifdef _MPI
          if(MpiStatus)then
             call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
          else
             call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
          endif
#else
          call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
#endif
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin)
          !
          call delete_Hv_sector()
          !
          deallocate(vvinit,alfa_,beta_)
          if(spH0%status)call sp_delete_matrix(spH0)
       endif
       !
       !
       !
       !EVALUATE (c^+_iorb + i*c^+_jorb)|gs>
       jsector = getCDGsector(ispin,isector)
       if(jsector/=0)then 
          if(ed_verbose==3)write(LOGfile,"(A,2I3,I15)")' add particle:',getnup(jsector),getndw(jsector),jdim
          jdim   = getdim(jsector)
          jdimUp = getDimUp(jsector)
          jdimDw = getDimDw(jsector)
          allocate(vvinit(jdim)) ; vvinit=zero
          !
          call build_sector(jsector,HJ)
          do iup=1,idimUP
             iud(1)   = HI(1)%map(iup)
             nud(1,:) = bdecomp(iud(1),Ns)
             !
             do idw=1,idimDW
                iud(2)   = HI(2)%map(idw)
                nud(2,:) = bdecomp(iud(2),Ns)
                !
                i = iup + (idw-1)*idimUP
                !
                !
                if(nud(ispin,iorb)/=0)cycle
                call cdg(iorb,iud(ispin),r,sgn)
                !
                jud = [iup,idw]
                jud(ispin) = binary_search(HJ(ispin)%map,r)
                !
                j = jud(1) + (jud(2)-1)*jdimUP
                vvinit(j) = sgn*state_cvec(i)
                !
             enddo
          enddo
          do iup=1,idimUP
             iud(1)   = HI(1)%map(iup)
             nud(1,:) = bdecomp(iud(1),Ns)
             !
             do idw=1,idimDW
                iud(2)   = HI(2)%map(idw)
                nud(2,:) = bdecomp(iud(2),Ns)
                !
                i = iup + (idw-1)*idimUP
                !
                !
                if(nud(ispin,jorb)/=0)cycle
                call cdg(jorb,iud(ispin),r,sgn)
                !
                jud = [iup,idw]
                jud(ispin) = binary_search(HJ(ispin)%map,r)
                !
                j = jud(1) + (jud(2)-1)*jdimUP
                vvinit(j) = vvinit(j) + xi*sgn*state_cvec(i)
                !
             enddo
          enddo
          call delete_sector(jsector,HJ)
          !
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)
          !
          nlanc=min(jdim,lanc_nGFiter)
          allocate(alfa_(nlanc),beta_(nlanc))
          !
          call setup_Hv_sector(jsector)
          if(ed_sparse_H)call ed_buildH_c()
#ifdef _MPI
          if(MpiStatus)then
             call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
          else
             call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
          endif
#else
          call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
#endif
          call add_to_lanczos_gf_normal(-xi*norm2,state_e,alfa_,beta_,1,iorb,jorb,ispin)
          !
          call delete_Hv_sector()
          !
          deallocate(vvinit,alfa_,beta_)
          if(spH0%status)call sp_delete_matrix(spH0)
       endif
       !
       !EVALUATE (c_iorb - xi*c_jorb)|gs>
       jsector = getCsector(ispin,isector)
       if(jsector/=0)then
          if(ed_verbose==3)write(LOGfile,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
          !
          jdim   = getdim(jsector)
          jdimUp = getDimUp(jsector)
          jdimDw = getDimDw(jsector)
          allocate(vvinit(jdim)); vvinit=zero
          !
          call build_sector(jsector,HJ)
          do iup=1,idimUP
             iud(1)   = HI(1)%map(iup)
             nud(1,:) = bdecomp(iud(1),Ns)
             !
             do idw=1,idimDW
                iud(2)   = HI(2)%map(idw)
                nud(2,:) = bdecomp(iud(2),Ns)
                !
                i = iup + (idw-1)*idimUP
                !
                !
                if(nud(ispin,iorb)/=1)cycle
                call c(iorb,iud(ispin),r,sgn)
                !
                jud = [iup,idw]
                jud(ispin) = binary_search(HJ(ispin)%map,r)
                !
                j = jud(1) + (jud(2)-1)*jdimUP
                vvinit(j) = sgn*state_cvec(i)
                !
             enddo
          enddo
          do iup=1,idimUP
             iud(1)   = HI(1)%map(iup)
             nud(1,:) = bdecomp(iud(1),Ns)
             !
             do idw=1,idimDW
                iud(2)   = HI(2)%map(idw)
                nud(2,:) = bdecomp(iud(2),Ns)
                !
                i = iup + (idw-1)*idimUP
                !
                !
                if(nud(ispin,jorb)/=1)cycle
                call c(jorb,iud(ispin),r,sgn)
                !
                jud = [iup,idw]
                jud(ispin) = binary_search(HJ(ispin)%map,r)
                !
                j = jud(1) + (jud(2)-1)*jdimUP
                vvinit(j) = vvinit(j) - xi*sgn*state_cvec(i)
                !
             enddo
          enddo
          call delete_sector(jsector,HJ)
          !
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)
          !
          nlanc=min(jdim,lanc_nGFiter)
          allocate(alfa_(nlanc),beta_(nlanc))
          !
          call setup_Hv_sector(jsector)
          if(ed_sparse_H)call ed_buildH_c()
#ifdef _MPI
          if(MpiStatus)then
             call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
          else
             call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
          endif
#else
          call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
#endif
          call add_to_lanczos_gf_normal(-xi*norm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin)
          !
          call delete_Hv_sector()
          !
          deallocate(vvinit,alfa_,beta_)
          if(spH0%status)call sp_delete_matrix(spH0)
       endif
       !
       !
       nullify(state_cvec)
       call delete_sector(isector,HI)
       !
    enddo
    return
  end subroutine lanc_build_gf_normal_mix_c





  !################################################################
  !################################################################
  !################################################################
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
    diag             = 0.d0
    subdiag          = 0.d0
    Z                = eye(Nlanc)
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call tql2(Nlanc,diag,subdiag,Z,ierr)
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





  !###################################################################################################
  !###################################################################################################
  !###################################################################################################
  !###################################################################################################
  !###################################################################################################
  !###################################################################################################






  subroutine build_sigma_normal
    integer                                           :: i,ispin,iorb
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: invG0mats,invGmats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: invG0real,invGreal
    complex(8),dimension(Norb,Norb)                   :: invGimp
    !
    ! if(.not.allocated(wm))allocate(wm(Lmats))
    ! if(.not.allocated(wr))allocate(wr(Lreal))
    ! wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    ! wr     = linspace(wini,wfin,Lreal)
    !
    invG0mats = zero
    invGmats  = zero
    invG0real = zero
    invGreal  = zero
    !
    !Get G0^-1
    invG0mats(:,:,:,:,:) = invg0_bath_mats(dcmplx(0d0,wm(:)),dmft_bath)
    invG0real(:,:,:,:,:) = invg0_bath_real(dcmplx(wr(:),eps),dmft_bath)
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
    impG0mats(:,:,:,:,:) = g0and_bath_mats(dcmplx(0d0,wm(:)),dmft_bath)
    impG0real(:,:,:,:,:) = g0and_bath_real(dcmplx(wr(:),eps),dmft_bath)
    !!
    !
    ! if(allocated(wm))deallocate(wm)
    ! if(allocated(wr))deallocate(wr)
    !
  end subroutine build_sigma_normal


END MODULE ED_GF_NORMAL
