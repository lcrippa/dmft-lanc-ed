!###################################################################
!PURPOSE  : Build the impurity Green's function using spectral sum 
!NOTE: in the MPI implementation we may require all the nodes to 
!evaluate the GF, this is safer, simpler and works for both Lanc &
!Ed. For Lanc we can indeed assign the contribution from each state 
!to different node and accumulate the result at the end.
!AUTHORS  : Adriano Amaricci
!###################################################################
MODULE ED_GREENS_FUNCTIONS
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,free_unit,reg,free_units,txtfy,splot
  USE SF_ARRAYS,  only: arange,linspace
  USE SF_LINALG,  only: inv,inv_sym,inv_her,eye
  USE SF_SP_LINALG, only: sp_lanc_tridiag
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_IO                     !< this contains the routine to print GF,Sigma and G0
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



  public :: buildGf_impurity

  public :: ed_greens_functions_set_MPI

  public :: ed_greens_functions_del_MPI



  !Lanczos shared variables
  !=========================================================
  real(8),dimension(:),pointer                :: state_vec
  complex(8),dimension(:),pointer             :: state_cvec
  real(8)                                     :: state_e

  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable            :: wm,tau,wr,vm

  !Auxiliary functions GF
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:,:) :: impDeltamats,impDeltareal
  complex(8),allocatable,dimension(:,:,:,:,:) :: invimpG0mats,invimpG0real
  complex(8),allocatable,dimension(:,:,:,:,:) :: invimpGmats,invimpGreal

  !AUX GF
  !=========================================================
  complex(8),allocatable,dimension(:,:)       :: auxGmats,auxGreal


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


  subroutine ed_greens_functions_set_MPI(comm)
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
  end subroutine ed_greens_functions_set_MPI


  subroutine ed_greens_functions_del_MPI()
#ifdef _MPI
    MpiComm  = MPI_UNDEFINED
    MpiStatus = .false.
    MPI_RANK=0
    MPI_SIZE=1
    MPI_MASTER=.true.
#endif
  end subroutine ed_greens_functions_del_MPI






  subroutine buildgf_impurity()
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(wr))allocate(wr(Lreal))
    wm     = pi/beta*dble(2*arange(1,Lmats)-1)
    wr     = linspace(wini,wfin,Lreal)
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
    ! GFpoles=zero
    ! GFweights=zero
    !
    if(.not.allocated(impDeltamats)) allocate(impDeltamats(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(invimpG0mats)) allocate(invimpG0mats(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(invimpGmats))  allocate( invimpGmats(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(impDeltareal)) allocate(impDeltareal(Nspin,Nspin,Norb,Norb,Lreal))
    if(.not.allocated(invimpG0real)) allocate(invimpG0real(Nspin,Nspin,Norb,Norb,Lreal))
    if(.not.allocated(invimpGreal))  allocate( invimpGreal(Nspin,Nspin,Norb,Norb,Lreal))
    impDeltamats=zero
    invimpGmats=zero
    invimpG0mats=zero
    impDeltareal=zero
    invimpGreal=zero
    invimpG0real=zero
    !
    write(LOGfile,"(A)")"Get impurity Greens functions:"
    call build_gf_normal()
    call build_sigma_normal()
    !
    if(MPI_MASTER)then
       if(ed_print_Sigma)call ed_print_impSigma()
       if(ed_print_G)call ed_print_impG()
       if(ed_print_G0)call ed_print_impG0()
    endif
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    if(allocated(invimpG0mats))deallocate(invimpG0mats)
    if(allocated(invimpGmats))deallocate(invimpGmats)
    if(allocated(impDeltamats))deallocate(impDeltamats)
    if(allocated(invimpG0real))deallocate(invimpG0real)
    if(allocated(invimpGreal))deallocate(invimpGreal)
    if(allocated(impDeltareal))deallocate(impDeltareal)
  end subroutine buildgf_impurity





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
          write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_s"//str(ispin)
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
                write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_m"//str(jorb)//"_s"//str(ispin)
                !if(hybrid)always T; if(replica)T iff following condition is T
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
    if(ed_para.and.Nspin==2)then
       impGmats(1,1,:,:,:) = (impGmats(1,1,:,:,:)+impGmats(2,2,:,:,:))/2.d0
       impGmats(2,2,:,:,:) =  impGmats(1,1,:,:,:)
       impGreal(1,1,:,:,:) = (impGreal(1,1,:,:,:)+impGreal(2,2,:,:,:))/2.d0
       impGreal(2,2,:,:,:) =  impGreal(1,1,:,:,:)
    endif
    !
  end subroutine build_gf_normal


  subroutine lanc_build_gf_normal_c(iorb,ispin)
    complex(8),allocatable           :: vvinit(:)
    real(8),allocatable              :: alfa_(:),beta_(:)
    integer                          :: iorb,ispin,isite,isector,istate
    integer                          :: idim,jsector
    integer                          :: jdim
    integer                          :: ib(Nlevels)
    integer                          :: m,i,j,r,numstates
    real(8)                          :: sgn,norm2,norm0
    complex(8)                       :: cnorm2
    integer                          :: Nitermax,Nlanc
    type(sector_map)                 :: HI,HJ
    !
    isite=impIndex(iorb,ispin)
    !

    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
       state_cvec => es_return_cvector(state_list,istate)
       norm0=sqrt(dot_product(state_cvec,state_cvec))
       if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
       idim  = getdim(isector)
       call build_sector(isector,HI)
       !
       !ADD ONE PARTICLE:
       jsector = getCDGsector(ispin,isector)
       if(jsector/=0)then 
          if(ed_verbose==3)write(LOGfile,"(A,2I3)")' add particle:',getnup(jsector),getndw(jsector)
          jdim  = getdim(jsector)
          allocate(vvinit(jdim))
          call build_sector(jsector,HJ)
          vvinit=zero
          do m=1,idim
             i=HI%map(m)
             ib = bdecomp(i,2*Ns)
             if(ib(isite)==0)then
                call cdg(isite,i,r,sgn)
                j=binary_search(HJ%map,r)
                vvinit(j) = sgn*state_cvec(m)
             endif
          enddo
          deallocate(HJ%map)
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)
          !
          call setup_Hv_sector(jsector)
          if(ed_sparse_H)call ed_buildH_c()
          !
          nlanc=min(jdim,lanc_nGFiter)
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
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,1,iorb,iorb,ispin)
          !
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
          jdim  = getdim(jsector)        
          allocate(vvinit(jdim))
          call build_sector(jsector,HJ)
          vvinit=zero
          do m=1,idim
             i=HI%map(m)
             ib = bdecomp(i,2*Ns)
             if(ib(isite)==1)then
                call c(isite,i,r,sgn)
                j=binary_search(HJ%map,r)
                vvinit(j) = sgn*state_cvec(m)
             endif
          enddo
          deallocate(HJ%map)
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)
          !
          call setup_Hv_sector(jsector)
          if(ed_sparse_H)call ed_buildH_c()
          !
          nlanc=min(jdim,lanc_nGFiter)
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
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,-1,iorb,iorb,ispin)
          !
          call delete_Hv_sector()
          !
          deallocate(vvinit,alfa_,beta_)
          if(spH0%status)call sp_delete_matrix(spH0)
       endif
       !
       nullify(state_cvec)
       deallocate(HI%map)
       !
    enddo

  end subroutine lanc_build_gf_normal_c

  subroutine lanc_build_gf_normal_mix_c(iorb,jorb,ispin)
    integer                          :: iorb,jorb,ispin,isite,jsite,isector,istate
    integer                          :: idim,jsector
    integer                          :: jdim
    integer                          :: ib(Nlevels)
    integer                          :: m,i,j,r,numstates
    real(8)                          :: sgn,norm2,norm0
    complex(8)                       :: cnorm2
    complex(8),allocatable           :: vvinit(:)
    complex(8),allocatable           :: cvinit(:)
    real(8),allocatable              :: alfa_(:),beta_(:)
    integer                          :: Nitermax,Nlanc
    type(sector_map)                 :: HI,HJ
    !
    isite=impIndex(iorb,ispin)  !orbital 1
    jsite=impIndex(jorb,ispin)  !orbital 2
    !

    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
       state_cvec => es_return_cvector(state_list,istate)
       norm0=sqrt(dot_product(state_cvec,state_cvec))
       if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
       !
       idim  = getdim(isector)
       call build_sector(isector,HI)
       !
       !EVALUATE (c^+_iorb + c^+_jorb)|gs>
       jsector = getCDGsector(ispin,isector)
       if(jsector/=0)then 
          jdim  = getdim(jsector)
          if(ed_verbose==3)write(LOGfile,"(A,2I3,I15)")' add particle:',getnup(jsector),getndw(jsector),jdim
          allocate(vvinit(jdim))
          call build_sector(jsector,HJ)
          vvinit=zero
          do m=1,idim
             i=HI%map(m)
             ib = bdecomp(i,2*Ns)
             if(ib(isite)==0)then
                call cdg(isite,i,r,sgn)
                j=binary_search(HJ%map,r)
                vvinit(j) = sgn*state_cvec(m)
             endif
          enddo
          do m=1,idim
             i=HI%map(m)
             ib = bdecomp(i,2*Ns)
             if(ib(jsite)==0)then
                call cdg(jsite,i,r,sgn)
                j=binary_search(HJ%map,r)
                vvinit(j) = vvinit(j) + sgn*state_cvec(m)
             endif
          enddo
          deallocate(HJ%map)
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)
          !
          call setup_Hv_sector(jsector)
          if(ed_sparse_H)call ed_buildH_c()
          !
          nlanc=min(jdim,lanc_nGFiter)
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
          jdim   = getdim(jsector)
          if(ed_verbose==3)write(LOGfile,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
          allocate(vvinit(jdim))
          call build_sector(jsector,HJ)
          vvinit=zero
          do m=1,idim
             i=HI%map(m)
             ib = bdecomp(i,2*Ns)
             if(ib(isite)==1)then
                call c(isite,i,r,sgn)
                j=binary_search(HJ%map,r)
                vvinit(j) = sgn*state_cvec(m)
             endif
          enddo
          do m=1,idim
             i=HI%map(m)
             ib = bdecomp(i,2*Ns)
             if(ib(jsite)==1)then
                call c(jsite,i,r,sgn)
                j=binary_search(HJ%map,r)
                vvinit(j) = vvinit(j) + sgn*state_cvec(m)
             endif
          enddo
          deallocate(HJ%map)
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)
          !
          call setup_Hv_sector(jsector)
          if(ed_sparse_H)call ed_buildH_c()
          !
          nlanc=min(jdim,lanc_nGFiter)
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
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin)
          !
          call delete_Hv_sector()
          !
          deallocate(vvinit,alfa_,beta_)
          if(spH0%status)call sp_delete_matrix(spH0)
       endif
       !
       !EVALUATE (c^+_iorb + i*c^+_jorb)|gs>
       jsector = getCDGsector(ispin,isector)
       if(jsector/=0)then 
          jdim  = getdim(jsector)
          if(ed_verbose==3)write(LOGfile,"(A,2I3,I15)")' add particle:',getnup(jsector),getndw(jsector),jdim
          allocate(cvinit(jdim))
          call build_sector(jsector,HJ)
          cvinit=zero
          do m=1,idim
             i=HI%map(m)
             ib = bdecomp(i,2*Ns)
             if(ib(isite)==0)then
                call cdg(isite,i,r,sgn)
                j=binary_search(HJ%map,r)
                cvinit(j) = sgn*state_cvec(m)
             endif
          enddo
          do m=1,idim
             i=HI%map(m)
             ib = bdecomp(i,2*Ns)
             if(ib(jsite)==0)then
                call cdg(jsite,i,r,sgn)
                j=binary_search(HJ%map,r)
                cvinit(j) = cvinit(j) + xi*sgn*state_cvec(m)
             endif
          enddo
          deallocate(HJ%map)
          norm2=dot_product(cvinit,cvinit)
          cvinit=cvinit/sqrt(norm2)
          !
          call setup_Hv_sector(jsector)
          if(ed_sparse_H)call ed_buildH_c()
          !
          nlanc=min(jdim,lanc_nGFiter)
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
          call add_to_lanczos_gf_normal(-xi*norm2,state_e,alfa_,beta_,1,iorb,jorb,ispin)
          !
          call delete_Hv_sector()
          !
          deallocate(cvinit,alfa_,beta_)
          if(spH0%status)call sp_delete_matrix(spH0)
       endif
       !
       !EVALUATE (c_iorb - xi*c_jorb)|gs>
       jsector = getCsector(ispin,isector)
       if(jsector/=0)then
          jdim   = getdim(jsector)
          if(ed_verbose==3)write(LOGfile,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
          allocate(cvinit(jdim))
          call build_sector(jsector,HJ)
          cvinit=zero
          do m=1,idim
             i=HI%map(m)
             ib = bdecomp(i,2*Ns)
             if(ib(isite)==1)then
                call c(isite,i,r,sgn)
                j=binary_search(HJ%map,r)
                cvinit(j) = sgn*state_cvec(m)
             endif
          enddo
          do m=1,idim
             i=HI%map(m)
             ib = bdecomp(i,2*Ns)
             if(ib(jsite)==1)then
                call c(jsite,i,r,sgn)
                j=binary_search(HJ%map,r)
                cvinit(j) = cvinit(j) - xi*sgn*state_cvec(m)
             endif
          enddo
          deallocate(HJ%map)
          norm2=dot_product(cvinit,cvinit)
          cvinit=cvinit/sqrt(norm2)
          !
          call setup_Hv_sector(jsector)
          if(ed_sparse_H)call ed_buildH_c()
          !
          nlanc=min(jdim,lanc_nGFiter)
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
          call add_to_lanczos_gf_normal(-xi*norm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin)
          !
          call delete_Hv_sector()
          !
          deallocate(cvinit,alfa_,beta_)
          if(spH0%status)call sp_delete_matrix(spH0)
       endif
       !
       nullify(state_cvec)
       deallocate(HI%map)
       !
    enddo
    !

    !
  end subroutine lanc_build_gf_normal_mix_c

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






  !+------------------------------------------------------------------+
  !                    SELF-ENERGY FUNCTIONS 
  !+------------------------------------------------------------------+
  subroutine build_sigma_normal
    integer                                           :: i,ispin,iorb
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: invG0mats,invGmats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: invG0real,invGreal
    complex(8),dimension(Norb,Norb)                   :: invGimp
    !
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(wr))allocate(wr(Lreal))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr     = linspace(wini,wfin,Lreal)
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
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    !
  end subroutine build_sigma_normal




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



end MODULE ED_GREENS_FUNCTIONS
