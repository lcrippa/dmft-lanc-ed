!########################################################################
!PURPOSE  : Build the impurity Hamiltonian
!|ImpUP,(2ImpUP),BathUP;,ImpDW,(2ImpDW),BathDW >
! |1,2;3...Ns>_UP * |Ns+1,Ns+2;Ns+3,...,2*Ns>_DOWN
!########################################################################
MODULE ED_HAMILTONIAN_MATVEC
  USE SF_MISC,    only: assert_shape
  USE SF_LINALG,  only: kronecker_product,eye
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_SETUP
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none
  private


  !>Build sparse hamiltonian of the sector
  public  :: vecDim_Hv_sector
  public  :: build_Hv_sector
  public  :: delete_Hv_sector

  !
  !>Sparse Mat-Vec product using stored sparse matrix
  !LL format
  public  :: spMatVec_main
#ifdef _MPI
  public  :: spMatVec_MPI_main
#endif
  !
  !>Sparse Mat-Vec direct on-the-fly product 
  public  :: directMatVec_main
#ifdef _MPI
  public  :: directMatVec_MPI_main
#endif
  !
  !>Related auxiliary routines:
  public  :: ed_hamiltonian_matvec_set_MPI
  public  :: ed_hamiltonian_matvec_del_MPI


  !> MPI local variables (shared)
#ifdef _MPI
  integer                                   :: MpiComm=MPI_UNDEFINED
#else
  integer                                   :: MpiComm=0
#endif
  logical                                   :: MpiStatus=.false.
  logical                                   :: MpiMaster=.true.
  integer                                   :: MpiIerr
  integer                                   :: MpiRank=0
  integer                                   :: MpiSize=1
  integer                                   :: MpiQ=1
  integer                                   :: MpiQup=1
  integer                                   :: MpiQdw=1
  integer                                   :: MpiR=0
  integer                                   :: MpiRup=0
  integer                                   :: MpiRdw=0
  integer                                   :: MpiIstart
  integer                                   :: MpiIend
  integer                                   :: MpiIshift
  !
  integer                                   :: Dim
  integer                                   :: DimUp
  integer                                   :: DimDw
  integer,allocatable,dimension(:)          :: DimUps
  integer,allocatable,dimension(:)          :: DimDws
  !
  integer                                   :: Hsector=0
  logical                                   :: Hstatus=.false.
  type(sector_map),dimension(:),allocatable :: Hs



contains




  !####################################################################
  !                      AUXILIARY MPI ROUTINES
  !####################################################################
  subroutine ed_hamiltonian_matvec_set_MPI(comm_)
#ifdef _MPI
    integer :: comm_
    MpiComm = comm_
    MpiStatus=.true.
    MpiRank = get_Rank_MPI(MpiComm)
    MpiSize = get_Size_MPI(MpiComm)
    MpiMaster= get_Master_MPI(MpiComm)
#else
    integer,optional :: comm_
#endif
  end subroutine ed_hamiltonian_matvec_set_MPI


  subroutine ed_hamiltonian_matvec_del_MPI()
#ifdef _MPI
    MpiComm = MPI_UNDEFINED
#else
    MpiComm = 0
#endif
    MpiStatus=.false.
    MpiRank=0
    MpiSize=1
  end subroutine ed_hamiltonian_matvec_del_MPI






  !####################################################################
  !                 MAIN ROUTINES: BUILD/DELETE SECTOR
  !####################################################################
  function vecDim_Hv_sector(isector) result(vecDim)
    integer :: isector
    integer :: vecDim
    integer :: mpiQdw
    integer :: DimUps(Ns_Ud),DimUp
    integer :: DimDws(Ns_Ud),DimDw
    !
    call get_DimUp(isector,DimUps) ; DimUp = product(DimUps)
    call get_DimDw(isector,DimDws) ; DimDw = product(DimDws)
    !
#ifdef _MPI
    if(MpiStatus)then
       !Dw split:
       mpiQdw = DimDw/MpiSize
       if(MpiRank < mod(DimDw,MpiSize) ) MpiQdw = MpiQdw+1
    else
       mpiQdw = DimDw
    endif
#else
    mpiQdw = DimDw
#endif
    !
    vecDim=DimUp*mpiQdw
    !
  end function vecDim_Hv_sector






  subroutine build_Hv_sector(isector,Hmat)
    integer                         :: isector,SectorDim
    real(8),dimension(:,:),optional :: Hmat   
    integer                         :: irank
    integer                         :: i,iup,idw
    integer                         :: j,jup,jdw
    !
    Hsector=isector
    Hstatus=.true.
    !
    allocate(Hs(2*Ns_Ud))
    allocate(DimUps(Ns_Ud))
    allocate(DimDws(Ns_Ud))
    call build_sector(isector,Hs)
    !
    call get_DimUp(isector,DimUps)
    call get_DimDw(isector,DimDws)
    DimUp = product(DimUps)
    DimDw = product(DimDws)
    Dim = getDim(isector)
    !
    !Dw split:
    mpiQdw = DimDw/MpiSize
    mpiRdw = mod(DimDw,MpiSize)
    if(MpiRank < mod(DimDw,MpiSize) ) then
       mpiRdw = 0
       MpiQdw = MpiQdw+1
    endif
    !Total split: split DW \times UP
    mpiQ = DimUp*mpiQdw
    mpiR = DimUp*mpiRdw
    mpiIstart = 1 + MpiRank*mpiQ+mpiR
    mpiIend   = (MpiRank+1)*mpiQ+mpiR
    mpiIshift = MpiRank*mpiQ+mpiR
#ifdef _MPI    
    if(MpiStatus.AND.ed_verbose>4)then
       if(MpiMaster)write(LOGfile,*)&
            "         mpiRank,   mpi_Q,   mpi_R,   mpi_Istart,   mpi_Iend,   mpi_Iend-mpi_Istart"
       do irank=0,MpiSize-1
          call Barrier_MPI(MpiComm)
          if(MpiRank==irank)then
             write(LOGfile,*)MpiRank,MpiQ,MpiR,MpiIstart,MpiIend,MpiIstart-MpiIend+1
          endif
       enddo
       call Barrier_MPI(MpiComm)
    endif
#endif
    !
    !
    if(present(Hmat))then
       call ed_buildh_main(isector,Hmat)
       return
    endif
    !
    select case (ed_sparse_H)
    case (.true.)
       call ed_buildh_main(isector)
    case (.false.)
       !nothing to be done: direct matrix-vector product
    end select
    !
  end subroutine build_Hv_sector


  subroutine delete_Hv_sector()
    integer :: iud
    call delete_sector(Hsector,Hs)
    deallocate(Hs)
    deallocate(DimUps)
    deallocate(DimDws)    
    Hsector=0
    Hstatus=.false.
    !
    !There is no difference here between Mpi and serial version, as Local part was removed.
#ifdef _MPI
    if(MpiStatus)then
       call sp_delete_matrix(MpiComm,spH0d)
    else
       call sp_delete_matrix(spH0d)
    endif
#else
    call sp_delete_matrix(spH0d)
#endif
    call sp_delete_matrix(spH0up)
    call sp_delete_matrix(spH0dw)
    !
  end subroutine delete_Hv_sector








  !####################################################################
  !             BUILD SPARSE HAMILTONIAN of the SECTOR
  !####################################################################
  subroutine ed_buildh_main(isector,Hmat)
    integer                             :: isector   
    real(8),dimension(:,:),optional     :: Hmat
    real(8),dimension(:,:),allocatable  :: Hredux,Htmp_up,Htmp_dw
    integer,dimension(Ns)               :: nup,ndw
    integer                             :: i,iup,idw
    integer                             :: j,jup,jdw
    integer                             :: m,mup,mdw
    integer                             :: ms
    integer                             :: impi,impi_up,impi_dw
    integer                             :: iorb,jorb,ispin,jspin,ibath
    integer                             :: kp,k1,k2,k3,k4
    integer                             :: alfa,beta
    real(8)                             :: sg1,sg2,sg3,sg4
    real(8)                             :: htmp,htmpup,htmpdw
    real(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    logical                             :: Jcondition,dryrun_
    !
    if(.not.Hstatus)stop "ed_buildh_main ERROR: Hsector NOT set"
    isector=Hsector
    !
    if(present(Hmat))call assert_shape(Hmat,[getdim(isector), getdim(isector)],"ed_buildh_main","Hmat")
    !
    !Get diagonal hybridization
    diag_hybr=0d0
    do ibath=1,Nbath
       do ispin=1,Nspin
          do iorb=1,Norb
             if(bath_type/="replica")then
                diag_hybr(ispin,iorb,ibath)=dmft_bath%v(ispin,iorb,ibath)
             else
                diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
             end if
          enddo
       enddo
    enddo
    !
    !
#ifdef _MPI
    if(MpiStatus)then
       call sp_set_mpi_matrix(MpiComm,spH0d,mpiIstart,mpiIend,mpiIshift)
       call sp_init_matrix(MpiComm,spH0d,Dim)
    else
       call sp_init_matrix(spH0d,Dim)
    endif
#else
    call sp_init_matrix(spH0d,Dim)
#endif
    call sp_init_matrix(spH0dw,DimDw)
    call sp_init_matrix(spH0up,DimUp)
    !
    !
    !-----------------------------------------------!
    !IMPURITY  HAMILTONIAN
    include "ED_HAMILTONIAN_MATVEC/stored/Himp.f90"
    !
    !LOCAL INTERACTION
    include "ED_HAMILTONIAN_MATVEC/stored/Hint.f90"
    !
    !BATH HAMILTONIAN
    include "ED_HAMILTONIAN_MATVEC/stored/Hbath.f90"
    !
    !IMPURITY- BATH HYBRIDIZATION
    include "ED_HAMILTONIAN_MATVEC/stored/Hhyb.f90"
    !-----------------------------------------------!
    !
    if(present(Hmat))then
       Hmat = 0d0
       allocate(Htmp_up(DimUp,DimUp));Htmp_up=0d0
       allocate(Htmp_dw(DimDw,DimDw));Htmp_dw=0d0
       !
#ifdef _MPI
       if(MpiStatus)then
          call sp_dump_matrix(MpiComm,spH0d,Hmat)
       else
          call sp_dump_matrix(spH0d,Hmat)
       endif
#else
       call sp_dump_matrix(spH0d,Hmat)
#endif
       call sp_dump_matrix(spH0up,Htmp_up)
       call sp_dump_matrix(spH0dw,Htmp_dw)
       !
       Hmat = Hmat + kronecker_product(eye(DimUp),Htmp_dw)
       Hmat = Hmat + kronecker_product(Htmp_up,eye(DimDw))
       !
       ! NEW:
       ! Hmat = Hmat + kronecker_product(eye(DimDw),Htmp_up)
       ! Hmat = Hmat + kronecker_product(Htmp_dw,eye(DimUp))
       !
       deallocate(Htmp_up,Htmp_dw)
    endif
    !
  end subroutine ed_buildh_main












  !####################################################################
  !        SPARSE MAT-VEC PRODUCT USING STORED SPARSE MATRIX 
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: Perform the matrix-vector product H*v used in the
  ! - serial
  ! - MPI
  !                     USE LL sparse format
  !+------------------------------------------------------------------+
  subroutine spMatVec_main(Nloc,v,Hv)
    integer                         :: Nloc
    real(8),dimension(Nloc)         :: v
    real(8),dimension(Nloc)         :: Hv
    real(8)                         :: val
    integer                         :: i,iup,idw,j,jup,jdw,jj
    !
    !
    Hv=0d0
    !
    do i = 1,Nloc
       do j=1,spH0d%row(i)%Size
          Hv(i) = Hv(i) + spH0d%row(i)%vals(j)*v(spH0d%row(i)%cols(j))
       enddo
    enddo
    !
    !DW:
    do iup=1,DimUp
       !
       do idw=1,DimDw
          i = iup + (idw-1)*DimUp
          do jj=1,spH0dw%row(idw)%Size
             jup = iup
             jdw = spH0dw%row(idw)%cols(jj)
             val = spH0dw%row(idw)%vals(jj)
             j     = jup +  (jdw-1)*DimUp
             Hv(i) = Hv(i) + val*V(j)
          enddo
       enddo
       !
    enddo
    !
    !UP:
    do idw=1,DimDw
       !
       do iup=1,DimUp
          i = iup + (idw-1)*DimUp
          do jj=1,spH0up%row(iup)%Size
             jup = spH0up%row(iup)%cols(jj)
             jdw = idw
             val = spH0up%row(iup)%vals(jj)
             j =  jup + (jdw-1)*DimUp
             Hv(i) = Hv(i) + val*V(j)
          enddo
       enddo
       !
    enddo
    !
  end subroutine spMatVec_main



#ifdef _MPI
  subroutine spMatVec_mpi_main(Nloc,v,Hv)
    integer                          :: Nloc
    real(8),dimension(Nloc)          :: v
    real(8),dimension(Nloc)          :: Hv
    !
    integer                          :: N
    real(8),dimension(:),allocatable :: vt,Hvt
    real(8)                          :: val
    integer                          :: i,iup,idw,j,jup,jdw,jj
    !local MPI
    integer                          :: irank
    integer                          :: MpiQup,MpiQdw
    !
    if(MpiComm==MPI_UNDEFINED)stop "spMatVec_mpi_cc ERROR: MpiComm = MPI_UNDEFINED"
    if(.not.MpiStatus)stop "spMatVec_mpi_cc ERROR: MpiStatus = F"
    !
    !Evaluate the local contribution: Hv_loc = Hloc*v
    Hv=0d0
    do i=1,Nloc                 !==spH0%Nrow
       do j=1,spH0d%row(i)%Size
          Hv(i) = Hv(i) + spH0d%row(i)%vals(j)*v(i)
       end do
    end do
    !
    mpiQdw=DimDw/MpiSize
    if(MpiRank<mod(DimDw,MpiSize))MpiQdw=MpiQdw+1
    !
    mpiQup=DimUp/MpiSize
    if(MpiRank<mod(DimUp,MpiSize))MpiQup=MpiQup+1
    !
    !
    !Non-local terms.
    !UP part: contiguous in memory.
    do idw=1,MpiQdw
       do iup=1,DimUp
          i = iup + (idw-1)*DimUp
          hxv_up: do jj=1,spH0up%row(iup)%Size
             jup = spH0up%row(iup)%cols(jj)
             jdw = idw
             val = spH0up%row(iup)%vals(jj)
             j   = jup + (idw-1)*DimUp
             Hv(i) = Hv(i) + val*v(j)
          end do hxv_up
       enddo
    end do
    !
    !DW part: non-contiguous in memory -> MPI transposition
    !Transpose the input vector as a whole:
    allocate(vt(mpiQup*DimDw)) ;vt=0d0
    allocate(Hvt(mpiQup*DimDw));Hvt=0d0
    call vector_transpose_MPI(DimUp,MpiQdw,v,DimDw,MpiQup,vt)
    Hvt=0d0    
    do iup=1,MpiQup
       do idw=1,DimDw
          i = idw + (iup-1)*DimDw
          hxv_dw: do jj=1,spH0dw%row(idw)%Size
             jdw = spH0dw%row(idw)%cols(jj)
             jup = iup
             val = spH0dw%row(idw)%vals(jj)
             j   = jdw + (jup-1)*DimDw
             Hvt(i) = Hvt(i) + val*vt(j)
          end do hxv_dw
       enddo
    end do
    deallocate(vt)
    allocate(vt(DimUp*mpiQdw));vt=0d0
    call vector_transpose_MPI(DimDw,mpiQup,Hvt,DimUp,mpiQdw,vt)
    Hv = Hv + vt
  end subroutine spMatVec_mpi_main
#endif
















     !####################################################################
     !            SPARSE MAT-VEC DIRECT ON-THE-FLY PRODUCT 
     !####################################################################
  subroutine directMatVec_main(Nloc,vin,Hv)
    integer                             :: Nloc
    real(8),dimension(Nloc)             :: vin
    real(8),dimension(Nloc)             :: Hv
    real(8),dimension(:),allocatable    :: vt,Hvt
    integer                             :: isector
    integer,dimension(Ns)               :: nup,ndw
    integer                             :: i,iup,idw
    integer                             :: j,jup,jdw
    integer                             :: m,mup,mdw
    integer                             :: ishift
    integer                             :: ms
    integer                             :: impi
    integer                             :: iorb,jorb,ispin,jspin,ibath
    integer                             :: kp,k1,k2,k3,k4
    integer                             :: alfa,beta
    real(8)                             :: sg1,sg2,sg3,sg4
    real(8)                             :: htmp,htmpup,htmpdw
    real(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    logical                             :: Jcondition
    integer                             :: MpiQup,MpiQdw
    !
    if(.not.Hstatus)stop "directMatVec_cc ERROR: Hsector NOT set"
    isector=Hsector
    !
    if(Nloc/=getdim(isector))stop "directMatVec_cc ERROR: Nloc != dim(isector)"
    !
    !Get diagonal hybridization
    diag_hybr=0d0
    do ibath=1,Nbath
       do ispin=1,Nspin
          do iorb=1,Norb
             if(bath_type/="replica")then
                diag_hybr(ispin,iorb,ibath)=dmft_bath%v(ispin,iorb,ibath)
             else
                diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
             end if
          enddo
       enddo
    enddo
    !
    !
    Hv=0d0
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "ED_HAMILTONIAN_MATVEC/direct/HxV_local.f90"
    !
    !NON-LOCAL TERMS:
    mpiQdw=DimDw
    !
    mpiQup=DimUp
    ! 
    !UP HAMILTONIAN TERMS: MEMORY CONTIGUOUS
    include "ED_HAMILTONIAN_MATVEC/direct/HxV_up.f90"
    !
    !DW HAMILTONIAN TERMS: MEMORY NON-CONTIGUOUS
    allocate(vt(mpiQup*DimDw)) ;vt=0d0
    allocate(Hvt(mpiQup*DimDw));Hvt=0d0
    include "ED_HAMILTONIAN_MATVEC/direct/HxV_dw.f90"
    deallocate(vt) ; allocate(vt(DimUp*mpiQdw)) ; vt=0d0
    Hv = Hv + Vt
    !-----------------------------------------------!
    !
    return
  end subroutine directMatVec_main





#ifdef _MPI
     subroutine directMatVec_MPI_main(Nloc,vin,Hv)
       integer                             :: Nloc
       real(8),dimension(Nloc)             :: Vin
       real(8),dimension(Nloc)             :: Hv
       real(8),dimension(:),allocatable    :: vt,Hvt
       integer,allocatable,dimension(:)    :: Counts,Displs
       integer                             :: isector
       integer,dimension(Ns)               :: nup,ndw
       integer                             :: i,iup,idw
       integer                             :: j,jup,jdw
       integer                             :: m,mup,mdw
       integer                             :: ms
       integer                             :: impi
       integer                             :: iorb,jorb,ispin,jspin,ibath
       integer                             :: kp,k1,k2,k3,k4
       integer                             :: alfa,beta
       real(8)                             :: sg1,sg2,sg3,sg4
       real(8)                             :: htmp,htmpup,htmpdw
       real(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
       logical                             :: Jcondition
       integer                             :: MpiQup,MpiQdw
       !
       if(.not.Hstatus)stop "directMatVec_cc ERROR: Hsector NOT set"
       isector=Hsector
       !
       !Get diagonal hybridization
       diag_hybr=0d0
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                if(bath_type/="replica")then
                   diag_hybr(ispin,iorb,ibath)=dmft_bath%v(ispin,iorb,ibath)
                else
                   diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
                end if
             enddo
          enddo
       enddo
       !    
       if(MpiComm==MPI_UNDEFINED)stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
       if(.not.MpiStatus)stop "directMatVec_MPI_cc ERROR: MpiStatus = F"
       !
       !
       Hv=0d0
       !
       !-----------------------------------------------!
       !LOCAL HAMILTONIAN PART: H_loc*vin = vout
       include "ED_HAMILTONIAN_MATVEC/direct/HxV_local.f90"
       !
       !NON-LOCAL TERMS:
       mpiQdw=DimDw/MpiSize
       if(MpiRank<mod(DimDw,MpiSize))MpiQdw=MpiQdw+1
       !
       mpiQup=DimUp/MpiSize
       if(MpiRank<mod(DimUp,MpiSize))MpiQup=MpiQup+1
       ! 
       !UP HAMILTONIAN TERMS: MEMORY CONTIGUOUS
       include "ED_HAMILTONIAN_MATVEC/direct/HxV_up.f90"
       !
       !DW HAMILTONIAN TERMS: MEMORY NON-CONTIGUOUS
       allocate(vt(mpiQup*DimDw)) ;vt=0d0
       allocate(Hvt(mpiQup*DimDw));Hvt=0d0
       call vector_transpose_MPI(DimUp,MpiQdw,Vin,DimDw,MpiQup,vt)
       include "ED_HAMILTONIAN_MATVEC/direct/HxV_dw.f90"
       deallocate(vt) ; allocate(vt(DimUp*mpiQdw)) ;vt=0d0
       call vector_transpose_MPI(DimDw,mpiQup,Hvt,DimUp,mpiQdw,vt)
       Hv = Hv + Vt
       !-----------------------------------------------!
       !
       return
     end subroutine directMatVec_MPI_main
#endif













     !####################################################################
     !               ALL-2-ALL-V VECTOR MPI TRANSPOSITION 
     !####################################################################
#ifdef _MPI
     subroutine vector_transpose_MPI(nrow,qcol,a,ncol,qrow,b)    
       integer                            :: nrow,ncol,qrow,qcol
       real(8)                            :: a(nrow,qcol)
       real(8)                            :: b(ncol,qrow)
       integer,allocatable,dimension(:,:) :: send_counts,send_offset
       integer,allocatable,dimension(:,:) :: recv_counts,recv_offset
       integer                            :: counts,Ntot
       integer :: i,j,irank,ierr
       !
       counts = Nrow/MpiSize
       Ntot   = Ncol/MpiSize
       if(mod(Ncol,MpiSize)/=0)Ntot=Ntot+1
       !
       allocate(send_counts(0:MpiSize-1,Ntot));send_counts=0
       allocate(send_offset(0:MpiSize-1,Ntot));send_offset=0
       allocate(recv_counts(0:MpiSize-1,Ntot));recv_counts=0
       allocate(recv_offset(0:MpiSize-1,Ntot));recv_offset=0
       !
       do i=1,qcol
          do irank=0,MpiSize-1
             if(irank < mod(Nrow,MpiSize))then
                send_counts(irank,i) = counts+1
             else
                send_counts(irank,i) = counts
             endif
          enddo
       enddo
       !
       do i=1,Ntot
          call MPI_AllToAll(&
               send_counts(:,i),1,MPI_INTEGER,&
               recv_counts(:,i),1,MPI_INTEGER,&
               MpiComm,ierr)
       enddo
       !
       do i=1,Ntot
          do irank=1,MpiSize-1
             send_offset(irank,i) = send_counts(irank-1,i) + send_offset(irank-1,i)
          enddo
       enddo
       !
       !Get the irank=0 elements, i.e. first entries:
       recv_offset(0,1) = 0
       do i=2,Ntot
          recv_offset(0,i) = sum(recv_counts(0,:i-1))
       enddo
       !the rest of the entries:
       do i=1,Ntot
          do irank=1,MpiSize-1
             recv_offset(irank,i) = recv_offset(irank-1,i) + sum(recv_counts(irank-1,:))
          enddo
       enddo
       !
       !
       do j=1,Ntot
          call MPI_AllToAllV(&
               A(:,j),send_counts(:,j),send_offset(:,j),MPI_DOUBLE_PRECISION,&
               B(:,:),recv_counts(:,j),recv_offset(:,j),MPI_DOUBLE_PRECISION,&
               MpiComm,ierr)
       enddo
       !
       call local_transpose(b,ncol,qrow)
       !
       return
     end subroutine vector_transpose_MPI


     subroutine local_transpose(mat,nrow,ncol)
       integer                      :: nrow,ncol
       real(8),dimension(Nrow,Ncol) :: mat
       mat = transpose(reshape(mat,[Ncol,Nrow]))
     end subroutine local_transpose
#endif



   end MODULE ED_HAMILTONIAN_MATVEC













! subroutine check_first_last(nup,ndw)
!   integer                          :: isector,nup,ndw
!   integer                          :: dim,dimUp,dimDw
!   integer                          :: i,iup,idw
!   integer                          :: j,jup,jdw
!   integer                          :: m,mup,mdw
!   integer                          :: ishift,irank
!   integer                          :: first_state   ,last_state
!   integer                          :: first_state_dw,last_state_dw
!   integer                          :: first_state_up,last_state_up
!   integer,allocatable,dimension(:) :: map_first_state_up,map_last_state_up
!   integer,allocatable,dimension(:) :: map_first_state_dw,map_last_state_dw
!   !
!   isector=getSector(nup,ndw)
!   !
!   Dim  = getDim(isector)
!   DimDw  = getDimDw(isector)
!   DimUp  = getDimUp(isector)
!   !
!   mpiQdw = dimdw/MpiSize
!   mpiRdw = 0
!   if(MpiRank==(MpiSize-1))mpiRdw=mod(dimdw,MpiSize)
!   first_state_dw = MpiRank*mpiQdw + 1
!   last_state_dw  = (MpiRank+1)*mpiQdw + mpiRdw
!   if(MpiMaster)print*,"         mpiRank,   mpi_Q,   mpi_R,   mpi_CHunk,   mpi_Istart,   mpi_Iend,   mpi_Iend-mpi_Istart"
!   do irank=0,MpiSize-1
!      if(MpiRank==irank)then
!         print*,MpiRank,MpiQdw,MpiRdw,mpiQdw+mpiRdw,first_state_dw,last_state_dw,last_state_dw-first_state_dw+1
!         call MPI_Barrier(MpiComm,MpiIerr)
!      endif
!      call MPI_Barrier(MpiComm,MpiIerr)
!   enddo
!   call MPI_Barrier(MpiComm,MpiIerr)
!   !
!   call sleep(1)
!   call MPI_Barrier(MpiComm,MpiIerr)
!   !
!   !
!   mpiQup = dimup/MpiSize
!   mpiRup = 0
!   if(MpiRank==(MpiSize-1))mpiRup=mod(dimup,MpiSize)
!   first_state_up = MpiRank*mpiQup + 1
!   last_state_up  = (MpiRank+1)*mpiQup + mpiRup
!   if(MpiMaster)print*,"         mpiRank,   mpi_Q,   mpi_R,   mpi_CHunk,   mpi_Istart,   mpi_Iend,   mpi_Iend-mpi_Istart"
!   do irank=0,MpiSize-1
!      if(MpiRank==irank)then
!         print*,MpiRank,MpiQup,MpiRup,mpiQup+mpiRup,first_state_up,last_state_up,last_state_up-first_state_up+1
!         call MPI_Barrier(MpiComm,MpiIerr)
!      endif
!      call MPI_Barrier(MpiComm,MpiIerr)
!   enddo
!   call MPI_Barrier(MpiComm,MpiIerr)
!   !
!   call sleep(1)
!   call MPI_Barrier(MpiComm,MpiIerr)
!   !
!   mpiQ = dim/MpiSize
!   mpiR = 0
!   if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
!   ishift      = MpiRank*mpiQ
!   first_state = MpiRank*mpiQ + 1
!   last_state  = (MpiRank+1)*mpiQ + mpiR
!   if(MpiMaster)print*,"         mpiRank,   mpi_Q,   mpi_R,   mpi_CHunk,   mpi_Istart,   mpi_Iend,   mpi_Iend-mpi_Istart"
!   do irank=0,MpiSize-1
!      if(MpiRank==irank)then
!         print*,MpiRank,MpiQ,MpiR,mpiQ+mpiR,first_state,last_state,last_state-first_state+1
!         call MPI_Barrier(MpiComm,MpiIerr)
!      endif
!      call MPI_Barrier(MpiComm,MpiIerr)
!   enddo
!   call MPI_Barrier(MpiComm,MpiIerr)
!   !
!   !
!   !The threads generates all the states of the sector, however in a highly non-local way
!   do iup=first_state_up,last_state_up
!      do idw=1,DimDw
!         i = iup + (idw-1)*DimUp
!         write(250+MpiRank,*)i,idw,iup
!      enddo
!   end do
!   !
!   !This instead is local.
!   do idw=first_state_dw,last_state_dw
!      do iup=1,DimUp
!         i = iup + (idw-1)*DimUp
!         write(300+MpiRank,*)i,idw,iup
!      enddo
!   end do
!   !
!   !
!   allocate(map_first_state_up(DimDw),map_last_state_up(DimDw))
!   first_state_dw = huge(1)
!   last_state_dw  = -huge(1)
!   !
!   map_first_state_up = huge(1)
!   map_last_state_up  = -huge(1)
!   !
!   do i=first_state,last_state
!      idw=idw_index(i,DimUp)
!      !
!      if(idw < first_state_dw) first_state_dw = idw
!      if(idw > last_state_dw ) last_state_dw  = idw
!      !
!      iup=iup_index(i,DimUp)
!      if(iup < map_first_state_up(idw)) map_first_state_up(idw) = iup
!      if(iup > map_last_state_up(idw) ) map_last_state_up(idw)  = iup
!      !
!      write(200+MpiRank,*)i,idw,iup
!      write(100+MpiRank,*)i,first_state_dw,last_state_dw
!   enddo
!   !
!   stop
! end subroutine check_first_last




!   subroutine build_up_MPI(isector)
!     integer :: isector
!     integer :: Dim,DimUp,DimDw
!     integer :: irank
!     integer :: i,iup,idw
!     integer :: j,jup,jdw
!     !
!     Dim  = getDim(isector)
!     DimDw  = getDimDw(isector)
!     DimUp  = getDimUp(isector)
!     !
!     !
!     !DOWN part:
!     MpiQdw = DimDw/MpiSize
!     MpiRdw = 0
!     if(MpiRank==(MpiSize-1))MpiRdw=mod(DimDw,MpiSize)
!     ishiftDw = MpiRank*MpiQdw
!     first_state_dw = MpiRank*mpiQdw + 1
!     last_state_dw  = (MpiRank+1)*mpiQdw + mpiRdw
!     !
! #ifdef _MPI    
!     if(MpiStatus.AND.ed_verbose>3)then
!        if(MpiMaster)print*,"         mpiRank,   mpi_Q,   mpi_R,   mpi_CHunk,   mpi_Istart,   mpi_Iend,   mpi_Iend-mpi_Istart"
!        do irank=0,MpiSize-1
!           if(MpiRank==irank)then
!              write(LOGfile,*)MpiRank,MpiQdw,MpiRdw,mpiQdw+mpiRdw,first_state_dw,last_state_dw,last_state_dw-first_state_dw+1
!              call MPI_Barrier(MpiComm,MpiIerr)
!           endif
!           call MPI_Barrier(MpiComm,MpiIerr)
!        enddo
!        call MPI_Barrier(MpiComm,MpiIerr)
!     endif
! #endif
!     !
!     !
!     !UP part:
!     MpiQup = DimUp/MpiSize
!     MpiRup = 0
!     if(MpiRank==(MpiSize-1))MpiRup=mod(DimUp,MpiSize)
!     ishiftUp = MpiRank*MpiQup
!     first_state_up = MpiRank*mpiQup + 1
!     last_state_up  = (MpiRank+1)*mpiQup + mpiRup
!     !
! #ifdef _MPI    
!     if(MpiStatus.AND.ed_verbose>3)then
!        if(MpiMaster)print*,"         mpiRank,   mpi_Q,   mpi_R,   mpi_CHunk,   mpi_Istart,   mpi_Iend,   mpi_Iend-mpi_Istart"
!        do irank=0,MpiSize-1
!           if(MpiRank==irank)then
!              write(LOGfile,*)MpiRank,MpiQup,MpiRup,mpiQup+mpiRup,first_state_up,last_state_up,last_state_up-first_state_up+1
!              call MPI_Barrier(MpiComm,MpiIerr)
!           endif
!           call MPI_Barrier(MpiComm,MpiIerr)
!        enddo
!        call MPI_Barrier(MpiComm,MpiIerr)
!     endif
! #endif
!     !
!     !
!     !FULL:
!     MpiQ = Dim/MpiSize
!     MpiR = 0
!     if(MpiRank==(MpiSize-1))MpiR=mod(Dim,MpiSize)
!     ishift = MpiRank*MpiQ
!     first_state = MpiRank*mpiQ + 1
!     last_state  = (MpiRank+1)*mpiQ + mpiR
!     !
! #ifdef _MPI    
!     if(MpiStatus.AND.ed_verbose>3)then
!        if(MpiMaster)print*,"         mpiRank,   mpi_Q,   mpi_R,   mpi_CHunk,   mpi_Istart,   mpi_Iend,   mpi_Iend-mpi_Istart"
!        do irank=0,MpiSize-1
!           if(MpiRank==irank)then
!              write(LOGfile,*)MpiRank,MpiQ,MpiR,mpiQ+mpiR,first_state,last_state,last_state-first_state+1
!              call MPI_Barrier(MpiComm,MpiIerr)
!           endif
!           call MPI_Barrier(MpiComm,MpiIerr)
!        enddo
!        call MPI_Barrier(MpiComm,MpiIerr)
!     endif
! #endif
!     !
!     !
!     if(allocated(map_first_state_dw))deallocate(map_first_state_dw)
!     if(allocated(map_last_state_dw))deallocate(map_last_state_dw)
!     if(allocated(map_first_state_up))deallocate(map_first_state_up)
!     if(allocated(map_last_state_up))deallocate(map_last_state_up)
!     !
!     allocate(map_first_state_dw(1),map_last_state_dw(1))
!     allocate(map_first_state_up(DimDw),map_last_state_up(DimDw))
!     !
!     map_first_state_dw =  huge(1)
!     map_last_state_dw  = -huge(1)
!     !
!     map_first_state_up =  huge(1)
!     map_last_state_up  = -huge(1)
!     !
!     do i=first_state,last_state
!        idw=idw_index(i,DimUp)
!        iup=iup_index(i,DimUp)
!        !
!        if(idw < map_first_state_dw(1)) map_first_state_dw(1) = idw
!        if(idw > map_last_state_dw(1) ) map_last_state_dw(1)  = idw
!        !
!        if(iup < map_first_state_up(idw)) map_first_state_up(idw) = iup
!        if(iup > map_last_state_up(idw) ) map_last_state_up(idw)  = iup
!     enddo
!   end subroutine build_up_MPI









! subroutine ed_buildh_orbs(isector,Hmat)
!   integer                             :: isector   
!   real(8),dimension(:,:),optional     :: Hmat
!   real(8),dimension(:,:),allocatable  :: Htmp_up,Htmp_dw,Hrdx
!   integer,dimension(2*Ns_Ud)          :: Indices
!   integer,dimension(Ns_Ud,Ns_Orb)     :: Nups,Ndws
!   integer,dimension(Ns_Ud)            :: Nup,Ndw
!   integer,dimension(Ns_Orb)           :: Ibup,Ibdw
!   integer                             :: i,iup,idw
!   integer                             :: j,jup,jdw
!   integer                             :: m,mup,mdw
!   integer                             :: ms,iud
!   integer                             :: iorb,jorb,ispin,jspin,ibath
!   integer                             :: kp,k1,k2,k3,k4
!   integer                             :: alfa,beta
!   real(8)                             :: sg1,sg2,sg3,sg4
!   real(8)                             :: htmp,htmpup,htmpdw
!   real(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
!   logical                             :: Jcondition
!   !
!   if(.not.Hstatus)stop "ed_buildh_orbs ERROR: Hsector NOT set"
!   isector=Hsector
!   !
!   if(present(Hmat))call assert_shape(Hmat,[getdim(isector), getdim(isector)],"ed_buildh_orbs","Hmat")
!   !
!   !Get diagonal hybridization
!   diag_hybr=0d0
!   do ibath=1,Nbath
!      do ispin=1,Nspin
!         do iorb=1,Norb
!            if(bath_type/="replica")then
!               diag_hybr(ispin,iorb,ibath)=dmft_bath%v(ispin,iorb,ibath)
!            else
!               diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
!            end if
!         enddo
!      enddo
!   enddo
!   !
!   !
!   call sp_init_matrix(spH0d,Dim)
!   do iud=1,Ns_Ud
!      call sp_init_matrix(spH0dws(iud),DimDws(iud))
!      call sp_init_matrix(spH0ups(iud),DimUps(iud))
!   enddo
!   !
!   !-----------------------------------------------!
!   !IMPURITY  HAMILTONIAN
!   include "ED_HAMILTONIAN_MATVEC/stored/Himp_orbs.f90"
!   !
!   !LOCAL INTERACTION
!   include "ED_HAMILTONIAN_MATVEC/stored/Hint_orbs.f90"
!   !
!   !BATH HAMILTONIAN
!   include "ED_HAMILTONIAN_MATVEC/stored/Hbath_orbs.f90"
!   !
!   !IMPURITY- BATH HYBRIDIZATION
!   include "ED_HAMILTONIAN_MATVEC/stored/Hhyb_orbs.f90"
!   !-----------------------------------------------!
!   !
!   if(present(Hmat))then
!      Hmat = 0d0
!      call sp_dump_matrix(spH0d,Hmat)
!      allocate(Htmp_up(DimUp,DimUp));Htmp_up=0d0
!      allocate(Htmp_dw(DimDw,DimDw));Htmp_dw=0d0
!      !
!      do iud=1,Ns_Ud
!         allocate(Hrdx(DimUps(iud),DimUps(iud)));Hrdx=0d0
!         call sp_dump_matrix(spH0ups(iud),Hrdx)
!         call build_Htmp_up(iud,Hrdx,DimUps(iud),Htmp_up)
!         Hmat = Hmat + kronecker_product(Htmp_up,eye(DimDw))          
!         deallocate(Hrdx)
!         allocate(Hrdx(DimDws(iud),DimDws(iud)));Hrdx=0d0
!         call sp_dump_matrix(spH0dws(iud),Hrdx)
!         call build_Htmp_dw(iud,Hrdx,DimDws(iud),Htmp_dw)
!         Hmat = Hmat + kronecker_product(eye(DimUp),Htmp_dw)
!         deallocate(Hrdx)
!      enddo
!      deallocate(Htmp_up,Htmp_dw)
!   endif
!   !
! contains

!   subroutine build_Htmp_up(iud,H,Dim,Hup)
!     integer                        :: iud,Dim,i
!     real(8),dimension(Dim,Dim)     :: H
!     real(8),dimension(DimUp,DimUp) :: Hup
!     if(dim/=DimUps(iud))stop "error in build_Htmp_up"
!     if(iud==1)then
!        Hup= kronecker_product(H,eye(product(DimUps(2:))))
!     else if(iud==Ns_Ud)then
!        Hup= kronecker_product(eye(product(DimUps(1:Ns_Ud-1))),H)
!     else
!        Hup= kronecker_product( &
!             kronecker_product( &
!             eye(product(DimUps(1:iud-1))), H) , eye(product(DimUps(iud+1:Ns_Ud))) )
!     end if
!   end subroutine build_Htmp_up

!   subroutine build_Htmp_dw(iud,H,Dim,Hdw)
!     integer                        :: iud,Dim,i
!     real(8),dimension(Dim,Dim)     :: H
!     real(8),dimension(DimDw,DimDw) :: Hdw
!     if(dim/=DimDws(iud))stop "error in build_Htmp_dw"
!     if(iud==1)then
!        Hdw= kronecker_product(H,eye(product(DimDws(2:))))
!     else if(iud==Ns_Ud)then
!        Hdw= kronecker_product(eye(product(DimDws(1:Ns_Ud-1))),H)
!     else
!        Hdw= kronecker_product( &
!             kronecker_product( &
!             eye(product(DimDws(1:iud-1))), H) , eye(product(DimDws(iud+1:Ns_Ud))) )
!     end if
!   end subroutine build_Htmp_dw
! end subroutine ed_buildh_orbs









! subroutine ed_dryrunH_main(isector)
!   integer                             :: isector   
!   integer,dimension(Ns)               :: nup,ndw
!   integer                             :: i,iup,idw
!   integer                             :: j,jup,jdw
!   integer                             :: m,mup,mdw
!   integer                             :: ms
!   integer                             :: impi,impi_up,impi_dw
!   integer                             :: iorb,jorb,ispin,jspin,ibath
!   integer                             :: kp,k1,k2,k3,k4
!   integer                             :: alfa,beta
!   real(8)                             :: sg1,sg2,sg3,sg4
!   real(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
!   logical                             :: Jcondition
!   !
!   if(.not.Hstatus)stop "ed_dryrun_main ERROR: Hsector NOT set"
!   isector=Hsector
!   !
!   !Get diagonal hybridization
!   diag_hybr=0d0
!   do ibath=1,Nbath
!      do ispin=1,Nspin
!         do iorb=1,Norb
!            if(bath_type/="replica")then
!               diag_hybr(ispin,iorb,ibath)=dmft_bath%v(ispin,iorb,ibath)
!            else
!               diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
!            end if
!         enddo
!      enddo
!   enddo
!   !
!   allocate(vecNnz_d(MpiQ+mpiR))
!   allocate(vecNnz_dw(MpiQdw+MpiRdw))
!   allocate(vecNnz_up(MpiQup+MpiRup))
!   if(Jhflag)allocate(vecNnz_nd(MpiQ+MpiR))
!   vecNnz_d  = 1               !this stores the diagonal, as such it is identically 1 for any row
!   vecNnz_dw = 0
!   vecNnz_up = 0
!   if(Jhflag)vecNnz_nd = 0
!   !
!   !-----------------------------------------------!
!   !IMPURITY  HAMILTONIAN
!   include "ED_HAMILTONIAN_MATVEC/dryrun/Himp_dryrun.f90"
!   !
!   !LOCAL INTERACTION
!   include "ED_HAMILTONIAN_MATVEC/dryrun/Hint_dryrun.f90"
!   !
!   !BATH HAMILTONIAN
!   include "ED_HAMILTONIAN_MATVEC/dryrun/Hbath_dryrun.f90"
!   !
!   !IMPURITY- BATH HYBRIDIZATION
!   include "ED_HAMILTONIAN_MATVEC/dryrun/Hhyb_dryrun.f90"
!   !-----------------------------------------------!
! end subroutine ed_dryrunH_main








! subroutine spMatVec_orbs(Nloc,v,Hv)
!   integer                         :: Nloc
!   real(8),dimension(Nloc)         :: v
!   real(8),dimension(Nloc)         :: Hv
!   integer                         :: i,iup,idw,j,iorb
!   integer,dimension(2*Norb)       :: Indices,Jndices
!   type(sparse_element_ll),pointer :: c
!   !
!   !
!   Hv=0d0
!   !
!   do i = 1,Nloc
!      c => spH0d%row(i)%root%next
!      do while(associated(c))
!         Hv(i) = Hv(i) + c%cval*V(c%col)    !<== V
!         c => c%next
!      enddo
!      nullify(c)
!   enddo
!   !
!   do i=1,Dim
!      call state2indices(i,[DimUps,DimDws],Indices)
!      do iorb=1,Norb
!         idw = Indices(iorb+Norb)
!         !
!         c => spH0dws(iorb)%row(idw)%root%next
!         do while(associated(c))
!            Jndices = Indices ; Jndices(iorb+Norb) = c%col
!            call indices2state(Jndices,[DimUps,DimDws],j)
!            Hv(i) = Hv(i) + c%cval*V(j)
!            c => c%next
!         enddo
!         nullify(c)
!         !
!         !
!         iup = Indices(iorb)
!         c => spH0ups(iorb)%row(iup)%root%next
!         do while(associated(c))
!            Jndices = Indices ; Jndices(iorb) = c%col
!            call indices2state(Jndices,[DimUps,DimDws],j)
!            Hv(i) = Hv(i) + c%cval*V(j)
!            c => c%next
!         enddo
!         nullify(c)
!         !
!      enddo
!   enddo
!   !
! end subroutine spMatVec_orbs




!   !+------------------------------------------------------------------+
!   !PURPOSE: Perform the matrix-vector product H*v used in the
!   ! - serial
!   ! - MPI
!   !                     USE ELL sparse format
!   !+------------------------------------------------------------------+
!   subroutine dpMatVec_main(Nloc,v,Hv)
!     integer                 :: Nloc
!     real(8),dimension(Nloc) :: v
!     real(8),dimension(Nloc) :: Hv
!     integer                 :: i,iup,idw,j,ij
!     !
!     Hv=0d0
!     !
!     do i=1,Nloc
!        do j=1,dpH0d%row(i)%Size
!           Hv(i) = Hv(i) + dpH0d%row(i)%vals(j)*v(dpH0d%row(i)%cols(j))
!        end do
!     end do
!     !
!     if(jhflag)then
!        do i=1,Nloc
!           do j=1,dpH0nd%row(i)%Size
!              Hv(i) = Hv(i) + dpH0nd%row(i)%vals(j)*v(dpH0nd%row(i)%cols(j))
!           end do
!        end do
!     endif
!     !
!     do idw=1,DimDw
!        do iup=1,DimUp
!           i = iup + (idw-1)*DimUp
!           do ij=1,dpH0dw%row(idw)%Size             
!              j = iup +  (dpH0dw%row(idw)%cols(ij)-1)*DimUp
!              Hv(i) = Hv(i) + dpH0dw%row(idw)%vals(ij)*v(j)
!           end do
!        enddo
!     enddo
!     !
!     do idw=1,DimDw
!        do iup=1,DimUp
!           i = iup + (idw-1)*DimUp
!           do ij=1,dpH0up%row(iup)%Size
!              j = dpH0up%row(iup)%cols(ij) + (idw-1)*DimUp
!              Hv(i) = Hv(i) + dpH0up%row(iup)%vals(ij)*V(j)
!           end do
!        enddo
!     enddo
!     !
!   end subroutine dpMatVec_main



! #ifdef _MPI
!   subroutine dpMatVec_mpi_main(Nloc,v,Hv)
!     integer                          :: Nloc
!     real(8),dimension(Nloc)          :: v
!     real(8),dimension(Nloc)          :: Hv
!     integer                          :: i
!     integer                          :: N
!     real(8),dimension(:),allocatable :: Vin,Vout
!     integer,allocatable,dimension(:) :: Counts,Displs
!     integer                          :: iup,idw,j,ij
!     integer                          :: impi_up,impi_dw,impi
!     !
!     if(MpiComm==MPI_UNDEFINED)stop "spMatVec_mpi_cc ERROR: MpiComm = MPI_UNDEFINED"
!     MpiSize = get_Size_MPI(MpiComm)
!     MpiRank = get_Rank_MPI(MpiComm)
!     !
!     N=0
!     call MPI_AllReduce(Nloc,N,1,MPI_Integer,MPI_Sum,MpiComm,MpiIerr)
!     if(N/=Dim)stop "spMatVec_mpi_cc ERROR: N != Dim"
!     !
!     !
!     Hv=0d0
!     !
!     !DIAGONAL PART:
!     do i=1,Nloc
!        do j=1,dpH0d%row(i)%Size
!           Hv(i) = Hv(i) + dpH0d%row(i)%vals(j)*V(dpH0d%row(i)%cols(j)-Ishift)
!        end do
!     end do
!     !
!     allocate(Counts(0:MpiSize-1))
!     Counts(0:)        = mpiQ
!     Counts(MpiSize-1) = mpiQ+mod(N,MpiSize)
!     !
!     allocate(Displs(0:MpiSize-1))
!     forall(i=0:MpiSize-1)Displs(i)=i*mpiQ
!     !
!     allocate(Vin(N)) ; Vin = 0d0
!     call MPI_Allgatherv(V(1:Nloc),Nloc,MPI_Double_Precision,Vin,Counts,Displs,MPI_Double_Precision,MpiComm,MpiIerr)
!     !
!     !NON-DIAGONAL PART: including S-E and P-H terms.
!     if(jhflag)then
!        do i=1,Nloc
!           do j=1,dpH0nd%row(i)%Size
!              Hv(i) = Hv(i) + dpH0nd%row(i)%vals(j)*Vin(dpH0nd%row(i)%cols(j))
!           end do
!        end do
!     endif
!     !
!     !
!     !DW PART: this is local and contiguous in memory, i.e. index i corresponds to consecutive states.
!     allocate(Vout(N)) ; Vout = 0d0
!     !
!     do idw=first_state_dw,last_state_dw
!        do iup=1,DimUp
!           i = iup + (idw-1)*DimUp
!           do ij=1,dpH0dw%row(idw-IshiftDw)%Size
!              j = iup +  (dpH0dw%row(idw-IshiftDw)%cols(ij)-1)*DimUp
!              Vout(i) = Vout(i) + dpH0dw%row(idw-IshiftDw)%vals(ij)*Vin(j)
!           end do
!        enddo
!     enddo
!     !
!     !UP PART: this is non-local and non-contiguous in memory, i.e. index i corresponds to non-consecutive states
!     do iup=first_state_up,last_state_up
!        do idw=1,DimDw
!           i = iup + (idw-1)*DimUp
!           do ij=1,dpH0up%row(iup-IshiftUp)%Size
!              j = dpH0up%row(iup-IshiftUp)%cols(ij) + (idw-1)*DimUp
!              Vout(i) = Vout(i) + dpH0up%row(iup-IshiftUp)%vals(ij)*Vin(j)
!           end do
!        enddo
!     enddo
!     !
!     !Now we pack back the Vout content to Hv vectors:
!     Vin=0d0
!     call AllReduce_MPI(MpiComm,Vout,Vin)
!     do i=first_state,last_state
!        Hv(i-Ishift) = Hv(i-Ishift) + Vin(i)
!     end do
!   end subroutine dpMatVec_mpi_main
! #endif
