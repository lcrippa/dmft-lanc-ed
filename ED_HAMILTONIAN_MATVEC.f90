!########################################################################
!PURPOSE  : Build the impurity Hamiltonian
!|ImpUP,(2ImpUP),BathUP;,ImpDW,(2ImpDW),BathDW >
! |1,2;3...Ns>_UP * |Ns+1,Ns+2;Ns+3,...,2*Ns>_DOWN
!########################################################################
MODULE ED_HAMILTONIAN_MATVEC
  USE SF_CONSTANTS,only:zero
  USE SF_MISC,    only: assert_shape
  USE SF_LINALG, only: kronecker_product,zeye
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
  public  :: build_Hv_sector
  public  :: delete_Hv_sector

  !
  !>Sparse Mat-Vec product using stored sparse matrix 
  public  :: spMatVec_cc
#ifdef _MPI
  public  :: spMatVec_MPI_cc
#endif
  !
  !
  !>Sparse Mat-Vec direct on-the-fly product 
  public  :: directMatVec_cc
#ifdef _MPI
  public  :: directMatVec_MPI_cc
#endif
  !
  !>Related auxiliary routines:
  public  :: ed_hamiltonian_matvec_set_MPI
  public  :: ed_hamiltonian_matvec_del_MPI


  !> MPI local variables (shared)
#ifdef _MPI
  integer                          :: MpiComm=MPI_UNDEFINED
#else
  integer                          :: MpiComm=0
#endif
  logical                          :: MpiStatus=.false.
  logical                          :: MpiMaster=.true.
  integer                          :: MpiIerr
  integer                          :: MpiRank=0
  integer                          :: MpiSize=1
  integer                          :: MpiQ=1
  integer                          :: MpiQup=1
  integer                          :: MpiQdw=1
  integer                          :: MpiR=0
  integer                          :: MpiRup=0
  integer                          :: MpiRdw=0
  integer                          :: Ishift
  integer                          :: IshiftUp
  integer                          :: IshiftDw
  !
  integer                          :: first_state   ,last_state
  integer                          :: first_state_dw,last_state_dw
  integer                          :: first_state_up,last_state_up

  integer,allocatable,dimension(:) :: map_first_state_dw,map_last_state_dw
  integer,allocatable,dimension(:) :: map_first_state_up,map_last_state_up

  integer                          :: Dim
  integer                          :: DimUp
  integer                          :: DimDw

  integer                          :: Hsector=0
  logical                          :: Hstatus=.false.
  type(sector_map)                 :: H,Hs(2)





contains


  !####################################################################
  !                        AUXILIARY ROUTINES
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


  subroutine build_Hv_sector(isector,Hmat)
    integer                            :: isector,SectorDim
    complex(8),dimension(:,:),optional :: Hmat   
    integer                            :: irank
    integer                            :: i,iup,idw
    integer                            :: j,jup,jdw
    !
    Hsector=isector
    Hstatus=.true.
    !
    call build_sector(isector,H)
    call build_sector(isector,Hs)
    !
    !
    Dim  = getDim(isector)
    DimDw  = getDimDw(isector)
    DimUp  = getDimUp(isector)
    !
    !
    !DOWN part:
    MpiQdw = DimDw/MpiSize
    MpiRdw = 0
    if(MpiRank==(MpiSize-1))MpiRdw=mod(DimDw,MpiSize)
    ishiftDw = MpiRank*MpiQdw
    first_state_dw = MpiRank*mpiQdw + 1
    last_state_dw  = (MpiRank+1)*mpiQdw + mpiRdw
    !
    !
    !UP part:
    MpiQup = DimUp/MpiSize
    MpiRup = 0
    if(MpiRank==(MpiSize-1))MpiRup=mod(DimUp,MpiSize)
    ishiftUp = MpiRank*MpiQup
    first_state_up = MpiRank*mpiQup + 1
    last_state_up  = (MpiRank+1)*mpiQup + mpiRup
    !
    !FULL:
    MpiQ = Dim/MpiSize
    MpiR = 0
    if(MpiRank==(MpiSize-1))MpiR=mod(Dim,MpiSize)
    ishift = MpiRank*MpiQ
    first_state = MpiRank*mpiQ + 1
    last_state  = (MpiRank+1)*mpiQ + mpiR
    !
#ifdef _MPI    
    if(MpiStatus.AND.ed_verbose>3)then
       if(MpiMaster)print*,"         mpiRank,   mpi_Q,   mpi_R,   mpi_CHunk,   mpi_Istart,   mpi_Iend,   mpi_Iend-mpi_Istart"
       call MPI_Barrier(MpiComm,MpiIerr)
       if(MpiMaster)print*,"DW:"
       call MPI_Barrier(MpiComm,MpiIerr)
       do irank=0,MpiSize-1
          if(MpiRank==irank)then
             write(*,*)MpiRank,MpiQdw,MpiRdw,mpiQdw+mpiRdw,first_state_dw,last_state_dw,last_state_dw-first_state_dw+1
             call MPI_Barrier(MpiComm,MpiIerr)
          endif
          call MPI_Barrier(MpiComm,MpiIerr)
       enddo
       call MPI_Barrier(MpiComm,MpiIerr)
       !
       if(MpiMaster)print*,"UP:"
       call MPI_Barrier(MpiComm,MpiIerr)
       do irank=0,MpiSize-1
          if(MpiRank==irank)then
             write(*,*)MpiRank,MpiQup,MpiRup,mpiQup+mpiRup,first_state_up,last_state_up,last_state_up-first_state_up+1
             call MPI_Barrier(MpiComm,MpiIerr)
          endif
          call MPI_Barrier(MpiComm,MpiIerr)
       enddo
       call MPI_Barrier(MpiComm,MpiIerr)
       !
       if(MpiMaster)print*,"Full:"
       call MPI_Barrier(MpiComm,MpiIerr)
       do irank=0,MpiSize-1
          if(MpiRank==irank)then
             write(*,*)MpiRank,MpiQ,MpiR,mpiQ+mpiR,first_state,last_state,last_state-first_state+1
             call MPI_Barrier(MpiComm,MpiIerr)
          endif
          call MPI_Barrier(MpiComm,MpiIerr)
       enddo
       call MPI_Barrier(MpiComm,MpiIerr)
    endif
#endif
    !
    !
    if(allocated(map_first_state_dw))deallocate(map_first_state_dw)
    if(allocated(map_last_state_dw))deallocate(map_last_state_dw)
    if(allocated(map_first_state_up))deallocate(map_first_state_up)
    if(allocated(map_last_state_up))deallocate(map_last_state_up)
    !
    allocate(map_first_state_dw(1),map_last_state_dw(1))
    allocate(map_first_state_up(DimDw),map_last_state_up(DimDw))
    !
    map_first_state_dw =  huge(1)
    map_last_state_dw  = -huge(1)
    !
    map_first_state_up =  huge(1)
    map_last_state_up  = -huge(1)
    !
    do i=first_state,last_state
       idw=idw_index(i,DimUp)
       iup=iup_index(i,DimUp)
       !
       if(idw < map_first_state_dw(1)) map_first_state_dw(1) = idw
       if(idw > map_last_state_dw(1) ) map_last_state_dw(1)  = idw
       !
       if(iup < map_first_state_up(idw)) map_first_state_up(idw) = iup
       if(iup > map_last_state_up(idw) ) map_last_state_up(idw)  = iup
    enddo
    !



    if(present(Hmat))then
       if(any( shape(Hmat) /= [Dim,Dim]))stop "setup_Hv_sector ERROR: size(Hmat) != SectorDim**2"
       !
       call ed_buildH_c(isector,Hmat)
       return
    endif
    !
    select case (ed_sparse_H)
    case (.true.)
       !
       !
       call ed_buildH_c(isector)

    case (.false.)
       !nothing to be done
    end select
    !
  end subroutine build_Hv_sector


  subroutine delete_Hv_sector()
    call delete_sector(Hsector,H)
    call delete_sector(Hsector,Hs)
    Hsector=0
    Hstatus=.false.
    if(spH0%status)call sp_delete_matrix(spH0)
    if(spH0up%status)call sp_delete_matrix(spH0up)
    if(spH0dw%status)call sp_delete_matrix(spH0dw)
    if(spH0nl%status)call sp_delete_matrix(spH0nl)
    !
  end subroutine delete_Hv_sector











  !####################################################################
  !             BUILD SPARSE HAMILTONIAN of the SECTOR
  !####################################################################
  subroutine ed_buildH_c(isector,Hmat,dryrun)
    integer                                :: isector   
    complex(8),dimension(:,:),optional     :: Hmat
    logical,optional                       :: dryrun
    !
    complex(8),dimension(:,:),allocatable  :: Hredux,Htmp_up,Htmp_dw
    integer,dimension(Nlevels)             :: ib
    integer,dimension(Ns)                  :: nup,ndw
    integer                                :: i,iup,idw
    integer                                :: j,jup,jdw
    integer                                :: m,mup,mdw
    integer                                :: ms
    integer                                :: impi
    integer                                :: iorb,jorb,ispin,jspin,ibath
    integer                                :: kp,k1,k2,k3,k4
    integer                                :: alfa,beta
    real(8)                                :: sg1,sg2,sg3,sg4
    complex(8)                             :: htmp,htmpup,htmpdw
    complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    logical                                :: Jcondition
    !
    if(.not.Hstatus)stop "ed_buildH_c ERROR: Hsector NOT set"
    isector=Hsector
    !
    if(present(Hmat))call assert_shape(Hmat,[getdim(isector), getdim(isector)],"ed_buildH_c","Hmat")
    !
    !Get diagonal hybridization
    diag_hybr=zero
    if(bath_type/="replica")then
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dcmplx(dmft_bath%v(ispin,iorb,ibath),00d0)
             enddo
          enddo
       enddo
    else
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
             enddo
          enddo
       enddo
    endif
    !
    !
    call sp_init_matrix(spH0,mpiQ + mpiR)
    call sp_init_matrix(spH0dw,DimDw)
    call sp_init_matrix(spH0up,DimUp)
    if(Jhflag)call sp_init_matrix(spH0nl,mpiQ + mpiR)
    !
    !-----------------------------------------------!
    !IMPURITY  HAMILTONIAN
    include "ED_HAMILTONIAN_MATVEC/Himp.f90"
    !
    !LOCAL INTERACTION
    include "ED_HAMILTONIAN_MATVEC/Hint.f90"
    !
    !BATH HAMILTONIAN
    include "ED_HAMILTONIAN_MATVEC/Hbath.f90"
    !
    !IMPURITY- BATH HYBRIDIZATION
    include "ED_HAMILTONIAN_MATVEC/Hhyb.f90"
    !-----------------------------------------------!
    !
    if(present(Hmat))then
       Hmat = zero
       allocate(Htmp_up(DimUp,DimUp));Htmp_up=zero
       allocate(Htmp_dw(DimDw,DimDw));Htmp_dw=zero
       !
#ifdef _MPI       
       if(MpiStatus)then
          !Dump the diagonal and global-non-local part:
          allocate(Hredux(dim,dim));Hredux=zero          
          call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
          if(Jhflag)call sp_dump_matrix(spH0nl,Hredux(first_state:last_state,:))
          call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Complex,MPI_Sum,MpiComm,MpiIerr)
          deallocate(Hredux)
          !
          !
          !Dump the DW part:
          allocate(Hredux(DimDw,DimDw));Hredux=zero
          call sp_dump_matrix(spH0dw,Hredux(:,:))
          call MPI_AllReduce(Hredux,Htmp_dw,DimDw*DimDw,MPI_Double_Complex,MPI_Sum,MpiComm,MpiIerr)
          deallocate(Hredux)
          Hmat = Hmat + kronecker_product(zeye(DimUp),Htmp_dw)
          !
          !Dump the UP part:
          call sp_dump_matrix(spH0up,Htmp_up)
          Hmat = Hmat + kronecker_product(Htmp_up,zeye(DimDw))
          !
       else
          !Dump the diagonal and global-non-local part:
          call sp_dump_matrix(spH0,Hmat)
          if(Jhflag)call sp_dump_matrix(spH0nl,Hmat)
          !Dump the UP part:
          call sp_dump_matrix(spH0up,Htmp_up)
          Hmat = Hmat + kronecker_product(Htmp_up,zeye(DimDw))
          !
          !Dump the DW part:
          call sp_dump_matrix(spH0dw,Htmp_dw)
          Hmat = Hmat + kronecker_product(zeye(DimUp),Htmp_dw)
       endif
#else
       !
       !Dump the diagonal and global-non-local part:
       call sp_dump_matrix(spH0,Hmat)
       if(Jhflag)call sp_dump_matrix(spH0nl,Hmat)
       !Dump the UP part:
       call sp_dump_matrix(spH0up,Htmp_up)
       Hmat = Hmat + kronecker_product(Htmp_up,zeye(DimDw))
       !
       !Dump the DW part:
       call sp_dump_matrix(spH0dw,Htmp_dw)
       Hmat = Hmat + kronecker_product(zeye(DimUp),Htmp_dw)
#endif

       deallocate(Htmp_up,Htmp_dw)
    endif
    !
  end subroutine ed_buildH_c














  !####################################################################
  !        SPARSE MAT-VEC PRODUCT USING STORED SPARSE MATRIX 
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: Perform the matrix-vector product H*v used in the
  ! - serial cmplx(H)*cmplx(V)
  ! - MPI cmplx(H)*cmplx(V)
  !+------------------------------------------------------------------+
  subroutine spMatVec_cc(Nloc,v,Hv)
    integer                      :: Nloc
    complex(8),dimension(Nloc)   :: v
    complex(8),dimension(Nloc)   :: Hv
    integer                      :: i,iup,idw,j
    type(sparse_element),pointer :: c
    !
    ! !
    Hv=zero
    !
    do i = 1,Nloc
       c => spH0%row(i)%root%next
       do while(associated(c))
          Hv(i) = Hv(i) + c%cval*V(c%col)    !<== V
          c => c%next
       enddo
       nullify(c)
    enddo
    !
    if(jhflag)then
       do i = 1,Nloc
          c => spH0nl%row(i)%root%next
          do while(associated(c))
             Hv(i) = Hv(i) + c%cval*V(c%col)    !<== V
             c => c%next
          enddo
          nullify(c)
       enddo
    endif
    !
    do idw=1,DimDw
       do iup=1,DimUp
          i = iup + (idw-1)*DimUp          
          c => spH0up%row(iup)%root%next
          do while(associated(c))
             j = c%col + (idw-1)*DimUp
             Hv(i) = Hv(i) + c%cval*V(j)
             c => c%next
          enddo
          nullify(c)
       enddo
    enddo
    !
    do idw=1,DimDw
       do iup=1,DimUp
          i = iup + (idw-1)*DimUp          
          c => spH0dw%row(idw)%root%next
          do while(associated(c))
             j = iup +  (c%col-1)*DimUp
             Hv(i) = Hv(i) + c%cval*V(j)
             c => c%next
          enddo
          nullify(c)
       enddo
    enddo
    !
  end subroutine spMatVec_cc



#ifdef _MPI
  subroutine spMatVec_mpi_cc(Nloc,v,Hv)
    integer                             :: Nloc
    complex(8),dimension(Nloc)          :: v
    complex(8),dimension(Nloc)          :: Hv
    integer                             :: i
    integer                             :: N
    complex(8),dimension(:),allocatable :: vin
    integer,allocatable,dimension(:)    :: Counts,Displs
    type(sparse_element),pointer        :: c
    integer                             :: iup,idw,j
    integer                             :: impi_up,impi_dw,impi
    !
    if(MpiComm==MPI_UNDEFINED)stop "spMatVec_mpi_cc ERROR: MpiComm = MPI_UNDEFINED"
    MpiSize = get_Size_MPI(MpiComm)
    MpiRank = get_Rank_MPI(MpiComm)
    !
    N=0
    call MPI_AllReduce(Nloc,N,1,MPI_Integer,MPI_Sum,MpiComm,MpiIerr)
    if(N/=Dim)stop "spMatVec_mpi_cc ERROR: N != Dim"
    !
    !
    Hv=zero
    !
    !LOCAL PART:
    do i=1,Nloc
       c => spH0%row(i)%root%next
       do while(associated(c))
          Hv(i) = Hv(i) + c%cval*V(c%col-ishift)    !<== V
          c => c%next
       enddo
       nullify(c)
    enddo
    !
    !
    allocate(Counts(0:MpiSize-1))
    Counts(0:)        = mpiQ
    Counts(MpiSize-1) = mpiQ+mod(N,MpiSize)
    !
    allocate(Displs(0:MpiSize-1))
    forall(i=0:MpiSize-1)Displs(i)=i*mpiQ
    !
    allocate(vin(N)) ; vin = zero
    call MPI_Allgatherv(v(1:Nloc),Nloc,MPI_Double_Complex,Vin,Counts,Displs,MPI_Double_Complex,MpiComm,MpiIerr)
    !
    !NON-LOCAL PART: including S-E and P-H terms.
    if(jhflag)then
       do i=1,Nloc
          c => spH0nl%row(i)%root%next
          do while(associated(c))
             Hv(i) = Hv(i) + c%cval*Vin(c%col)    !<== Vin
             c => c%next
          enddo
          nullify(c)
       enddo
    endif
    !
    !
    !NON-LOCAL PART:
    do idw=map_first_state_dw(1),map_last_state_dw(1)
       do iup=map_first_state_up(idw),map_last_state_up(idw)
          i    = iup + (idw-1)*DimUp
          !
          impi    = i   - ishift
          !
          c => spH0up%row(iup)%root%next
          do while(associated(c))
             j = c%col + (idw-1)*DimUp
             Hv(impi) = Hv(impi) + c%cval*Vin(j) !<== Vin
             c => c%next
          enddo
          nullify(c)
          !
          c => spH0dw%row(idw)%root%next
          do while(associated(c))
             j = iup +  (c%col-1)*DimUp
             Hv(impi) = Hv(impi) + c%cval*Vin(j)
             c => c%next
          enddo
          nullify(c)
       enddo
    enddo
  end subroutine spMatVec_mpi_cc
#endif












  !####################################################################
  !            SPARSE MAT-VEC DIRECT ON-THE-FLY PRODUCT 
  !####################################################################
  subroutine directMatVec_cc(Nloc,vin,Hv)
    integer                                :: Nloc
    complex(8),dimension(Nloc)             :: vin
    complex(8),dimension(Nloc)             :: Hv
    integer                                :: isector
    integer,dimension(Ns)                  :: nup,ndw
    integer                                :: i,iup,idw
    integer                                :: j,jup,jdw
    integer                                :: m,mup,mdw
    integer                                :: ishift
    integer                                :: ms
    integer                                :: impi
    integer                                :: iorb,jorb,ispin,jspin,ibath
    integer                                :: kp,k1,k2,k3,k4
    integer                                :: alfa,beta
    real(8)                                :: sg1,sg2,sg3,sg4
    complex(8)                             :: htmp,htmpup,htmpdw
    complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    logical                                :: Jcondition
    !
    if(.not.Hstatus)stop "directMatVec_cc ERROR: Hsector NOT set"
    isector=Hsector
    !
    if(Nloc/=getdim(isector))stop "directMatVec_cc ERROR: Nloc != dim(isector)"
    !
    !Get diagonal hybridization
    diag_hybr=zero
    if(bath_type/="replica")then
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dcmplx(dmft_bath%v(ispin,iorb,ibath),00d0)
             enddo
          enddo
       enddo
    else
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
             enddo
          enddo
       enddo
    endif
    !
    !
    Hv=zero
    !
    !-----------------------------------------------!
    !IMPURITY  HAMILTONIAN
    include "ED_HAMILTONIAN_MATVEC/HxVimp.f90"
    !
    !LOCAL INTERACTION
    include "ED_HAMILTONIAN_MATVEC/HxVint.f90"
    !
    !BATH HAMILTONIAN
    include "ED_HAMILTONIAN_MATVEC/HxVbath.f90"
    !
    !IMPURITY- BATH HYBRIDIZATION
    include "ED_HAMILTONIAN_MATVEC/HxVhyb.f90"
    !-----------------------------------------------!
    !
  end subroutine directMatVec_cc





#ifdef _MPI
  subroutine directMatVec_MPI_cc(Nloc,v,Hv)
    integer                                :: Nloc
    complex(8),dimension(Nloc)             :: v
    complex(8),dimension(Nloc)             :: Hv
    integer                                :: N
    complex(8),dimension(:),allocatable    :: vin
    integer,allocatable,dimension(:)       :: Counts,Displs
    integer                                :: isector
    integer,dimension(Ns)                  :: nup,ndw
    integer                                :: i,iup,idw
    integer                                :: j,jup,jdw
    integer                                :: m,mup,mdw
    integer                                :: ms
    integer                                :: impi
    integer                                :: iorb,jorb,ispin,jspin,ibath
    integer                                :: kp,k1,k2,k3,k4
    integer                                :: alfa,beta
    real(8)                                :: sg1,sg2,sg3,sg4
    complex(8)                             :: htmp,htmpup,htmpdw
    complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    logical                                :: Jcondition
    !
    if(.not.Hstatus)stop "directMatVec_cc ERROR: Hsector NOT set"
    isector=Hsector
    !
    !Get diagonal hybridization
    diag_hybr=zero
    if(bath_type/="replica")then
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dcmplx(dmft_bath%v(ispin,iorb,ibath),00d0)
             enddo
          enddo
       enddo
    else
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
             enddo
          enddo
       enddo
    endif
    !
    if(MpiComm==MPI_UNDEFINED)stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
    !
    N=0
    call MPI_AllReduce(Nloc,N,1,MPI_Integer,MPI_Sum,MpiComm,MpiIerr)
    if(N/=Dim)stop "directMatVec_MPI_cc ERROR: N != dim(isector)"
    !
    !
    allocate(Counts(0:MpiSize-1))
    Counts(0:)        = mpiQ
    Counts(MpiSize-1) = mpiQ+mod(N,MpiSize)
    !
    allocate(Displs(0:MpiSize-1))
    forall(i=0:MpiSize-1)Displs(i)=i*mpiQ
    !
    allocate(vin(N)) ; vin = zero
    call MPI_Allgatherv(v(1:Nloc),Nloc,MPI_Double_Complex,Vin,Counts,Displs,MPI_Double_Complex,MpiComm,MpiIerr)
    !
    Hv=zero
    !
    !-----------------------------------------------!
    !IMPURITY  HAMILTONIAN
    include "ED_HAMILTONIAN_MATVEC/HxVimp.f90"
    !
    !LOCAL INTERACTION
    include "ED_HAMILTONIAN_MATVEC/HxVint.f90"
    !
    !BATH HAMILTONIAN
    include "ED_HAMILTONIAN_MATVEC/HxVbath.f90"
    !
    !IMPURITY- BATH HYBRIDIZATION
    include "ED_HAMILTONIAN_MATVEC/HxVhyb.f90"
    !-----------------------------------------------!
    !
  end subroutine directMatVec_MPI_cc
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
