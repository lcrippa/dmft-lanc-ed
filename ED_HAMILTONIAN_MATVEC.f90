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
  public  :: ed_buildH_c
  !
  !>Sparse Mat-Vec product using stored sparse matrix 
  public  :: spMatVec_cc
#ifdef _MPI
  public  :: spMatVec_MPI_cc
#endif
  !
  !>Related auxiliary routines:
  public  :: ed_hamiltonian_matvec_set_MPI
  public  :: ed_hamiltonian_matvec_del_MPI
  public  :: setup_Hv_sector
  public  :: delete_Hv_sector



  !> MPI local variables (shared)
#ifdef _MPI
  integer                      :: MpiComm=MPI_UNDEFINED
#else
  integer                      :: MpiComm=0
#endif
  logical                      :: MpiStatus=.false.
  integer                      :: MpiIerr
  integer                      :: MpiRank=0
  integer                      :: MpiSize=1
  integer                      :: mpiQ=1
  integer                      :: mpiQup=1
  integer                      :: mpiQdw=1
  integer                      :: mpiR=0
  integer                      :: mpiRup=0
  integer                      :: mpiRdw=0
  !
  integer                      :: Hsector=0
  logical                      :: Hstatus=.false.
  type(sector_map)             :: H,Hs(2)





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
    MpiQ=1
    MpiR=0
  end subroutine ed_hamiltonian_matvec_del_MPI


  subroutine setup_Hv_sector(isector)
    integer                   :: isector
    Hsector=isector
    Hstatus=.true.
    call build_sector(isector,H)
    call build_sector(isector,Hs)
  end subroutine setup_Hv_sector


  subroutine delete_Hv_sector()
    call delete_sector(Hsector,H)
    call delete_sector(Hsector,Hs)
    Hsector=0
    Hstatus=.false.
  end subroutine delete_Hv_sector








  !####################################################################
  !             BUILD SPARSE HAMILTONIAN of the SECTOR
  !####################################################################
  subroutine ed_buildH_c(Hmat)
    complex(8),dimension(:,:),optional     :: Hmat
    complex(8),dimension(:,:),allocatable  :: Hredux,Htmp_up,Htmp_dw
    integer                                :: isector
    integer,dimension(Nlevels)             :: ib
    integer,dimension(Ns)                  :: nup,ndw
    integer                                :: dim,dimUp,dimDw
    integer                                :: i,iup,idw
    integer                                :: j,jup,jdw
    integer                                :: m,mup,mdw
    integer                                :: ishift,ishift_up,ishift_dw
    integer                                :: ms
    integer                                :: impi,impi_up,impi_dw
    integer                                :: iorb,jorb,ispin,jspin,ibath
    integer                                :: kp,k1,k2,k3,k4
    integer                                :: alfa,beta
    real(8)                                :: sg1,sg2,sg3,sg4
    complex(8)                             :: htmp,htmpup,htmpdw
    complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    logical                                :: Jcondition
    integer                                :: first_state,last_state
    integer                                :: first_state_up,last_state_up
    integer                                :: first_state_dw,last_state_dw
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
    if(spH0%status)call sp_delete_matrix(spH0)
    if(spH0up%status)call sp_delete_matrix(spH0up)
    if(spH0dw%status)call sp_delete_matrix(spH0dw)
    !
    !FULL:
    dim  = getdim(isector)
    mpiQ = dim/MpiSize
    mpiR = 0
    if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
    !
    call sp_init_matrix(spH0,mpiQ + mpiR)
    ishift      = MpiRank*mpiQ
    first_state = MpiRank*mpiQ + 1
    last_state  = (MpiRank+1)*mpiQ + mpiR
    !
    !UP:
    dimUp  = getDimUp(isector)
    mpiQup = dimUp/MpiSize
    mpiRup = 0
    if(MpiRank==(MpiSize-1))mpiRup=mod(dimUp,MpiSize)
    call sp_init_matrix(spH0up,mpiQup + mpiRup)
    ishift_up      = MpiRank*mpiQup
    first_state_up = MpiRank*mpiQup + 1
    last_state_up  = (MpiRank+1)*mpiQup + mpiRup
    !
    !DW
    dimDw  = getDimDw(isector)
    mpiQdw = dimDw/MpiSize
    mpiRdw = 0
    if(MpiRank==(MpiSize-1))mpiRdw=mod(dimDw,MpiSize)
    call sp_init_matrix(spH0dw,mpiQdw + mpiRdw)
    ishift_dw      = MpiRank*mpiQdw
    first_state_dw = MpiRank*mpiQdw + 1
    last_state_dw  = (MpiRank+1)*mpiQdw + mpiRdw
    !
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
       if(MpiStatus)then
          allocate(Hredux(dim,dim));Hredux=zero
          !Dump the diagonal and global-non-local part:
          call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
#ifdef _MPI
          call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Complex,MPI_Sum,MpiComm,MpiIerr)
#endif
          deallocate(Hredux)
          !
          !Dump the UP part:
          allocate(Hredux(DimUp,DimUp));Hredux=zero
          call sp_dump_matrix(spH0up,Hredux(first_state_up:last_state_up,:))
#ifdef _MPI
          call MPI_AllReduce(Hredux,Htmp_up,DimUp*DimUp,MPI_Double_Complex,MPI_Sum,MpiComm,MpiIerr)
#endif
          Hmat = Hmat + kronecker_product(Htmp_up,zeye(DimDw))
          deallocate(Hredux)
          !
          !Dump the DW part:
          allocate(Hredux(DimDw,DimDw));Hredux=zero
          call sp_dump_matrix(spH0dw,Hredux(first_state_dw:last_state_dw,:))
#ifdef _MPI
          call MPI_AllReduce(Hredux,Htmp_dw,DimDw*DimDw,MPI_Double_Complex,MPI_Sum,MpiComm,MpiIerr)
#endif
          Hmat = Hmat + kronecker_product(zeye(DimUp),Htmp_dw)
          deallocate(Hredux)
          !
       else
          !
          !Dump the diagonal and global-non-local part:
          call sp_dump_matrix(spH0,Hmat)
          !Dump the UP part:
          call sp_dump_matrix(spH0up,Htmp_up)
          Hmat = Hmat + kronecker_product(Htmp_up,zeye(DimDw))
          !
          !Dump the DW part:
          call sp_dump_matrix(spH0dw,Htmp_dw)
          Hmat = Hmat + kronecker_product(zeye(DimUp),Htmp_dw)
       endif
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
    integer                      :: Dim,DimUp,DimDw
    type(sparse_element),pointer :: c
    !
    Dim = spH0%Nrow
    DimUp = spH0up%Nrow
    DimDw = spH0dw%Nrow
    !
    Hv=zero
    !
    do idw=1,DimDw
       do iup=1,DimUp
          i = iup + (idw-1)*DimUp
          !
          c => spH0%row(i)%root%next
          do while(associated(c))
             Hv(i) = Hv(i) + c%cval*V(c%col)    !<== V
             c => c%next
          enddo
          nullify(c)

          c => spH0up%row(iup)%root%next
          do while(associated(c))
             j = c%col + (idw-1)*DimUp
             Hv(i) = Hv(i) + c%cval*V(j)
             c => c%next
          enddo
          nullify(c)

          c => spH0dw%row(idw)%root%next
          do while(associated(c))
             j = iup +  (c%col-1)*DimUp
             Hv(i) = Hv(i) + c%cval*V(j)
             c => c%next
          enddo

       enddo
    enddo
  end subroutine spMatVec_cc



#ifdef _MPI
  subroutine spMatVec_mpi_cc(Nloc,v,Hv)
    integer                             :: Nloc
    complex(8),dimension(Nloc)          :: v
    complex(8),dimension(Nloc)          :: Hv
    integer                             :: i
    integer                             :: N
    complex(8),dimension(:),allocatable :: vin
    integer,allocatable,dimension(:)    :: SendCounts,Displs
    type(sparse_element),pointer        :: c
    integer                             :: iup,idw,j
    integer                             :: Dim,DimUp,DimDw,DimUp_,DimDw_
    N=0
    !
    if(MpiComm==MPI_UNDEFINED)stop "ud_spHtimesV_cc ERRROR: MpiComm = MPI_UNDEFINED"    
    call MPI_AllReduce(Nloc,N,1,MPI_Integer,MPI_Sum,MpiComm,MpiIerr)

    DimUp_ = spH0up%Nrow
    DimUp  = 0
    call MPI_AllReduce(DimUp_,DimUp,1,MPI_Integer,MPI_Sum,MpiComm,MpiIerr)

    DimDw_ = spH0dw%Nrow
    DimDw  = 0
    call MPI_AllReduce(DimDw_,DimDw,1,MPI_Integer,MPI_Sum,MpiComm,MpiIerr)
    !
    MpiSize = get_Size_MPI(MpiComm)
    mpiQ = get_Q_MPI(MpiComm,N)
    mpiR = get_R_MPI(MpiComm,N)
    !
    allocate(vin(N))
    allocate(SendCounts(0:MpiSize-1),displs(0:MpiSize-1))
    vin                   = zero
    SendCounts(0:)        = mpiQ
    SendCounts(MpiSize-1) = mpiQ+mod(N,MpiSize)
    forall(i=0:MpiSize-1)Displs(i)=i*mpiQ
    call MPI_Allgatherv(v(1:Nloc),Nloc,MPI_Double_Complex,vin,SendCounts,Displs,MPI_Double_Complex,MpiComm,MpiIerr)
    !
    Hv=zero
    !    
    ! print*,DimUp,DimDw,getDimUp(Hsector),getDimDw(Hsector)
    do idw=1,DimDw
       do iup=1,DimUp
          i = iup + (idw-1)*DimUp
          !
          c => spH0%row(i)%root%next
          do while(associated(c))
             Hv(i) = Hv(i) + c%cval*V(c%col)    !<== V
             c => c%next
          enddo
          nullify(c)

          c => spH0up%row(iup)%root%next
          do while(associated(c))
             j = c%col + (idw-1)*DimUp
             Hv(i) = Hv(i) + c%cval*V(j)
             c => c%next
          enddo
          nullify(c)

          c => spH0dw%row(idw)%root%next
          do while(associated(c))
             j = iup +  (c%col-1)*DimUp
             Hv(i) = Hv(i) + c%cval*V(j)
             c => c%next
          enddo

       enddo
    enddo
    nullify(c)
    !
  end subroutine spMatVec_mpi_cc
#endif












!   !####################################################################
!   !            SPARSE MAT-VEC DIRECT ON-THE-FLY PRODUCT 
!   !####################################################################
!   subroutine directMatVec_cc(Nloc,vin,Hv)
!     integer                                :: Nloc
!     complex(8),dimension(Nloc)             :: vin
!     complex(8),dimension(Nloc)             :: Hv
!     integer                                :: isector
!     integer,dimension(Nlevels)             :: ib
!     integer,dimension(Ns)                  :: ibup,ibdw
!     integer                                :: dim,dimUp,dimDw
!     integer                                :: i,iup,idw
!     integer                                :: m,mup,mdw
!     integer                                :: ishift,ishift_up,ishift_dw
!     integer                                :: j,ms,impi
!     integer                                :: iorb,jorb,ispin,jspin,ibath
!     integer                                :: kp,k1,k2,k3,k4
!     integer                                :: alfa,beta
!     real(8)                                :: sg1,sg2,sg3,sg4
!     real(8),dimension(Norb)                :: nup,ndw
!     complex(8)                             :: htmp,htmpup,htmpdw
!     complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
!     logical                                :: Jcondition
!     integer                                :: first_state,last_state
!     integer                                :: first_state_up,last_state_up
!     integer                                :: first_state_dw,last_state_dw
!     !
!     if(.not.Hstatus)stop "directMatVec_cc ERROR: Hsector NOT set"
!     isector=Hsector
!     !
!     dim=getdim(isector)
!     if(Nloc/=dim)stop "directMatVec_cc ERROR: Nloc != dim(isector)"
!     !
!     mpiQ = dim/MpiSize
!     mpiR = 0
!     if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
!     ishift      = MpiRank*mpiQ
!     first_state = MpiRank*mpiQ + 1
!     last_state  = (MpiRank+1)*mpiQ + mpiR
!     !
!     !Get diagonal hybridization
!     diag_hybr=zero
!     if(bath_type/="replica")then
!        do ibath=1,Nbath
!           do ispin=1,Nspin
!              do iorb=1,Norb
!                 diag_hybr(ispin,iorb,ibath)=dcmplx(dmft_bath%v(ispin,iorb,ibath),00d0)
!              enddo
!           enddo
!        enddo
!     else
!        do ibath=1,Nbath
!           do ispin=1,Nspin
!              do iorb=1,Norb
!                 diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
!              enddo
!           enddo
!        enddo
!     endif
!     !
!     Hv=zero
!     !-----------------------------------------------!
!     !IMPURITY  HAMILTONIAN
!     include "ED_HAMILTONIAN_MATVEC/HxVimp.f90"
!     !
!     !LOCAL INTERACTION
!     include "ED_HAMILTONIAN_MATVEC/HxVint.f90"
!     !
!     !BATH HAMILTONIAN
!     include "ED_HAMILTONIAN_MATVEC/HxVbath.f90"
!     !
!     !IMPURITY- BATH HYBRIDIZATION
!     include "ED_HAMILTONIAN_MATVEC/HxVhyb.f90"
!     !-----------------------------------------------!
!     !
!   end subroutine directMatVec_cc



! #ifdef _MPI
!   subroutine directMatVec_MPI_cc(Nloc,v,Hv)
!     integer                                :: Nloc
!     complex(8),dimension(Nloc)             :: v
!     complex(8),dimension(Nloc)             :: Hv
!     integer                                :: N
!     complex(8),dimension(:),allocatable    :: vin
!     integer,allocatable,dimension(:)       :: SendCounts,Displs
!     integer                                :: isector
!     integer,dimension(Nlevels)             :: ib
!     integer,dimension(Ns)                  :: ibup,ibdw
!     integer                                :: dim,dimUp,dimDw
!     integer                                :: i,iup,idw
!     integer                                :: m,mup,mdw
!     integer                                :: ishift,ishift_up,ishift_dw
!     integer                                :: j,ms,impi
!     integer                                :: iorb,jorb,ispin,jspin,ibath
!     integer                                :: kp,k1,k2,k3,k4
!     integer                                :: alfa,beta
!     real(8)                                :: sg1,sg2,sg3,sg4
!     real(8),dimension(Norb)                :: nup,ndw
!     complex(8)                             :: htmp,htmpup,htmpdw
!     complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
!     logical                                :: Jcondition
!     integer                                :: first_state,last_state
!     integer                                :: first_state_up,last_state_up
!     integer                                :: first_state_dw,last_state_dw
!     !
!     if(.not.Hstatus)stop "directMatVec_MPI_cc ERROR: Hsector NOT set"
!     isector=Hsector
!     !
!     dim=getdim(isector)
!     !
!     !
!     !Get diagonal hybridization
!     diag_hybr=zero
!     if(bath_type/="replica")then
!        do ibath=1,Nbath
!           do ispin=1,Nspin
!              do iorb=1,Norb
!                 diag_hybr(ispin,iorb,ibath)=dcmplx(dmft_bath%v(ispin,iorb,ibath),00d0)
!              enddo
!           enddo
!        enddo
!     else
!        do ibath=1,Nbath
!           do ispin=1,Nspin
!              do iorb=1,Norb
!                 diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
!              enddo
!           enddo
!        enddo
!     endif
!     !
!     if(MpiComm==MPI_UNDEFINED)stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
!     !
!     mpiQ = dim/MpiSize
!     mpiR = 0
!     if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
!     ishift      = MpiRank*mpiQ
!     first_state = MpiRank*mpiQ + 1
!     last_state  = (MpiRank+1)*mpiQ + mpiR
!     !
!     N=0
!     if(MpiComm==MPI_UNDEFINED)stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
!     call MPI_AllReduce(Nloc,N,1,MPI_Integer,MPI_Sum,MpiComm,MpiIerr)
!     if(N/=dim)stop "directMatVec_MPI_cc ERROR: N != dim(isector)"
!     !
!     allocate(vin(N))
!     allocate(SendCounts(0:MpiSize-1),displs(0:MpiSize-1))
!     vin                   = zero
!     SendCounts(0:)        = mpiQ
!     SendCounts(MpiSize-1) = mpiQ+mod(N,MpiSize)
!     forall(i=0:MpiSize-1)Displs(i)=i*mpiQ
!     call MPI_Allgatherv(v(1:Nloc),Nloc,MPI_Double_Complex,vin,SendCounts,Displs,MPI_Double_Complex,MpiComm,MpiIerr)
!     call MPI_Bcast(vin,N,MPI_Double_Complex,0,MpiComm,MpiIerr)
!     !
!     Hv=zero
!     !
!     !-----------------------------------------------!
!     !IMPURITY  HAMILTONIAN
!     include "ED_HAMILTONIAN_MATVEC/HxVimp.f90"
!     !
!     !LOCAL INTERACTION
!     include "ED_HAMILTONIAN_MATVEC/HxVint.f90"
!     !
!     !BATH HAMILTONIAN
!     include "ED_HAMILTONIAN_MATVEC/HxVbath.f90"
!     !
!     !IMPURITY- BATH HYBRIDIZATION
!     include "ED_HAMILTONIAN_MATVEC/HxVhyb.f90"
!     !-----------------------------------------------!
!     !
!   end subroutine directMatVec_MPI_cc
! #endif








end MODULE ED_HAMILTONIAN_MATVEC






