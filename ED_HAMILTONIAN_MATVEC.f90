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
  integer                          :: MpiIerr
  integer                          :: MpiRank=0
  integer                          :: MpiSize=1
  integer                          :: mpiQ=1
  integer                          :: mpiQup=1
  integer                          :: mpiQdw=1
  integer                          :: mpiR=0
  integer                          :: mpiRup=0
  integer                          :: mpiRdw=0
  !
  integer                          :: first_state,last_state
  integer                          :: first_state_dw,last_state_dw
  integer,allocatable,dimension(:) :: first_state_up,last_state_up
  !  
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


  subroutine build_Hv_sector(isector,Hmat)
    integer                            :: isector,SectorDim,irank
    complex(8),dimension(:,:),optional :: Hmat
    !
    Hsector=isector
    Hstatus=.true.
    !
    call build_sector(isector,H)
    call build_sector(isector,Hs)
    !
    SectorDim=getDim(isector)
    if(present(Hmat))then
       if(any( shape(Hmat) /= [SectorDim,SectorDim]))&
            stop "setup_Hv_sector ERROR: size(Hmat) != SectorDim**2"
       call ed_buildH_c(isector,Hmat)
       return
    endif
    !
    select case (ed_sparse_H)
    case (.true.)
       !if(ed_sparse_type="ell")call ed_buildH_c(isector,dryrun=.true.)
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
    if(spH0%status)then
       ! if(MpiStatus)then          
       !    call sp_delete_matrix(MpiComm,spH0)
       ! else
       call sp_delete_matrix(spH0)
       ! endif
    endif
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
    integer                                :: dim,dimUp,dimDw
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
    !
    Dim  = getDim(isector)
    DimDw  = getDimDw(isector)
    DimUp  = getDimUp(isector)


    !FULL:
    mpiQ = dim/MpiSize
    mpiR = 0
    if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
    call sp_init_matrix(spH0,mpiQ + mpiR)
    ishift      = MpiRank*mpiQ
    first_state = MpiRank*mpiQ + 1
    last_state  = (MpiRank+1)*mpiQ + mpiR
    !
    !DW
    call sp_init_matrix(spH0dw,DimDw)
    !
    !UP:
    call sp_init_matrix(spH0up,DimUp)

    !Select the MPI states from the tensor product:
    include "ED_HAMILTONIAN_MATVEC/Build_Up_MPI_Range.f90"
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
    integer                      :: Dim,DimUp,DimDw
    type(sparse_element),pointer :: c
    !
    Dim   = getDim(Hsector)
    DimDw = getDimDw(Hsector)
    DimUp = getDimUp(Hsector)
    !
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
    integer                             :: Dim,DimUp,DimDw,DimUp_,DimDw_
    integer                             :: first_state,last_state
    integer                             :: first_up,last_up
    integer                             :: first_dw,last_dw
    integer                             :: ishift,ishift_up,ishift_dw
    integer                             :: impi_up,impi_dw,impi
    !
    Dim   = getDim(Hsector)
    DimDw = getDimDw(Hsector)
    DimUp = getDimUp(Hsector)
    !
    if(MpiComm==MPI_UNDEFINED)stop "ud_spHtimesV_cc ERRROR: MpiComm = MPI_UNDEFINED"
    !
    N=0
    call MPI_AllReduce(Nloc,N,1,MPI_Integer,MPI_Sum,MpiComm,MpiIerr)
    !
    MpiSize = get_Size_MPI(MpiComm)
    MpiRank = get_Rank_MPI(MpiComm)
    !
    if(N/=Dim)stop "ud_spHtimesV_cc ERRRO: N != Dim"
    !FULL:
    mpiQ = Dim/MpiSize
    mpiR = 0
    if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
    ishift      = MpiRank*mpiQ
    first_state = MpiRank*mpiQ + 1
    last_state  = (MpiRank+1)*mpiQ + mpiR
    !
    !Select the MPI states from the tensor product:
    include "ED_HAMILTONIAN_MATVEC/Build_Up_MPI_Range.f90"
    !
    allocate(vin(N))
    vin = zero
    !
    allocate(Counts(0:MpiSize-1),Displs(0:MpiSize-1))
    Counts(0:)        = mpiQ
    Counts(MpiSize-1) = mpiQ+mod(N,MpiSize)
    forall(i=0:MpiSize-1)Displs(i)=i*mpiQ
    call MPI_Allgatherv(v(1:Nloc),Nloc,MPI_Double_Complex,Vin,Counts,Displs,MPI_Double_Complex,MpiComm,MpiIerr)
    !
    Hv=zero
    !
    !LOCAL PART:
    !note that this is not completely local because it includes possible S-E and P-H terms.
    !otherwise it should be just a vector of dim=DimUp*DimDw
    do i=1,Nloc
       c => spH0%row(i)%root%next
       do while(associated(c))
          Hv(i) = Hv(i) + c%cval*Vin(c%col)    !<== V
          c => c%next
       enddo
       nullify(c)
    enddo
    !
    !
    !NON-LOCAL PART:
    do idw=first_state_dw,last_state_dw
       do iup=first_state_up(idw),last_state_up(idw)
          i    = iup + (idw-1)*DimUp
          !
          impi    = i   - ishift
          !
          c => spH0up%row(iup)%root%next
          do while(associated(c))
             j = c%col + (idw-1)*DimUp
             Hv(impi) = Hv(impi) + c%cval*Vin(j)
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
    integer                                :: dim,dimUp,dimDw
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
    Dim  = getDim(isector)
    DimDw  = getDimDw(isector)
    DimUp  = getDimUp(isector)
    !
    !FULL:
    mpiQ = dim/MpiSize
    mpiR = 0
    if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
    ishift      = MpiRank*mpiQ
    first_state = MpiRank*mpiQ + 1
    last_state  = (MpiRank+1)*mpiQ + mpiR
    !
    !Select the MPI states from the tensor product:
    include "ED_HAMILTONIAN_MATVEC/Build_Up_MPI_Range.f90"
    !
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
    integer                                :: dim,dimUp,dimDw
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
    Dim   = getDim(Hsector)
    DimDw = getDimDw(Hsector)
    DimUp = getDimUp(Hsector)
    !
    N=0
    call MPI_AllReduce(Nloc,N,1,MPI_Integer,MPI_Sum,MpiComm,MpiIerr)
    !
    if(N/=Dim)stop "directMatVec_MPI_cc ERROR: N != dim(isector)"
    !
    !FULL:
    mpiQ = Dim/MpiSize
    mpiR = 0
    if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
    ishift      = MpiRank*mpiQ
    first_state = MpiRank*mpiQ + 1
    last_state  = (MpiRank+1)*mpiQ + mpiR
    !
    !Select the MPI states from the tensor product:
    include "ED_HAMILTONIAN_MATVEC/Build_Up_MPI_Range.f90"
    !
    allocate(vin(N))
    vin = zero
    !
    allocate(Counts(0:MpiSize-1),Displs(0:MpiSize-1))
    Counts(0:)        = mpiQ
    Counts(MpiSize-1) = mpiQ+mod(N,MpiSize)
    forall(i=0:MpiSize-1)Displs(i)=i*mpiQ
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






