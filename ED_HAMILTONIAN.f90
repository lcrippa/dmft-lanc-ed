MODULE ED_HAMILTONIAN
  USE ED_HAMILTONIAN_COMMON
  USE ED_HAMILTONIAN_SPARSE_HxV
  USE ED_HAMILTONIAN_DIRECT_HxV
  !
  implicit none
  private


  !>Related auxiliary routines:
  public  :: ed_hamiltonian_matvec_set_MPI
  public  :: ed_hamiltonian_matvec_del_MPI


  !>Build sparse hamiltonian of the sector
  public  :: vecDim_Hv_sector
  public  :: build_Hv_sector
  public  :: delete_Hv_sector


  !>Sparse Mat-Vec product using stored sparse matrix
  public  :: spMatVec_main
  public  :: spMatVec_orbs
#ifdef _MPI
  public  :: spMatVec_MPI_main
  public  :: spMatVec_MPI_orbs
#endif


  !>Sparse Mat-Vec direct on-the-fly product 
  public  :: directMatVec_main
  public  :: directMatVec_orbs
#ifdef _MPI
  public  :: directMatVec_MPI_main
  public  :: directMatVec_MPI_orbs
#endif

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
       if(ed_total_ud)then
          spHtimesV_p => null()
          call ed_buildh_main(isector,Hmat)          
       else
          spHtimesV_p => null()
          call ed_buildh_orbs(isector,Hmat)
       end if
       return
    endif
    !
    select case (ed_sparse_H)
    case (.true.)
       if(ed_total_ud)then
          spHtimesV_p => spMatVec_main
#ifdef _MPI
          if(MpiStatus)spHtimesV_p => spMatVec_MPI_main
#endif
          call ed_buildh_main(isector)
       else
          spHtimesV_p => spMatVec_orbs
#ifdef _MPI
          if(MpiStatus)spHtimesV_p => spMatVec_MPI_orbs
#endif
          call ed_buildh_orbs(isector)
       endif
    case (.false.)
       if(ed_total_ud)then
          spHtimesV_p => directMatVec_main
#ifdef _MPI
          if(MpiStatus)spHtimesV_p => directMatVec_MPI_main
#endif
       else
          spHtimesV_p => directMatVec_orbs
#ifdef _MPI
          if(MpiStatus)spHtimesV_p => directMatVec_MPI_orbs
#endif
       endif
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
       if(spH0d%status)call sp_delete_matrix(MpiComm,spH0d)
    else
       if(spH0d%status)call sp_delete_matrix(spH0d)
    endif
#else
    if(spH0d%status)call sp_delete_matrix(spH0d)
#endif
    do iud=1,Ns_Ud
       if(spH0ups(iud)%status)call sp_delete_matrix(spH0ups(iud))
       if(spH0dws(iud)%status)call sp_delete_matrix(spH0dws(iud))
    enddo
    !
    spHtimesV_p => null()
    !
  end subroutine delete_Hv_sector




end MODULE ED_HAMILTONIAN
