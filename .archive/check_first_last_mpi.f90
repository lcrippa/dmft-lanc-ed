

subroutine check_first_last(nup,ndw)
  integer :: nup,ndw
  integer :: DimUp,DimDw,Dim
  integer :: isector
  integer :: i,iup,idw
  integer :: j,jup,jdw
  integer :: m,mup,mdw
  integer :: ishift,ishift_up,ishift_dw
  integer :: impi,impi_up,impi_dw
  integer :: first_state,last_state
  integer :: first_state_up,last_state_up
  integer :: first_state_dw,last_state_dw
  integer,allocatable :: first_up(:),last_up(:)


  isector = getSector(nup,ndw)
  print*,isector
  DimUp=getDimUp(isector)
  DimDw=getDimDw(isector)
  print*,DimUp,DimDw
  Dim=getDim(isector)
  print*,Dim,DimUp*DimDw

  allocate(first_up(DimDw),last_up(DimDw))

  !FULL:
  mpiQ = dim/MpiSize
  mpiR = 0
  if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
  ishift      = MpiRank*mpiQ
  first_state = MpiRank*mpiQ + 1
  last_state  = (MpiRank+1)*mpiQ + mpiR
  !
  !DW
  mpiQdw = dimDw/MpiSize
  mpiRdw = 0
  if(MpiRank==(MpiSize-1))mpiRdw=mod(dimDw,MpiSize)
  ishift_dw      = MpiRank*mpiQdw
  first_state_dw = MpiRank*mpiQdw + 1
  last_state_dw  = (MpiRank+1)*mpiQdw + mpiRdw
  !
  !UP:
  mpiQup = dimUp/MpiSize
  mpiRup = 0
  if(MpiRank==(MpiSize-1))mpiRup=mod(dimUp,MpiSize)
  ishift_up      = 0
  first_state_up = 1
  last_state_up  = DimUp

  do i=0,mpiSize-1
     if(mpiRank==i)then
        print*,"First-Last   :",first_state,last_state
        print*,"First-Last DW:",first_state_dw,last_state_dw
        print*,"First-Last UP:",first_state_up,last_state_up
        print*,"First-Last UP:",MpiRank*mpiQup + 1,(MpiRank+1)*mpiQup + mpiRup
        print*,""
        call MPI_Barrier(MpiComm,MpiIerr)
     endif
  enddo

  first_state_dw=1000
  last_state_dw=0
  !
  first_up=1000
  last_up=0
  do i=first_state,last_state
     idw=idw_index(i,DimUp)
     iup=iup_index(i,DimUp)
     if(idw < first_state_dw)first_state_dw=idw
     if(idw > last_state_dw) last_state_dw =idw
     !
     if(iup < first_up(idw))first_up(idw)=iup
     if(iup > last_up(idw)) last_up(idw) =iup
     !
     write(200+mpiRank,*)i,iup,idw,iup + (idw-1)*DimUp
  enddo
  print*,"First-Last DW:",first_state_dw,last_state_dw
  call MPI_Barrier(MpiComm,MpiIerr)    
  do idw=first_state_dw,last_state_dw
     print*,idw,first_up(idw),last_up(idw)
  enddo
  call MPI_Barrier(MpiComm,MpiIerr)    
  do idw=first_state_dw,last_state_dw
     do iup=first_up(idw),last_up(idw)
        i = iup + (idw-1)*DimUp
        write(100+mpiRank,*)i,i-ishift,ishift
     enddo
  enddo
  call MPI_Barrier(MpiComm,MpiIerr)


  stop


end subroutine check_first_last


