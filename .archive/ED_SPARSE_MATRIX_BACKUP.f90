!>>> ACTHUNG: IN THIS VERSION THE LOCAL (in the memory) PART IS DISABLED <<<<
MODULE ED_SPARSE_MATRIX  !THIS VERSION CONTAINS ONLY DBLE ELEMENT: (SYMMETRIC MATRIX) 
  USE SF_IOTOOLS, only: str,free_unit
#ifdef _MPI
  USE SF_MPI
  USE MPI
#endif
  implicit none
  private


  !SPARSE MATRIX: LINKED LIST FORMAT:
  type sparse_element_ll
     real(8)                                  :: cval !value of the entry: double precision
     integer                                  :: col  !col connected to this compress value
     type(sparse_element_ll),pointer          :: next !link to next entry in the row
  end type sparse_element_ll

  type sparse_row_ll
     integer                                  :: size=0    !size of the list
     integer                                  :: min_column= huge(1)
     integer                                  :: max_column=-huge(1)
     type(sparse_element_ll),pointer          :: root    !head/root of the list\== list itself
  end type sparse_row_ll

  type sparse_matrix_ll
     integer                                  :: Nrow
     integer                                  :: Ncol
     logical                                  :: status=.false.
     type(sparse_row_ll),dimension(:),pointer :: row
#ifdef _MPI
     logical                                  :: mpi=.false.
     integer                                  :: istart=0 !global start index for MPI storage
     integer                                  :: iend=0
     integer                                  :: ishift=0
#endif
  end type sparse_matrix_ll




  !INIT SPARSE MATRICES 
  interface sp_init_matrix
     module procedure :: sp_init_matrix_ll
#ifdef _MPI
     module procedure :: mpi_sp_init_matrix_ll
#endif
  end interface sp_init_matrix



  !DELETE SPARSE MATRIX 
  interface sp_delete_matrix
     module procedure :: sp_delete_matrix_ll
#ifdef _MPI
     module procedure :: mpi_sp_delete_matrix_ll
#endif
  end interface sp_delete_matrix


  !INSERT ELEMENTS
  interface sp_insert_element
     module procedure :: sp_insert_element_ll
#ifdef _MPI
     module procedure :: mpi_sp_insert_element_ll
#endif
  end interface sp_insert_element



  !DUMP SPARSE MATRIX INTO STANDARD MATRIX
  interface sp_dump_matrix
     module procedure :: sp_dump_matrix_ll
#ifdef _MPI
     module procedure :: mpi_sp_dump_matrix_ll
#endif
  end interface sp_dump_matrix


  !PRETTY PRINTING
  interface sp_print_matrix
     module procedure :: sp_print_matrix_ll
#ifdef _MPI
     module procedure :: mpi_sp_print_matrix_ll
#endif
  end interface sp_print_matrix


  !SPY PRINT SPARSE MATRIX
  interface sp_spy_matrix
     module procedure :: sp_spy_matrix_ll
#ifdef _MPI
     module procedure :: mpi_sp_spy_matrix_ll
#endif
  end interface sp_spy_matrix


  !Linked-List Sparse Matrix
  public :: sparse_element_ll
  public :: sparse_matrix_ll

  public :: sp_init_matrix      !init the sparse matrix   !checked
  public :: sp_delete_matrix    !delete the sparse matrix !checked
  public :: sp_insert_element   !insert an element        !checked
  public :: sp_dump_matrix      !dump sparse into array   !checked
  public :: sp_print_matrix     !print sparse             !checked
  public :: sp_spy_matrix       !
#ifdef _MPI
  public :: sp_set_mpi_ll
#endif



  integer :: MpiRank=0
  integer :: MpiSize=1
  integer :: MpiIerr
  logical :: MpiMaster=.true.






contains       


  !+------------------------------------------------------------------+
  !PURPOSE:  initialize the sparse matrix list
  !+------------------------------------------------------------------+
  subroutine sp_init_matrix_ll(sparse,N,N1)
    type(sparse_matrix_ll),intent(inout) :: sparse
    integer                           :: N
    integer,optional                  :: N1
    integer                           :: i
    !put here a delete statement to avoid problems
    if(sparse%status)stop "sp_init_matrix LL: alreay allocate can not init"
    sparse%Nrow=N
    sparse%Ncol=N
    if(present(N1))sparse%Ncol=N1
    !
    allocate(sparse%row(N))
    do i=1,N
       allocate(sparse%row(i)%root)
       sparse%row(i)%root%next => null()
       sparse%row(i)%size=0
       sparse%row(i)%min_column=huge(1)
       sparse%row(i)%max_column=-huge(1)
    end do
    !
    sparse%status=.true.
    !
  end subroutine sp_init_matrix_ll



#ifdef _MPI
  subroutine mpi_sp_init_matrix_ll(MpiComm,sparse,N,N1)
    integer                              :: MpiComm
    type(sparse_matrix_ll),intent(inout) :: sparse
    integer                              :: N
    integer,optional                     :: N1
    integer                              :: i,Ncol,Nloc
    !
    call sp_test_matrix_mpi(MpiComm,sparse,"mpi_sp_init_matrix_ll")
    !
    Ncol = N
    if(present(N1))Ncol=N1
    Nloc = sparse%iend-sparse%istart+1
    call sp_init_matrix_ll(sparse,Nloc,Ncol)
  end subroutine mpi_sp_init_matrix_ll
#endif







  !+------------------------------------------------------------------+
  !PURPOSE: delete an entire sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_delete_matrix_ll(sparse)    
    type(sparse_matrix_ll),intent(inout) :: sparse
    integer                           :: i
    type(sparse_row_ll),pointer          :: row
    type(sparse_element_ll),pointer      :: p,c
    if(.not.sparse%status)stop "Warning SPARSE/sp_delete_matrix: sparse not allocated already."
    do i=1,sparse%Nrow
       row=>sparse%row(i)
       do
          p => row%root
          c => p%next
          if(.not.associated(c))exit  !empty list
          p%next => c%next !
          c%next=>null()
          deallocate(c)
       end do
       sparse%row(i)%size=0
       sparse%row(i)%min_column=huge(1)
       sparse%row(i)%max_column=-huge(1)
       deallocate(sparse%row(i)%root)
    enddo
    deallocate(sparse%row)
    !
    sparse%Nrow=0
    sparse%Ncol=0
    sparse%status=.false.
#ifdef _MPI
    sparse%istart = 0
    sparse%iend   = 0
    sparse%ishift = 0
    sparse%mpi    = .false.
#endif 
  end subroutine sp_delete_matrix_ll


#ifdef _MPI
  subroutine mpi_sp_delete_matrix_ll(MpiComm,sparse)
    integer                              :: MpiComm
    type(sparse_matrix_ll),intent(inout) :: sparse
    integer                              :: i
    type(sparse_row_ll),pointer          :: row
    type(sparse_element_ll),pointer      :: p,c
    if(.not.sparse%status)stop "Error SPARSE/mpi_sp_delete_matrix: sparse is not allocated."
    do i=1,sparse%Nrow
       row=>sparse%row(i)
       do
          p => row%root
          c => p%next
          if(.not.associated(c))exit  !empty list
          p%next => c%next !
          c%next=>null()
          deallocate(c)
       end do
       sparse%row(i)%size=0
       sparse%row(i)%min_column=huge(1)
       sparse%row(i)%max_column=-huge(1)
       deallocate(sparse%row(i)%root)
    enddo
    deallocate(sparse%row)
    !
    sparse%istart=0
    sparse%iend=0
    sparse%ishift=0
    sparse%mpi=.false.
    !
    sparse%Nrow=0
    sparse%Ncol=0
    sparse%status=.false.
  end subroutine mpi_sp_delete_matrix_ll
#endif    










  !+------------------------------------------------------------------+
  !PURPOSE: insert an element value at position (i,j) in the sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_insert_element_ll(sparse,value,i,j)
    type(sparse_matrix_ll),intent(inout) :: sparse
    real(8),intent(in)                   :: value
    integer,intent(in)                   :: i,j
    type(sparse_row_ll),pointer          :: row
    integer                              :: column
    type(sparse_element_ll),pointer      :: p,c
    logical                              :: iadd
    !
    column = j
    !
    row => sparse%row(i)
    !
    p => row%root
    c => p%next
    iadd = .false.                !check if column already exist
    do                            !traverse the list
       if(.not.associated(c))exit !empty list or end of the list
       if(c%col == column)then
          iadd=.true.
          exit
       endif
       !if(c%col > column)exit
       if(column <= c%col)exit
       p => c
       c => c%next
    end do
    if(iadd)then
       c%cval=c%cval + value
    else
       allocate(p%next)                !Create a new element in the list
       p%next%cval= value
       p%next%col = column
       row%size   = row%size+1
       if(column<row%min_column)row%min_column=column
       if(column>row%max_column)row%max_column=column
       if(.not.associated(c))then !end of the list special case (current=>current%next)
          p%next%next  => null()
       else
          p%next%next  => c      !the %next of the new node come to current
       end if
    endif
  end subroutine sp_insert_element_ll



#ifdef _MPI
  subroutine mpi_sp_insert_element_ll(MpiComm,sparse,value,i,j)
    integer                              :: MpiComm
    type(sparse_matrix_ll),intent(inout) :: sparse
    real(8),intent(in)                   :: value
    integer,intent(in)                   :: i,j
    type(sparse_row_ll),pointer          :: row
    integer                              :: column
    type(sparse_element_ll),pointer      :: p,c
    logical                              :: iadd
    !
    call sp_test_matrix_mpi(MpiComm,sparse," mpi_sp_insert_element_ll")
    !
    column = j
    !
    row => sparse%row(i-sparse%Ishift)
    !
    p => row%root
    c => p%next
    iadd = .false.                !check if column already exist
    do                            !traverse the list
       if(.not.associated(c))exit !empty list or end of the list
       if(c%col == column)then
          iadd=.true.
          exit
       endif
       !if(c%col > column)exit
       if(column <= c%col)exit
       p => c
       c => c%next
    end do
    if(iadd)then
       c%cval=c%cval + value
    else
       allocate(p%next)                !Create a new element in the list
       p%next%cval= value
       p%next%col = column
       row%size   = row%size+1
       if(column<row%min_column)row%min_column=column
       if(column>row%max_column)row%max_column=column
       if(.not.associated(c))then !end of the list special case (current=>current%next)
          p%next%next  => null()
       else
          p%next%next  => c      !the %next of the new node come to current
       end if
    endif
  end subroutine mpi_sp_insert_element_ll
#endif

    !
    !
    !












  !+------------------------------------------------------------------+
  !PURPOSE: dump a sparse matrix into a regular 2dim array
  !+------------------------------------------------------------------+
  subroutine sp_dump_matrix_ll(sparse,matrix)
    type(sparse_matrix_ll),intent(in)    :: sparse
    real(8),dimension(:,:),intent(inout) :: matrix
    type(sparse_element_ll),pointer      :: c
    integer                              :: i,Ndim1,Ndim2
    !
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !
    if(sparse%Nrow/=Ndim1 .OR. sparse%Ncol/=Ndim2)stop "Warning SPARSE/dump_matrix: dimensions error"
    !
    matrix=0.d0
    do i=1,Ndim1
       c => sparse%row(i)%root%next
       do while(associated(c))
          matrix(i,c%col) = matrix(i,c%col) + c%cval
          c => c%next  !traverse list
       enddo
    enddo
  end subroutine sp_dump_matrix_ll


  !
#ifdef _MPI
  subroutine mpi_sp_dump_matrix_ll(MpiComm,sparse,matrix)
    integer                              :: MpiComm
    type(sparse_matrix_ll),intent(in)    :: sparse
    real(8),dimension(:,:),intent(inout) :: matrix
    type(sparse_element_ll),pointer      :: c
    integer                              :: i,N1_,N2_,Ndim1,Ndim2,Nrow,Ncol
    !
    call sp_test_matrix_mpi(MpiComm,sparse," mpi_sp_dump_matrix_ll")
    !
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !
    N1_=sparse%Nrow
    N2_=sparse%Ncol
    call MPI_AllReduce(N1_,Nrow,1,MPI_Integer,MPI_SUM,MpiComm,MpiIerr)
    call MPI_AllReduce(N2_,Ncol,1,MPI_Integer,MPI_MAX,MpiComm,MpiIerr)
    !
    if(Nrow>Ndim1 .OR. Ncol>Ndim2)stop "Warning SPARSE/mpi_dump_matrix: dimensions error"
    ! !
    matrix=0d0
    do i=sparse%Istart,sparse%Iend
       c => sparse%row(i-sparse%Ishift)%root%next
       do while(associated(c))
          matrix(i,c%col) = matrix(i,c%col) + c%cval
          c => c%next  !traverse list
       enddo
    enddo
  end subroutine mpi_sp_dump_matrix_ll
#endif





  !+------------------------------------------------------------------+
  !PURPOSE: pretty print a sparse matrix on a given unit using format fmt
  !+------------------------------------------------------------------+
  subroutine sp_print_matrix_ll(sparse,unit,fmt)
    type(sparse_matrix_ll)          :: sparse
    integer,optional             :: unit
    integer                      :: unit_
    integer                      :: i,j,Ns
    character(len=*),optional    :: fmt
    character(len=64)            :: fmt_
    type(sparse_row_ll),pointer     :: row
    type(sparse_element_ll),pointer :: c
    integer                      :: count=0
    unit_=6;if(present(unit))unit_=unit
    fmt_='F15.9';if(present(fmt))fmt_=fmt
    write(*,*)"Print sparse matrix (compact mode) ->",unit_
    do i=1,sparse%Nrow
       row => sparse%row(i)
       c => row%root%next   !assume is associated,ie list exists
       do
          if(.not.associated(c))exit
          count=count+1
          write(unit_,"("//trim(fmt_)//",A1,I0,3X)",advance='no')c%cval,',',c%col
          c => c%next  !traverse list
       end do
       write(unit_,*)
    enddo
    write(unit_,*)
  end subroutine sp_print_matrix_ll




#ifdef _MPI
  subroutine mpi_sp_print_matrix_ll(MpiComm,sparse,unit,fmt)
    integer                         :: MpiComm
    type(sparse_matrix_ll)          :: sparse
    integer,optional                :: unit
    integer                         :: unit_
    integer                         :: i,j,Ns
    character(len=*),optional       :: fmt
    character(len=64)               :: fmt_
    type(sparse_row_ll),pointer     :: row
    type(sparse_element_ll),pointer :: c
    integer                         :: count=0
    unit_=6;if(present(unit))unit_=unit
    fmt_='F15.9';if(present(fmt))fmt_=fmt
    write(*,*)"Print sparse matrix (compact mode) ->",unit_
    do i=1,sparse%Nrow
       row => sparse%row(i)
       c => row%root%next   !assume is associated,ie list exists
       do
          if(.not.associated(c))exit
          count=count+1
          write(unit_,"("//trim(fmt_)//",A1,I0,3X)",advance='no')c%cval,',',c%col
          c => c%next  !traverse list
       end do
       write(unit_,*)
    enddo
    write(unit_,*)
  end subroutine mpi_sp_print_matrix_ll
#endif








  !+------------------------------------------------------------------+
  !PURPOSE: pretty print a sparse matrix on a given unit using format fmt
  !+------------------------------------------------------------------+  
  subroutine sp_spy_matrix_ll(sparse,header)
    type(sparse_matrix_ll)          :: sparse
    character ( len = * )           :: header
    integer                         :: N1,N2
    type(sparse_element_ll),pointer :: c
    character ( len = 255 )         :: command_filename
    integer                         :: command_unit
    character ( len = 255 )         :: data_filename
    integer                         :: data_unit
    integer                         :: i, j
    character ( len = 6 )           :: n1_s,n2_s,n1_i,n2_i
    integer                         :: nz_num
    character ( len = 255 )         :: png_filename
    !
    !  Create data file.
    !
    !
    N1 = sparse%Nrow
    N2 = sparse%Ncol
    data_filename = trim ( header ) // '_data.dat'
    open (unit=free_unit(data_unit), file = data_filename, status = 'replace' )
    nz_num = 0
    do i=1,N1
       c => sparse%row(i)%root%next
       do while(associated(c))
          write(data_unit,'(2x,i6,2x,i6)') c%col,i
          nz_num = nz_num + 1
          c => c%next  !traverse list
       enddo
    enddo
    close(data_unit)
    !
    !  Create command file.
    !
    command_filename = "plot_"//str(header)//'_commands.gp'
    open(unit = free_unit(command_unit), file = command_filename, status = 'replace' )
    write(command_unit,'(a)') '#unset key'
    write(command_unit,'(a)') 'set terminal postscript eps enhanced color font "Times-Roman,16"'
    write(command_unit,'(a)') 'set output "|ps2pdf -sEPSCrop - '//str(header)//".pdf"//'"'
    write(command_unit,'(a)') 'set size ratio -1'
    write(command_unit,'(a)') 'set xlabel "<--- J --->"'
    write(command_unit,'(a)') 'set ylabel "<--- I --->"'
    write(command_unit,'(a,i6,a)')'set title "',nz_num,' nonzeros for '//str(header)//'"'
    write(command_unit,'(a)') 'set timestamp'
    write(command_unit,'(a)' )'plot [x=1:'//str(N1)//'] [y='//str(N2)//':1] "'//&
         str(data_filename)//'" w p pt 5 ps 0.4 lc rgb "red"'
    close ( unit = command_unit )
    return
  end subroutine sp_spy_matrix_ll



#ifdef _MPI
  subroutine mpi_sp_spy_matrix_ll(MpiComm,sparse,header)
    integer                         :: MpiComm
    type(sparse_matrix_ll)          :: sparse
    character ( len = * )           :: header
    integer                         :: N1,N2,N1_,N2_
    type(sparse_element_ll),pointer :: c
    character ( len = 255 )         :: command_filename
    integer                         :: command_unit
    character ( len = 255 )         :: data_filename
    integer                         :: data_unit
    integer                         :: i, j
    character ( len = 6 )           :: n1_s,n2_s,n1_i,n2_i
    integer                         :: nz_num
    character ( len = 255 )         :: png_filename
    !
    call sp_test_matrix_mpi(MpiComm,sparse," mpi_sp_spy_matrix_ll")
    !
    MpiRank  = get_Rank_MPI(MpiComm)
    !
    !  Create data file.
    !
    N1_=sparse%Nrow
    N2_=sparse%Ncol
    N1=0
    N2=0
    call MPI_AllReduce(N1_,N1,1,MPI_Integer,MPI_SUM,MpiComm,MpiIerr)
    call MPI_AllReduce(N2_,N2,1,MPI_Integer,MPI_MAX,MpiComm,MpiIerr)
    !
    nz_num = 0
    !
    data_filename = trim(header)//"_rank"//str(MpiRank,4)//'_matrix.dat'
    open(unit=free_unit(data_unit),file=data_filename, status = 'replace' )
    do i=1,sparse%Nrow
       c => sparse%row(i)%root%next
       do while(associated(c))
          write(data_unit,'(2x,i6,2x,i6)') c%col,i+sparse%Ishift
          nz_num = nz_num + 1
          c => c%next  !traverse list
       enddo
    enddo
    write(data_unit,'(2x,i6,2x,i6)')
    close(data_unit)
    !
    call MPI_Barrier(MpiComm,MpiIerr)
    !
    !  Create command file.
    !
    command_filename = "plot_"//trim(header)//"_rank"//str(MpiRank,4)//'_commands.gp'
    open(unit = free_unit(command_unit), file = command_filename, status = 'replace' )
    write(command_unit,'(a)') '#unset key'
    write(command_unit,'(a)') 'set terminal postscript eps enhanced color font "Times-Roman,16"'
    write(command_unit,'(a)') 'set output "|ps2pdf -sEPSCrop - '//str(header)//"_rank"//str(MpiRank,4)//".pdf"//'"'
    write(command_unit,'(a)') 'set size ratio -1'
    write(command_unit,'(a)') 'set xlabel "<--- J --->"'
    write(command_unit,'(a)') 'set ylabel "<--- I --->"'
    write(command_unit,'(a,i6,a)' ) &
         'set title "',nz_num,' nonzeros for '//str(header)//"_rank"//str(MpiRank,4)//'"'
    write(command_unit,'(a)') 'set timestamp'
    write(command_unit,'(a)' )'plot [x=1:'//str(N1)//'] [y='//str(N2)//':1] "'//&
         str(data_filename)//'" w p pt 5 ps 0.4 lc rgb "red"'
    close ( unit = command_unit )
    return
  end subroutine mpi_sp_spy_matrix_ll
#endif







#ifdef _MPI
  subroutine sp_set_mpi_ll(MpiComm,sparse,istart,iend,ishift)
    integer                              :: MpiComm
    type(sparse_matrix_ll),intent(inout) :: sparse
    integer                              :: istart,iend,ishift
    sparse%istart = istart
    sparse%iend   = iend
    sparse%ishift = ishift
    sparse%mpi    = .true.
  end subroutine sp_set_mpi_ll

  subroutine sp_test_matrix_mpi(MpiComm,sparse,text)
    integer                              :: MpiComm
    type(sparse_matrix_ll),intent(in)    :: sparse
    character(len=*)                     :: text
    integer                              :: MpiRank
    MpiRank = get_Rank_MPI(MpiComm)
    if(.not.sparse%mpi)then
       print*,"Rank, Error in "//trim(text)//": mpi no set"
       stop
    endif
  end subroutine sp_test_matrix_mpi
#endif


end module ED_SPARSE_MATRIX
