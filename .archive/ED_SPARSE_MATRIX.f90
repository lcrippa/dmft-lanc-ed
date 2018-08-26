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
     integer                                  :: istart=0 !global start index for MPI storage
     integer                                  :: iend=0
     integer                                  :: ishift=0
     integer                                  :: Q
     integer                                  :: R
     integer                                  :: Chunk
#endif
  end type sparse_matrix_ll


  !SPARSE MATRIX: ELL STORAGE
  type sparse_row_ell
     integer                                   :: Size
     integer                                   :: Count
     real(8),dimension(:),allocatable          :: vals
     integer,dimension(:),allocatable          :: cols
  end type sparse_row_ell

  type sparse_matrix_ell
     logical                                   :: status=.false.
     integer                                   :: Nrow
     integer                                   :: Ncol
     type(sparse_row_ell),dimension(:),pointer :: Row
#ifdef _MPI
     integer                                   :: istart=0 !global start index for MPI storage
     integer                                   :: iend=0
     integer                                   :: ishift=0
     integer                                   :: Q
     integer                                   :: R
     integer                                   :: Chunk
#endif     
  end type sparse_matrix_ell



  !INIT SPARSE MATRICES 
  interface sp_init_matrix
     module procedure :: sp_init_matrix_ll
     module procedure :: sp_init_matrix_ell
#ifdef _MPI
     module procedure :: mpi_sp_init_matrix_ll
     module procedure :: mpi_sp_init_matrix_ell
#endif
  end interface sp_init_matrix



  !DELETE SPARSE MATRIX 
  interface sp_delete_matrix
     module procedure :: sp_delete_matrix_ll
     module procedure :: sp_delete_matrix_ell
#ifdef _MPI
     module procedure :: mpi_sp_delete_matrix_ll
     module procedure :: mpi_sp_delete_matrix_ell
#endif
  end interface sp_delete_matrix


  !INSERT ELEMENTS
  interface sp_insert_element
     module procedure :: sp_insert_element_ll
     module procedure :: sp_insert_element_ell
#ifdef _MPI
     module procedure :: mpi_sp_insert_element_ll
     module procedure :: mpi_sp_insert_element_ell
#endif
  end interface sp_insert_element



  !DUMP SPARSE MATRIX INTO STANDARD MATRIX
  interface sp_dump_matrix
     module procedure :: sp_dump_matrix_ll
     module procedure :: sp_dump_matrix_ell
#ifdef _MPI
     module procedure :: mpi_sp_dump_matrix_ll
     module procedure :: mpi_sp_dump_matrix_ell
#endif
  end interface sp_dump_matrix


  !PRETTY PRINTING
  interface sp_print_matrix
     module procedure :: sp_print_matrix_ll
     module procedure :: sp_print_matrix_ell
#ifdef _MPI
     module procedure :: mpi_sp_print_matrix_ll
     module procedure :: mpi_sp_print_matrix_ell
#endif
  end interface sp_print_matrix


  !SPY PRINT SPARSE MATRIX
  interface sp_spy_matrix
     module procedure :: sp_spy_matrix_ll
     module procedure :: sp_spy_matrix_ell
#ifdef _MPI
     module procedure :: mpi_sp_spy_matrix_ll
     module procedure :: mpi_sp_spy_matrix_ell
#endif
  end interface sp_spy_matrix


  !Linked-List Sparse Matrix
  public  :: sparse_element_ll
  public  :: sparse_matrix_ll
  public  :: sparse_matrix_ell

  public  :: sp_init_matrix      !init the sparse matrix   !checked
  public  :: sp_delete_matrix    !delete the sparse matrix !checked
  public  :: sp_insert_element   !insert an element        !checked
  public  :: sp_dump_matrix      !dump sparse into array   !checked
  public  :: sp_print_matrix     !print sparse             !checked
  public  :: sp_spy_matrix       !


  integer :: MpiRank=0
  integer :: MpiSize=1
  integer :: MpiIerr
  logical :: MpiMaster=.true.
  integer :: MpiQ
  integer :: MpiR
  integer :: MpiChunk
  integer :: MpiIstart,MpiIend,MpiIshift






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
#ifdef _MPI
    sparse%Q      = N
    sparse%R      = 0
    sparse%Chunk  = sparse%Q + sparse%R
    sparse%istart = 1
    sparse%iend   = N
    sparse%ishift = 0
#endif
    !
  end subroutine sp_init_matrix_ll

  subroutine sp_init_matrix_ell(sparse,vecNnz,N,N1)
    type(sparse_matrix_ell) :: sparse
    integer                 :: N
    integer,dimension(N)    :: vecNnz
    integer,optional        :: N1
    integer                 :: i,Nnz
    if(sparse%status)stop "sp_init_matrix ELL: alreay allocate can not init"
    sparse%Nrow    = N
    sparse%Ncol    = N
    if(present(N1))sparse%Ncol=N1
    !
    allocate(sparse%row(N))
    do i=1,N
       Nnz = vecNnz(i)
       sparse%row(i)%Size  = Nnz
       sparse%row(i)%Count = 0
       allocate(sparse%row(i)%vals(Nnz))  ;sparse%row(i)%vals=0d0
       allocate(sparse%row(i)%cols(Nnz))  ;sparse%row(i)%cols=0
    end do
    !
    sparse%status=.true.
    !
#ifdef _MPI
    sparse%Q      = N
    sparse%R      = 0
    sparse%Chunk  = sparse%Q + sparse%R
    sparse%istart = 1
    sparse%iend   = N
    sparse%ishift = 0
#endif 
  end subroutine sp_init_matrix_ell


#ifdef _MPI
  subroutine mpi_sp_init_matrix_ll(MpiComm,sparse,N,N1)
    integer                              :: MpiComm
    type(sparse_matrix_ll),intent(inout) :: sparse
    integer                              :: N
    integer,optional                     :: N1
    integer                              :: i,Ncol
    !
    if(sparse%status)stop "mpi_sp_init_matrix LL: alreay allocate can not init"
    !
    call sp_MPI_setup(MpiComm,N)    
    !
    Ncol = N
    if(present(N1))Ncol=N1
    call sp_init_matrix_ll(sparse,MpiChunk,Ncol)
    !
    sparse%Q      = MpiQ
    sparse%R      = MpiR
    sparse%Chunk  = sparse%Q + sparse%R
    sparse%istart = MpiIstart
    sparse%iend   = MpiIend
    sparse%ishift = MpiIshift
  end subroutine mpi_sp_init_matrix_ll

  subroutine mpi_sp_init_matrix_ell(MpiComm,sparse,vecNnz,N,N1)
    integer                 :: MpiComm
    type(sparse_matrix_ell) :: sparse
    integer                 :: N
    integer                 :: vecNnz(N)
    integer,optional        :: N1
    integer                 :: i,Ncol,Nnz
    !
    if(sparse%status)stop "mpi_sp_init_matrix ELL: alreay allocate can not init"
    !
    call sp_MPI_setup(MpiComm,N)    
    !
    sparse%Nrow    = MpiChunk
    sparse%Ncol    = N
    if(present(N1))sparse%Ncol=N1
    !
    allocate(sparse%row(MpiChunk))
    do i=MpiIstart,MpiIend
       Nnz = vecNnz(i-MpiIshift)
       sparse%row(i-MpiIshift)%Size  = Nnz
       sparse%row(i-MpiIshift)%Count = 0
       allocate(sparse%row(i-MpiIshift)%vals(Nnz))  ;sparse%row(i-MpiIshift)%vals=0d0
       allocate(sparse%row(i-MpiIshift)%cols(Nnz))  ;sparse%row(i-MpiIshift)%cols=0
    end do
    !
    sparse%Q      = MpiQ
    sparse%R      = MpiR
    sparse%Chunk  = sparse%Q + sparse%R
    sparse%istart = MpiIstart
    sparse%iend   = MpiIend
    sparse%ishift = MpiIshift
    !
    sparse%status=.true.
  end subroutine mpi_sp_init_matrix_ell
#endif










  !+------------------------------------------------------------------+
  !PURPOSE: delete an entire sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_delete_matrix_ll(sparse)    
    type(sparse_matrix_ll),intent(inout) :: sparse
    integer                           :: i
    type(sparse_row_ll),pointer          :: row
    type(sparse_element_ll),pointer      :: p,c
    if(.not.sparse%status)return
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
    sparse%Q      = 0
    sparse%R      = 0
    sparse%Chunk  = 0
    sparse%istart = 0
    sparse%iend   = 0
    sparse%ishift = 0
#endif 
  end subroutine sp_delete_matrix_ll

  subroutine sp_delete_matrix_ell(sparse)    
    type(sparse_matrix_ell),intent(inout) :: sparse
    integer                               :: i
    if(.not.sparse%status)return
    !
    do i=1,sparse%Nrow
       deallocate(sparse%row(i)%vals)
       deallocate(sparse%row(i)%cols)
       sparse%row(i)%Size  = 0
       sparse%row(i)%Count = 0
    enddo
    deallocate(sparse%row)
    !
    sparse%Nrow=0
    sparse%Ncol=0
    sparse%status=.false.
    !
#ifdef _MPI
    sparse%Q      = 0
    sparse%R      = 0
    sparse%Chunk  = 0
    sparse%istart = 0
    sparse%iend   = 0
    sparse%ishift = 0
#endif
  end subroutine sp_delete_matrix_ell

#ifdef _MPI
  subroutine mpi_sp_delete_matrix_ll(MpiComm,sparse)
    integer                           :: MpiComm
    type(sparse_matrix_ll),intent(inout) :: sparse
    integer                           :: i
    type(sparse_row_ll),pointer          :: row
    type(sparse_element_ll),pointer      :: p,c
    if(.not.sparse%status)return
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
    !
    sparse%Q      = 0
    sparse%R      = 0
    sparse%Chunk  = 0
    sparse%istart = 0
    sparse%iend   = 0
    sparse%ishift = 0
  end subroutine mpi_sp_delete_matrix_ll

  subroutine mpi_sp_delete_matrix_ell(MpiComm,sparse)
    integer                               :: MpiComm
    type(sparse_matrix_ell),intent(inout) :: sparse
    integer                               :: i
    if(.not.sparse%status)return
    !
    do i=1,sparse%Nrow
       deallocate(sparse%row(i)%vals)
       deallocate(sparse%row(i)%cols)
       sparse%row(i)%Size  = 0
       sparse%row(i)%count = 0
    enddo
    deallocate(sparse%row)
    !
    sparse%Nrow=0
    sparse%Ncol=0
    sparse%status=.false.
    !
    sparse%Q      = 0
    sparse%R      = 0
    sparse%Chunk  = 0
    sparse%istart = 0
    sparse%iend   = 0
    sparse%ishift = 0
  end subroutine mpi_sp_delete_matrix_ell
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

  subroutine sp_insert_element_ell(sparse,value,i,j)
    type(sparse_matrix_ell),intent(inout) :: sparse
    real(8),intent(in)                    :: value
    integer,intent(in)                    :: i,j
    type(sparse_row_ell),pointer          :: row
    integer                               :: column,pos
    logical                               :: iadd
    !
    column = j
    !
    row => sparse%row(i)
    iadd = .false.                !check if column already exist
    if(any(row%cols == column))iadd=.true.
    !
    if(iadd)then                !this column exists so just sum it up
       !find the position  column in %cols 
       do pos=1,size(row%cols)
          if(row%cols(pos)==column)exit
       end do
       !add up value to the current one in %vals
       row%vals(pos)=row%vals(pos) + value
    else                        !this column is new. increase counter and store it 
       row%Count=row%Count+1
       pos = row%Count
       if(pos > row%Size) stop "sp_insert_element_ell ERROR: pos > row%Size"
       row%vals(pos) = value
       row%cols(pos) = column
    endif
  end subroutine sp_insert_element_ell


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
    column = j
    !
    ! if(column>=sparse%Istart.AND.column<=sparse%Iend)then
    !    ! row => sparse%loc(i)
    !    row => sparse%loc(i-sparse%Ishift)
    ! else
    ! row => sparse%row(i)
    row => sparse%row(i-sparse%Ishift)
    ! endif
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

  subroutine mpi_sp_insert_element_ell(MpiComm,sparse,value,i,j)
    integer                               :: MpiComm
    type(sparse_matrix_ell),intent(inout) :: sparse
    real(8),intent(in)                    :: value
    integer,intent(in)                    :: i,j
    type(sparse_row_ell),pointer          :: row
    integer                               :: column,pos
    logical                               :: iadd
    !
    column = j
    !
    ! if(column>=sparse%Istart.AND.column<=sparse%Iend)then
    !    ! row => sparse%loc(i)
    !    row => sparse%loc(i-sparse%Ishift)
    ! else
    ! row => sparse%row(i)
    row => sparse%row(i-sparse%Ishift)
    ! endif
    iadd = .false.                !check if column already exist
    if(any(row%cols == column))iadd=.true.
    !
    if(iadd)then                !this column exists so just sum it up
       !find the position  column in %cols 
       do pos=1,size(row%cols)
          if(row%cols(pos)==column)exit
       end do
       !add up value to the current one in %vals
       row%vals(pos)=row%vals(pos) + value
    else                        !this column is new. increase counter and store it 
       row%Count=row%Count+1
       pos = row%Count
       if(pos > row%Size) stop "sp_insert_element_ell ERROR: pos > row%Size"
       row%vals(pos) = value
       row%cols(pos) = column
    endif
  end subroutine mpi_sp_insert_element_ell
#endif












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

  subroutine sp_dump_matrix_ell(sparse,matrix)
    type(sparse_matrix_ell),intent(in)   :: sparse
    real(8),dimension(:,:),intent(inout) :: matrix
    integer                              :: i,j,Ndim1,Ndim2
    !
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !
    if(sparse%Nrow/=Ndim1 .OR. sparse%Ncol/=Ndim2)stop "Warning SPARSE/dump_matrix: dimensions error"
    !
    matrix=0.d0
    do i=1,Ndim1
       do j=1,sparse%row(i)%Size
          matrix(i,sparse%row(i)%cols(j)) = matrix(i,sparse%row(i)%cols(j)) + sparse%row(i)%vals(j)
       enddo
    enddo
  end subroutine sp_dump_matrix_ell


  !
#ifdef _MPI
  subroutine mpi_sp_dump_matrix_ll(MpiComm,sparse,matrix)
    integer                                          :: MpiComm
    type(sparse_matrix_ll),intent(in)                :: sparse
    real(8),dimension(:,:),intent(inout)             :: matrix
    real(8),dimension(size(matrix,1),size(matrix,2)) :: mtmp
    type(sparse_element_ll),pointer                  :: c
    integer                                          :: i,N1_,N2_,Ndim1,Ndim2,Nrow,Ncol
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
    !
    call sp_MPI_setup(MpiComm,Ndim1)
    !
    mtmp=0d0
    do i=MpiIstart,MpiIend
       c => sparse%row(i-MpiIshift)%root%next
       do while(associated(c))
          mtmp(i,c%col) = mtmp(i,c%col) + c%cval
          c => c%next  !traverse list
       enddo
    enddo
    !
    call MPI_AllReduce(Mtmp,Matrix,Ndim1*Ndim2,MPI_Double_Precision,MPI_Sum,MpiComm,MpiIerr)
  end subroutine mpi_sp_dump_matrix_ll

  subroutine mpi_sp_dump_matrix_ell(MpiComm,sparse,matrix)
    integer                                          :: MpiComm
    type(sparse_matrix_ell),intent(in)               :: sparse
    real(8),dimension(:,:),intent(inout)             :: matrix
    real(8),dimension(size(matrix,1),size(matrix,2)) :: mtmp
    integer                                          :: i,j,N1_,N2_,Ndim1,Ndim2,Nrow,Ncol
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
    !
    call sp_MPI_setup(MpiComm,Ndim1)
    !
    mtmp=0d0
    do i=MpiIstart,MpiIend
       do j=1,sparse%row(i-mpiIshift)%Size
          mtmp(i,sparse%row(i-mpiIshift)%cols(j))=mtmp(i,sparse%row(i-mpiIshift)%cols(j))+sparse%row(i-mpiIshift)%vals(j)
       enddo
    enddo
    !
    call MPI_AllReduce(Mtmp,Matrix,Ndim1*Ndim2,MPI_Double_Precision,MPI_Sum,MpiComm,MpiIerr)
  end subroutine mpi_sp_dump_matrix_ell
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
    fmt_='2F15.9';if(present(fmt))fmt_=fmt
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

  subroutine sp_print_matrix_ell(sparse,unit,fmt)
    type(sparse_matrix_ell)      :: sparse
    integer,optional             :: unit
    integer                      :: unit_
    integer                      :: i,j,Ns
    character(len=*),optional    :: fmt
    character(len=64)            :: fmt_
    type(sparse_row_ell),pointer :: row
    integer                      :: count=0
    unit_=6;if(present(unit))unit_=unit
    fmt_='2F15.9';if(present(fmt))fmt_=fmt
    write(*,*)"Print sparse matrix (compact mode) ->",unit_
    do i=1,sparse%Nrow
       row => sparse%row(i)
       do j=1,row%Size
          write(unit_,"("//trim(fmt_)//",A1,I0,3X)",advance='no')row%vals(j),',',row%cols(j)
       end do
       write(unit_,*)
    enddo
    write(unit_,*)
  end subroutine sp_print_matrix_ell



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
    fmt_='2F15.9';if(present(fmt))fmt_=fmt
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

  subroutine mpi_sp_print_matrix_ell(MpiComm,sparse,unit,fmt)
    integer                      :: MpiComm
    type(sparse_matrix_ell)      :: sparse
    integer,optional             :: unit
    integer                      :: unit_
    integer                      :: i,j,Ns
    character(len=*),optional    :: fmt
    character(len=64)            :: fmt_
    type(sparse_row_ell),pointer :: row
    integer                      :: count=0
    unit_=6;if(present(unit))unit_=unit
    fmt_='2F15.9';if(present(fmt))fmt_=fmt
    write(*,*)"Print sparse matrix (compact mode) ->",unit_
    do i=1,sparse%Nrow
       row => sparse%row(i)
       do j=1,row%Size
          write(unit_,"("//trim(fmt_)//",A1,I0,3X)",advance='no')row%vals(j),',',row%cols(j)
       end do
       write(unit_,*)
    enddo
    write(unit_,*)
  end subroutine mpi_sp_print_matrix_ell
#endif








  !+------------------------------------------------------------------+
  !PURPOSE: pretty print a sparse matrix on a given unit using format fmt
  !+------------------------------------------------------------------+  
  subroutine sp_spy_matrix_ll(sparse,header)
    type(sparse_matrix_ll)          :: sparse
    character ( len = * )        :: header
    integer                      :: N1,N2
    type(sparse_element_ll),pointer :: c
    character ( len = 255 )      :: command_filename
    integer                      :: command_unit
    character ( len = 255 )      :: data_filename
    integer                      :: data_unit
    integer                      :: i, j
    character ( len = 6 )        :: n1_s,n2_s,n1_i,n2_i
    integer                      :: nz_num
    character ( len = 255 )      :: png_filename
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

  subroutine sp_spy_matrix_ell(sparse,header)
    type(sparse_matrix_ell)         :: sparse
    character ( len = * )           :: header
    integer                         :: N1,N2
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
       do j=1,sparse%row(i)%size
          write(data_unit,'(2x,i6,2x,i6)') sparse%row(i)%cols(j),i
          nz_num = nz_num + 1
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
  end subroutine sp_spy_matrix_ell


#ifdef _MPI
  subroutine mpi_sp_spy_matrix_ll(MpiComm,sparse,header)
    integer                      :: MpiComm
    type(sparse_matrix_ll)          :: sparse
    character ( len = * )        :: header
    integer                      :: N1,N2,N1_,N2_
    type(sparse_element_ll),pointer :: c
    character ( len = 255 )      :: command_filename
    integer                      :: command_unit
    character ( len = 255 )      :: data_filename(2)
    integer                      :: data_unit
    integer                      :: i, j
    character ( len = 6 )        :: n1_s,n2_s,n1_i,n2_i
    integer                      :: nz_num
    character ( len = 255 )      :: png_filename
    !
    !  Create data file.
    !
    N1_=sparse%Nrow
    N2_=sparse%Ncol
    call MPI_AllReduce(N1_,N1,1,MPI_Integer,MPI_SUM,MpiComm,MpiIerr)
    call MPI_AllReduce(N2_,N2,1,MPI_Integer,MPI_MAX,MpiComm,MpiIerr)
    !
    call sp_MPI_setup(MpiComm,N1)
    !
    nz_num = 0
    !
    data_filename(1) = trim(header)//"_rank"//str(MpiRank,4)//'_matrix.dat'
    open(unit=free_unit(data_unit),file=data_filename(1), status = 'replace' )
    do i=1,MpiChunk
       c => sparse%row(i)%root%next
       do while(associated(c))
          write(data_unit,'(2x,i6,2x,i6)') c%col,i+MpiIshift
          nz_num = nz_num + 1
          c => c%next  !traverse list
       enddo
    enddo
    write(data_unit,'(2x,i6,2x,i6)')
    close(data_unit)
    !
    ! data_filename(2) = trim(header)//"_rank"//str(MpiRank,4)//'_local.dat'
    ! open(unit=free_unit(data_unit),file=data_filename(2), status = 'replace' )
    ! ! do i=1,MpiChunk
    ! !    c => sparse%loc(i)%root%next
    ! !    do while(associated(c))
    ! !       write(data_unit,'(2x,i6,2x,i6)') c%col,i+MpiIshift
    ! !       nz_num = nz_num + 1
    ! !       c => c%next  !traverse list
    ! !    enddo
    ! ! enddo
    ! write(data_unit,'(2x,i6,2x,i6)')
    ! close(data_unit)
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
         str(data_filename(1))//'" w p pt 5 ps 0.4 lc rgb "red" title "Matrix", "'
    close ( unit = command_unit )
    return
  end subroutine mpi_sp_spy_matrix_ll


  subroutine mpi_sp_spy_matrix_ell(MpiComm,sparse,header)
    integer                         :: MpiComm
    type(sparse_matrix_ell)         :: sparse
    character ( len = * )           :: header
    integer                         :: N1,N2,N1_,N2_
    character ( len = 255 )         :: command_filename
    integer                         :: command_unit
    character ( len = 255 )         :: data_filename(2)
    integer                         :: data_unit
    integer                         :: i, j
    character ( len = 6 )           :: n1_s,n2_s,n1_i,n2_i
    integer                         :: nz_num
    character ( len = 255 )         :: png_filename
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
    call sp_MPI_setup(MpiComm,N1)
    !
    nz_num = 0
    !
    data_filename(1) = trim(header)//"_rank"//str(MpiRank,4)//'_matrix.dat'
    open(unit=free_unit(data_unit),file=data_filename(1), status = 'replace' )
    do i=1,MpiChunk
       do j=1,sparse%row(i)%Size
          write(data_unit,'(2x,i6,2x,i6)') sparse%row(i)%cols(j),i+MpiIshift
          nz_num = nz_num + 1
       enddo
    enddo
    write(data_unit,'(2x,i6,2x,i6)')
    close(data_unit)
    !
    ! data_filename(2) = trim(header)//"_rank"//str(MpiRank,4)//'_local.dat'
    ! open(unit=free_unit(data_unit),file=data_filename(2), status = 'replace' )
    ! ! do i=1,MpiChunk
    ! !    do j=1,sparse%loc(i)%Size
    ! !       write(data_unit,'(2x,i6,2x,i6)') sparse%loc(i)%cols(j),i+MpiIshift
    ! !       nz_num = nz_num + 1
    ! !    enddo
    ! ! enddo
    ! write(data_unit,'(2x,i6,2x,i6)')
    ! close(data_unit)
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
         str(data_filename(1))//'" w p pt 5 ps 0.4 lc rgb "red" title "Matrix", "'! //&
    ! str(data_filename(2))//'" w p pt 5 ps 0.4 lc rgb "blue" title "Local"'
    close ( unit = command_unit )
    return
  end subroutine mpi_sp_spy_matrix_ell
#endif













#ifdef _MPI
  subroutine sp_MPI_setup(MpiComm,N,bool)
    integer :: MpiComm,N,irank
    logical,optional :: bool
    !
    MpiRank   = get_rank_MPI(MpiComm)
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    MpiQ      = get_Q_MPI(MpiComm,N)
    MpiR      = get_R_MPI(MpiComm,N)
    MpiChunk  = MpiQ+MpiR
    MpiIstart = 1+MpiRank*MpiQ
    MpiIend   = (MpiRank+1)*MpiQ+MpiR
    MpiIshift = MpiRank*MpiQ
    if(present(bool).AND.bool)then
       if(MpiMaster)print*,"         mpiRank,   mpi_Q,   mpi_R,   mpi_CHunk,   mpi_Istart,   mpi_Iend,   mpi_Iend-mpi_Istart"
       do irank=0,MpiSize-1
          if(MpiRank==irank)then
             print*,MpiRank,MpiQ,MpiR,MpiChunk,MpiIstart,MpiIend,MpiIend-MpiIstart+1
             call MPI_Barrier(MpiComm,MpiIerr)
          endif
          call MPI_Barrier(MpiComm,MpiIerr)
       enddo
       call MPI_Barrier(MpiComm,MpiIerr)
    endif
  end subroutine sp_MPI_setup
#endif




end module ED_SPARSE_MATRIX
