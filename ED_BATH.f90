MODULE ED_BATH
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,txtfy
  USE SF_LINALG, only: eye,inv
  USE SF_ARRAYS, only:linspace
  USE SF_MISC, only: assert_shape
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  implicit none

  private



  !##################################################################
  !
  !     USER BATH ROUTINES:
  !
  !##################################################################

  !Interface for user bath I/O operations: get,set,copy
  interface get_component_bath
     module procedure get_full_component_bath
     module procedure get_spin_component_bath
     module procedure get_spin_orb_component_bath
  end interface get_component_bath

  interface set_component_bath
     module procedure set_full_component_bath
     module procedure set_spin_component_bath
     module procedure set_spin_orb_component_bath
  end interface set_component_bath

  interface copy_component_bath
     module procedure copy_full_component_bath
     module procedure copy_spin_component_bath
     module procedure copy_spin_orb_component_bath
  end interface copy_component_bath

  interface get_bath_dimension
     module procedure ::  get_bath_dimension_direct
     module procedure ::  get_bath_dimension_symmetries
  end interface get_bath_dimension


  !explicit symmetries:
  interface break_symmetry_bath
     module procedure break_symmetry_bath_site
     module procedure break_symmetry_bath_lattice
  end interface break_symmetry_bath

  interface spin_symmetrize_bath
     module procedure spin_symmetrize_bath_site
     module procedure spin_symmetrize_bath_lattice
  end interface spin_symmetrize_bath

  interface orb_symmetrize_bath
     module procedure orb_symmetrize_bath_site
     module procedure orb_symmetrize_bath_lattice
  end interface orb_symmetrize_bath


  interface orb_equality_bath
     module procedure orb_equality_bath_site
     module procedure orb_equality_bath_lattice
  end interface orb_equality_bath


  interface ph_trans_bath
     module procedure ph_trans_bath_site
     module procedure ph_trans_bath_lattice
  end interface ph_trans_bath


  interface ph_symmetrize_bath
     module procedure ph_symmetrize_bath_site
     module procedure ph_symmetrize_bath_lattice
  end interface ph_symmetrize_bath


  interface is_identity
     module procedure ::  is_identity_so
     module procedure ::  is_identity_nn
  end interface is_identity

  interface is_diagonal
     module procedure ::  is_diagonal_so
     module procedure ::  is_diagonal_nn
  end interface is_diagonal



  public  :: get_bath_dimension
  public  :: check_bath_dimension
  !
  public  :: get_component_bath_dimension
  public  :: get_spin_component_bath_dimension
  public  :: get_orb_component_bath_dimension
  public  :: get_spin_orb_component_bath_dimension
  !
  public  :: get_component_bath
  public  :: set_component_bath
  public  :: copy_component_bath
  !
  public  :: break_symmetry_bath              
  public  :: spin_symmetrize_bath
  public  :: orb_symmetrize_bath
  public  :: orb_equality_bath
  public  :: ph_trans_bath
  public  :: ph_symmetrize_bath


  !##################################################################
  !
  !     DMFT BATH ROUTINES:
  !
  !##################################################################
  !PUBLIC   = transparent to the final user
  !INTERNAL = opaque to the user but available for internal use in the code.
  !
  !DMFT BATH procedures:
  public  :: allocate_dmft_bath               !INTERNAL (for effective_bath)
  public  :: deallocate_dmft_bath             !INTERNAL (for effective_bath)
  public  :: init_dmft_bath                   !INTERNAL (for effective_bath)
  public  :: write_dmft_bath                  !INTERNAL (for effective_bath)
  public  :: save_dmft_bath                   !INTERNAL (for effective_bath)
  public  :: set_dmft_bath                    !INTERNAL (for effective_bath)
  public  :: get_dmft_bath                    !INTERNAL (for effective_bath)
  public  :: bath_from_sym                    !INTERNAL (for effective_bath)
  public  :: mask_hloc

  integer :: ibath,ilat,iorb



contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if a matrix is the identity
  !+-------------------------------------------------------------------+
  function is_identity_nn(mnnn) result(flag)
    real(8),dimension(nspin,nspin,norb,norb) :: mnnn
    real(8),dimension(nspin*norb,nspin*norb) :: mtmp
    integer                                  :: i,j
    logical                                  :: flag
    !
    flag=.true.
    !
    mtmp=nn2so_reshape(mnnn,nspin,norb)
    !
    do i=1,Nspin*Norb-1
       if((mtmp(i,i).ne.mtmp(i+1,i+1)).or.(mtmp(i,i).lt.1.d-6))flag=.false.
    enddo
    !
    do i=1,Nspin*Norb
       do j=1,Nspin*Norb
          if((i.ne.j).and.(mtmp(i,j).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_identity_nn

  function is_identity_so(mlso) result(flag)
    real(8),dimension(nspin*norb,nspin*norb) :: mlso
    real(8),dimension(nspin*norb,nspin*norb) :: mtmp
    integer                                  :: i,j
    logical                                  :: flag
    !
    flag=.true.
    !
    mtmp=mlso
    !
    do i=1,Nspin*Norb-1
       if((mtmp(i,i).ne.mtmp(i+1,i+1)).or.(mtmp(i,i).lt.1.d-6))flag=.false.
    enddo
    !
    do i=1,Nspin*Norb
       do j=1,Nspin*Norb
          if((i.ne.j).and.(mtmp(i,j).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_identity_so



  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if a matrix is diagonal
  !+-------------------------------------------------------------------+
  function is_diagonal_nn(mnnn) result(flag)
    real(8),dimension(nspin,nspin,norb,norb) :: mnnn
    real(8),dimension(nspin*norb,nspin*norb) :: mtmp
    integer                                  :: i,j
    logical                                  :: flag
    !
    flag=.true.
    !
    mtmp=abs((nn2so_reshape(mnnn,nspin,norb)))
    !
    do i=1,Nspin*Norb
       do j=1,Nspin*Norb
          if((i.ne.j).and.(mtmp(i,j).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_diagonal_nn

  function is_diagonal_so(mlso) result(flag)
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: mlso
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: mtmp
    integer                                  :: i,j
    logical                                  :: flag
    !
    flag=.true.
    !
    mtmp=abs((mlso))
    !
    do i=1,Nspin*Norb
       do j=1,Nspin*Norb
          if((i.ne.j).and.(mtmp(i,j).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_diagonal_so



  !+-------------------------------------------------------------------+
  !PURPOSE  : Create bath mask
  !+-------------------------------------------------------------------+
  function mask_hloc(hloc,wdiag,uplo) result(Hmask)
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hloc
    logical,optional                         :: wdiag,uplo
    logical                                  :: wdiag_,uplo_
    logical,dimension(Nspin,Nspin,Norb,Norb) :: Hmask
    integer                                  :: iorb,jorb,ispin,io,jo
    !
    wdiag_=.false.;if(present(wdiag))wdiag_=wdiag
    uplo_ =.false.;if(present(uplo))  uplo_=uplo
    !
    Hmask=.false.
    where(abs(Hloc)>1d-6)Hmask=.true.
    !
    !
    if(wdiag_)then
       do ispin=1,Nspin
          do iorb=1,Norb
             Hmask(ispin,ispin,iorb,iorb)=.true.
          enddo
       enddo
    endif
    !
    if(uplo_)then
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = index_stride_so(ispin,iorb)
                jo = index_stride_so(ispin,jorb)
                if(io>jo)Hmask(ispin,ispin,iorb,jorb)=.false.
             enddo
          enddo
       enddo
    endif
    !
  end function mask_hloc



  !##################################################################
  !
  !     USER BATH ROUTINES:
  !
  !##################################################################
  include 'ED_BATH/user_aux.f90'

  include 'ED_BATH/user_ctrl.f90'

  !##################################################################
  !
  !     DMFT BATH ROUTINES:
  !
  !##################################################################
  include 'ED_BATH/dmft_aux.f90'



END MODULE ED_BATH
