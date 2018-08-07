!+-------------------------------------------------------------------+
!PURPOSE  : Inquire the correct bath size to allocate the 
! the bath array in the calling program.
!
! Get size of each dimension of the component array. 
! The Result is an rank 1 integer array Ndim with dimension:
! 3 for get_component_size_bath
! 2 for get_spin_component_size_bath & get_orb_component_size_bath
! 1 for get_spin_orb_component_size_bath
!+-------------------------------------------------------------------+
function get_bath_dimension(Hloc_nn,ispin_) result(bath_size)
  complex(8),optional,intent(in) :: Hloc_nn(:,:,:,:)
  integer,optional               :: ispin_
  integer                        :: bath_size,ndx,ispin,iorb,jspin,jorb,io,jo,Maxspin
  real(8),allocatable            :: Hloc(:,:,:,:)

  select case(bath_type)
  case default
     !( e [Nspin][Norb][Nbath] + v [Nspin][Norb][Nbath] )
     bath_size = Norb*Nbath + Norb*Nbath
     if(.not.present(ispin_))bath_size=Nspin*bath_size
  case('hybrid')
     !(e [Nspin][1][Nbath] + v [Nspin][Norb][Nbath] )
     bath_size = Nbath + Norb*Nbath
     if(.not.present(ispin_))bath_size=Nspin*bath_size
  case('replica')
     !
     !off-diagonal non-vanishing elements
     if(present(Hloc_nn))then
        allocate(Hloc(Nspin,Nspin,Norb,Norb));Hloc=dreal(Hloc_nn)
     elseif(allocated(impHloc))then
        allocate(Hloc(Nspin,Nspin,Norb,Norb));Hloc=impHloc
     else
        stop "ERROR: get_bath_dimension: bath_type=replica neither Hloc_nn present nor impHloc allocated"
     endif
     !
     ndx=0
     do ispin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              io = iorb + (ispin-1)*Norb
              jo = jorb + (ispin-1)*Norb
              if(io.lt.jo)then
                 if(abs(Hloc(ispin,ispin,iorb,jorb)).gt.1d-6)ndx=ndx+1
              endif
           enddo
        enddo
     enddo
     !Real diagonal elements (always assumed)
     ndx= ndx + Nspin * Norb
     !number of non vanishing elements for each replica
     ndx = ndx * Nbath
     !diagonal hybridizations
     ndx = ndx + Nbath
     !
     bath_size = ndx
     !
  end select
end function get_bath_dimension


function get_component_bath_dimension(itype) result(Ndim)
  integer                  :: itype
  integer                  :: Ndim(3)
  call check_type_bath(itype)
  select case(bath_type)
  case default
     Ndim=[Nspin,Norb,Nbath]
  case('hybrid')
     if(itype==1)then
        Ndim=[Nspin,1,Nbath]
     else
        Ndim=[Nspin,Norb,Nbath]
     endif
  end select
end function get_component_bath_dimension
!
function get_spin_component_bath_dimension(itype) result(Ndim) 
  integer                  :: itype
  integer                  :: Ndim(2)
  call check_type_bath(itype)
  select case(bath_type)
  case default
     Ndim=[Norb,Nbath]
  case('hybrid')
     if(itype==1)then
        Ndim=[1,Nbath]
     else
        Ndim=[Norb,Nbath]
     endif
  end select
end function get_spin_component_bath_dimension
!
function get_orb_component_bath_dimension(itype) result(Ndim)
  integer                  :: itype
  integer                  :: Ndim(2)
  Ndim=[Nspin,Nbath]
end function get_orb_component_bath_dimension
!
function get_spin_orb_component_bath_dimension(itype) result(Ndim)
  integer                  :: itype
  integer                  :: Ndim
  Ndim=Nbath
end function get_spin_orb_component_bath_dimension









!+-----------------------------------------------------------------------------+!
!PURPOSE: Get a specified itype,ispin,iorb component of the user bath.
! The component is returned into an Array of rank D
! get_full_component_bath    : return the entire itype component (D=3)
! get_spin_component_bath    : return the itype component for the select ispin (D=2)
! get_spin_orb_component_bath: return the itype component for the select ispin & iorb (D=1)
!+-----------------------------------------------------------------------------+!
subroutine get_full_component_bath(array,bath_,itype)
  real(8),dimension(:,:,:) :: array
  real(8),dimension(:)     :: bath_
  integer                  :: itype
  logical                  :: check
  type(effective_bath)     :: dmft_bath_
  !
  check= check_bath_dimension(bath_)
  if(.not.check)stop "get_component_bath error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  call assert_component_size_bath(array,itype,"get_component_bath","Array")
  call check_type_bath(itype)
  if(itype==1)then
     Array = dmft_bath_%e(:,:,:)
  else
     Array = dmft_bath_%v(:,:,:)
  endif
  call deallocate_dmft_bath(dmft_bath_)
end subroutine get_full_component_bath

subroutine get_spin_component_bath(array,bath_,itype,ispin)
  real(8),dimension(:,:) :: array
  real(8),dimension(:)   :: bath_
  integer                :: itype,ispin
  logical                :: check
  type(effective_bath)   :: dmft_bath_
  !
  check= check_bath_dimension(bath_)
  if(.not.check)stop "get_spin_component_bath error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  call assert_spin_component_size_bath(array,itype,"get_spin_component_bath","Array")
  call check_type_bath(itype)
  if(ispin>Nspin)stop "get_spin_component_bath error: ispin > Nspin"
  if(itype==1)then
     Array = dmft_bath_%e(ispin,:,:)
  else
     Array = dmft_bath_%v(ispin,:,:)
  endif
  call deallocate_dmft_bath(dmft_bath_)
end subroutine get_spin_component_bath

subroutine get_spin_orb_component_bath(array,bath_,itype,ispin,iorb)
  real(8),dimension(:) :: array
  real(8),dimension(:) :: bath_
  integer              :: itype,ispin,iorb
  logical              :: check
  type(effective_bath) :: dmft_bath_
  !
  check= check_bath_dimension(bath_)
  if(.not.check)stop "get_spin_orb_component_bath error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  call assert_spin_orb_component_size_bath(array,itype,"get_spin_orb_component_bath","Array")
  call check_type_bath(itype)
  if(ispin>Nspin)stop "get_spin_orb_component_bath error: ispin > Nspin"
  if(iorb>Norb)stop "get_spin_orb_component_bath error: iorb > Norb"
  select case(bath_type)
  case default
     if(itype==1)then
        Array = dmft_bath_%e(ispin,iorb,:)
     else
        Array = dmft_bath_%v(ispin,iorb,:)
     endif
  case('hybrid')
     if(itype==1)then
        Array = dmft_bath_%e(ispin,1,:)
     else
        Array = dmft_bath_%v(ispin,iorb,:)
     endif
  end select
  call deallocate_dmft_bath(dmft_bath_)
end subroutine get_spin_orb_component_bath





!+-----------------------------------------------------------------------------+!
!PURPOSE: Set a specified itype,ispin,iorb component of the user bath.
! The component is set from an Array of rank D
! set_full_component_bath    : return the entire itype component (D=3)
! set_spin_component_bath    : return the itype component for the select ispin (D=2)
! set_spin_orb_component_bath: return the itype component for the select ispin & iorb (D=1)
!+-----------------------------------------------------------------------------+!
subroutine set_full_component_bath(array,bath_,itype)
  real(8),dimension(:,:,:) :: array
  real(8),dimension(:)     :: bath_
  integer                  :: itype
  logical                  :: check
  type(effective_bath)     :: dmft_bath_
  !
  check= check_bath_dimension(bath_)
  if(.not.check)stop "set_component_bath error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  call assert_component_size_bath(array,itype,"set_component_bath","Array")
  call check_type_bath(itype)
  if(itype==1)then
     dmft_bath_%e(:,:,:)  = Array
  else
     dmft_bath_%v(:,:,:)  = Array
  endif
  call get_dmft_bath(dmft_bath_,bath_)
  call deallocate_dmft_bath(dmft_bath_)
end subroutine set_full_component_bath

subroutine set_spin_component_bath(array,bath_,itype,ispin)
  real(8),dimension(:,:) :: array
  real(8),dimension(:)   :: bath_
  integer                :: itype,ispin
  logical                :: check
  type(effective_bath)   :: dmft_bath_
  !
  check= check_bath_dimension(bath_)
  if(.not.check)stop "set_spin_component_bath error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  call assert_spin_component_size_bath(array,itype,"set_spin_component_bath","Array")
  call check_type_bath(itype)
  if(ispin>Nspin)stop "set_spin_component_bath error: ispin > Nspin"
  if(itype==1)then
     dmft_bath_%e(ispin,:,:)  = Array
  else
     dmft_bath_%v(ispin,:,:)  = Array
  endif
  call get_dmft_bath(dmft_bath_,bath_)
  call deallocate_dmft_bath(dmft_bath_)
end subroutine set_spin_component_bath

subroutine set_spin_orb_component_bath(array,bath_,itype,ispin,iorb)
  real(8),dimension(:) :: array
  real(8),dimension(:) :: bath_
  integer              :: itype,ispin,iorb
  logical              :: check
  type(effective_bath) :: dmft_bath_
  !
  check= check_bath_dimension(bath_)
  if(.not.check)stop "set_spin_orb_component_bath error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  call assert_spin_orb_component_size_bath(array,itype,"set_spin_orb_component_bath","Array")
  call check_type_bath(itype)
  if(ispin>Nspin)stop "set_spin_orb_component_bath error: ispin > Nspin"
  if(iorb>Norb)stop "set_spin_orb_component_bath error: iorb > Norb"
  select case(bath_type)
  case default
     if(itype==1)then
        dmft_bath_%e(ispin,iorb,:)  = Array
     else
        dmft_bath_%v(ispin,iorb,:)  = Array
     endif
  case('hybrid')
     if(itype==1)then
        dmft_bath_%e(ispin,1,:)  = Array
     else
        dmft_bath_%v(ispin,iorb,:)  = Array
     endif
  end select
  call get_dmft_bath(dmft_bath_,bath_)
  call deallocate_dmft_bath(dmft_bath_)
end subroutine set_spin_orb_component_bath
!



!+-----------------------------------------------------------------------------+!
!PURPOSE: Copy a specified component of IN bath to the OUT bath.
! copy_full_component_bath    : copy the entire itype component
! copy_spin_component_bath    : copy ispin to jspin component
! copy_spin_orb_component_bath: copy ispin.iorb to jspin.jorb components
!+-----------------------------------------------------------------------------+!
subroutine copy_full_component_bath(bathIN,bathOUT,itype)
  real(8),dimension(:)     :: bathIN,bathOUT
  integer                  :: itype
  logical                  :: check
  type(effective_bath)     :: dIN,dOUT
  !
  check= check_bath_dimension(bathIN)
  if(.not.check)stop "copy_component_bath error: wrong bath dimensions IN"
  check= check_bath_dimension(bathOUT)
  if(.not.check)stop "copy_component_bath error: wrong bath dimensions OUT"
  call allocate_dmft_bath(dIN)
  call allocate_dmft_bath(dOUT)
  call set_dmft_bath(bathIN,dIN)
  call set_dmft_bath(bathOUT,dOUT)
  call check_type_bath(itype)
  if(itype==1)then
     dOUT%e(:,:,:)  = dIN%e(:,:,:)
  else
     dOUT%v(:,:,:)  = dIN%v(:,:,:)
  endif
  call get_dmft_bath(dOUT,bathOUT)
  call deallocate_dmft_bath(dIN)
  call deallocate_dmft_bath(dOUT)
end subroutine copy_full_component_bath

subroutine copy_spin_component_bath(bathIN,ispin,bathOUT,jspin,itype)
  real(8),dimension(:)     :: bathIN,bathOUT
  integer                  :: ispin,jspin
  integer,optional         :: itype
  logical                  :: check
  type(effective_bath)     :: dIN,dOUT
  !
  check= check_bath_dimension(bathIN)
  if(.not.check)stop "copy_component_bath error: wrong bath dimensions IN"
  check= check_bath_dimension(bathOUT)
  if(.not.check)stop "copy_component_bath error: wrong bath dimensions OUT"
  call allocate_dmft_bath(dIN)
  call allocate_dmft_bath(dOUT)
  call set_dmft_bath(bathIN,dIN)
  call set_dmft_bath(bathOUT,dOUT)
  if(present(itype))call check_type_bath(itype)
  if(ispin>Norb.OR.jspin>Nspin)stop "copy_spin_component_bath error: ispin/jspin > Nspin"
  if(present(itype))then          
     if(itype==1)then
        dOUT%e(jspin,:,:)  = dIN%e(ispin,:,:)
     else
        dOUT%v(jspin,:,:)  = dIN%v(ispin,:,:)
     endif
  else
     dOUT%e(jspin,:,:)  = dIN%e(ispin,:,:)
     dOUT%v(jspin,:,:)  = dIN%v(ispin,:,:)
  endif
  call get_dmft_bath(dOUT,bathOUT)
  call deallocate_dmft_bath(dIN)
  call deallocate_dmft_bath(dOUT)
end subroutine copy_spin_component_bath

subroutine copy_spin_orb_component_bath(bathIN,ispin,iorb,bathOUT,jspin,jorb,itype)
  real(8),dimension(:)     :: bathIN,bathOUT
  integer                  :: ispin,jspin
  integer                  :: iorb,jorb
  integer,optional         :: itype
  logical                  :: check
  type(effective_bath)     :: dIN,dOUT
  !
  check= check_bath_dimension(bathIN)
  if(.not.check)stop "copy_spin_orb_component_bath error: wrong bath dimensions IN"
  check= check_bath_dimension(bathOUT)
  if(.not.check)stop "copy_spin_orb_component_bath error: wrong bath dimensions OUT"
  call allocate_dmft_bath(dIN)
  call allocate_dmft_bath(dOUT)
  call set_dmft_bath(bathIN,dIN)
  call set_dmft_bath(bathOUT,dOUT)
  if(present(itype))call check_type_bath(itype)
  if(ispin>Norb.OR.jspin>Nspin)stop "copy_spin_orb_component_bath error: ispin/jspin > Nspin"
  if(iorb>Norb.OR.jorb>Norb)stop "copy_spin_orb_component_bath error: iorb/jorb > Norb"
  !
  select case(bath_type)      
  case default
     if(present(itype))then
        if(itype==1)then
           dOUT%e(jspin,jorb,:)  = dIN%e(ispin,iorb,:)
        else
           dOUT%v(jspin,jorb,:)  = dIN%v(ispin,iorb,:)
        endif
     else
        dOUT%e(jspin,jorb,:)  = dIN%e(ispin,iorb,:)
        dOUT%v(jspin,jorb,:)  = dIN%v(ispin,iorb,:)
     endif
  case('hybrid')
     if(present(itype))then
        if(itype==1)then
           dOUT%e(jspin,1,:)    = dIN%e(ispin,1,:)
        else
           dOUT%v(jspin,jorb,:) = dIN%v(ispin,iorb,:)
        endif
     else
        dOUT%e(jspin,1,:)     = dIN%e(ispin,1,:)
        dOUT%v(jspin,jorb,:)  = dIN%v(ispin,iorb,:)
     endif
  end select
  call get_dmft_bath(dOUT,bathOUT)
  call deallocate_dmft_bath(dIN)
  call deallocate_dmft_bath(dOUT)
end subroutine copy_spin_orb_component_bath
!















!##################################################################
!
!     USER BATH PREDEFINED SYMMETRIES:
!
!##################################################################

!+-------------------------------------------------------------------+
!PURPOSE  : given a bath array apply a specific transformation or 
! impose a given symmetry:
! - break spin symmetry by applying a symmetry breaking field
! - given a bath array set both spin components to have 
!    the same bath, i.e. impose non-magnetic solution
! - given a bath array enforces the particle-hole symmetry 
!    by setting the positive energies in modulo identical to the negative
!    ones.
! - given a bath enforce normal (i.e. non superconducting) solution
! - given a dmft bath pull/push the components W^{ss'}_\a(l) of the Hybridization 
!    matrix
! - given a dmft bath pull/push the nonsu2 components
!+-------------------------------------------------------------------+
subroutine break_symmetry_bath_site(bath_,field,sign,save)
  real(8),dimension(:)   :: bath_
  type(effective_bath)   :: dmft_bath_
  real(8)                :: field
  real(8)                :: sign
  logical,optional       :: save
  logical                :: save_
  save_=.true.;if(present(save))save_=save
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  dmft_bath_%e(1,:,:)    =dmft_bath_%e(1,:,:)      + sign*field
  dmft_bath_%e(Nspin,:,:)=dmft_bath_%e(Nspin,:,:)  - sign*field
  if(save_)call save_dmft_bath(dmft_bath_)
  call get_dmft_bath(dmft_bath_,bath_)
  call deallocate_dmft_bath(dmft_bath_)
end subroutine break_symmetry_bath_site
subroutine break_symmetry_bath_lattice(bath_,field,sign,save)
  real(8),dimension(:,:) :: bath_
  real(8)                :: field
  real(8)                :: sign
  logical,optional       :: save
  logical                :: save_
  integer                :: Nsites,ilat
  save_=.true.;if(present(save))save_=save
  Nsites=size(bath_,1)
  do ilat=1,Nsites
     ed_file_suffix=reg(ineq_site_suffix)//reg(txtfy(ilat,site_indx_padding))
     call break_symmetry_bath_site(bath_(ilat,:),field,sign,save_)
  enddo
  ed_file_suffix=""
end subroutine break_symmetry_bath_lattice

!---------------------------------------------------------!

subroutine spin_symmetrize_bath_site(bath_,save)
  real(8),dimension(:)   :: bath_
  type(effective_bath)   :: dmft_bath_
  logical,optional       :: save
  logical                :: save_
  save_=.true.;if(present(save))save_=save
  if(Nspin==1)then
     write(LOGfile,"(A)")"spin_symmetrize_bath: Nspin=1 nothing to symmetrize"
     return
  endif
  !
  call allocate_dmft_bath(dmft_bath_)
  if (bath_type=="replica")call init_dmft_bath_mask(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  !
  if (bath_type/="replica") then
     dmft_bath_%e(Nspin,:,:)=dmft_bath_%e(1,:,:)
     dmft_bath_%v(Nspin,:,:)=dmft_bath_%v(1,:,:)
  else
     stop "spin symmetrize not implemented for replica"
  endif
  if(save_)call save_dmft_bath(dmft_bath_)
  call get_dmft_bath(dmft_bath_,bath_)
  call deallocate_dmft_bath(dmft_bath_)
end subroutine spin_symmetrize_bath_site
subroutine spin_symmetrize_bath_lattice(bath_,save)
  real(8),dimension(:,:) :: bath_
  logical,optional       :: save
  logical                :: save_
  integer                :: Nsites,ilat
  save_=.true.;if(present(save))save_=save
  Nsites=size(bath_,1)
  do ilat=1,Nsites
     ed_file_suffix=reg(ineq_site_suffix)//reg(txtfy(ilat,site_indx_padding))
     call spin_symmetrize_bath_site(bath_(ilat,:),save_)
  enddo
  ed_file_suffix=""
end subroutine spin_symmetrize_bath_lattice

!---------------------------------------------------------!

subroutine orb_symmetrize_bath_site(bath_,save)
  real(8),dimension(:)   :: bath_
  type(effective_bath)   :: dmft_bath_
  logical,optional       :: save
  logical                :: save_
  integer                :: iorb
  real(8),allocatable    :: lvl(:,:),hyb(:,:)
  save_=.true.;if(present(save))save_=save
  if(Norb==1)then
     write(LOGfile,"(A)")"orb_symmetrize_bath: Norb=1 nothing to symmetrize"
     return
  endif
  !
  call allocate_dmft_bath(dmft_bath_)
  if (bath_type=="replica")call init_dmft_bath_mask(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  !
  if (bath_type/="replica") then
     if(allocated(lvl))deallocate(lvl);allocate(lvl(Nspin,Nbath));lvl=0d0;lvl=sum(dmft_bath_%e,dim=2)/Norb
     if(allocated(hyb))deallocate(hyb);allocate(hyb(Nspin,Nbath));hyb=0d0;hyb=sum(dmft_bath_%v,dim=2)/Norb
     do iorb=1,Norb
        dmft_bath_%e(:,iorb,:)=lvl
        dmft_bath_%v(:,iorb,:)=hyb
     enddo
  else
     stop "orb symmetrize not implemented for replica"
  endif
  if(save_)call save_dmft_bath(dmft_bath_)
  call get_dmft_bath(dmft_bath_,bath_)
  call deallocate_dmft_bath(dmft_bath_)
end subroutine orb_symmetrize_bath_site
subroutine orb_symmetrize_bath_lattice(bath_,save)
  real(8),dimension(:,:) :: bath_
  logical,optional       :: save
  logical                :: save_
  integer                :: Nsites,ilat
  save_=.true.;if(present(save))save_=save
  Nsites=size(bath_,1)
  do ilat=1,Nsites
     ed_file_suffix=reg(ineq_site_suffix)//reg(txtfy(ilat,site_indx_padding))
     call orb_symmetrize_bath_site(bath_(ilat,:),save_)
  enddo
  ed_file_suffix=""
end subroutine orb_symmetrize_bath_lattice

!---------------------------------------------------------!

subroutine ph_symmetrize_bath_site(bath_,save)
  real(8),dimension(:)   :: bath_
  type(effective_bath)   :: dmft_bath_
  integer                :: i
  logical,optional       :: save
  logical                :: save_
  save_=.true.;if(present(save))save_=save
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  if(Nbath==1)return
  if(mod(Nbath,2)==0)then
     do i=1,Nbath/2
        dmft_bath_%e(:,:,Nbath+1-i)=-dmft_bath_%e(:,:,i)
        dmft_bath_%v(:,:,Nbath+1-i)= dmft_bath_%v(:,:,i)
     enddo
  else
     do i=1,(Nbath-1)/2
        dmft_bath_%e(:,:,Nbath+1-i)=-dmft_bath_%e(:,:,i)
        dmft_bath_%v(:,:,Nbath+1-i)= dmft_bath_%v(:,:,i)
     enddo
     dmft_bath_%e(:,:,(Nbath-1)/2+1)=0.d0
  endif
  if(save_)call save_dmft_bath(dmft_bath_)
  call get_dmft_bath(dmft_bath_,bath_)
  call deallocate_dmft_bath(dmft_bath_)
end subroutine ph_symmetrize_bath_site
subroutine ph_symmetrize_bath_lattice(bath_,save)
  real(8),dimension(:,:) :: bath_
  logical,optional       :: save
  logical                :: save_
  integer                :: Nsites,ilat
  save_=.true.;if(present(save))save_=save
  Nsites=size(bath_,1)
  do ilat=1,Nsites
     ed_file_suffix=reg(ineq_site_suffix)//reg(txtfy(ilat,site_indx_padding))
     call ph_symmetrize_bath_site(bath_(ilat,:),save_)
  enddo
  ed_file_suffix=""
end subroutine ph_symmetrize_bath_lattice

!---------------------------------------------------------!

subroutine ph_trans_bath_site(bath_,save)
  real(8),dimension(:)   :: bath_
  type(effective_bath)   :: dmft_bath_
  type(effective_bath)   :: tmp_dmft_bath
  integer                :: i
  logical,optional       :: save
  logical                :: save_
  save_=.true.;if(present(save))save_=save
  call allocate_dmft_bath(dmft_bath_)
  call allocate_dmft_bath(tmp_dmft_bath)
  call set_dmft_bath(bath_,dmft_bath_)
  if(Nbath==1)return
  do i=1,Nbath
     select case(Norb)
     case default
        ! do nothing
        dmft_bath_%e(:,:,i)= dmft_bath_%e(:,:,i)
        dmft_bath_%v(:,:,i)= dmft_bath_%v(:,:,i)
     case(1)
        dmft_bath_%e(:,:,i)= -dmft_bath_%e(:,:,i)
        dmft_bath_%v(:,:,i)=  dmft_bath_%v(:,:,i)
     case(2)
        tmp_dmft_bath%e(:,1,i) = -dmft_bath_%e(:,2,i)
        tmp_dmft_bath%e(:,2,i) = -dmft_bath_%e(:,1,i)
        dmft_bath_%e(:,:,i)    = tmp_dmft_bath%e(:,:,i)
        tmp_dmft_bath%v(:,1,i) = dmft_bath_%v(:,2,i)
        tmp_dmft_bath%v(:,2,i) = dmft_bath_%v(:,1,i)
        dmft_bath_%v(:,:,i)    = tmp_dmft_bath%v(:,:,i)
     end select
  end do
  if(save_)call save_dmft_bath(dmft_bath_)
  call get_dmft_bath(dmft_bath_,bath_)
  call deallocate_dmft_bath(dmft_bath_)
end subroutine ph_trans_bath_site
subroutine ph_trans_bath_lattice(bath_,save)
  real(8),dimension(:,:) :: bath_
  logical,optional       :: save
  logical                :: save_
  integer                :: Nsites,ilat
  save_=.true.;if(present(save))save_=save
  Nsites=size(bath_,1)
  do ilat=1,Nsites
     ed_file_suffix=reg(ineq_site_suffix)//reg(txtfy(ilat,site_indx_padding))
     call ph_trans_bath_site(bath_(ilat,:),save_)
  enddo
  ed_file_suffix=""
end subroutine ph_trans_bath_lattice
