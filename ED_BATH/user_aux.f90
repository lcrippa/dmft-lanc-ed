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
function get_bath_dimension_direct(Hloc_nn,ispin_) result(bath_size)
  complex(8),optional,intent(in) :: Hloc_nn(:,:,:,:)
  integer,optional               :: ispin_
  integer                        :: bath_size,ndx,ispin,iorb,jspin,jorb,io,jo,counter
  real(8),allocatable            :: Hloc(:,:,:,:)
  select case(bath_type)
  case default
     !( e [Nspin][Norb][Nbath] + v [Nspin][Norb][Nbath] )
     bath_size = Norb*Nbath + Norb*Nbath
     if(.not.present(ispin_))bath_size=Nspin*bath_size
     !
     !
  case('hybrid')
     !(e [Nspin][1][Nbath] + v [Nspin][Norb][Nbath] )
     bath_size = Nbath + Norb*Nbath
     if(.not.present(ispin_))bath_size=Nspin*bath_size
     !
     !
  case('replica')
     !
     if(present(Hloc_nn))then
        allocate(Hloc(Nspin,Nspin,Norb,Norb));Hloc=dreal(Hloc_nn)
     elseif(allocated(impHloc))then
        allocate(Hloc(Nspin,Nspin,Norb,Norb));Hloc=impHloc
     else
        stop "ERROR: get_bath_dimension: bath_type=replica neither Hloc_nn present nor impHloc allocated"
     endif
     counter=0
     !
     !Real part of nonzero elements
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 io=index_stride_so(ispin,iorb)
                 jo=index_stride_so(jspin,jorb)
                 if((Hloc(ispin,jspin,iorb,jorb).ne.0d0).and.(io.le.jo))then
                    if(Hloc(ispin,jspin,iorb,jorb).ne.0.d0)counter=counter+1
                 endif
              enddo
           enddo
        enddo
     enddo
     ndx   = counter         !all elements
     ndx   = ndx + 1         !we also print n_Dec
     !
     !number of non vanishing elements for each replica
     ndx = ndx * Nbath
     !diagonal hybridizations: Vs
     ndx = ndx + Nbath
     !
     bath_size = ndx
  end select
end function get_bath_dimension_direct

function get_bath_dimension_symmetries(Hloc_nn) result(bath_size)
  real(8),dimension(:,:,:,:,:),intent(in) :: Hloc_nn
  integer                                 :: bath_size,ndx,isym,Nsym
  !
  !number of symmetries
  Nsym=size(Hloc_nn,5)
  !
  !add identity
  ndx=Nsym
  !
  !for each replica we also print N_dec
  ndx=ndx+1
  !
  !number of replicas
  ndx = ndx * Nbath
  !diagonal hybridizations: Vs
  ndx = ndx + Nbath
  !
  bath_size = ndx
  !
end function get_bath_dimension_symmetries





!+-------------------------------------------------------------------+
!PURPOSE  : Check if the dimension of the bath array are consistent
!+-------------------------------------------------------------------+
function check_bath_dimension(bath_,Hloc_nn) result(bool)
  real(8),dimension(:)        :: bath_
  integer                     :: Ntrue,i
  logical                     :: bool
  real(8),optional,intent(in) :: Hloc_nn(:,:,:,:)
  real(8),allocatable         :: Hbasis_rebuild(:,:,:,:,:)![Nspin][:][Norb][:][Nsym]
  select case (bath_type)
  case default
     if (present(Hloc_nn))then
        Ntrue = get_bath_dimension(one*Hloc_nn)
     else
        Ntrue = get_bath_dimension()
     endif
  case ('replica')
     if(.not.allocated(H_basis))STOP "check_bath_dimension: Hbasis not allocated"
     if(.not.allocated(Hbasis_rebuild))allocate(Hbasis_rebuild(Nspin,Nspin,Norb,Norb,size(H_basis)))
     do i=1,size(H_basis)
        Hbasis_rebuild(:,:,:,:,i)=H_basis(i)%O
     enddo
     Ntrue   = get_bath_dimension(Hbasis_rebuild)
  end select
  bool  = ( size(bath_) == Ntrue )
end function check_bath_dimension









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
  if(bath_type=="replica")stop "break_symmetry_bath_site ERROR: can not be used with bath_type=replica"
  save_=.true.;if(present(save))save_=save
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  dmft_bath_%e(1,:,:)    =dmft_bath_%e(1,:,:)      + sign*field
  dmft_bath_%e(Nspin,:,:)=dmft_bath_%e(Nspin,:,:)  - sign*field
  if(save_)call save_dmft_bath(dmft_bath_)
  call get_dmft_bath(dmft_bath_,bath_)
  call deallocate_dmft_bath(dmft_bath_)
end subroutine break_symmetry_bath_site
!
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
  if(bath_type=="replica")stop "spin_symmetry_bath_site ERROR: can not be used with bath_type=replica"
  save_=.true.;if(present(save))save_=save
  if(Nspin==1)then
     write(LOGfile,"(A)")"spin_symmetrize_bath: Nspin=1 nothing to symmetrize"
     return
  endif
  !
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  !
  dmft_bath_%e(Nspin,:,:)=dmft_bath_%e(1,:,:)
  dmft_bath_%v(Nspin,:,:)=dmft_bath_%v(1,:,:)
  !
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
  if(bath_type=="replica")stop "orb_symmetry_bath_site ERROR: can not be used with bath_type=replica"
  save_=.true.;if(present(save))save_=save
  if(Norb==1)then
     write(LOGfile,"(A)")"orb_symmetrize_bath: Norb=1 nothing to symmetrize"
     return
  endif
  !
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  !
  if(allocated(lvl))deallocate(lvl);allocate(lvl(Nspin,Nbath));lvl=0d0;lvl=sum(dmft_bath_%e,dim=2)/Norb
  if(allocated(hyb))deallocate(hyb);allocate(hyb(Nspin,Nbath));hyb=0d0;hyb=sum(dmft_bath_%v,dim=2)/Norb
  do iorb=1,Norb
     dmft_bath_%e(:,iorb,:)=lvl
     dmft_bath_%v(:,iorb,:)=hyb
  enddo
  !
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


subroutine orb_equality_bath_site(bath_,indx,save)
  real(8),dimension(:)   :: bath_
  type(effective_bath)   :: dmft_bath_
  integer,optional       :: indx
  logical,optional       :: save
  integer                :: indx_
  logical                :: save_
  integer                :: iorb
  real(8),allocatable    :: lvl(:,:),hyb(:,:)
  if(bath_type=="replica")stop "orb_equality_bath_site ERROR: can not be used with bath_type=replica"
  indx_=1     ;if(present(indx))indx_=indx
  save_=.true.;if(present(save))save_=save
  if(Norb==1)then
     write(LOGfile,"(A)")"orb_symmetrize_bath: Norb=1 nothing to symmetrize"
     return
  endif
  !
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  !
  if(allocated(lvl))deallocate(lvl);allocate(lvl(Nspin,Nbath));lvl=0d0;lvl=dmft_bath_%e(:,indx_,:)
  if(allocated(hyb))deallocate(hyb);allocate(hyb(Nspin,Nbath));hyb=0d0;hyb=dmft_bath_%v(:,indx_,:)
  do iorb=1,Norb
     if(iorb==indx_)cycle
     dmft_bath_%e(:,iorb,:)=lvl
     dmft_bath_%v(:,iorb,:)=hyb
  enddo
  !
  if(save_)call save_dmft_bath(dmft_bath_)
  call get_dmft_bath(dmft_bath_,bath_)
  call deallocate_dmft_bath(dmft_bath_)
end subroutine orb_equality_bath_site
subroutine orb_equality_bath_lattice(bath_,indx,save)
  real(8),dimension(:,:) :: bath_
  integer,optional       :: indx
  logical,optional       :: save
  integer                :: indx_
  logical                :: save_
  integer                :: iorb
  integer                :: Nsites,ilat
  indx_=1     ;if(present(indx))indx_=indx
  save_=.true.;if(present(save))save_=save
  Nsites=size(bath_,1)
  do ilat=1,Nsites
     ed_file_suffix=reg(ineq_site_suffix)//reg(txtfy(ilat,site_indx_padding))
     call orb_equality_bath_site(bath_(ilat,:),indx_,save_)
  enddo
  ed_file_suffix=""
end subroutine orb_equality_bath_lattice



subroutine ph_symmetrize_bath_site(bath_,save)
  real(8),dimension(:)   :: bath_
  type(effective_bath)   :: dmft_bath_
  integer                :: i
  logical,optional       :: save
  logical                :: save_
  if(bath_type=="replica")stop "ph_symmetry_bath_site ERROR: can not be used with bath_type=replica"
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
  if(bath_type=="replica")stop "ph_trans_bath_site ERROR: can not be used with bath_type=replica"
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
