MODULE ED_SETUP
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE SF_TIMER
  USE SF_IOTOOLS, only:free_unit,reg,file_length
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none
  private


  interface build_sector
     ! module procedure :: build_sector_full
     module procedure :: build_sector_sresolved
  end interface build_sector

  interface delete_sector
     ! module procedure :: delete_sector_full
     module procedure :: delete_sector_sresolved
  end interface delete_sector

  interface get_Nup
     module procedure :: get_Nup_Main
     module procedure :: get_Nup_Orbs
  end interface get_Nup

  interface get_Ndw
     module procedure :: get_Ndw_Main
     module procedure :: get_Ndw_Orbs
  end interface get_Ndw

  interface get_DimUp
     module procedure :: get_DimUp_Main
     module procedure :: get_DimUp_Orbs
  end interface get_DimUp

  interface get_DimDw
     module procedure :: get_DimDw_Main
     module procedure :: get_DimDw_Orbs
  end interface get_DimDw



  public :: setup_ed_dimensions
  public :: init_ed_structure
  public :: setup_global
  !
  public :: build_sector
  public :: delete_sector
  !  
  public :: get_Sector
  public :: get_Indices
  public :: get_Nup
  public :: get_Ndw
  public :: get_DimUp
  public :: get_DimDw
  public :: get_Dim
  !
  public :: indices2state
  public :: state2indices
  public :: iup_index
  public :: idw_index
  !
  public :: bdecomp
  public :: bjoin
  !
  public :: c,cdg
  !    
  public :: twin_sector_order
  public :: get_twin_sector
  public :: flip_state
  !
  public :: binary_search



contains



  !+------------------------------------------------------------------+
  !PURPOSE  : Setup Dimensions of the problem
  ! Norb    = # of impurity orbitals
  ! Nbath   = # of bath levels (depending on bath_type)
  ! Ns      = # of levels (per spin)
  ! Nlevels = 2*Ns = Total # of levels (counting spin degeneracy 2) 
  !+------------------------------------------------------------------+
  subroutine setup_ed_dimensions()
    integer                                           :: maxtwoJz,inJz,dimJz
    integer                                           :: isector,in,shift
    select case(bath_type)
    case default
       Ns = (Nbath+1)*Norb
    case ('hybrid')
       Ns = Nbath+Norb
       if(.not.ed_total_ud)stop "setup_ed_dimension: bath_type==hybrid AND .NOT.ed_total_ud"
    case ('replica')
       Ns = Norb*(Nbath+1)
    end select
    Nlevels  = 2*Ns
    !
    Ns_Orb   = Ns/Norb
    !
    select case (ed_total_ud)
    case (.true.)
       Nsectors = (Ns+1)*(Ns+1) !nup=0:Ns;ndw=0:Ns
    case (.false.)
       Nsectors = ((Ns_Orb+1)*(Ns_Orb+1))**Norb !prod_{a=1,Norb}nup(a)=0:Ns;ndw(a)=0:Ns
    end select
  end subroutine setup_ed_dimensions



  !+------------------------------------------------------------------+
  !PURPOSE  : Init ED structure and calculation
  !+------------------------------------------------------------------+
  subroutine init_ed_structure(MpiComm)
    integer,optional                                  :: MpiComm
    logical                                           :: control
    real(8),dimension(Nspin,Nspin,Norb,Norb)          :: reHloc         !local hamiltonian, real part 
    real(8),dimension(Nspin,Nspin,Norb,Norb)          :: imHloc         !local hamiltonian, imag part
    integer                                           :: i,iorb,jorb,ispin,jspin
    integer,allocatable,dimension(:,:)                :: dim_sector_max
    logical                                           :: MPI_MASTER=.true.
    !
#ifdef _MPI
    if(present(MpiComm))MPI_MASTER=get_Master_MPI(MpiComm)
#endif
    !
    call setup_ed_dimensions()
    !
    select case(ed_total_ud)
    case (.true.)
       allocate(dim_sector_max(2,1))
       dim_sector_max(1,1)=get_sector_dimension(Ns,Ns/2)
       dim_sector_max(2,1)=get_sector_dimension(Ns,Ns-Ns/2)
    case (.false.)
       allocate(dim_sector_max(2,Norb))
       do iorb=1,Norb
          dim_sector_max(1,iorb) = get_sector_dimension(Ns_Orb,Ns_Orb/2)
          dim_sector_max(2,iorb) = get_sector_dimension(Ns_Orb,Ns_Orb-Ns_Orb/2)
       enddo
    end select
    !
    write(LOGfile,"(A)")"Summary:"
    write(LOGfile,"(A)")"--------------------------------------------"
    write(LOGfile,"(A,I15)")'# of levels/spin      = ',Ns
    write(LOGfile,"(A,I15)")'Total size            = ',Nlevels
    write(LOGfile,"(A,I15)")'# of impurities       = ',Norb
    write(LOGfile,"(A,I15)")'# of bath/impurity    = ',Nbath
    write(LOGfile,"(A,I15)")'# of Bath levels/spin = ',Ns-Norb
    write(LOGfile,"(A,I15)")'Number of sectors     = ',Nsectors
    write(LOGfile,"(A,10I15)")'Largest Sector(s)     = ',dim_sector_max,product(dim_sector_max)
    write(LOGfile,"(A)")"--------------------------------------------"
    !
    allocate(impHloc(Nspin,Nspin,Norb,Norb))
    impHloc=zero
    !
    !Allocate indexing arrays
    select case (ed_total_ud)
    case (.true.)
       allocate(getCsector(1,2,Nsectors));getCsector=0
       allocate(getCDGsector(1,2,Nsectors));getCDGsector=0
    case (.false.)
       allocate(getCsector(Norb,2,Nsectors));getCsector=0
       allocate(getCDGsector(Norb,2,Nsectors));getCDGsector=0
    end select
    !
    allocate(impIndex(Norb,2));impIndex=0
    !
    allocate(getDim(Nsectors));getDim=0
    !
    allocate(getBathStride(Norb,Nbath));getBathStride=0
    allocate(twin_mask(Nsectors))
    allocate(neigen_sector(Nsectors))
    !
    !
    !check finiteT
    finiteT=.true.              !assume doing finite T per default
    if(lanc_nstates_total==1)then     !is you only want to keep 1 state
       finiteT=.false.          !set to do zero temperature calculations
       write(LOGfile,"(A)")"Required Lanc_nstates_total=1 => set T=0 calculation"
    endif
    !
    !
    !check whether lanc_nstates_sector and lanc_states are even (we do want to keep doublet among states)
    if(finiteT)then
       if(mod(lanc_nstates_sector,2)/=0)then
          lanc_nstates_sector=lanc_nstates_sector+1
          write(LOGfile,"(A,I10)")"Increased Lanc_nstates_sector:",lanc_nstates_sector
       endif
       if(mod(lanc_nstates_total,2)/=0)then
          lanc_nstates_total=lanc_nstates_total+1
          write(LOGfile,"(A,I10)")"Increased Lanc_nstates_total:",lanc_nstates_total
       endif
    endif
    !
    !
    if(finiteT)then
       write(LOGfile,"(A)")"Lanczos FINITE temperature calculation:"
       write(LOGfile,"(A,I3)")"Nstates x Sector = ", lanc_nstates_sector
       write(LOGfile,"(A,I3)")"Nstates   Total  = ", lanc_nstates_total
       call sleep(1)
    else
       write(LOGfile,"(A)")"Lanczos ZERO temperature calculation:"
       call sleep(1)
    endif
    !
    !
    !CHECKS:
    if(Lfit>Lmats)Lfit=Lmats
    if(Nspin>2)stop "ED ERROR: Nspin > 2 is currently not supported"
    if(Norb>3)stop "ED ERROR: Norb > 3 is currently not supported"
    !
    Jhflag=.FALSE.
    if(Norb>1.AND.(Jx/=0d0.OR.Jp/=0d0))Jhflag=.TRUE.
    !
    if(.not.ed_total_ud)then
       if(bath_type=="hybrid")stop "ED ERROR: ed_total_ud=F can not be used with bath_type=hybrid" 
       if(Jhflag)stop "ED ERROR: ed_total_ud=F can not be used with Jx!=0 OR Jp!=0"
       if(ed_sparse_H.AND.ed_sparse_format=="ELL")stop "ED ERROR: ed_total_ud=F can not be used with ed_sparse_H=T AND ed_sparse_format=ELL"
       if(.not.ed_sparse_H)stop "ED ERROR: ed_total_ud=F can not be used with ed_sparse_H=F"
    endif
    !
    if(nread/=0.d0)then
       i=abs(floor(log10(abs(nerr)))) !modulus of the order of magnitude of nerror
       niter=nloop/3
    endif
    if(Nspin>1.AND.ed_twin.eqv..true.)then
       write(LOGfile,"(A)")"WARNING: using twin_sector with Nspin>1"
       call sleep(1)
    end if
    !
    !
    !allocate functions
    allocate(impSmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impSreal(Nspin,Nspin,Norb,Norb,Lreal))
    impSmats=zero
    impSreal=zero
    !
    allocate(impGmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impGreal(Nspin,Nspin,Norb,Norb,Lreal))
    impGmats=zero
    impGreal=zero
    !
    allocate(impG0mats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impG0real(Nspin,Nspin,Norb,Norb,Lreal))
    impG0mats=zero
    impG0real=zero
    !
    !allocate observables
    allocate(ed_dens(Norb),ed_docc(Norb),ed_dens_up(Norb),ed_dens_dw(Norb))
    ed_dens=0d0
    ed_docc=0d0
    ed_dens_up=0d0
    ed_dens_dw=0d0
    !
    if(chiflag)then
       allocate(spinChi_tau(Norb+1,0:Ltau))
       allocate(spinChi_w(Norb+1,Lreal))
       allocate(spinChi_iv(Norb+1,0:Lmats))
       !
       ! allocate(densChi_tau(Norb,Norb,0:Ltau))
       ! allocate(densChi_w(Norb,Norb,Lreal))
       ! allocate(densChi_iv(Norb,Norb,0:Lmats))
       ! allocate(densChi_mix_tau(Norb,Norb,0:Ltau))
       ! allocate(densChi_mix_w(Norb,Norb,Lreal))
       ! allocate(densChi_mix_iv(Norb,Norb,0:Lmats))
       ! allocate(densChi_tot_tau(0:Ltau))
       ! allocate(densChi_tot_w(Lreal))
       ! allocate(densChi_tot_iv(0:Lmats))
       ! !
       ! allocate(pairChi_tau(Norb,0:Ltau))
       ! allocate(pairChi_w(Norb,Lreal))
       ! allocate(pairChi_iv(Norb,0:Lmats))
    endif
    !
  end subroutine init_ed_structure





  !+------------------------------------------------------------------+
  !PURPOSE: SETUP THE GLOBAL POINTERS FOR THE ED CALCULAIONS.
  !+------------------------------------------------------------------+
  subroutine setup_global
    integer                          :: DimUp,DimDw
    integer                          :: DimUps(Norb),DimDws(Norb)
    integer                          :: Nup,Ndw
    integer                          :: Nups(Norb),Ndws(Norb)
    integer                          :: Jups(Norb),Jdws(Norb)
    integer                          :: Indices(2*Norb)
    integer                          :: Iup,Idw
    integer                          :: Jup,Jdw
    integer                          :: i,iorb
    integer                          :: isector,jsector
    integer                          :: unit,status,istate
    logical                          :: IOfile
    integer                          :: list_len
    integer,dimension(:),allocatable :: list_sector
    !
    !Store full dimension of the sectors:
    select case(ed_total_ud)
    case (.true.)
       do isector=1,Nsectors
          call get_DimUp(isector,DimUp)
          call get_DimDw(isector,DimDw)
          getDim(isector)  = DimUp*DimDw
       enddo
       !
    case (.false.)
       do isector=1,Nsectors
          call get_DimUp(isector,DimUps)
          call get_DimDw(isector,DimDws)
          DimUp = product(DimUps)
          DimDw = product(DimDws)       
          getDim(isector)  = DimUp*DimDw
       enddo
       !
    end select
    !
    inquire(file="state_list"//reg(ed_file_suffix)//".restart",exist=IOfile)
    if(IOfile)then
       write(LOGfile,"(A)")"Restarting from a state_list file:"
       list_len=file_length("state_list"//reg(ed_file_suffix)//".restart")
       allocate(list_sector(list_len))
       !
       open(free_unit(unit),file="state_list"//reg(ed_file_suffix)//".restart",status="old")
       status=0
       select case(ed_total_ud)
       case (.true.)
          do while(status>=0)
             read(unit,*,iostat=status)istate,isector,nup,ndw
             list_sector(istate)=isector
             call get_Nup(isector,Iup)
             call get_Ndw(isector,Idw)
             if(nup/=iup.OR.ndw/=idw)&
                  stop "setup_global error: nup!=nup(isector).OR.ndw!=ndw(isector)"
          enddo
          !
       case (.false.)
          do while(status>=0)
             read(unit,*,iostat=status)istate,isector,indices
             list_sector(istate)=isector
             call get_Nup(isector,Nups)
             call get_Ndw(isector,Ndws)
             if(any(Indices /= [Nups,Ndws]))&
                  stop "setup_global error: nups!=nups(isector).OR.ndws!=ndws(isector)"
          enddo
          close(unit)
          !
       end select
       close(unit)
       !
       lanc_nstates_total = list_len
       do isector=1,Nsectors
          neigen_sector(isector) = max(1,count(list_sector==isector))
       enddo
    else
       do isector=1,Nsectors
          !init every sector to required eigenstates
          neigen_sector(isector) = min(getDim(isector),lanc_nstates_sector)
       enddo
    endif
    !
    twin_mask=.true.
    if(ed_twin)then
       select case(ed_total_ud)
       case (.true.)
          do isector=1,Nsectors
             call get_Nup(isector,Nup)
             call get_Ndw(isector,Ndw)
             if(nup<ndw)twin_mask(isector)=.false.
          enddo
          write(LOGfile,"(A,I4,A,I4)")"Looking into ",count(twin_mask)," sectors out of ",Nsectors
          !
       case (.false.)
          write(LOGfile,"(A)")"Warning: ed_twin=T AND ed_total_ud=F are incompatible at the moment. set to F"
          ed_twin=.false.
          call sleep(2)
          !
       end select
    endif
    !
    do iorb=1,Norb
       impIndex(iorb,1)=iorb
       impIndex(iorb,2)=iorb+Ns
    enddo
    !
    select case(bath_type)
    case default
       do i=1,Nbath
          do iorb=1,Norb
             getBathStride(iorb,i) = Norb + (iorb-1)*Nbath + i
          enddo
       enddo
    case ('hybrid')
       do i=1,Nbath
          getBathStride(:,i)       = Norb + i
       enddo
    case ('replica')
       do i=1,Nbath
          do iorb=1,Norb
             getBathStride(iorb,i) = iorb + i*Norb !Norb + (i-1)*Norb + iorb
          enddo
       enddo
    end select
    !
    getCsector  = 0
    getCDGsector= 0
    select case (ed_total_ud)
    case (.true.)
       do isector=1,Nsectors
          call get_Nup(isector,Nup)
          call get_Ndw(isector,Ndw)
          !
          jup=nup-1; jdw=ndw; if(jup < 0)cycle
          call get_Sector([jup,jdw],Ns,jsector)
          getCsector(1,1,isector)=jsector
          !
          jup=nup; jdw=ndw-1; if(jdw < 0)cycle
          call get_Sector([jup,jdw],Ns,jsector)
          getCsector(1,2,isector)=jsector
          !
          jup=nup+1; jdw=ndw; if(jup > Ns)cycle
          call get_Sector([jup,jdw],Ns,jsector)
          getCDGsector(1,1,isector)=jsector
          !
          jup=nup; jdw=ndw+1; if(jdw > Ns)cycle
          call get_Sector([jup,jdw],Ns,jsector)
          getCDGsector(1,2,isector)=jsector
       enddo
       !
    case (.false.)
       do isector=1,Nsectors
          call get_Nup(isector,Nups)
          call get_Ndw(isector,Ndws)
          !
          do iorb=1,Norb
             !UPs:
             Jups=Nups
             Jdws=Ndws 
             Jups(iorb)=Jups(iorb)-1; if(Jups(iorb) < 0)cycle
             call get_Sector([Jups,Jdws],Ns_Orb,jsector)
             getCsector(iorb,1,isector)=jsector
             !
             Jups=Nups
             Jdws=Ndws 
             Jups(iorb)=Jups(iorb)+1; if(Jups(iorb) > Ns)cycle
             call get_Sector([Jups,Jdws],Ns_Orb,jsector)
             getCDGsector(iorb,1,isector)=jsector
             !
             !
             !DWs:
             Jups=Nups
             Jdws=Ndws 
             Jdws(iorb)=Jdws(iorb)-1; if(Jdws(iorb) < 0)cycle
             call get_Sector([Jups,Jdws],Ns_Orb,jsector)
             getCsector(iorb,2,isector)=jsector
             !
             Jups=Nups
             Jdws=Ndws 
             Jdws(iorb)=Jdws(iorb)+1; if(Jdws(iorb) > Ns)cycle
             call get_Sector([Jups,Jdws],Ns_Orb,jsector)
             getCDGsector(iorb,2,isector)=jsector
          enddo
       enddo
       !
    end select
    return
  end subroutine setup_global










  !##################################################################
  !##################################################################
  !AUXILIARY PROCEDURES - Sectors,Nup,Ndw,DimUp,DimDw,...
  !##################################################################
  !##################################################################
  elemental function get_sector_dimension(n,np) result(dim)
    integer,intent(in) :: n,np
    integer            :: dim
    dim = binomial(n,np)
  end function get_sector_dimension


  subroutine get_Sector(indices,Ns,isector)
    integer,dimension(:) :: indices
    integer              :: isector
    integer              :: i,Ns,N,factor
    N = size(indices)
    Factor = Ns+1
    isector = 1
    do i=N,1,-1
       isector = isector + indices(i)*(Factor)**(N-i)
    enddo
  end subroutine get_Sector


  subroutine get_Indices(isector,Ns,indices)
    integer                          :: isector,Ns
    integer,dimension(:)             :: indices
    integer                          :: i,count,Dim
    integer,dimension(size(indices)) :: indices_
    !
    Dim = size(indices)
    if(mod(Dim,2)/=0)stop "get_Indices_main error: Dim%2 != 0"
    count=isector-1
    do i=1,Dim
       indices_(i) = mod(count,Ns+1)
       count      = count/(Ns+1)
    enddo
    indices = indices_(Dim:1:-1)
  end subroutine get_Indices


  subroutine get_Nup_Main(isector,Nup)
    integer              :: isector,Nup
    integer              :: i,count
    integer,dimension(2) :: indices_
    count=isector-1
    do i=1,2
       indices_(i) = mod(count,Ns+1)
       count      = count/(Ns+1)
    enddo
    Nup     = indices_(2)
  end subroutine get_Nup_Main
  subroutine get_Nup_Orbs(isector,Nup)
    integer                   :: isector,Nup(Norb)
    integer                   :: i,count
    integer,dimension(2*Norb) :: indices_
    count=isector-1
    do i=1,2*Norb
       indices_(i) = mod(count,Ns_Orb+1)
       count      = count/(Ns_Orb+1)
    enddo
    Nup = indices_(2*Norb:Norb+1:-1)
  end subroutine get_Nup_Orbs


  subroutine get_Ndw_Main(isector,Ndw)
    integer              :: isector,Ndw
    integer              :: i,count
    integer,dimension(2) :: indices_
    count=isector-1
    do i=1,2
       indices_(i) = mod(count,Ns+1)
       count      = count/(Ns+1)
    enddo
    Ndw     = indices_(1)
  end subroutine get_Ndw_Main
  subroutine get_Ndw_Orbs(isector,Ndw)
    integer                   :: isector,Ndw(Norb)
    integer                   :: i,count
    integer,dimension(2*Norb) :: indices_
    count=isector-1
    do i=1,2*Norb
       indices_(i) = mod(count,Ns_Orb+1)
       count      = count/(Ns_Orb+1)
    enddo
    Ndw = indices_(Norb:1:-1)
  end subroutine get_Ndw_Orbs



  subroutine get_DimUp_main(isector,DimUp)
    integer                :: isector,DimUp
    integer                :: Nup
    call get_Nup_Main(isector,Nup)
    DimUp = binomial(Ns,Nup)
  end subroutine get_DimUp_main
  subroutine  get_DimUp_Orbs(isector,DimUps)
    integer                :: isector,DimUps(Norb)
    integer                :: Nups(Norb),iorb
    call get_Nup_Orbs(isector,Nups)
    do iorb=1,Norb
       DimUps(iorb) = binomial(Ns_Orb,Nups(iorb))
    enddo
  end subroutine get_DimUp_Orbs


  subroutine get_DimDw_main(isector,DimDw)
    integer                :: isector,DimDw
    integer                :: Ndw
    call get_Ndw_Main(isector,Ndw)
    DimDw = binomial(Ns,Ndw)
  end subroutine get_DimDw_main
  subroutine get_DimDw_Orbs(isector,DimDws)
    integer                :: isector,DimDws(Norb)
    integer                :: Ndws(Norb),iorb
    call get_Ndw_Orbs(isector,Ndws)
    do iorb=1,Norb
       DimDws(iorb) = binomial(Ns_Orb,Ndws(iorb))
    enddo
  end subroutine get_DimDw_Orbs


  subroutine get_Dim(isector,Dim)
    integer                  :: isector,Dim
    integer,dimension(Norb)  :: Nups,Ndws,DimUps,DimDws
    integer                  :: Nup,Ndw,DimUp,DimDw
    integer                  :: iorb
    select case(ed_total_ud)
    case (.true.)
       call get_Indices(isector,Ns,[Nup,Ndw])
       Dim = binomial(Ns,Nup)*binomial(Ns,Ndw)
    case (.false.)
       call get_Indices(isector,Ns_Orb,[Nups,Ndws])
       Dim=1
       do iorb=1,Norb
          Dim = Dim*binomial(Ns_Orb,Nups(iorb))*binomial(Ns_Orb,Ndws(iorb))
       end do
    end select
  end subroutine get_Dim




  subroutine indices2state(ivec,Nvec,istate)
    integer,dimension(:)          :: ivec
    integer,dimension(size(ivec)) :: Nvec
    integer                       :: istate,i
    istate=ivec(1)
    do i=2,size(ivec)
       istate = istate + (ivec(i)-1)*product(Nvec(1:i-1))
    enddo
  end subroutine indices2state

  subroutine state2indices(istate,Nvec,ivec)
    integer                       :: istate
    integer,dimension(:)          :: Nvec
    integer,dimension(size(Nvec)) :: Ivec
    integer                       :: i,count,N
    count = istate-1
    N     = size(Nvec)
    do i=1,N
       Ivec(i) = mod(count,Nvec(i))+1
       count   = count/Nvec(i)
    enddo
  end subroutine state2indices


  function iup_index(i,DimUp) result(iup)
    integer :: i
    integer :: DimUp
    integer :: iup
    iup = mod(i,DimUp);if(iup==0)iup=DimUp
  end function iup_index


  function idw_index(i,DimUp) result(idw)
    integer :: i
    integer :: DimUp
    integer :: idw
    idw = (i-1)/DimUp+1
  end function idw_index



  !##################################################################
  !##################################################################
  !BUILD SECTORS
  !##################################################################
  !##################################################################
  subroutine build_sector_sresolved(isector,H)
    integer                       :: isector
    type(sector_map),dimension(:) :: H
    integer                       :: Nup,Ndw
    integer                       :: Nup_,Ndw_
    integer                       :: DimUp,DimDw
    integer                       :: Nups(Norb),Ndws(Norb)
    integer                       :: Nups_(Norb),Ndws_(Norb)
    integer                       :: DimUps(Norb),DimDws(Norb)
    integer                       :: iup,idw
    integer                       :: dim,N,Nh,iorb
    !
    N = size(H)
    !
    select case(ed_total_ud)
    case (.true.)
       if(N/=2)stop "build_sector_sresolved error: size(H) != 2"
       call get_Nup(isector,Nup)
       call get_Ndw(isector,Ndw)
       call get_DimUp(isector,DimUp)
       call get_DimDw(isector,DimDw)
       !
       call map_allocate(H,[DimUp,DimDw])
       !
       !UP    
       dim=0
       do iup=0,2**Ns-1
          nup_ = popcnt(iup)
          if(nup_ /= nup)cycle
          dim  = dim+1
          H(1)%map(dim) = iup
       enddo
       !DW
       dim=0
       do idw=0,2**Ns-1
          ndw_ = popcnt(idw)
          if(ndw_ /= ndw)cycle
          dim  = dim+1
          H(2)%map(dim) = idw
       enddo
       !
    case (.false.)
       if( N/=2*Norb )stop "build_sector_sresolved error: size(H) != 2*Norb"
       !
       call get_Nup(isector,Nups)
       call get_Ndw(isector,Ndws)
       call get_DimUp(isector,DimUps)
       call get_DimDw(isector,DimDws)
       !
       call map_allocate(H,[DimUps,DimDws])
       do iorb=1,Norb
          !UP    
          dim=0
          do iup=0,2**Ns-1
             nup_ = popcnt(iup)
             if(nup_ /= Nups(iorb))cycle
             dim  = dim+1
             H(iorb)%map(dim) = iup
          enddo
          !DW
          dim=0
          do idw=0,2**Ns-1
             ndw_= popcnt(idw)
             if(ndw_ /= Ndws(iorb))cycle
             dim = dim+1
             H(iorb+Norb)%map(dim) = idw
          enddo
       enddo
    end select
    !
  end subroutine build_sector_sresolved

  subroutine delete_sector_sresolved(isector,H)
    integer                   :: isector
    type(sector_map)          :: H(:)
    call map_deallocate(H)
  end subroutine delete_sector_sresolved








  !##################################################################
  !##################################################################
  !CREATION / DESTRUCTION OPERATORS
  !##################################################################
  !##################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE: input state |in> of the basis and calculates 
  !   |out>=C_pos|in>  OR  |out>=C^+_pos|in> ; 
  !   the sign of |out> has the phase convention, pos labels the sites
  !+-------------------------------------------------------------------+
  subroutine c(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn    
    integer               :: l
    if(.not.btest(in,pos-1))stop "C error: C_i|...0_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibclr(in,pos-1)
  end subroutine c

  subroutine cdg(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn    
    integer               :: l
    if(btest(in,pos-1))stop "C^+ error: C^+_i|...1_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibset(in,pos-1)
  end subroutine cdg






  !##################################################################
  !##################################################################
  !TWIN SECTORS ROUTINES:
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE  : Build the re-ordering map to go from sector A(nup,ndw)
  ! to its twin sector B(ndw,nup), with nup!=ndw.
  !
  !- build the map from the A-sector to \HHH
  !- get the list of states in \HHH corresponding to sector B twin of A
  !- return the ordering of B-states in \HHH with respect to those of A
  !+------------------------------------------------------------------+
  subroutine twin_sector_order(isector,order)
    integer                   :: isector
    integer,dimension(:)      :: order
    integer                   :: Dim
    integer                   :: DimUp,DimDw
    integer,dimension(Norb)   :: DimUps,DimDws
    type(sector_map)          :: H(2),HI(2*Norb)
    integer                   :: i,iorb
    integer,dimension(2*Norb) :: Indices,Istates
    integer                   :: iup,idw
    !
    dim = getdim(isector)
    if(size(Order)/=dim)stop "twin_sector_order error: wrong dimensions of *order* array"
    select case(ed_total_ud)
    case (.true.)
       call get_DimUp(isector,DimUp)
       call get_DimDw(isector,DimDw)
       call build_sector(isector,H)
       do idw=1,DimDw
          do iup=1,DimUp
             i = iup + (idw-1)*DimUp
             Order(i) = flip_state( [H(1)%map(iup), H(2)%map(idw)] )
          enddo
       enddo
       call delete_sector(isector,H)
       !
    case (.false.)
       call get_DimUp(isector,DimUps)
       call get_DimDw(isector,DimDws)
       call build_sector(isector,HI)
       do i=1,Dim
          call state2indices(i,[DimUps,DimDws],Indices)
          forall(iorb=1:2*Norb)Istates(iorb) = HI(iorb)%map(Indices(iorb))
          Order(i) = flip_state( Istates )
       enddo
       call delete_sector(isector,HI)
    end select
    !
    call sort_array(Order)
    !
  end subroutine twin_sector_order



  !+------------------------------------------------------------------+
  !PURPOSE  : Flip an Hilbert space state m=|{up}>|{dw}> into:
  !
  ! normal: j=|{dw}>|{up}>  , nup --> ndw
  !+------------------------------------------------------------------+
  function flip_state(istate) result(j)
    integer,dimension(:)      :: istate
    integer                   :: j
    integer,dimension(Norb)   :: jups,jdws
    integer,dimension(2*Norb) :: dims
    integer                   :: iup,idw
    integer                   :: jup,jdw
    !
    select case (ed_total_ud)
    case (.true.)
       if(size(istate)/=2)stop "flip_state error: size(istate)!=2"
       jup = istate(2)!idw
       jdw = istate(1)!iup
       j   = jup + (jdw-1)*2**Ns
    case (.false.)
       if(size(istate)/=2*Norb)stop "flip_state error: size(istate)!= 2*Norb"
       jups = istate(Norb+1:2*Norb)
       jdws = istate(1:Norb)
       dims = 2**Ns_Orb
       call indices2state([jups,jdws],Dims,j)
    end select
    !
  end function flip_state


  !+------------------------------------------------------------------+
  !PURPOSE  : get the twin of a given sector (the one with opposite 
  ! quantum numbers): 
  ! nup,ndw ==> ndw,nup (spin-exchange)
  !+------------------------------------------------------------------+
  function get_twin_sector(isector) result(jsector)
    integer,intent(in) :: isector
    integer            :: jsector
    integer            :: iup,idw
    integer,dimension(Norb) :: Iups,Idws
    select case(ed_total_ud)
    case (.true.)
       call get_Nup(isector,iup)
       call get_Ndw(isector,idw)
       call get_Sector([idw,iup],Ns,jsector)
    case (.false.)
       call get_Nup(isector,iups)
       call get_Ndw(isector,idws)
       call get_Sector([idws,iups],Ns_Orb,jsector)
    end select
  end function get_twin_sector










  !##################################################################
  !##################################################################
  !AUXILIARY COMPUTATIONAL ROUTINES ARE HERE BELOW:
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE  : input a state |i> and output a vector ivec(Nlevels)
  !with its binary decomposition
  !(corresponds to the decomposition of the number i-1)
  !+------------------------------------------------------------------+
  function bdecomp(i,Ntot) result(ivec)
    integer :: Ntot,ivec(Ntot),l,i
    logical :: busy
    !this is the configuration vector |1,..,Ns,Ns+1,...,Ntot>
    !obtained from binary decomposition of the state/number i\in 2^Ntot
    do l=0,Ntot-1
       busy=btest(i,l)
       ivec(l+1)=0
       if(busy)ivec(l+1)=1
    enddo
  end function bdecomp



  !+------------------------------------------------------------------+
  !PURPOSE  : input a vector ib(Nlevels) with the binary sequence 
  ! and output the corresponding state |i>
  !(corresponds to the recomposition of the number i-1)
  !+------------------------------------------------------------------+
  function bjoin(ib,Ntot) result(i)
    integer                 :: Ntot
    integer,dimension(Ntot) :: ib
    integer                 :: i,j
    i=0
    do j=0,Ntot-1
       i=i+ib(j+1)*2**j
    enddo
  end function bjoin



  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the factorial of an integer N!=1.2.3...(N-1).N
  !+------------------------------------------------------------------+
  recursive function factorial(n) result(f)
    integer            :: f
    integer,intent(in) :: n
    if(n<=0)then
       f=1
    else
       f=n*factorial(n-1)
    end if
  end function factorial



  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the binomial factor n1 over n2
  !+------------------------------------------------------------------+
  elemental function binomial(n1,n2) result(nchoos)
    integer,intent(in) :: n1,n2
    real(8)            :: xh
    integer            :: i
    integer nchoos
    xh = 1.d0
    if(n2<0) then
       nchoos = 0
       return
    endif
    if(n2==0) then
       nchoos = 1
       return
    endif
    do i = 1,n2
       xh = xh*dble(n1+1-i)/dble(i)
    enddo
    nchoos = int(xh + 0.5d0)
  end function binomial



  !+------------------------------------------------------------------+
  !PURPOSE : binary search of a value in an array
  !+------------------------------------------------------------------+
  recursive function binary_search(a,value) result(bsresult)
    integer,intent(in) :: a(:), value
    integer            :: bsresult, mid
    mid = size(a)/2 + 1
    if (size(a) == 0) then
       bsresult = 0        ! not found
       !stop "binary_search error: value not found"
    else if (a(mid) > value) then
       bsresult= binary_search(a(:mid-1), value)
    else if (a(mid) < value) then
       bsresult = binary_search(a(mid+1:), value)
       if (bsresult /= 0) then
          bsresult = mid + bsresult
       end if
    else
       bsresult = mid      ! SUCCESS!!
    end if
  end function binary_search







  !+------------------------------------------------------------------+
  !PURPOSE : sort array of integer using random algorithm
  !+------------------------------------------------------------------+
  subroutine sort_array(array)
    integer,dimension(:),intent(inout)      :: array
    integer,dimension(size(array))          :: order
    integer                                 :: i
    forall(i=1:size(array))order(i)=i
    call qsort_sort( array, order, 1, size(array) )
    array=order
  contains
    recursive subroutine qsort_sort( array, order, left, right )
      integer, dimension(:)                 :: array
      integer, dimension(:)                 :: order
      integer                               :: left
      integer                               :: right
      integer                               :: i
      integer                               :: last
      if ( left .ge. right ) return
      call qsort_swap( order, left, qsort_rand(left,right) )
      last = left
      do i = left+1, right
         if ( compare(array(order(i)), array(order(left)) ) .lt. 0 ) then
            last = last + 1
            call qsort_swap( order, last, i )
         endif
      enddo
      call qsort_swap( order, left, last )
      call qsort_sort( array, order, left, last-1 )
      call qsort_sort( array, order, last+1, right )
    end subroutine qsort_sort
    !---------------------------------------------!
    subroutine qsort_swap( order, first, second )
      integer, dimension(:)                 :: order
      integer                               :: first, second
      integer                               :: tmp
      tmp           = order(first)
      order(first)  = order(second)
      order(second) = tmp
    end subroutine qsort_swap
    !---------------------------------------------!
    function qsort_rand( lower, upper )
      implicit none
      integer                               :: lower, upper
      real(8)                               :: r
      integer                               :: qsort_rand
      call random_number(r)
      qsort_rand =  lower + nint(r * (upper-lower))
    end function qsort_rand
    function compare(f,g)
      integer                               :: f,g
      integer                               :: compare
      compare=1
      if(f<g)compare=-1
    end function compare
  end subroutine sort_array




end MODULE ED_SETUP







! interface print_state_vector
!    module procedure print_state_vector_ivec
!    module procedure print_state_vector_ivec_ud
!    module procedure print_state_vector_int
! end interface print_state_vector


! !+------------------------------------------------------------------+
! !PURPOSE  : print a state vector |{up}>|{dw}>
! !+------------------------------------------------------------------+
! subroutine print_state_vector_ivec(ivec,unit)
!   integer,intent(in) :: ivec(:)
!   integer,optional   :: unit
!   integer            :: unit_
!   integer            :: i,j,Ntot
!   character(len=2)   :: fbt
!   character(len=16)  :: fmt
!   unit_=6;if(present(unit))unit_=unit
!   Ntot = size(ivec)
!   write(fbt,'(I2.2)')Ntot
!   fmt="(B"//adjustl(trim(fbt))//"."//adjustl(trim(fbt))//")"
!   i= bjoin(ivec,Ntot)
!   write(unit_,"(I9,1x,A1)",advance="no")i,"|"
!   write(unit_,"(10I1)",advance="no")(ivec(j),j=1,Ntot)
!   write(unit_,"(A4)",advance="no")"> - "
!   write(unit_,fmt,advance="yes")i
! end subroutine print_state_vector_ivec
! !
! subroutine  print_state_vector_ivec_ud(ivec,jvec,unit)
!   integer,intent(in) :: ivec(:),jvec(size(ivec))
!   integer,optional   :: unit
!   integer            :: unit_
!   integer            :: i,j,iup,idw,Ntot
!   character(len=2)   :: fbt
!   character(len=20)  :: fmt
!   unit_=6;if(present(unit))unit_=unit
!   Ntot = size(ivec)
!   write(fbt,'(I2.2)')Ntot
!   fmt="(B"//adjustl(trim(fbt))//"."//adjustl(trim(fbt))//",1x,B"//adjustl(trim(fbt))//"."//adjustl(trim(fbt))//")"
!   iup = bjoin(ivec,Ntot)
!   idw = bjoin(jvec,Ntot)
!   i = bjoin([ivec,jvec],2*Ntot)
!   write(unit_,"(I9,1x,I4,1x,A1)",advance="no")i,iup,"|"
!   write(unit_,"(10I1)",advance="no")(ivec(j),j=1,Ntot)
!   write(unit_,"(A1,I4,A2)",advance="no")">",idw," |"
!   write(unit_,"(10I1)",advance="no")(jvec(j),j=1,Ntot)
!   write(unit_,"(A4)",advance="no")"> - "
!   write(unit_,fmt,advance="yes")ibits(i,0,Ntot),ibits(i,Ntot,2*Ntot)
! end subroutine print_state_vector_ivec_ud
! !
! subroutine print_state_vector_int(i,Ntot,unit)
!   integer,intent(in) :: i
!   integer,intent(in) :: Ntot
!   integer,optional   :: unit
!   integer            :: unit_
!   integer            :: j
!   integer            :: ivec(Ntot)
!   character(len=2)   :: fbt
!   character(len=16)  :: fmt
!   unit_=6;if(present(unit))unit_=unit
!   write(fbt,'(I2.2)')Ntot
!   fmt="(B"//adjustl(trim(fbt))//"."//adjustl(trim(fbt))//")"
!   ivec = bdecomp(i,Ntot)
!   write(unit_,"(I9,1x,A1)",advance="no")i,"|"
!   write(unit_,"(10I1)",advance="no")(ivec(j),j=1,Ntot)
!   write(unit_,"(A4)",advance="no")"> - "
!   write(unit_,fmt,advance="yes")i
! end subroutine print_state_vector_int







! getCsector=0
! do isector=1,Nsectors
!    call get_Nup(isector,Nup)
!    call get_Ndw(isector,Ndw)
!    !
!    jup=nup-1; jdw=ndw; if(jup < 0)cycle
!    !
!    call get_Sector([jup,jdw],Ns,jsector)
!    getCsector(1,1,isector)=jsector
! enddo
! !
! !
! !
! do isector=1,Nsectors
!    call get_Nup(isector,Nup)
!    call get_Ndw(isector,Ndw)
!    !
!    jup=nup;jdw=ndw-1;if(jdw < 0)cycle
!    !
!    call get_Sector([jup,jdw],Ns,jsector)
!    getCsector(1,2,isector)=jsector
! enddo
! !
! !
! !
! getCDGsector=0
! do isector=1,Nsectors
!    call get_Nup(isector,Nup)
!    call get_Ndw(isector,Ndw)
!    !
!    jup=nup+1;jdw=ndw;if(jup > Ns)cycle
!    !
!    call get_Sector([jup,jdw],Ns,jsector)
!    getCDGsector(1,1,isector)=jsector
! enddo
! !
! !
! !
! do isector=1,Nsectors
!    call get_Nup(isector,Nup)
!    call get_Ndw(isector,Ndw)
!    !
!    jup=nup;jdw=ndw+1;if(jdw > Ns)cycle
!    !
!    call get_Sector([jup,jdw],Ns,jsector)
!    getCDGsector(1,2,isector)=jsector
! enddo
