program ed_hm_daghofer
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  !
  USE MPI
  !
  implicit none

  integer                                     :: iloop,Lk,Nso
  logical                                     :: converged
  integer                                     :: ispin,ilat!,i,j

  !Bath:
  integer                                     :: Nb
  real(8),allocatable,dimension(:)            :: Bath,Bath_prev

  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:) :: Weiss
  complex(8),allocatable,dimension(:,:,:,:,:) :: Smats,Sreal
  complex(8),allocatable,dimension(:,:,:,:,:) :: Gmats,Greal

  !hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)     :: Hk
  complex(8),allocatable,dimension(:,:)       :: modelHloc
  complex(8),allocatable,dimension(:,:,:,:)   :: Hloc
  real(8),allocatable,dimension(:)            :: Wtk

  integer,allocatable,dimension(:)            :: ik2ix,ik2iy
  real(8),dimension(2)                        :: bk1,bk2

  !variables for the model:
  integer                                     :: Nk,Nkpath
  real(8)                                     :: ts,wmixing
  character(len=32)                           :: finput
  character(len=32)                           :: hkfile
  logical                                     :: spinsym,bathsym,iget_akw
  !
  real(8),dimension(2)                        :: Eout
  real(8),allocatable,dimension(:)            :: dens
  !
  real(8),dimension(:,:),allocatable          :: Zmats
  complex(8),dimension(:,:),allocatable       :: Zfoo
  complex(8),allocatable,dimension(:,:,:,:)   :: S0
  !
  integer                                     :: comm,rank
  logical                                     :: master
  !
  !modify daghofer hamiltonian
  real(8)                                     :: alpha,theta,etanm


  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)



  !parse input file
  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  !
  call parse_input_variable(alpha,"ALPHA",finput,default=1.d0)
  call parse_input_variable(theta,"THETA",finput,default=0.d0)
  call parse_input_variable(etanm,"ETANM",finput,default=0.d0)
  !Parse additional variables && read Input && read H(k)^2x2
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(ts,"TS","inputED.conf",default=1d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  !
  call parse_input_variable(bathsym,"BATHSYM",finput,default=.false.)
  !
  call parse_input_variable(iget_akw,"IGET_AKW",finput,default=.false.)
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")

  call ed_read_input(trim(finput),comm)

  if(Norb/=3)stop "Wrong setup from input file: Norb=3"
  Nso=Nspin*Norb



  !Allocate Weiss Field:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats));Weiss=zero
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats));Smats=zero
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats));Gmats=zero
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal));Sreal=zero
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal));Greal=zero
  allocate(Hloc(Nspin,Nspin,Norb,Norb));Hloc=zero
  allocate(S0(Nspin,Nspin,Norb,Norb));S0=zero
  allocate(Zmats(Nso,Nso));Zmats=eye(Nso)
  allocate(Zfoo(Nso,Nso));Zfoo=0d0
  allocate(dens(Norb));dens=0d0


  if(iget_akw)then
     call read_sigma_real(Sreal)
     call get_Akw(Sreal)
     stop
  endif






  !Build the Hamiltonian on a grid or on a path
  call build_hk(trim(hkfile))
  Hloc = so2nn_reshape(modelHloc,Nspin,Norb)

  !Setup solver
  Nb=get_bath_dimension()
  allocate(Bath(Nb))
  allocate(Bath_prev(Nb))
  call ed_init_solver(Bath,Hloc)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(comm,Bath,Hloc)
     call ed_get_sigma_matsubara(Smats)


     ! compute the local gf:
     call dmft_gloc_matsubara(comm,Hk,Wtk,Gmats,Smats)
     if(master)call dmft_print_gf_matsubara(Gmats,"LG",iprint=1)


     ! compute the Weiss field (only the Nineq ones)
     call dmft_self_consistency(comm,Gmats,Smats,Weiss,Hloc,SCtype=cg_scheme)

     !Fit the new bath, starting from the old bath + the supplied Weiss
     call ed_chi2_fitgf(comm,Weiss,Bath,ispin=1,iorb=1)
     call ed_chi2_fitgf(comm,Weiss,Bath,ispin=1,iorb=3)

     !copy iorb=1 component in iorb=2 to enforce 1=2 symmetry
     if(bathsym)call copy_component_bath(Bath,1,1,Bath,1,2)

     !MIXING:
     if(iloop>1)Bath=wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath

     if(master)then
        converged = check_convergence(Weiss(1,1,1,1,:),dmft_error,nsuccess,nloop)
        if(nread/=0.d0) then
           call ed_get_dens(dens)
           call ed_search_variable(xmu,sum(dens),converged)
        endif
     endif
     call Bcast_MPI(comm,converged)
     call Bcast_MPI(comm,xmu)
     !
     if(master)call end_loop
  enddo

  call ed_get_sigma_real(Sreal)
  call dmft_gloc_realaxis(comm,Hk,Wtk,Greal,Sreal)
  if(master)call dmft_print_gf_realaxis(Greal,"LG",iprint=1)
  if(master)call dmft_print_gf_matsubara(Smats,"LSigma",iprint=1)
  if(master)call dmft_print_gf_realaxis(Sreal,"LSigma",iprint=1)



contains




  !--------------------------------------------------------------------!
  !Lattice Hamitonian:
  !--------------------------------------------------------------------!
  function hk_model(kpoint,Nso) result(hk)
    real(8),dimension(:)            :: kpoint
    integer                         :: Nso
    complex(8),dimension(Nso,Nso)   :: hk
    real(8)                         :: kx,ky
    real(8)                         :: t1,t2,t3,t4,t5,t6,t7,t8,dxy,xmu_tb
    !
    kx=kpoint(1)
    ky=kpoint(2)
    !
    hk(:,:) = zero
    !
    !daghofer model 
    t1  =  0.02d0
    t2  =  0.06d0
    t3  =  0.03d0
    t4  = -0.01d0
    t5  =  0.2d0*alpha
    t6  =  0.3d0*alpha
    t7  = -0.2d0*alpha
    t8  = -t7/2.d0
    dxy =  0.4d0-theta
    !
    xmu_tb = 0.212d0
    !
    hk(1,1) = 2.d0*t2*cos(kx) + 2.d0*t1*cos(ky) + 4.d0*t3*cos(kx)*cos(ky) - xmu_tb + etanm
    hk(2,2) = 2.d0*t1*cos(kx) + 2.d0*t2*cos(ky) + 4.d0*t3*cos(kx)*cos(ky) - xmu_tb - etanm
    hk(3,3) = 2.d0*t5*(cos(kx)+cos(ky)) + 4.d0*t6*cos(kx)*cos(ky) + dxy   - xmu_tb
    hk(1,2) = 4.d0*t4*sin(kx)*sin(ky)
    hk(1,3) = 2.d0*t7*sin(kx)*xi + 4.d0*t8*sin(kx)*cos(ky)*xi
    hk(2,3) = 2.d0*t7*sin(ky)*xi + 4.d0*t8*sin(ky)*cos(kx)*xi 
    !
    hk(2,1) = hk(1,2)
    hk(3,1) = dconjg(hk(1,3))
    hk(3,2) = dconjg(hk(2,3))
    !
  end function hk_model






  !---------------------------------------------------------------------
  !PURPOSE: get model Hamiltonian
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional             :: file
    integer                               :: i,j,ik
    integer                               :: ix,iy
    real(8)                               :: kx,ky  
    integer                               :: iorb,jorb
    integer                               :: isporb,jsporb
    integer                               :: ispin,jspin
    integer                               :: unit
    complex(8),dimension(Nso,Nso,Lmats)   :: Gmats
    complex(8),dimension(Nso,Nso,Lreal)   :: Greal
    real(8),dimension(2)                  :: kvec
    real(8)                               :: blen,area_hex,area_rect,points_in,points_tot
    real(8),allocatable,dimension(:)      :: kxgrid,kygrid
    real(8),dimension(:,:),allocatable    :: kpath

    Lk= Nk*Nk

    write(LOGfile,*)"Build H(k)    :",Lk
    write(LOGfile,*)"# of SO-bands :",Nso

    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nso,Nso,Lk));Hk=zero
    allocate(wtk(Lk));Wtk=0d0

    call TB_set_bk([pi2,0d0,0d0],[0d0,pi2,0d0],[0d0,0d0,pi2])

    call TB_build_model(Hk, hk_model, Nso, [Nk,Nk])

    Wtk = 1d0/Lk

    if(present(file))&
         call TB_write_hk(Hk, trim(file), &
         No = Nso,Nd = Norb,Np = 0,Nineq = 1,&
         Nkvec=[Nk,Nk])
    !
    allocate(modelHloc(Nso,Nso))
    modelHloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(modelHloc))<1.d-4)modelHloc=0d0

    if(master)then
       !path: G X M G
       allocate(kpath(4,3))
       kpath(1,:)=[0d0,0d0,0d0]
       kpath(2,:)=[ pi,0d0,0d0]
       kpath(3,:)=[ pi, pi,0d0]
       kpath(4,:)=[0d0,0d0,0d0]
       call TB_solve_model(hk_model,Nso,kpath,Nkpath,&
            colors_name=[red1,green1,blue1],&
            points_name=[character(len=10) :: "G","X","M", "G"],&
            file="Eigenbands.nint")

       !draw bands along both X and Y high-symmetry points (nematic)
       if(allocated(kpath))deallocate(kpath)
       !path: G M Y G X M G
       allocate(kpath(7,3))
       kpath(1,:)=[0d0,0d0,0d0]
       kpath(2,:)=[ pi, pi,0d0]
       kpath(3,:)=[0d0, pi,0d0]
       kpath(4,:)=[0d0,0d0,0d0]
       kpath(5,:)=[ pi,0d0,0d0]
       kpath(6,:)=[ pi, pi,0d0]
       kpath(7,:)=[0d0,0d0,0d0]
       call TB_solve_model(hk_model,Nso,kpath,Nkpath,&
            colors_name=[red1,green1,blue1],&
            points_name=[character(len=10) :: "G", "M", "Y", "G","X","M", "G"],&
            file="Eigenbands_XY.nint")

    endif
  end subroutine build_hk






  !----------------------------------------------------------------------------------------!
  ! purpose: read the real local self-energy from disk
  !----------------------------------------------------------------------------------------!
  subroutine read_sigma_real(Sreal)
    complex(8),intent(inout)         :: Sreal(:,:,:,:,:)
    character(len=30)                :: suffix
    integer                          :: ilat,ispin,iorb
    real(8),dimension(:),allocatable :: wm,wr

    if(size(Sreal,1)/=Nspin) stop "save_sigma: error in dim 2. Nspin"
    if(size(Sreal,2)/=Nspin) stop "save_sigma: error in dim 3. Nspin"
    if(size(Sreal,3)/=Norb) stop "save_sigma: error in dim 4. Norb"
    if(size(Sreal,4)/=Norb) stop "save_sigma: error in dim 5. Norb"

    allocate(wr(Lreal))

    wr = linspace(wini,wfin,Lreal)
    write(LOGfile,*)"write spin-orbital diagonal elements:"
    do ispin=1,Nspin
       do iorb=1,Norb
          suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
          call sread("LSigma"//trim(suffix),wr,Sreal(ispin,ispin,iorb,iorb,:))
       enddo
    enddo

  end subroutine read_sigma_real






  subroutine get_Akw(Sreal)
    integer                                       :: ik,ispin,iorb,ilat
    integer                                       :: Npts,Nktot
    !
    complex(8),allocatable,dimension(:,:,:,:,:)   :: Sreal
    complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gkreal
    real(8),allocatable,dimension(:,:)            :: Akreal
    real(8),dimension(:),allocatable              :: wr
    real(8),dimension(:),allocatable              :: kgrid
    real(8),dimension(:,:),allocatable            :: kpath
    !
    Nspin = size(Sreal,1)
    Norb  = size(Sreal,3)
    Lreal = size(Sreal,5)
    Nso=Nspin*Norb

    allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    !path: G X M G
    allocate(kpath(4,3))
    kpath(1,:)=[0d0,0d0,0d0]
    kpath(2,:)=[ pi,0d0,0d0]
    kpath(3,:)=[ pi, pi,0d0]
    kpath(4,:)=[0d0,0d0,0d0]
    !
    !get kpoints
    Npts  = size(kpath,1)
    Nktot = (Npts-1)*Nkpath
    !
    allocate(kgrid(Nktot))
    !
    !allocate Hamiltonian and build model along path
    allocate(Hk(Nso,Nso,Nktot));Hk=zero
    call TB_build_model(hk,hk_model,Nso,kpath,Nkpath)
    !
    !allocate and compute Gkw
    allocate(Gkreal(Nktot,Nspin,Nspin,Norb,Norb,Lreal))
    do ik=1,Nktot
       kgrid(ik)=ik
       call dmft_gk_realaxis(Hk(:,:,ik),1d0,Gkreal(ik,:,:,:,:,:),Sreal) 
    enddo
    !
    !get Akw
    allocate(Akreal(Nktot,Lreal))
    Akreal=zero
    do ispin=1,Nspin
       do iorb=1,Norb
          Akreal = Akreal - dimag(Gkreal(:,ispin,ispin,iorb,iorb,:))/pi/Nspin/Norb
       enddo
    enddo
    call splot3d("Akw.dat",kgrid,wr,Akreal(:,:))
    ! 
  end subroutine get_Akw



end program ed_hm_daghofer



