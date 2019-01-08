program hm_Nbands_bethe
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  integer                                     :: iloop,Nb,Le,Nso,iorb
  logical                                     :: converged
  real(8),dimension(5)                        :: Wbethe,Dbethe
  !Bath:
  real(8),allocatable                         :: Bath(:),Bath_(:)
  !  
  complex(8),allocatable                      :: Hloc(:,:,:,:)
  real(8),dimension(:,:),allocatable          :: Dbands
  real(8),dimension(:,:),allocatable          :: Ebands
  real(8),dimension(:),allocatable            :: H0
  real(8),dimension(:),allocatable            :: de,dens
  !
  real(8),dimension(:),allocatable            :: Wband
  !
  complex(8),allocatable,dimension(:,:,:,:,:) :: Weiss,Smats,Sreal,Gmats,Greal
  character(len=16)                           :: finput
  real(8)                                     :: wmixing
  !
  integer                                     :: comm,rank
  logical                                     :: master


  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)


  call parse_cmd_variable(finput,"FINPUT",default='inputED.in')
  call parse_input_variable(Le,"LE",finput,default=500)
  call parse_input_variable(Wbethe,"WBETHE",finput,default=[1d0,1d0,1d0,1d0,1d0])
  call parse_input_variable(Dbethe,"DBETHE",finput,default=[0d0,0d0,0d0,0d0,0d0])
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  !
  call ed_read_input(trim(finput),comm)

  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  if(Nspin/=1.OR.Norb>5)stop "Wrong setup from input file: Nspin=1 OR Norb>5"
  Nso=Nspin*Norb


  allocate(Ebands(Nso,Le))
  allocate(Dbands(Nso,Le))
  allocate(Wband(Nso))
  allocate(H0(Nso))
  allocate(de(Nso))
  !
  Wband = Wbethe(:Norb)
  H0    = Dbethe(:Norb)
  do iorb=1,Norb
     Ebands(iorb,:) = linspace(-Wband(iorb),Wband(iorb),Le,mesh=de(iorb))
     Dbands(iorb,:) = dens_bethe(Ebands(iorb,:),Wband(iorb))*de(iorb)
  enddo
  if(master)call TB_write_Hloc(one*diag(H0))


  !Allocate Weiss Field:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  allocate(dens(Norb))
  Hloc(1,1,:,:)=diag(H0)


  !setup solver
  Nb=get_bath_dimension()
  allocate(bath(Nb))
  allocate(bath_(Nb))
  call ed_init_solver(comm,bath,Hloc)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")
     !
     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(comm,bath)
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_real(Sreal)
     !
     ! compute the local gf:
     call dmft_gloc_matsubara(comm,Ebands,Dbands,H0,Gmats,Smats)
     if(master)call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)
     !
     call dmft_gloc_realaxis(Comm,Ebands,Dbands,H0,Greal,Sreal)
     if(master)call dmft_print_gf_realaxis(Greal,"Greal",iprint=1)
     !
     !Get the Weiss field/Delta function to be fitted
     call dmft_self_consistency(comm,Gmats,Smats,Weiss,Hloc,SCtype=cg_scheme)
     !
     !Perform the SELF-CONSISTENCY by fitting the new bath
     call ed_chi2_fitgf(comm,Weiss,bath,ispin=1)
     !
     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
     Bath_=Bath

     !Check convergence (if required change chemical potential)
     if(master)then
        converged = check_convergence(Weiss(1,1,1,1,:)+Weiss(1,1,2,2,:),dmft_error,nsuccess,nloop,reset=.false.)
        call ed_get_dens(dens)
        if(nread/=0d0)call search_chemical_potential(xmu,sum(dens),converged)
     endif
     call Bcast_MPI(comm,converged)
     call Bcast_MPI(comm,xmu)
     !
     if(master)call end_loop
  enddo


  call dmft_gloc_realaxis(Comm,Ebands,Dbands,H0,Greal,Sreal)
  if(master)call dmft_print_gf_realaxis(Greal,"Greal",iprint=1)
  call dmft_kinetic_energy(comm,Ebands,Dbands,H0,Smats(1,1,:,:,:))

  call finalize_MPI()


end program hm_Nbands_bethe



