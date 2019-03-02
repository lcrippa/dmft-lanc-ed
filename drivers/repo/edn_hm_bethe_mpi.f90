program lancED
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  integer                                     :: iloop,Nb,Le
  logical                                     :: converged
  real(8)                                     :: wband
  !Bath:
  real(8),allocatable                         :: Bath(:),Bath_(:)
  !The local hybridization function:
  complex(8),allocatable                      :: Hloc(:,:,:,:)
  complex(8),allocatable,dimension(:,:,:,:,:) :: Weiss,Smats,Sreal,Gmats,Greal
  character(len=16)                           :: finput
  real(8)                                     :: wmixing,Eout(2),de,dens
  real(8),allocatable                         :: Gtau(:)
  real(8),dimension(:,:,:),allocatable        :: He
  real(8),dimension(:),allocatable            :: Wte
  integer                                     :: comm,rank
  logical                                     :: master


  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)


  call parse_cmd_variable(finput,"FINPUT",default='inputED.in')
  call parse_input_variable(wband,"wband",finput,default=1d0)
  call parse_input_variable(Le,"LE",finput,default=500)
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


  allocate(He(1,1,Le),Wte(Le))
  He(1,1,:) = linspace(-Wband,Wband,Le,mesh=de)
  Wte       = dens_bethe(He(1,1,:),wband)*de


  !Allocate Weiss Field:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  Hloc=zero

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

     ! Retrieve the Self-energy from ED
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_real(Sreal)
     !
     ! compute the local gf & print it
     call dmft_gloc_matsubara(Comm,one*He,Wte,Gmats,Smats)
     if(master)call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)
     !
     !Get the Weiss field/Delta function to be fitted
     if(master)call dmft_self_consistency(Gmats,Smats,Weiss,Hloc,cg_scheme)
     call Bcast_MPI(comm,Weiss)
     !
     !
     !Perform the SELF-CONSISTENCY by fitting the new bath
     if(master)then
        call ed_chi2_fitgf(Weiss,bath,ispin=1)
        !
        !MIXING:
        if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
        Bath_=Bath
        !
        !Check convergence (if required change chemical potential)
        converged = check_convergence(Weiss(1,1,1,1,:),dmft_error,nsuccess,nloop,reset=.false.)
     endif
     call Bcast_MPI(comm,bath)
     call Bcast_MPI(comm,converged)
     call Bcast_MPI(comm,xmu)
     !
     if(master)call end_loop
  enddo


  call dmft_gloc_realaxis(Comm,one*He,Wte,Greal,Sreal)
  if(master)call dmft_print_gf_realaxis(Greal,"Greal",iprint=1)
  call dmft_kinetic_energy(Comm,one*He,Wte,Smats)



  call finalize_MPI()



end program lancED

