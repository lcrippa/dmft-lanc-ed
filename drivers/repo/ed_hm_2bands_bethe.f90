program hm_2bands_dos
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  integer                                     :: iloop,Nb,Le,Nso
  logical                                     :: converged,converged1,converged2
  !Bath:
  real(8),allocatable                         :: Bath(:),Bath_prev(:)
  !
  real(8)                                     :: Wratio,Wband(2),dens(2)
  complex(8),allocatable                      :: Hloc(:,:,:,:),zeta(:,:)
  real(8),dimension(:,:),allocatable          :: Dbands
  real(8),dimension(:,:),allocatable          :: Ebands
  real(8),dimension(:),allocatable            :: H0,de,wmats
  real(8)                                     :: Delta
  complex(8),allocatable,dimension(:,:,:,:,:) :: Weiss,Smats,Sreal,Gmats,Greal,Gimp,Weiss_prev
  character(len=16)                           :: finput
  real(8)                                     :: wmixing,Eout(2)
  integer                                     :: sc_method,mix_method
  !
  !MPI VARIABLES (local use -> ED code has its own set of MPI variables)
  integer :: MpiComm
  integer :: MpiRank
  integer :: MpiSize
  logical :: MpiMaster
  !
  !Init MPI: use of MPI overloaded functions in SciFor
  call init_MPI(MpiComm,.true.) !init MPI and print check message on screen
  MpiRank   = get_Rank_MPI(MpiComm)
  MpiMaster = get_Master_MPI(MpiComm)
  !
  call parse_cmd_variable(finput,"FINPUT",default='inputED.in')
  call parse_input_variable(Le,"LE",finput,default=500)
  call parse_input_variable(Wratio,"Wratio",finput,default=1d0)
  call parse_input_variable(delta,"DELTA",finput,default=0d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  call parse_input_variable(mix_method,"MIX_METHOD",finput,default=1,comment="Pick the way to mix:1=mix bath;2=mix Weiss/Delta")  
  call parse_input_variable(sc_method,"SC_METHOD",finput,default=1,comment="Pick the way to update the Weiss/Delta function:1=normal self-consistency;2=bethe")
  !
  call ed_read_input(trim(finput),MpiComm)

  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  if(Nspin/=1.OR.Norb/=2)stop "Wrong setup from input file: Nspin=1; Norb=2"
  Nso=Nspin*Norb

  allocate(Ebands(Nso,Le))
  allocate(Dbands(Nso,Le))
  allocate(de(Nso))
  Wband=[1d0,Wratio]
  Ebands(1,:) = linspace(-Wband(1),Wband(1),Le,mesh=de(1))
  Ebands(2,:) = linspace(-Wband(2),Wband(2),Le,mesh=de(2))
  !
  Dbands(1,:) = dens_bethe(Ebands(1,:),Wband(1))*de(1)
  Dbands(2,:) = dens_bethe(Ebands(2,:),Wband(2))*de(2)

  !Allocate Local Functions
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Weiss_prev(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gimp(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  allocate(H0(Nso))
  H0=[-Delta/2,Delta/2]
  Hloc(1,1,:,:)=diag(H0)
  call TB_write_Hloc(Hloc(1,1,:,:))

  !Build Matsubara frequencies
  allocate(wmats(Lmats))
  allocate(zeta(2,Lmats))
  wmats     = pi/beta*(2*arange(1,Lmats)-1)
  zeta(1,:) = xi*wmats + xmu - Hloc(1,1,1,1)
  zeta(2,:) = xi*wmats + xmu - Hloc(1,1,2,2)


  !setup solver
  Nb=get_bath_dimension()
  allocate(bath(Nb))
  allocate(bath_prev(Nb))
  call ed_init_solver(MpiComm,bath,Hloc) !pass MPI communicator to the ED code:


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(MpiMaster)call start_loop(iloop,nloop,"DMFT-loop")
     !
     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(MpiComm,bath) !pass MPI communicator to the ED code
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_real(Sreal)
     call ed_get_gimp_matsubara(Gimp)
     call ed_get_dens(dens)
     !
     !
     !Mpi computation of the Gloc on matsubara
     call dmft_gloc_matsubara(MpiComm,Ebands,Dbands,H0,Gmats,Smats)
     if(MpiMaster)call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)
     !
     !
     !Get the Weiss field/Delta function to be fitted
     if(cg_scheme=='weiss')then
        if(sc_method==1)then
           call dmft_weiss(Gmats,Smats,Weiss,Hloc)
        elseif(sc_method==2)then
           Weiss=zero
           Weiss(1,1,1,1,:)=zeta(1,:)-0.25d0*Wband(1)**2*Gimp(1,1,1,1,:);Weiss(1,1,1,1,:)=one/Weiss(1,1,1,1,:)
           Weiss(1,1,2,2,:)=zeta(2,:)-0.25d0*Wband(2)**2*Gimp(1,1,2,2,:);Weiss(1,1,2,2,:)=one/Weiss(1,1,2,2,:)
        endif
        if(MpiMaster)call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=1)
     else
        if(sc_method==1)then
           call dmft_delta(Gmats,Smats,Weiss,Hloc)
        elseif(sc_method==2)then
           Weiss=zero
           Weiss(1,1,1,1,:)=0.25d0*Wband(1)**2*Gimp(1,1,1,1,:)
           Weiss(1,1,2,2,:)=0.25d0*Wband(2)**2*Gimp(1,1,2,2,:)
        endif
        if(MpiMaster)call dmft_print_gf_matsubara(Weiss,"Delta",iprint=1)
     endif
     !
     !
     !Only the master does the fit and mix;  bcast at the end
     if(MpiMaster)then
        if(mix_method==1)then     !fit and mix the bath:
           call ed_chi2_fitgf(Weiss,bath,ispin=1)
           if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_prev
           Bath_prev=Bath        
        elseif(mix_method==2)then  !mix the weiss/delta and fit:
           if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_prev        
           call ed_chi2_fitgf(Weiss,bath,ispin=1)
           Weiss_prev=Weiss
        endif
     endif
     call Bcast_MPI(MpiComm,Bath)
     call Bcast_MPI(MpiComm,Weiss)
     !
     !
     !Master checks convergence and if required changes mu; Bcast at the end
     if(MpiMaster)then
        converged1 = check_convergence(Weiss(1,1,1,1,:),dmft_error,nsuccess,nloop,reset=.false.,index=1,total=2)
        converged2 = check_convergence(Weiss(1,1,2,2,:),dmft_error,nsuccess,nloop,reset=.false.,index=2,total=2)
        converged  = converged1.AND.converged2
        !
        if(nread/=0.d0) call search_chemical_potential(xmu,sum(dens),converged)
     endif
     call Bcast_MPI(MpiComm,converged)
     call Bcast_MPI(MpiComm,xmu)
     !
     !
     if(MpiMaster)call end_loop
  enddo
  !
  !Mpi computation of the Gloc on real-axis
  call dmft_gloc_realaxis(MpiComm,Ebands,Dbands,H0,Greal,Sreal)
  if(MpiMaster)call dmft_print_gf_realaxis(Greal,"Gloc",iprint=1)

  !Mpi computation of Kinetic energy
  call dmft_kinetic_energy(MpiComm,Ebands,Dbands,H0,Smats)


  call finalize_MPI()


end program hm_2bands_dos


