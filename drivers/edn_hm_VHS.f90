program hm_VHS
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
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
  complex(8),allocatable,dimension(:,:,:,:,:) :: Weiss,Smats,Sreal,Gmats,Greal,Weiss_
  complex(8),allocatable,dimension(:) :: Gtest
  character(len=16)                           :: finput
  real(8)                                     :: wmixing
  !
  real(8)                                     :: valen,valdos
  integer                                     :: ic,lendos,i
  logical                                     :: master
  logical :: betheSC,wGimp,mixG0,symOrbs


  call parse_cmd_variable(finput,"FINPUT",default='inputED.in')
  call parse_input_variable(Le,"LE",finput,default=500)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  call parse_input_variable(mixG0,"mixG0",finput,default=.false.)
  !
  call ed_read_input(trim(finput))

  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  if(Nspin/=1.OR.Norb/=1)stop "Wrong setup from input file: Nspin=1 OR Norb=1"
  Nso=Nspin*Norb

  allocate(Wband(Nso))
  allocate(H0(Nso))
  allocate(de(Nso))
  !
  lendos=file_length("dos.dat")
  if(lendos>0)then
     allocate(Ebands(1,lendos))
     allocate(Dbands(1,lendos))
     call sread("dos.dat",Ebands(1,:),Dbands(1,:))
     de(1) = Ebands(1,2)-Ebands(1,1)
     print*,"Int(DOS)=",simps(Dbands(1,:),Ebands(1,2)-Ebands(1,1))
     Dbands(1,:)=Dbands(1,:)*de(1)

  else

     lendos=Le
     allocate(Ebands(1,lendos))
     allocate(Dbands(1,lendos))
     Ebands(1,:) = linspace(wini,wfin,lendos,mesh=de(1));print*,de(1)
     do i=1,lendos
        Dbands(1,i) = dens_2dsquare(Ebands(1,i),1d0)
     enddo
  endif

  print*,"Int(DOS)=",simps(Dbands(1,:),Ebands(1,2)-Ebands(1,1))
  Dbands(1,:)=Dbands(1,:)*de(1)



  H0    = 0.0d0
  call TB_write_Hloc(one*diag(H0))


  !Allocate Weiss Field:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Weiss_(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Gtest(Lmats))
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  allocate(dens(Norb))
  Hloc(1,1,:,:)=diag(H0)


  !setup solver
  Nb=get_bath_dimension()
  allocate(bath(Nb))
  allocate(bath_(Nb))
  call ed_init_solver(bath,Hloc)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")
     !
     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath)
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_real(Sreal)
     !
     ! compute the local gf:
     call dmft_gloc_matsubara(Ebands,Dbands,H0,Gmats,Smats)
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)
     !
     call dmft_gloc_realaxis(Ebands,Dbands,H0,Greal,Sreal)
     call dmft_print_gf_realaxis(Greal,"Greal",iprint=1)
     !
     !Get the Weiss field/Delta function to be fitted
     call dmft_self_consistency(Gmats,Smats,Weiss,Hloc,SCtype=cg_scheme)     
     call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=1)
     !
     !
     !
     if(mixG0)then
        if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_
        Weiss_=Weiss
     endif
     !
     !Perform the SELF-CONSISTENCY by fitting the new bath
     call ed_chi2_fitgf(Weiss,bath,ispin=1)
     !
     !MIXING:
     if(.not.mixG0)then
        if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
        Bath_=Bath
     endif
     !
     !Check convergence (if required change chemical potential)
     Gtest=zero
     do iorb=1,Norb
        Gtest=Gtest+Weiss(1,1,iorb,iorb,:)/Norb
     enddo
     converged = check_convergence(Gtest,dmft_error,nsuccess,nloop,reset=.false.)

     if(nread/=0d0)then
        call ed_get_dens(dens)
        call ed_search_variable(xmu,sum(dens),converged)
        !call search_chemical_potential(xmu,sum(dens),converged)
     endif


     call end_loop
  enddo


  call dmft_gloc_realaxis(Ebands,Dbands,H0,Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Greal",iprint=1)
  call dmft_kinetic_energy(Ebands,Dbands,H0,Smats(1,1,:,:,:))


end program hm_VHS



