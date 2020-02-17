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
  complex(8),allocatable,dimension(:,:,:,:,:) :: Weiss,Smats,Sreal,Gmats,Greal,Weiss_
  complex(8),allocatable,dimension(:)         :: Gtest
  character(len=16)                           :: finput
  real(8)                                     :: wmixing
  !
  integer                                     :: comm,rank
  logical                                     :: master
  logical                                     :: betheSC,wGimp,mixG0,symOrbs
  integer                                     :: isolve,info
  real(8)                                     :: impHloc
  real(8),dimension(1)                        :: aparams


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
  call parse_input_variable(betheSC,"BETHESC",finput,default=.false.)
  call parse_input_variable(wGimp,"wGimp",finput,default=.false.)
  call parse_input_variable(mixG0,"mixG0",finput,default=.false.)
  call parse_input_variable(symOrbs,"symOrbs",finput,default=.false.)
  call parse_input_variable(isolve,"isolve",finput,default=0,comment="isolve:")
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
  call ed_init_solver(comm,bath,Hloc)

  !LEGEND:
  !isolve==1: solve using Broyden within DMFT loop
  !isolve==2: solve using Fsolve within DMFT loop
  !isolve==3: regular DMFT loop using ed_search_variable
  if(isolve==0)then
     iloop=0;converged=.false.

     iloop=0;converged=.false.
     do while(.not.converged.AND.iloop<nloop)
        iloop=iloop+1
        if(master)call start_loop(iloop,nloop,"DMFT-loop")
        !
        aparams=[xmu]
        call broyden1(solve_loop_params,aparams,tol=nerr,tol1=nerr)
        ! call fsolve(solve_loop_params,aparams,tol=nerr,info=info)
        ! print*,"FSOLVE info=",info
        !
        !Check convergence (if required change chemical potential)
        if(master)then
           Gtest=zero
           do iorb=1,Norb
              Gtest=Gtest+Weiss(1,1,iorb,iorb,:)/Norb
           enddo
           converged = check_convergence(Gtest,dmft_error,nsuccess,nloop,reset=.false.)
           converged = converged.AND.(abs(sum(dens)-nread)<=nerr)
        endif
        call Bcast_MPI(comm,converged)
        !
        if(master)call end_loop
     enddo

  else
     !
     iloop=0;converged=.false.
     do while(.not.converged.AND.iloop<nloop)
        iloop=iloop+1
        if(master)call start_loop(iloop,nloop,"DMFT-loop")
        !
        !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
        if(isolve==1)then
           aparams=[xmu]
           call broyden1(solve_dmft_params,aparams,tol=nerr,tol1=nerr)
        elseif(isolve==2)then
           aparams=[xmu]
           call fsolve(solve_dmft_params,aparams,tol=nerr,info=info)
           print*,"FSOLVE info=",info
        endif
        !
        call ed_solve(comm,bath,sflag=.true.)
        call ed_get_sigma_matsubara(Smats)
        !
        ! compute the local gf:
        call dmft_gloc_matsubara(comm,Ebands,Dbands,H0,Gmats,Smats)
        !
        !Get the Weiss field/Delta function to be fitted
        if(.not.betheSC)then
           call dmft_self_consistency(comm,Gmats,Smats,Weiss,Hloc,SCtype=cg_scheme)
        else
           if(wGimp)call ed_get_gimp_matsubara(Gmats)
           call dmft_self_consistency(comm,Gmats,Weiss,Hloc,SCtype=cg_scheme,wbands=Wband)
        endif
        !
        !
        !
        if(mixG0)then
           if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_
           Weiss_=Weiss
        endif
        !
        !Perform the SELF-CONSISTENCY by fitting the new bath
        if(symOrbs)then
           call ed_chi2_fitgf(comm,Weiss,bath,ispin=1,iorb=1)
           call orb_equality_bath(bath,save=.true.)
        else
           call ed_chi2_fitgf(comm,Weiss,bath,ispin=1)
        endif
        !
        !MIXING:
        if(.not.mixG0)then
           if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
           Bath_=Bath
        endif
        !
        !Check convergence (if required change chemical potential)
        if(master)then
           Gtest=zero
           do iorb=1,Norb
              Gtest=Gtest+Weiss(1,1,iorb,iorb,:)/Norb
           enddo
           converged = check_convergence(Gtest,dmft_error,nsuccess,nloop,reset=.false.)
           if(isolve==3)then
              if(nread/=0d0)then
                 call ed_get_dens(dens)
                 call ed_search_variable(xmu,sum(dens),converged)
              endif
           endif
        endif
        call Bcast_MPI(comm,converged)
        call Bcast_MPI(comm,xmu)
        !
        if(master)call end_loop
     enddo

  endif


  call ed_get_sigma_realaxis(Sreal)
  call dmft_gloc_realaxis(Comm,Ebands,Dbands,H0,Greal,Sreal)
  call dmft_kinetic_energy(comm,Ebands,Dbands,H0,Smats(1,1,:,:,:))
  if(master)call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)
  if(master)call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=1)
  if(master)call dmft_print_gf_realaxis(Greal,"Greal",iprint=1)


  call finalize_MPI()


contains

  function solve_loop_params(params) result(fss)
    real(8),dimension(:),intent(in) :: params ![xmu]
    real(8),dimension(size(params)) :: fss    ![dens-nread]
    integer,save                    :: iter=0
    xmu = params(1)
    !
    iter=iter+1
    !
    call ed_solve(comm,bath)
    call ed_get_sigma_matsubara(Smats)
    !
    ! compute the local gf:
    call dmft_gloc_matsubara(comm,Ebands,Dbands,H0,Gmats,Smats)
    !
    !Get the Weiss field/Delta function to be fitted
    if(.not.betheSC)then
       call dmft_self_consistency(comm,Gmats,Smats,Weiss,Hloc,SCtype=cg_scheme)
    else
       if(wGimp)call ed_get_gimp_matsubara(Gmats)
       call dmft_self_consistency(comm,Gmats,Weiss,Hloc,SCtype=cg_scheme,wbands=Wband)
    endif
    !
    !
    if(mixG0)then
       if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_
       Weiss_=Weiss
    endif
    !
    !Perform the SELF-CONSISTENCY by fitting the new bath
    if(symOrbs)then
       call ed_chi2_fitgf(comm,Weiss,bath,ispin=1,iorb=1)
       call orb_equality_bath(bath,save=.true.)
    else
       call ed_chi2_fitgf(comm,Weiss,bath,ispin=1)
    endif
    !
    !MIXING:
    if(.not.mixG0)then
       if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
       Bath_=Bath
    endif
    !
    !
    call ed_get_dens(dens)
    fss(1) = sum(dens)-nread
    !
    open(1991,file="iter_optimize.ed",access='append')
    write(*,*)iloop,iter,xmu,fss
    write(1991,*)iloop,iter,xmu,fss,sum(dens),Nread
    close(1991)

  end function solve_loop_params





  function solve_dmft_params(params) result(fss)
    real(8),dimension(:),intent(in) :: params ![xmu,xh]
    real(8),dimension(size(params)) :: fss    ![dens-nread,Sz-Szread]==[dens-nread,Sz-0]
    integer,save                    :: iter=0
    xmu = params(1)
    !
    iter=iter+1
    !
    !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
    call ed_solve(comm,bath,sflag=.false.)
    !
    call ed_get_dens(dens)
    fss(1) = sum(dens) - nread

    open(1991,file="iter_optimize.ed",access='append')
    write(*,*)iloop,xmu,fss
    write(1991,*)iloop,xmu,fss,sum(dens),Nread
    close(1991)
  end function solve_dmft_params



end program hm_Nbands_bethe



