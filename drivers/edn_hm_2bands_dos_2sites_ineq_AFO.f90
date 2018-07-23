program hm_2bands_dos
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS,   only:free_unit,reg,txtfy
  USE SF_ARRAYS,    only:arange
  USE SF_TIMER
  USE SF_LINALG,    only:inv,eye,eigh,diag
  USE SF_MISC,      only:assert_shape
  USE SF_SPECIAL,   only:fermi
  USE DMFT_CTRL_VARS
  USE MPI
  implicit none
  integer                                         :: iloop,Nb,Le,Nso,unit1,ip,ispin
  integer                                         :: Nineq,Nlat,Nlso,unit,unit_tx,unit_ty
  logical                                         :: converged
  !Bath:
  real(8),allocatable                             :: Bath_ineq(:,:),Bath_prev(:,:)
  !
  real(8),dimension(2)                            :: Wband
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Weiss_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Smats,Smats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Sreal,Sreal_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Gmats,Gmats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Greal,Greal_ineq

  complex(8),allocatable,dimension(:,:,:,:,:)     :: Hloc
  complex(8),allocatable,dimension(:,:,:,:,:)     :: Hloc_ineq

  real(8),dimension(:,:),allocatable              :: Dbands
  real(8),dimension(:,:),allocatable              :: Ebands
  real(8),dimension(:),allocatable                :: de,dos_wt
  real(8),dimension(:,:,:,:),allocatable          :: H0_b
  real(8),dimension(:,:),allocatable              :: H0_a
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: He_b,Smatsnn
  complex(8),dimension(:,:,:),allocatable         :: He,He1,He2,Smats1,Smats2
  complex(8),dimension(:,:,:,:,:,:),allocatable   :: Smats11,Smats22
  real(8),dimension(:),allocatable                :: dens,wr,H00
  real(8)                                         :: Delta,hyb,ndelta_even,ndelta_odd,dr
  character(len=16)                               :: finput
  real(8)                                         :: wmixing,Eout1(2),Eout2(2),tx,ty
  character(len=10)                               :: dos_model
  character(len=4)                                :: dir_iter
  logical                                         :: spinsym,fullsym



  !Parse additional variables
  call parse_cmd_variable(finput,"FINPUT",default='inputED.in')
  call parse_input_variable(Le,"LE",finput,default=500)
  call parse_input_variable(Wband,"WBAND",finput,default=[1d0,0.5d0])
  call parse_input_variable(xmu,"XMU",finput,default=0.d0)
  call parse_input_variable(delta,"DELTA",finput,default=0d0)
  call parse_input_variable(dos_model,"DOS_MODEL",finput,"bethe")
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  call parse_input_variable(spinsym,"spinsym",finput,default=.true.)
  call parse_input_variable(fullsym,"fullsym",finput,default=.false.)
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

  Nso=Nspin*Norb
  Nlat=2                      !number of independent sites, 2 for AF0 ordering
  Nineq=Nlat
  Nlso=Nlat*Nso
  if(Norb/=2)stop "Norb != 2"
  if(Nspin/=2)stop "Nspin != 2"
  if(Nso/=4)stop "Nso != 4"
  if(Nlso/=8)stop "Nlso != 8"

  if(fullsym)then
     Nineq=1
     write(*,*)"Using Nineq sites=",Nineq
     open(free_unit(unit),file="symmetries.used")
     write(unit,*)"Symmetries used are:"
     write(unit,*)"(site=2,l,s)=(site=1,l,-s)"
     close(unit)
  endif

  if(spinsym)sb_field=0.d0

  allocate(Ebands(Nso,Le))
  allocate(Dbands(Nso,Le))
  allocate(de(Nso))
  Ebands(1,:) = linspace(-Wband(1),Wband(1),Le,mesh=de(1))
  Ebands(2,:) = linspace(-Wband(2),Wband(2),Le,mesh=de(2))
  Ebands(3,:) = Ebands(1,:)
  Ebands(4,:) = Ebands(2,:)
  !
  select case(dos_model)
  case ("bethe")
     Dbands(1,:) = dens_bethe(Ebands(1,:),Wband(1))*de(1)
     Dbands(2,:) = dens_bethe(Ebands(2,:),Wband(2))*de(2)
  case ("flat")
     Dbands(1,:) = dens_flat(Ebands(1,:),Wband(1))*de(1)
     Dbands(2,:) = dens_flat(Ebands(2,:),Wband(2))*de(2)
  case default
     stop "error: dos_model not in {bethe,flat}. Add your own if needed"
  end select
  Dbands(3,:) = Dbands(1,:)
  Dbands(4,:) = Dbands(2,:)

  allocate(H0_b(Nspin,Nspin,Norb,Norb))
  H0_b = 0.d0
  H0_b(1,1,1,1) = -Delta/2.d0
  H0_b(2,2,1,1) = H0_b(1,1,1,1)
  H0_b(1,1,2,2) = Delta/2.d0
  H0_b(2,2,2,2) = H0_b(1,1,2,2)

  allocate(H0_a(Nso,Nso))
  H0_a = d_nn2nso(H0_b,Nspin,Norb)
  allocate(H00(Nso))
  do ip=1,Nso
     H00(ip)=H0_a(ip,ip)
  enddo
  !  call TB_write_Hloc(H0)


  !Allocate Weiss Field:
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  !
  allocate(Weiss_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Gmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
  !
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))
  allocate(Hloc_ineq(Nineq,Nspin,Nspin,Norb,Norb))
  allocate(He(Nlso,Nlso,Le))
  allocate(He1(Nso,Nso,Le))
  allocate(He2(Nso,Nso,Le))

  call build_he_b()
  do ip=1,Le
     He(:,:,ip) = d_nnn2nlso(He_b(:,:,:,:,:,:,ip),Nlat,Nspin,Norb)
  enddo
  do ip=1,Le
     He1(:,:,ip) = d_nn2nso_c(He_b(1,1,:,:,:,:,ip),Nspin,Norb)
     He2(:,:,ip) = d_nn2nso_c(He_b(2,2,:,:,:,:,ip),Nspin,Norb)
  enddo
  Hloc = 0.d0
  Hloc(1,1,1,1,1)=H0_b(1,1,1,1)
  Hloc(1,2,2,1,1)=H0_b(1,1,1,1)
  Hloc(1,1,1,2,2)=H0_b(1,1,2,2)
  Hloc(1,2,2,2,2)=H0_b(1,1,2,2)
  Hloc(2,1,1,1,1)=H0_b(1,1,1,1)
  Hloc(2,2,2,1,1)=H0_b(1,1,1,1)
  Hloc(2,1,1,2,2)=H0_b(1,1,2,2)
  Hloc(2,2,2,2,2)=H0_b(1,1,2,2)
  do ip=1,Nineq
     Hloc_ineq(ip,:,:,:,:) = Hloc(ip,:,:,:,:)
  enddo
  allocate(dens(Norb))


  !setup solver
  Nb=get_bath_dimension()
  allocate(Bath_ineq(Nineq,Nb))
  allocate(Bath_prev(Nineq,Nb))

  call ed_init_solver(Bath_ineq,Hloc_ineq)
  do ip=1,Nineq
     call break_symmetry_bath(Bath_ineq(ip,:),sb_field,(-1d0)**(ip+1))
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")
     !
     !
     call ed_solve(Bath_ineq,Hloc_ineq)
     call ed_get_sigma_matsubara(Smats_ineq,Nineq)
     call ed_get_sigma_real(Sreal_ineq,Nineq)
     do ip=1,Nineq
        Smats(ip,:,:,:,:,:) = Smats_ineq(ip,:,:,:,:,:)
        Sreal(ip,:,:,:,:,:) = Sreal_ineq(ip,:,:,:,:,:)
     enddo
     if(fullsym)then
        do ispin=1,2
           Smats(2,ispin,ispin,:,:,:)=Smats(1,3-ispin,3-ispin,:,:,:)
           Sreal(2,ispin,ispin,:,:,:)=Sreal(1,3-ispin,3-ispin,:,:,:)
        enddo
     endif
     !
     !
     ! compute the local gf:
     call dmft_gloc_matsubara(He,dos_wt,Gmats,Smats) !tridiag option off
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=4)
     do ip=1,Nineq
        Gmats_ineq(ip,:,:,:,:,:) = Gmats(ip,:,:,:,:,:)
     enddo
     !
     !
     ! Compute the Weiss field
     ! compute the Weiss field (only the Nineq ones)
     if(cg_scheme=='weiss')then
        call dmft_weiss(Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq)
     else
        call dmft_delta(Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq)
     endif
     !
     !
     ! fit baths and mix result with old baths
     call ed_chi2_fitgf(Bath_ineq,Weiss_ineq,Hloc_ineq,ispin=1)
     if(spinsym)then
        call spin_symmetrize_bath(Bath_ineq,save=.true.)
     else
        call ed_chi2_fitgf(Bath_ineq,Weiss_ineq,Hloc_ineq,ispin=2)
     endif
     !
     !
     ! Mixing:
     if(iloop>1)Bath_ineq = wmixing*Bath_ineq + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath_ineq
     !
     ! Convergence
     converged = check_convergence(Weiss_ineq(:,1,1,1,1,:),dmft_error,nsuccess,nloop)
     !
     ! Chemical potential
     call ed_get_dens(dens)
     if(nread/=0d0)then
        call search_chemical_potential(xmu,sum(dens),converged)
        unit1=12
        open(unit=unit1,file="xmu.restart",form="FORMATTED",status="REPLACE",action="WRITE")
        write(unit=unit1,fmt=*) xmu
        close(unit=unit1)
     endif
     !
     call end_loop
  enddo


  call dmft_gloc_realaxis(He,dos_wt,Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=4)
  allocate(Smats1(Nso,Nso,Lmats))
  allocate(Smats2(Nso,Nso,Lmats))
  allocate(Smatsnn(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Le))

  do ip=1,Le
     Smats1(:,:,ip) = d_nn2nso_c(Smats(1,:,:,:,:,ip),Nspin,Norb)
     Smats2(:,:,ip) = d_nn2nso_c(Smats(2,:,:,:,:,ip),Nspin,Norb)
  enddo
  !Eout1 = dmft_kinetic_energy(He,dos_wt,H00,Smats1)
  !Eout2 = dmft_kinetic_energy(He,dos_wt,H00,Smats2)
  Eout1 = dmft_kinetic_energy_good(Ebands,Dbands,H00,Smats1)
  Eout2 = dmft_kinetic_energy_good(Ebands,Dbands,H00,Smats2)
  unit_tx=free_unit()
  open(unit_tx,file="ekin_site1.dat")
  write(unit_tx,*)Eout1
  close(unit_tx)
  unit_tx=free_unit()
  open(unit_tx,file="ekin_site2.dat")
  write(unit_tx,*)Eout2
  close(unit_tx)


  if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
  wr = linspace(wini,wfin,Lreal)
  dr=wr(2)-wr(1)
  tx=0.d0
  ty=0.d0
  do ip=1,Lreal-1
     if(wr(ip).le.0.d0)then
        tx=tx-2.d0*(dreal(Greal(1,1,1,1,2,ip))+dreal(Greal(1,1,1,1,2,ip+1)))*dr/(2.d0*3.14159265d0)
        ty=ty-2.d0*(dimag(Greal(1,1,1,1,2,ip))+dimag(Greal(1,1,1,1,2,ip+1)))*dr/(2.d0*3.14159265d0)
     endif
  enddo
  unit_tx=free_unit()
  open(unit_tx,file="tx.dat")
  write(unit_tx,*)tx
  close(unit_tx)
  unit_ty=free_unit()
  open(unit_ty,file="ty.dat")
  write(unit_ty,*)ty
  close(unit_ty)

  tx=0.d0
  ty=0.d0
  do ip=1,Lreal-1
     if(wr(ip).le.0.d0)then
        tx=tx-2.d0*(dreal(Greal(2,1,1,1,2,ip))+dreal(Greal(2,1,1,1,2,ip+1)))*dr/(2.d0*3.14159265d0)
        ty=ty-2.d0*(dimag(Greal(2,1,1,1,2,ip))+dimag(Greal(2,1,1,1,2,ip+1)))*dr/(2.d0*3.14159265d0)
     endif
  enddo
  unit_tx=free_unit()
  open(unit_tx,file="tx1.dat")
  write(unit_tx,*)tx
  close(unit_tx)
  unit_ty=free_unit()
  open(unit_ty,file="ty1.dat")
  write(unit_ty,*)ty
  close(unit_ty)


contains



  function dmft_kinetic_energy_good(Ebands,Dbands,Hloc,Sigma) result(Eout)
    real(8),dimension(:,:),intent(in)                           :: Ebands ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1),size(Ebands,2)),intent(in) :: Dbands ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1)),intent(in)                :: Hloc ![Nspin*Norb]
    complex(8),dimension(:,:,:)                                 :: Sigma ![Nspin*Norb][Nspin*Norb][L]
    !
    integer                                                     :: Lk,Nso,Liw
    integer                                                     :: i,ik,iso
    !
    integer                                                     :: Norb,Nporb
    integer                                                     :: Nspin
    real(8)                                                     :: beta
    real(8)                                                     :: xmu
    !
    real(8),dimension(size(Ebands,1),size(Ebands,1))            :: Sigma_HF
    !
    complex(8)                                                  :: Ak,Bk,Ck,Dk
    complex(8)                                                  :: Gk,Tk
    real(8)                                                     :: Tail0,Tail1,Lail0,Lail1,spin_degeneracy
    !
    real(8)                                                     :: H0,Hl
    real(8)                                                     :: Ekin,Eloc
    real(8)                                                     :: Eout(2)
    !
    real(8),dimension(:),allocatable                            :: wm
    complex(8),dimension(:),allocatable                         :: wmc
    !Retrieve parameters:
    call get_ctrl_var(Norb,"NORB")
    call get_ctrl_var(Nspin,"NSPIN")
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    Nso = size(Ebands,1)
    Lk  = size(Ebands,2)
    Liw = size(Sigma,3)
    !Testing:
    if(Nso/=Norb*Nspin)stop "dmft_kinetic_energy_normal_dos: Nso != Norb*Nspin [from Hk]"
    call assert_shape(Sigma,[Nso,Nso,Liw],"dmft_kinetic_energy_normal_main","Sigma")
    !
    !Allocate and setup the Matsubara freq.
    if(allocated(wm))deallocate(wm);allocate(wm(Liw))
    if(allocated(wmc))deallocate(wmc);allocate(wmc(Liw))
    wm = pi/beta*dble(2*arange(1,Liw)-1)
    wmc(:)=-wm(:)*(xi)**2.d0
    !write(*,*)'wm(1),wm(50)',wm(1),wm(50)
    !write(*,*)'xmu',xmu
    !
    !Get HF part of the self-energy
    Sigma_HF = dreal(Sigma(:,:,Liw))
    !
    !
    write(*,"(A)") "Kinetic energy computation"
    call start_timer()
    H0=0.d0
    Hl=0.d0
    !Get principal part: Tr[ Hk.(Gk-Tk) ]
    do ik=1,Lk
       do iso=1,2
          Ak=(0.d0,0.d0)
          Bk=(0.d0,0.d0)
          Ak = Ebands(iso,ik)
          Bk =-Ebands(iso,ik) - Sigma_HF(iso,iso)
          !if(abs(Bk).gt.100.d0)then
          !   write(*,*)'Bk',Bk
          !endif
          do i=1,Liw
             Tk=(0.d0,0.d0)
             Gk=(0.d0,0.d0)
             Ck=(0.d0,0.d0)
             Dk=(0.d0,0.d0)
             Gk = (xi*wmc(i)+xmu) - Sigma(iso,iso,i) - Ebands(iso,ik)
             !  if(abs(dreal(Gk)).lt.10.d0**(-5.d0).and.abs(dimag(Gk)).lt.10.d0**(-5.d0))then
             !     write(*,*)'Gk,1.d0/Gk',Gk,1.d0/Gk
             !  endif
             Gk = 1.d0/Gk
             Tk = -xi/(wmc(i)) + Bk/(wmc(i))**2.d0
             !  if(abs(Tk).gt.100.d0)then
             !     write(*,*)'-----------------------'
             !     write(*,*)'ik,i,iso,Tk',ik,i,iso,Tk
             !     write(*,*)'Bk,wm(i)',Bk,wmc(i)
             !     write(*,*)'Ak,Gk,Hloc(iso)',Ak,Gk,Hloc(iso)
             !     write(*,*)'Dbands(iso,ik)',Dbands(iso,ik)
             !     write(*,*)'-----------------------'
             !  endif
             Ck = Ak*(Gk - Tk)
             Dk = Hloc(iso)*(Gk - Tk)
             !  if(abs(Ck).gt.100.d0)then
             !     write(*,*)'Ck',Ck
             !  endif
             !if(abs(Dbands(iso,ik)).gt.100.d0)then
             !   write(*,*)'iso,ik,Dbands(iso,ik)',iso,ik,Dbands(iso,ik)
             !endif
             H0 = H0 + Dbands(iso,ik)*dreal(Ck)/beta
             Hl = Hl + Dbands(iso,ik)*dreal(Dk)/beta
          enddo
       enddo
       !if(ik.eq.30)then
       !   write(*,*)'ik,Gk,Tk,Ck,Dk,H0,Hl',ik,Gk,Tk,Ck,Dk,H0,Hl
       !   write(*,*)'ik,Ak,Bk',ik,Ak,Bk
       !endif
       !if(ik.eq.300)then
       !   write(*,*)'ik,Gk,Tk,Ck,Dk',ik,Gk,Tk,Ck,Dk
       !   write(*,*)'ik,H0,Hl',H0,Hl
       !   write(*,*)'ik,Ak,Bk',ik,Ak,Bk
       !endif
       call eta(ik,Lk)
    enddo
    call stop_timer()
    !write(*,*)'Sigma_HF(1,1)',Sigma_HF(1,1)
    !write(*,*)'Ak,Bk',Ak,Bk
    !write(*,*)'Dbands(1,60),Dbands(2,60)',Dbands(1,60),Dbands(2,60)
    !write(*,*)'Ebands(1,60),Ebands(2,60)',Ebands(1,60),Ebands(2,60)
    !write(*,*)'Gk,Tk,Ck,Dk',Gk,Tk,Ck,Dk
    !write(*,*)'H0,Hl',H0,Hl
    spin_degeneracy=2.d0     !2 if Nspin=1, 1 if Nspin=2
    !write(*,*)'spin_degeneracy,beta',spin_degeneracy,beta
    H0=H0*2.d0*spin_degeneracy
    Hl=Hl*2.d0*spin_degeneracy
    !
    !get tail subtracted contribution: Tr[ Hk.Tk ]
    Tail0=0.d0
    Tail1=0.d0
    Lail0=0.d0
    Lail1=0.d0
    do ik=1,Lk
       do iso=1,2
          Ak = Ebands(iso,ik)
          Bk =-Ebands(iso,ik) - Sigma_HF(iso,iso)
          Ck= Ak*Bk
          Dk= Hloc(iso)*Bk
          Tail0 = Tail0 + 0.5d0*Dbands(iso,ik)*dreal(Ak)
          Tail1 = Tail1 + 0.25d0*Dbands(iso,ik)*dreal(Ck)
          Lail0 = Lail0 + 0.5d0*Dbands(iso,ik)*Hloc(iso)
          Lail1 = Lail1 + 0.25d0*Dbands(iso,ik)*dreal(Dk)
       enddo
    enddo
    Tail0=Tail0*spin_degeneracy
    Tail1=Tail1*beta*spin_degeneracy
    Lail0=Lail0*spin_degeneracy
    Lail1=Lail1*beta*spin_degeneracy
    !
    Ekin=H0+Tail0+Tail1
    Eloc=Hl+Lail0+Lail1
    Eout = [Ekin,Eloc]
    write(*,*)'H0,Tail0,Tail1',H0,Tail0,Tail1
    write(*,*)'Hl,Lail0,Lail1',Hl,Lail0,Lail1
    write(*,*)'ekin,eloc',Ekin,Eloc
    !
    call write_kinetic_info()
    call write_kinetic_value(Eout)
    !
    deallocate(wm)
    deallocate(wmc)
  end function dmft_kinetic_energy_good


  !+-------------------------------------------------------------------+
  !PURPOSE  : write legend, i.e. info about columns 
  !+-------------------------------------------------------------------+
  subroutine write_kinetic_info()
    integer :: unit
    unit = free_unit()
    open(unit,file="dmft_kinetic_energy.info")
    write(unit,"(A1,90(A14,1X))")"#",reg(txtfy(1))//"<K>",reg(txtfy(2))//"<Eloc>"
    close(unit)
  end subroutine write_kinetic_info


  !+-------------------------------------------------------------------+
  !PURPOSE  : Write energies to file
  !+-------------------------------------------------------------------+
  subroutine write_kinetic_value(Eout)
    real(8) :: Eout(2),Ekin,Eloc
    integer :: unit
    unit = free_unit()
    Ekin=Eout(1)
    Eloc=Eout(2)
    open(unit,file="dmft_kinetic_energy.dat")
    write(unit,"(90F15.9)")Ekin,Eloc
    close(unit)
  end subroutine write_kinetic_value

  function dens_flat(ebands,wband) result(dens)
    real(8),dimension(:)            :: ebands
    real(8)                         :: wband
    real(8),dimension(size(ebands)) :: dens
    integer                         :: i
    real(8)                         :: e
    do i=1,size(ebands)
       e=ebands(i)
       dens(i)= step(wband-abs(e))/(2*wband)
    enddo
  end function dens_flat


  function dens_flat_one(ee,wband) result(dens) 
    real(8)                         :: ee,wband
    real(8)                         :: dens
    dens = step(wband-abs(ee))/(2*wband)
  end function dens_flat_one

  function dens_bethe_one(ee,wband) result(dens) 
    real(8)                         :: ee,wband
    real(8)                         :: dens
    dens = 2.d0*sqrt(wband**2.d0 - ee**2.d0)/(pi*wband**2.d0)*step(wband-abs(ee))
  end function dens_bethe_one


  subroutine build_he_b()
    integer                :: ie,ispin,unit_dos1,unit_dos2
    real(8)                :: de,e,alp
    allocate(He_b(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Le))
    allocate(dos_wt(Le))
    alp=Wband(2)/Wband(1)
    de=2.d0*Wband(1)/dble(Le)
    do ie=1,Le
       e = -Wband(1) + dble(ie-1)*de
       He_b(:,:,:,:,:,:,ie) = 0.d0
       He_b(1,1,1,1,1,1,ie) = H0_a(1,1)
       He_b(1,1,2,2,1,1,ie) = H0_a(1,1)
       He_b(1,1,1,1,2,2,ie) = H0_a(2,2)
       He_b(1,1,2,2,2,2,ie) = H0_a(2,2)
       He_b(2,2,1,1,1,1,ie) = H0_a(1,1)
       He_b(2,2,2,2,1,1,ie) = H0_a(1,1)
       He_b(2,2,1,1,2,2,ie) = H0_a(2,2)
       He_b(2,2,2,2,2,2,ie) = H0_a(2,2)
       He_b(1,2,1,1,1,1,ie) = e
       He_b(1,2,2,2,1,1,ie) = e
       He_b(2,1,1,1,1,1,ie) = He_b(1,2,1,1,1,1,ie)
       He_b(2,1,2,2,1,1,ie) = He_b(1,2,2,2,1,1,ie)
       He_b(1,2,1,1,2,2,ie) = alp*e
       He_b(1,2,2,2,2,2,ie) = alp*e
       He_b(2,1,1,1,2,2,ie) = He_b(1,2,1,1,2,2,ie)
       He_b(2,1,2,2,2,2,ie) = He_b(1,2,2,2,2,2,ie)
       select case(dos_model)
       case ("bethe")
          dos_wt(ie)=de*dens_bethe_one(e,Wband(1))
       case ("flat")
          dos_wt(ie)=de*dens_flat_one(e,Wband(1))
       end select
    enddo
  end subroutine build_he_b


  function d_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,jlat,js
    Hlso=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                      js = jorb + (jspin-1)*Norb + (jlat-1)*Norb*Nspin !lattice-spin-orbit stride
                      Hlso(is,js) = Hnnn(ilat,jlat,ispin,jspin,iorb,jorb)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nnn2nlso

  function d_nn2nso(Hnnn,Nspin,Norb) result(Hso)
    real(8),dimension(Nspin,Nspin,Norb,Norb)              :: Hnnn
    integer                                               :: Nspin,Norb
    real(8),dimension(Nspin*Norb,Nspin*Norb)              :: Hso
    integer                                               :: iorb,ispin,is
    integer                                               :: jorb,jspin,js
    Hso=zero

    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function d_nn2nso

  function d_nn2nso_c(Hnnn,Nspin,Norb) result(Hso)
    complex(8),dimension(Nspin,Nspin,Norb,Norb)           :: Hnnn
    integer                                               :: Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb)           :: Hso
    integer                                               :: iorb,ispin,is
    integer                                               :: jorb,jspin,js
    Hso=zero

    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function d_nn2nso_c


  !  subroutine get_gloc_mats
  !    integer                                     :: i,ie,unit_Gm1,unit_Gm2
  !    complex(8),dimension(:,:),allocatable       :: zeta,fg
  !    complex(8),dimension(:,:),allocatable       :: Hh
  !    real(8),dimension(:),allocatable            :: wm
  !    !
  !    allocate(wm(Lmats))
  !    allocate(zeta(Norb,Norb))
  !    allocate(Hh(Norb,Norb))
  !    allocate(fg(Norb,Norb))
  !    !
  !    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
  !    !
  !    unit_Gm1=free_unit()
  !    open(unit_Gm1,file="newGloc_l1_s1_iw.dat")
  !    !
  !    do i=1,Lmats
  !       zeta(1,1)= xi*wm(i) + xmu - Smats(1,1,1,1,i)
  !       zeta(1,2)= - Smats(1,1,1,2,i)
  !       zeta(2,1)= - Smats(1,1,2,1,i)
  !       zeta(2,2)= xi*wm(i) + xmu - Smats(1,1,2,2,i)
  !       fg=zero
  !       do ie=1,Le
  !          Hh(:,:)=He(:,:,ie)
  !          fg=fg+inverse_ge(zeta,Hh)*dos_wt(ie)
  !       enddo
  !       Gmats(1,1,:,:,i)  = fg(:,:)
  !       write(unit_Gm1,*) wm(i),dreal(Gmats(1,1,1,1,i)),dimag(Gmats(1,1,1,1,i))
  !       !
  !    enddo
  !    close(unit_Gm1)
  !    !
  !    unit_Gm2=free_unit()
  !    open(unit_Gm2,file="newGloc_l2_s1_iw.dat")
  !    do i=1,Lmats
  !       write(unit_Gm2,*) wm(i),dreal(Gmats(1,1,2,2,i)),dimag(Gmats(1,1,2,2,i))
  !    enddo
  !    close(unit_Gm2)
  !    !
  !  end subroutine get_gloc_mats



  !  subroutine get_gloc_real
  !    integer                                     :: i,ie,unit_Gr1,unit_Gr2
  !    complex(8),dimension(:,:),allocatable       :: zeta,fg
  !    complex(8),dimension(:,:),allocatable       :: Hh
  !    real(8),dimension(:),allocatable            :: wr
  !    !
  !    allocate(wr(Lreal))
  !    allocate(zeta(Norb,Norb))
  !    allocate(Hh(Norb,Norb))
  !    allocate(fg(Norb,Norb))
  !    !
  !    wr = linspace(wini,wfin,Lreal)
  !    !
  !    unit_Gr1=free_unit()
  !    open(unit_Gr1,file="newGloc_l1_s1_realw.dat")
  !    !
  !    do i=1,Lreal
  !       zeta(1,1)= dcmplx(wr(i),eps) + xmu - Sreal(1,1,1,1,i)
  !       zeta(1,2)= - Sreal(1,1,1,2,i)
  !       zeta(2,1)= - Sreal(1,1,2,1,i)
  !       zeta(2,2)= dcmplx(wr(i),eps) + xmu - Sreal(1,1,2,2,i)
  !       fg=zero
  !       do ie=1,Le
  !          Hh(:,:)=He(:,:,ie)
  !          fg=fg+inverse_ge(zeta,Hh)*dos_wt(ie)
  !       enddo
  !       Greal(1,1,:,:,i)  = fg(:,:)
  !       write(unit_Gr1,*) wr(i),dreal(Greal(1,1,1,1,i)),dimag(Greal(1,1,1,1,i))
  !       !
  !    enddo
  !    close(unit_Gr1)
  !    !
  !    unit_Gr2=free_unit()
  !    open(unit_Gr2,file="newGloc_l2_s1_realw.dat")
  !    do i=1,Lreal
  !       write(unit_Gr2,*) wr(i),dreal(Greal(1,1,2,2,i)),dimag(Greal(1,1,2,2,i))
  !    enddo
  !    close(unit_Gr2)
  !    !
  !  end subroutine get_gloc_real



  function inverse_ge(zeta,Hh) result(ge)
    complex(8),dimension(2,2)   :: Hh
    complex(8),dimension(2,2)   :: zeta
    complex(8),dimension(2,2)   :: ge
    complex(8)                  :: a11,a12,a22,det
    ge=zero
    a11 = zeta(1,1) - Hh(1,1)
    a22 = zeta(2,2) - Hh(2,2)
    a12 = 0.
    det = a11*a22
    ge(1,1) = a22/det
    ge(2,2) = a11/det
    ge(2,1) = 0.
    ge(1,2) = 0.
  end function inverse_ge



  function inverse_ge_diag(zeta,Hh) result(ge)
    complex(8),dimension(2,2)   :: Hh
    complex(8),dimension(2,2)   :: zeta
    complex(8),dimension(2,2)   :: ge
    complex(8)                  :: a11,a12,a22,det
    ge=zero
    a11 = zeta(1,1) - Hh(1,1)
    a22 = zeta(2,2) - Hh(2,2)
    a12 = zeta(1,2) - Hh(1,2)
    det = a11*a22 - a12*conjg(a12)
    ge(1,1) = ((a11 + a22)/2.d0 - sqrt(((a11-a22)/2.d0)**2.d0 + a12*conjg(a12)))/det
    ge(2,2) = ((a11 + a22)/2.d0 + sqrt(((a11-a22)/2.d0)**2.d0 + a12*conjg(a12)))/det
    ge(2,1) = 0.d0
    ge(1,2) = 0.d0
  end function inverse_ge_diag


  ! subroutine search_chemical_potential_even(var_even,ntmp_even,converged_even)
  !   real(8),intent(inout) :: var_even
  !   real(8),intent(in)    :: ntmp_even
  !   logical,intent(inout) :: converged_even
  !   logical               :: bool_even
  !   real(8)               :: ndiff_even
  !   integer,save          :: count_even=0,totcount_even=0,i
  !   integer,save          :: nindex_even=0
  !   integer               :: nindex_old_even(3)
  !   real(8)               :: ndelta_old_even,nratio_even
  !   integer,save          :: nth_magnitude_even=-2,nth_magnitude_old_even=-2
  !   real(8),save          :: nth_even=1.d-2
  !   logical,save          :: ireduce_even=.true.
  !   integer               :: unit_even
  !   !
  !   ndiff_even=ntmp_even-nread
  !   nratio_even = 0.5d0;!nratio = 1.d0/(6.d0/11.d0*pi)
  !   !
  !   !check actual value of the density *ntmp* with respect to goal value *nread*
  !   count_even=count_even+1
  !   totcount_even=totcount_even+1
  !   if(count_even>2)then
  !      do i=1,2
  !         nindex_old_even(i+1)=nindex_old_even(i)
  !      enddo
  !   endif
  !   nindex_old_even(1)=nindex_even
  !   !
  !   if(ndiff_even >= nth_even)then
  !      nindex_even=-1
  !   elseif(ndiff_even <= -nth_even)then
  !      nindex_even=1
  !   else
  !      nindex_even=0
  !   endif
  !   !
  !   ndelta_old_even=ndelta_even
  !   bool_even=nindex_even/=0.AND.( (nindex_even+nindex_old_even(1)==0).OR.(nindex_even+sum(nindex_old_even(:))==0) )
  !   !if(nindex_old(1)+nindex==0.AND.nindex/=0)then !avoid loop forth and back
  !   if(bool_even)then
  !      ndelta_even=ndelta_old_even*nratio_even !decreasing the step
  !   else
  !      ndelta_even=ndelta_old_even
  !   endif
  !   !
  !   if(ndelta_old_even<1.d-9)then
  !      ndelta_old_even=0.d0
  !      nindex_even=0
  !   endif
  !   !update chemical potential
  !   var_even=var_even+dble(nindex_even)*ndelta_even
  !   !xmu=xmu+dble(nindex)*ndelta
  !   !
  !   !Print information
  !   write(LOGfile,*)"EVEN"
  !   write(LOGfile,"(A,f16.9,A,f15.9)")"n    = ",ntmp_even," /",nread
  !   if(nindex_even>0)then
  !      write(LOGfile,"(A,es16.9,A)")"shift= ",nindex_even*ndelta_even," ==>"
  !   elseif(nindex_even<0)then
  !      write(LOGfile,"(A,es16.9,A)")"shift= ",nindex_even*ndelta_even," <=="
  !   else
  !      write(LOGfile,"(A,es16.9,A)")"shift= ",nindex_even*ndelta_even," == "
  !   endif
  !   write(LOGfile,"(A,f15.9)")"var  = ",var_even
  !   write(LOGfile,"(A,ES16.9,A,ES16.9)")"dn   = ",ndiff_even,"/",nth_even
  !   unit_even=free_unit()
  !   open(unit_even,file="search_mu_iteration_even"//reg(ed_file_suffix)//".ed",position="append")
  !   write(unit_even,*)var_even,ntmp_even,ndiff_even
  !   close(unit_even)
  !   !
  !   !check convergence within actual threshold
  !   !if reduce is activetd
  !   !if density is in the actual threshold
  !   !if DMFT is converged
  !   !if threshold is larger than nerror (i.e. this is not last loop)
  !   bool_even=ireduce_even.AND.(abs(ndiff_even)<nth_even).AND.converged_even.AND.(nth_even>nerr)
  !   if(bool_even)then
  !      nth_magnitude_old_even=nth_magnitude_even        !save old threshold magnitude
  !      nth_magnitude_even=nth_magnitude_old_even-1      !decrease threshold magnitude || floor(log10(abs(ntmp-nread)))
  !      nth_even=max(nerr,10.d0**(nth_magnitude_even))        !set the new threshold 
  !      count_even=0                                     !reset the counter
  !      converged_even=.false.                                !reset convergence
  !      ndelta_even=ndelta_old_even*nratio_even               !reduce the delta step
  !      !
  !   endif
  !   !
  !   !if density is not converged set convergence to .false.
  !   if(abs(ntmp_even-nread)>nth_even)converged_even=.false.
  !   !
  !   !check convergence for this threshold
  !   !!---if smallest threshold-- NO MORE
  !   !if reduce is active (you reduced the treshold at least once)
  !   !if # iterations > max number
  !   !if not yet converged
  !   !set threshold back to the previous larger one.
  !   !bool=(nth==nerr).AND.ireduce.AND.(count>niter).AND.(.not.converged)
  !   bool_even=ireduce_even.AND.(count_even>niter).AND.(.not.converged_even)
  !   if(bool_even)then
  !      ireduce_even=.false.
  !      nth_even=10.d0**(nth_magnitude_old_even)
  !   endif
  !   !
  !   write(LOGfile,"(A,I5)")"count= ",count_even
  !   write(LOGfile,"(A,L2)"),"Converged=",converged_even
  !   print*,""
  !   !
  ! end subroutine search_chemical_potential_even


  ! subroutine search_chemical_potential_odd(var_odd,ntmp_odd,converged_odd)
  !   real(8),intent(inout) :: var_odd
  !   real(8),intent(in)    :: ntmp_odd
  !   logical,intent(inout) :: converged_odd
  !   logical               :: bool_odd
  !   real(8)               :: ndiff_odd
  !   integer,save          :: count_odd=0,totcount_odd=0,j
  !   integer,save          :: nindex_odd=0
  !   integer               :: nindex_old_odd(3)
  !   real(8)               :: ndelta_old_odd,nratio_odd
  !   integer,save          :: nth_magnitude_odd=-2,nth_magnitude_old_odd=-2
  !   real(8),save          :: nth_odd=1.d-2
  !   logical,save          :: ireduce_odd=.true.
  !   integer               :: unit_odd
  !   !
  !   ndiff_odd=ntmp_odd-nread
  !   nratio_odd = 0.5d0;!nratio = 1.d0/(6.d0/11.d0*pi)
  !   !
  !   !check actual value of the density *ntmp* with respect to goal value *nread*
  !   count_odd=count_odd+1
  !   totcount_odd=totcount_odd+1
  !   if(count_odd>2)then
  !      do j=1,2
  !         nindex_old_odd(j+1)=nindex_old_odd(j)
  !      enddo
  !   endif
  !   nindex_old_odd(1)=nindex_odd
  !   !
  !   if(ndiff_odd >= nth_odd)then
  !      nindex_odd=-1
  !   elseif(ndiff_odd <= -nth_odd)then
  !      nindex_odd=1
  !   else
  !      nindex_odd=0
  !   endif
  !   !
  !   ndelta_old_odd=ndelta_odd
  !   bool_odd=nindex_odd/=0.AND.( (nindex_odd+nindex_old_odd(1)==0).OR.(nindex_odd+sum(nindex_old_odd(:))==0) )
  !   !if(nindex_old(1)+nindex==0.AND.nindex/=0)then !avoid loop forth and back
  !   if(bool_odd)then
  !      ndelta_odd=ndelta_old_odd*nratio_odd !decreasing the step
  !   else
  !      ndelta_odd=ndelta_old_odd
  !   endif
  !   !
  !   if(ndelta_old_odd<1.d-9)then
  !      ndelta_old_odd=0.d0
  !      nindex_odd=0
  !   endif
  !   !update chemical potential
  !   var_odd=var_odd+dble(nindex_odd)*ndelta_odd
  !   !xmu=xmu+dble(nindex)*ndelta
  !   !
  !   !Print information
  !   write(LOGfile,*)"ODD"
  !   write(LOGfile,"(A,f16.9,A,f15.9)")"n    = ",ntmp_odd," /",nread
  !   if(nindex_odd>0)then
  !      write(LOGfile,"(A,es16.9,A)")"shift= ",nindex_odd*ndelta_odd," ==>"
  !   elseif(nindex_odd<0)then
  !      write(LOGfile,"(A,es16.9,A)")"shift= ",nindex_odd*ndelta_odd," <=="
  !   else
  !      write(LOGfile,"(A,es16.9,A)")"shift= ",nindex_odd*ndelta_odd," == "
  !   endif
  !   write(LOGfile,"(A,f15.9)")"var  = ",var_odd
  !   write(LOGfile,"(A,ES16.9,A,ES16.9)")"dn   = ",ndiff_odd,"/",nth_odd
  !   unit_odd=free_unit()
  !   open(unit_odd,file="search_mu_iteration"//reg(ed_file_suffix)//".ed",position="append")
  !   write(unit_odd,*)var_odd,ntmp_odd,ndiff_odd
  !   close(unit_odd)
  !   !
  !   !check convergence within actual threshold
  !   !if reduce is activetd
  !   !if density is in the actual threshold
  !   !if DMFT is converged
  !   !if threshold is larger than nerror (i.e. this is not last loop)
  !   bool_odd=ireduce_odd.AND.(abs(ndiff_odd)<nth_odd).AND.converged_odd.AND.(nth_odd>nerr)
  !   if(bool_odd)then
  !      nth_magnitude_old_odd=nth_magnitude_odd        !save old threshold magnitude
  !      nth_magnitude_odd=nth_magnitude_old_odd-1      !decrease threshold magnitude || floor(log10(abs(ntmp-nread)))
  !      nth_odd=max(nerr,10.d0**(nth_magnitude_odd))   !set the new threshold 
  !      count_odd=0                                !reset the counter
  !      converged_odd=.false.                      !reset convergence
  !      ndelta_odd=ndelta_old_odd*nratio_odd                  !reduce the delta step
  !      !
  !   endif
  !   !
  !   !if density is not converged set convergence to .false.
  !   if(abs(ntmp_odd-nread)>nth_odd)converged_odd=.false.
  !   !
  !   !check convergence for this threshold
  !   !!---if smallest threshold-- NO MORE
  !   !if reduce is active (you reduced the treshold at least once)
  !   !if # iterations > max number
  !   !if not yet converged
  !   !set threshold back to the previous larger one.
  !   !bool=(nth==nerr).AND.ireduce.AND.(count>niter).AND.(.not.converged)
  !   bool_odd=ireduce_odd.AND.(count_odd>niter).AND.(.not.converged_odd)
  !   if(bool_odd)then
  !      ireduce_odd=.false.
  !      nth_odd=10.d0**(nth_magnitude_old_odd)
  !   endif
  !   !
  !   write(LOGfile,"(A,I5)")"count= ",count_odd
  !   write(LOGfile,"(A,L2)"),"Converged=",converged_odd
  !   print*,""
  !   !
  ! end subroutine search_chemical_potential_odd



end program hm_2bands_dos



