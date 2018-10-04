program hm_2bands_dos
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  integer                                     :: iloop,Nb,Le,Nso,unit1
  logical                                     :: converged
  !Bath:
  real(8),allocatable                         :: Bath(:),Bath_(:)
  !
  real(8),dimension(2)                        :: Wband
  complex(8),allocatable                      :: Hloc(:,:,:,:)
  real(8),dimension(:,:),allocatable          :: Dbands
  real(8),dimension(:,:),allocatable          :: Ebands
  real(8),dimension(:),allocatable            :: H0,de
  real(8),dimension(:),allocatable            :: dens
  real(8)                                     :: Delta
  complex(8),allocatable,dimension(:,:,:,:,:) :: Weiss,Smats,Sreal,Gmats,Greal
  character(len=16)                           :: finput
  real(8)                                     :: wmixing,Eout(2)
  real(8)                                     :: norm1,norm2,x1,x2,lambda,cfp,www,wlx
  character(len=10)                           :: dos_model
  !
  integer                                     :: comm,rank
  logical                                     :: master

  call init_MPI(comm,msg=.true.)
  rank   = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  call parse_cmd_variable(finput,"FINPUT",default='inputED.in')
  call parse_input_variable(Le,"LE",finput,default=500)
  call parse_input_variable(x1,"X1",finput,default=0.0d0)
  call parse_input_variable(x2,"X2",finput,default=0.0d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=1.0d0)
  call parse_input_variable(cfp,"CFP",finput,default=0.1d0)
  call parse_input_variable(Wband,"WBAND",finput,default=[1d0,0.5d0])
  call parse_input_variable(delta,"DELTA",finput,default=0.d0)
  call parse_input_variable(dos_model,"DOS_MODEL",finput,"bethe")
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  !
  call ed_read_input(trim(finput),comm)

  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(delta,"delta")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  if(Nspin/=1.OR.Norb/=2)stop "Wrong setup from input file: Nspin=1; Norb=2"
  Nso=Nspin*Norb

  delta = delta + cfp*x2**2.d0 !I include the phoninic crystal field in the crystal field of the metal

  allocate(Ebands(Nso,Le))
  allocate(Dbands(Nso,Le))
  allocate(de(Nso))

  wlx = abs(lambda*x1)
  www = sqrt(Wband(1)**2.d0 + wlx**2.d0)
  de(1) = (www - wlx)/(Le/2.d0 -1.d0)
  do iloop=1,Le/2
     Ebands(1,iloop) = -www + (iloop - 1)*de(1)
     Ebands(1,Le-iloop+1) = www - (iloop - 1)*de(1)
  enddo
  Ebands(2,:) = linspace(-Wband(2),Wband(2),Le,mesh=de(2))
  !
  select case(dos_model)
  case ("bethe")
     Dbands(2,:) = dens_bethe(Ebands(2,:),Wband(2))*de(2)
     Dbands(1,:) = dens_peaks_phon(Ebands(1,:),wlx,Wband(1),x1)*de(1)
  case ("flat")
     Dbands(1,:) = dens_flat(Ebands(1,:),Wband(1))*de(1)
     Dbands(2,:) = dens_flat(Ebands(2,:),Wband(2))*de(2)
  case default
     stop "error: dos_model not in {bethe,flat}. Add your own if needed"
  end select

  do iloop=1,Le
     if(Dbands(1,iloop)/de(1).gt.20.d0)then
        Dbands(1,iloop)=0.d0
     endif
  enddo

  norm1 = 0.d0
  norm2 = 0.d0
  do iloop=1,Le/2-1
     norm1 = norm1 + (Dbands(1,iloop) + Dbands(1,iloop+1))/2.d0
     norm2 = norm2 + (Dbands(1,Le-iloop) + Dbands(1,Le+1-iloop))/2.d0
  enddo

  do iloop=1,Le/2
     Dbands(1,iloop) = Dbands(1,iloop)/(2.d0*norm1)
     Dbands(1,Le+1-iloop) = Dbands(1,Le+1-iloop)/(2.d0*norm2)
  enddo

  allocate(H0(Nso))
  H0=[-Delta/2,Delta/2]
  call TB_write_Hloc(one*diag(H0))


  !Allocate Weiss Field:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  Hloc(1,1,:,:)=diag(H0)
  allocate(dens(Norb))


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
     if(master)call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=2)
     !
     !Get the Weiss field/Delta function to be fitted
     if(master)call dmft_self_consistency(Gmats,Smats,Weiss,Hloc,SCtype=cg_scheme)
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

        !Check convergence (if required change chemical potential)
        converged = check_convergence(Weiss(1,1,1,1,1:min(250,Lmats)),dmft_error,nsuccess,nloop,reset=.false.)
        !
        call ed_get_dens(dens)
        if(nread/=0d0) then
           call search_chemical_potential(xmu,sum(dens),converged)
           unit1=12
           open(unit=unit1,file="xmu.restart",form="FORMATTED",status="REPLACE",action="WRITE")
           write(unit=unit1,fmt=*) xmu
           close(unit=unit1)
        endif
     endif
     call Bcast_MPI(comm,bath)
     call Bcast_MPI(comm,converged)
     call Bcast_MPI(comm,xmu)
     !
     if(master)call end_loop
  enddo

  ! if(master)then
  !    unit1=12
  !    open(unit=unit1,file="xmu.restart",form="FORMATTED",status="REPLACE",action="WRITE")
  !    write(unit=unit1,fmt=*) xmu
  !    close(unit=unit1)
  ! endif

  call dmft_gloc_realaxis(comm,Ebands,Dbands,H0,Greal,Sreal)
  if(master)call dmft_print_gf_realaxis(Greal,"Greal",iprint=2)
  call dmft_kinetic_energy(comm,Ebands,Dbands,H0,Smats(1,1,:,:,:))


  call finalize_MPI()
  

contains



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

  function dens_peaks_phon(ebands,wlx,wband,xx) result(dens)
    real(8),dimension(:)             :: ebands
    real(8)                          :: wlx,e,wband,epss,xx
    integer                          :: i
    real(8),dimension(size(ebands))  :: dens
    epss = 0.0000001d0
    if(wlx.eq.0.d0)then
       do i=1,size(ebands)
          e=ebands(i)
          if(e.ge.0.d0)then
             dens(i) =  dens_peaks_one(sqrt(e**2.d0 ),wband)
          else
             dens(i) = dens_peaks_one(-sqrt(e**2.d0 ),wband)
          endif
       enddo
    else
       do i=1,size(ebands)
          e=ebands(i)
          if(e.ge.0.d0)then
             dens(i) =  (e/sqrt(e**2.d0 - wlx**2.d0 + epss ))*dens_peaks_one(sqrt(e**2.d0 - wlx**2.d0 ),wband)
          else
             dens(i) = -(e/sqrt(e**2.d0 - wlx**2.d0 + epss ))*dens_peaks_one(-sqrt(e**2.d0 - wlx**2.d0 ),wband)
          endif
       enddo
    endif
  end function dens_peaks_phon

  function dens_peaks_phon_one(e,wlx,wband,xx) result(dens)
    real(8)                          :: e
    real(8)                          :: wlx,wband,epss,xx
    integer                          :: i
    real(8)                          :: dens
    epss = 0.0000001d0
    if(e.ge.0.d0)then
       dens =  (e/sqrt(e**2.d0 - wlx**2.d0 + epss ))*dens_peaks_one(sqrt(e**2.d0 - wlx**2.d0 ),wband)
    else
       dens = -(e/sqrt(e**2.d0 - wlx**2.d0 + epss ))*dens_peaks_one(-sqrt(e**2.d0 - wlx**2.d0 ),wband)
    endif
  end function dens_peaks_phon_one

  function dens_peaks(ebands,wband) result(dens)
    real(8),dimension(:)             :: ebands
    real(8)                          :: wband,e,center
    integer                          :: i
    real(8),dimension(size(ebands))  :: dens
    real(8)                          :: a,b,norm
    center=0.d0
    a=1.9d0
    b=2.1d0
    norm=abs(2.d0*wband*(a**2.d0)/(15.d0*b)+4.d0*wband*a*sqrt((a/(2.d0*b))**2.d0+ &
         (wband**2.d0)*(b*(wband**2.d0)-a)/b)/15.d0+24.d0*(wband**3.d0)*(b*(wband**2.d0)-a)/15.d0)
    do i=1,size(ebands)
       e=ebands(i)
       if(abs(e).lt.abs(wband))then
          dens(i)=(a*(e-center)**2.d0-b*(e-center)**4.d0+(wband**2.d0)*(b*(wband**2.d0)-a))/norm
       else
          dens(i)=0.d0
       endif
    enddo
  end function dens_peaks

  function dens_peaks_one(e,wband) result(dens)
    real(8)                          :: e
    real(8)                          :: wband,center
    integer                          :: i
    real(8)                          :: dens
    real(8)                          :: a,b,norm
    center=0.d0
    a=1.9d0
    b=2.1d0
    norm=abs(2.d0*wband*(a**2.d0)/(15.d0*b)+4.d0*wband*a*sqrt((a/(2.d0*b))**2.d0+ &
         (wband**2.d0)*(b*(wband**2.d0)-a)/b)/15.d0+24.d0*(wband**3.d0)*(b*(wband**2.d0)-a)/15.d0)
    if(abs(e).lt.abs(wband))then
       dens = (a*(e-center)**2.d0-b*(e-center)**4.d0+(wband**2.d0)*(b*(wband**2.d0)-a))/norm
    else
       dens = 0.d0
    endif
  end function dens_peaks_one

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




end program hm_2bands_dos



