MODULE ED_OBSERVABLES
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  USE SF_ARRAYS, only: arange
  USE SF_LINALG
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_SETUP
  USE ED_HAMILTONIAN_MATVEC
  USE ED_BATH
  USE ED_AUX_FUNX
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif

  implicit none
  private
  !
  public :: observables_impurity
  public :: local_energy_impurity
  public :: ed_observables_set_MPI
  public :: ed_observables_del_MPI


  logical,save                       :: iolegend=.true.
  real(8),dimension(:),allocatable   :: dens,dens_up,dens_dw
  real(8),dimension(:),allocatable   :: docc
  real(8),dimension(:),allocatable   :: magz
  real(8),dimension(:,:),allocatable :: sz2,n2
  real(8),dimensioN(:,:),allocatable :: zimp,simp
  real(8)                            :: s2tot
  real(8)                            :: Egs
  !

#ifdef _MPI
  integer                            :: MpiComm=MPI_UNDEFINED
#endif
  logical                            :: MpiStatus=.false.
  logical                            :: MPI_MASTER=.true.  !


contains 

  subroutine ed_observables_set_MPI(comm)
#ifdef _MPI
    integer :: comm
    MpiComm  = comm
    MpiStatus = .true.
    MPI_MASTER= get_Master_MPI(MpiComm)
#else
    integer,optional :: comm
#endif
  end subroutine ed_observables_set_MPI


  subroutine ed_observables_del_MPI()
#ifdef _MPI
    MpiComm  = MPI_UNDEFINED
    MpiStatus = .false.
#endif
  end subroutine ed_observables_del_MPI



  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate and print out many interesting physical qties
  !+-------------------------------------------------------------------+
  subroutine observables_impurity()
    integer,dimension(Ns)           :: ibup,ibdw
    integer,dimension(Nlevels)      :: ib
    integer                         :: i,j
    integer                         :: istate,nud(2,Ns),iud(2),jud(2)
    integer                         :: isector,jsector
    integer                         :: idim,idimUP,idimDW
    integer                         :: jdim,jdimUP,jdimDW
    integer                         :: isz,jsz
    integer                         :: iorb,jorb,ispin,jspin,isite,jsite,ibath
    integer                         :: numstates
    integer                         :: r,m,k
    integer                         :: iup,idw,jup,jdw,mup,mdw
    real(8)                         :: sgn,sgn1,sgn2
    real(8)                         :: gs_weight
    real(8)                         :: Ei
    real(8)                         :: peso
    real(8)                         :: norm
    real(8),dimension(Norb)         :: nup,ndw,Sz,nt
    complex(8),dimension(:),pointer :: gscvec
    type(sector_map)                :: HI(2)
    complex(8),allocatable          :: vvinit(:)

    !
    !LOCAL OBSERVABLES:
    ! density, 
    ! double occupancy, 
    ! magnetization, 
    ! orbital//spin correlations  
    ! superconducting order parameter, etc..
    allocate(dens(Norb),dens_up(Norb),dens_dw(Norb))
    allocate(docc(Norb))
    allocate(magz(Norb),sz2(Norb,Norb),n2(Norb,Norb))
    allocate(simp(Norb,Nspin),zimp(Norb,Nspin))
    !
    Egs     = state_list%emin
    dens    = 0.d0
    dens_up = 0.d0
    dens_dw = 0.d0
    docc    = 0.d0
    magz    = 0.d0
    sz2     = 0.d0
    n2      = 0.d0
    s2tot   = 0.d0
    !
    numstates=state_list%size
    do istate=1,numstates
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
       !
       gscvec  => es_return_cvector(state_list,istate)       
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       idim  = getdim(isector)
       idimUp  = getDimUp(isector)
       idimDw  = getDimDw(isector)
       call build_sector(isector,HI)
       !
       !pdens=0d0
       do iup=1,idimUP
          mup  = HI(1)%map(iup)
          ibup = bdecomp(mup,Ns)
          !
          do idw=1,idimDW
             mdw  = HI(2)%map(idw)
             ibdw = bdecomp(mdw,Ns)
             !
             i = iup + (idw-1)*idimUP
             !
             gs_weight=peso*abs(gscvec(i))**2
             !
             !Get operators:
             do iorb=1,Norb
                nup(iorb)= ibup(iorb)
                ndw(iorb)= ibdw(iorb)
                sz(iorb) = (nup(iorb) - ndw(iorb))/2d0
                nt(iorb) =  nup(iorb) + ndw(iorb)
             enddo
             !
             !pdens     = pdens      +  nt(1)*gs_weight*zeta_function
             !Evaluate averages of observables:
             do iorb=1,Norb
                dens(iorb)     = dens(iorb)      +  nt(iorb)*gs_weight
                dens_up(iorb)  = dens_up(iorb)   +  nup(iorb)*gs_weight
                dens_dw(iorb)  = dens_dw(iorb)   +  ndw(iorb)*gs_weight
                docc(iorb)     = docc(iorb)      +  nup(iorb)*ndw(iorb)*gs_weight
                magz(iorb)     = magz(iorb)      +  (nup(iorb)-ndw(iorb))*gs_weight
                sz2(iorb,iorb) = sz2(iorb,iorb)  +  (sz(iorb)*sz(iorb))*gs_weight
                n2(iorb,iorb)  = n2(iorb,iorb)   +  (nt(iorb)*nt(iorb))*gs_weight
                do jorb=iorb+1,Norb
                   sz2(iorb,jorb) = sz2(iorb,jorb)  +  (sz(iorb)*sz(jorb))*gs_weight
                   sz2(jorb,iorb) = sz2(jorb,iorb)  +  (sz(jorb)*sz(iorb))*gs_weight
                   n2(iorb,jorb)  = n2(iorb,jorb)   +  (nt(iorb)*nt(jorb))*gs_weight
                   n2(jorb,iorb)  = n2(jorb,iorb)   +  (nt(jorb)*nt(iorb))*gs_weight
                enddo
             enddo
             s2tot = s2tot  + (sum(sz))**2*gs_weight
          enddo
       enddo
       if(associated(gscvec))nullify(gscvec)
       call delete_sector(isector,HI)
    enddo
    !
    !IMPURITY DENSITY MATRIX
    if(allocated(imp_density_matrix)) deallocate(imp_density_matrix)
    allocate(imp_density_matrix(Nspin,Nspin,Norb,Norb));imp_density_matrix=zero
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
       gscvec  => es_return_cvector(state_list,istate)
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       call build_sector(isector,HI)
       do ispin=1,Nspin
          !
          !Diagonal densities
          do iorb=1,Norb
             do iup=1,idimUP
                iud(1)   = HI(1)%map(iup)
                nud(1,:) = bdecomp(iud(1),Ns)
                do idw=1,idimDW
                   iud(2)   = HI(2)%map(idw)
                   nud(2,:) = bdecomp(iud(2),Ns)
                   !
                   i = iup + (idw-1)*idimUP
                   !
                   imp_density_matrix(ispin,ispin,iorb,iorb) = imp_density_matrix(ispin,ispin,iorb,iorb) + &
                        peso*nud(ispin,iorb)*conjg(gscvec(i))*gscvec(i)
                enddo
             enddo
          enddo
          !off-diagonal
          do iorb=1,Norb
             do jorb=1,Norb
                do iup=1,idimUP
                   iud(1)   = HI(1)%map(iup)
                   nud(1,:) = bdecomp(iud(1),Ns)
                   do idw=1,idimDW
                      iud(2)   = HI(2)%map(idw)
                      nud(2,:) = bdecomp(iud(2),Ns)
                      !
                      if((nud(ispin,jorb)==1).and.(nud(ispin,iorb)==0))then
                         call c(jorb,iud(ispin),r,sgn1)
                         call cdg(iorb,r,k,sgn2)
                         !
                         jud = [iup,idw]
                         jud(ispin) = binary_search(HI(ispin)%map,r)
                         !
                         i = iud(1) + (iud(2)-1)*idimUP
                         j = jud(1) + (jud(2)-1)*jdimUP
                         !
                         imp_density_matrix(ispin,ispin,iorb,jorb) = imp_density_matrix(ispin,ispin,iorb,jorb) + &
                              peso*sgn1*gscvec(i)*sgn2*conjg(gscvec(j))
                      endif
                   enddo
                enddo
             enddo
          enddo
          !
          !
       enddo
       !
       nullify(gscvec)
       call delete_sector(isector,HI)
    enddo
    !
    !
    !
    !
    !
    !
    if(MPI_MASTER)then
       call get_szr
       if(iolegend)call write_legend
       call write_observables()
    endif
    write(LOGfile,"(A,10f18.12,f18.12,A)")"dens"//reg(ed_file_suffix)//"=",(dens(iorb),iorb=1,Norb),sum(dens)
    write(LOGfile,"(A,10f18.12,A)")    "docc"//reg(ed_file_suffix)//"=",(docc(iorb),iorb=1,Norb)
    if(Nspin==2)then
       write(LOGfile,"(A,10f18.12,A)") "mag "//reg(ed_file_suffix)//"=",(magz(iorb),iorb=1,Norb)
    endif
    !
    do iorb=1,Norb
       ed_dens_up(iorb)=dens_up(iorb)
       ed_dens_dw(iorb)=dens_dw(iorb)
       ed_dens(iorb)   =dens(iorb)
       ed_docc(iorb)   =docc(iorb)
    enddo
    !
    deallocate(dens,docc,dens_up,dens_dw,magz,sz2,n2)
    deallocate(simp,zimp)    
  end subroutine observables_impurity








  !+-------------------------------------------------------------------+
  !PURPOSE  : Get internal energy from the Impurity problem.
  !+-------------------------------------------------------------------+
  subroutine local_energy_impurity(MpiComm)
    integer,optional                :: MpiComm
    integer                         :: i,j
    integer                         :: istate,nud(2,Ns),iud(2),jud(2)
    integer                         :: isector,jsector
    integer                         :: idim,idimUP,idimDW
    integer                         :: jdim,jdimUP,jdimDW
    integer                         :: isz,jsz
    integer                         :: iorb,jorb,ispin,jspin,isite,jsite,ibath
    integer                         :: r,m,k
    integer                         :: iup,idw,jup,jdw,mup,mdw
    integer                         :: numstates
    integer                         :: k1,k2,k3,k4
    real(8)                         :: sg1,sg2,sg3,sg4
    real(8)                         :: Egs,gs_weight
    real(8)                         :: Ei
    real(8)                         :: peso
    real(8)                         :: norm
    real(8),dimension(Norb)         :: nup,ndw
    real(8),dimension(Nspin,Norb)   :: eloc
    complex(8),dimension(:),pointer :: gscvec
    type(sector_map)                :: H(2)
    logical                         :: Jcondition
    !
    Egs     = state_list%emin
    ed_Ehartree= 0.d0
    ed_Eknot   = 0.d0
    ed_Epot    = 0.d0
    ed_Dust    = 0.d0
    ed_Dund    = 0.d0
    ed_Dse     = 0.d0
    ed_Dph     = 0.d0
    !
    !Get diagonal part of Hloc
    do ispin=1,Nspin
       do iorb=1,Norb
          Eloc(ispin,iorb)=impHloc(ispin,ispin,iorb,iorb)
       enddo
    enddo
    !
    numstates=state_list%size
    do istate=1,numstates
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
       gscvec  => es_return_cvector(state_list,istate)
       !
       idim  = getdim(isector)
       idimUp  = getDimUp(isector)
       idimDw  = getDimDw(isector)
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       call build_sector(isector,H)
       do iup=1,idimUP
          iud(1)   = H(1)%map(iup)
          nud(1,:) = bdecomp(iud(1),Ns)
          !
          do idw=1,idimDW
             iud(2)   = H(2)%map(idw)
             nud(2,:) = bdecomp(iud(2),Ns)
             !
             i = iup + (idw-1)*idimUP
             !
             gs_weight=peso*abs(gscvec(i))**2
             !
             !Get operators:
             do iorb=1,Norb
                nup(iorb)= nud(1,iorb)
                ndw(iorb)= nud(2,iorb)
             enddo
             !
             !start evaluating the Tr(H_loc) to estimate potential energy
             !
             !LOCAL ENERGY
             ed_Eknot = ed_Eknot + dot_product(eloc(1,:),nup)*gs_weight + dot_product(eloc(Nspin,:),ndw)*gs_weight
             !==> HYBRIDIZATION TERMS I: same or different orbitals, same spins.
             do iorb=1,Norb
                do jorb=1,Norb
                   !SPIN UP
                   ispin=1
                   if((nud(ispin,iorb)==0).AND.(nud(ispin,jorb)==1))then
                      call c(jorb,iud(ispin),k1,sg1)
                      call cdg(iorb,k1,k2,sg2)
                      !
                      jud        = [iup,idw]
                      jud(ispin) = binary_search(H(ispin)%map,r)
                      !                      
                      j = jud(1) + (jud(2)-1)*jdimUP
                      ed_Eknot = ed_Eknot + impHloc(ispin,ispin,iorb,jorb)*sg1*sg2*gscvec(i)*conjg(gscvec(j))
                   endif
                   !SPIN DW
                   ispin=2
                   if((nud(ispin,iorb)==0).AND.(nud(ispin,jorb)==1))then
                      call c(jorb,iud(ispin),k1,sg1)
                      call cdg(iorb,k1,k2,sg2)
                      !
                      jud        = [iup,idw]
                      jud(ispin) = binary_search(H(ispin)%map,r)
                      !                      
                      j = jud(1) + (jud(2)-1)*jdimUP
                      ed_Eknot = ed_Eknot + impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2*gscvec(i)*conjg(gscvec(j))
                   endif
                enddo
             enddo
             !
             !DENSITY-DENSITY INTERACTION: SAME ORBITAL, OPPOSITE SPINS
             !Euloc=\sum=i U_i*(n_u*n_d)_i
             !ed_Epot = ed_Epot + dot_product(uloc,nup*ndw)*gs_weight
             do iorb=1,Norb
                ed_Epot = ed_Epot + Uloc(iorb)*nup(iorb)*ndw(iorb)*gs_weight
             enddo
             !
             !DENSITY-DENSITY INTERACTION: DIFFERENT ORBITALS, OPPOSITE SPINS
             !Eust=\sum_ij Ust*(n_up_i*n_dn_j + n_up_j*n_dn_i)
             !    "="\sum_ij (Uloc - 2*Jh)*(n_up_i*n_dn_j + n_up_j*n_dn_i)
             if(Norb>1)then
                do iorb=1,Norb
                   do jorb=iorb+1,Norb
                      ed_Epot = ed_Epot + Ust*(nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))*gs_weight
                      ed_Dust = ed_Dust + (nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))*gs_weight
                   enddo
                enddo
             endif
             !
             !DENSITY-DENSITY INTERACTION: DIFFERENT ORBITALS, PARALLEL SPINS
             !Eund = \sum_ij Und*(n_up_i*n_up_j + n_dn_i*n_dn_j)
             !    "="\sum_ij (Ust-Jh)*(n_up_i*n_up_j + n_dn_i*n_dn_j)
             !    "="\sum_ij (Uloc-3*Jh)*(n_up_i*n_up_j + n_dn_i*n_dn_j)
             if(Norb>1)then
                do iorb=1,Norb
                   do jorb=iorb+1,Norb
                      ed_Epot = ed_Epot + (Ust-Jh)*(nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))*gs_weight
                      ed_Dund = ed_Dund + (nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))*gs_weight
                   enddo
                enddo
             endif
             !
             !SPIN-EXCHANGE (S-E) TERMS
             !S-E: Jh *( c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up )  (i.ne.j) 
             if(Norb>1.AND.Jhflag)then
                do iorb=1,Norb
                   do jorb=1,Norb
                      Jcondition=(&
                           (iorb/=jorb).AND.&
                           (nup(jorb)==1).AND.&
                           (ndw(iorb)==1).AND.&
                           (ndw(jorb)==0).AND.&
                           (nup(iorb)==0))
                      if(Jcondition)then
                         call c(iorb,iud(2),k1,sg1) !DW
                         call cdg(jorb,k1,k2,sg2) !DW
                         jdw=binary_search(H(2)%map,k2)
                         !
                         call c(jorb,iud(1),k1,sg3) !UP
                         call cdg(iorb,k1,k2,sg4)    !UP
                         jup=binary_search(H(1)%map,k2)
                         !
                         j = jup + (jdw-1)*idimUP
                         !
                         ed_Epot = ed_Epot + Jx*sg1*sg2*sg3*sg4*gscvec(i)*conjg(gscvec(j))
                         ed_Dse  = ed_Dse  + sg1*sg2*sg3*sg4*gscvec(i)*conjg(gscvec(j))
                      endif
                   enddo
                enddo
             endif
             !
             !PAIR-HOPPING (P-H) TERMS
             !P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
             !P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
             if(Norb>1.AND.Jhflag)then
                do iorb=1,Norb
                   do jorb=1,Norb
                      Jcondition=(&
                           (nup(jorb)==1).AND.&
                           (ndw(jorb)==1).AND.&
                           (ndw(iorb)==0).AND.&
                           (nup(iorb)==0))
                      if(Jcondition)then
                         call c(jorb,iud(2),k1,sg1)       !c_jorb_dw
                         call cdg(iorb,k1,k2,sg2)       !c^+_iorb_dw
                         jdw = binary_search(H(2)%map,k2)
                         !
                         call c(jorb,iud(1),k1,sg1)       !c_jorb_up
                         call cdg(iorb,k1,k2,sg4)       !c^+_iorb_up
                         jup = binary_search(H(1)%map,k2)
                         !
                         j = jup + (jdw-1)*idimup
                         !
                         ed_Epot = ed_Epot + Jp*sg1*sg2*sg3*sg4*gscvec(i)*conjg(gscvec(j))!gs_weight
                         ed_Dph  = ed_Dph  + sg1*sg2*sg3*sg4*gscvec(i)*conjg(gscvec(j))!gs_weight
                      endif
                   enddo
                enddo
             endif
             !
             !HARTREE-TERMS CONTRIBUTION:
             if(hfmode)then
                !ed_Ehartree=ed_Ehartree - 0.5d0*dot_product(uloc,nup+ndw)*gs_weight + 0.25d0*sum(uloc)*gs_weight
                do iorb=1,Norb
                   ed_Ehartree=ed_Ehartree - 0.5d0*uloc(iorb)*(nup(iorb)+ndw(iorb))*gs_weight + 0.25d0*uloc(iorb)*gs_weight
                enddo
                if(Norb>1)then
                   do iorb=1,Norb
                      do jorb=iorb+1,Norb
                         ed_Ehartree=ed_Ehartree - 0.5d0*Ust*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))*gs_weight + 0.25d0*Ust*gs_weight
                         ed_Ehartree=ed_Ehartree - 0.5d0*(Ust-Jh)*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))*gs_weight + 0.25d0*(Ust-Jh)*gs_weight
                      enddo
                   enddo
                endif
             endif
          enddo
       enddo
       if(associated(gscvec))nullify(gscvec)
       call delete_sector(isector,H)
    enddo
    ed_Epot = ed_Epot + ed_Ehartree
    !
    if(ed_verbose==3)then
       write(LOGfile,"(A,10f18.12)")"<Hint>  =",ed_Epot
       write(LOGfile,"(A,10f18.12)")"<V>     =",ed_Epot-ed_Ehartree
       write(LOGfile,"(A,10f18.12)")"<E0>    =",ed_Eknot
       write(LOGfile,"(A,10f18.12)")"<Ehf>   =",ed_Ehartree    
       write(LOGfile,"(A,10f18.12)")"Dust    =",ed_Dust
       write(LOGfile,"(A,10f18.12)")"Dund    =",ed_Dund
       write(LOGfile,"(A,10f18.12)")"Dse     =",ed_Dse
       write(LOGfile,"(A,10f18.12)")"Dph     =",ed_Dph
    endif
    if(MPI_MASTER)then
       call write_energy_info()
       call write_energy()
    endif
    !
    !
  end subroutine local_energy_impurity



  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : get scattering rate and renormalization constant Z
  !+-------------------------------------------------------------------+
  subroutine get_szr()
    integer                  :: ispin,iorb
    real(8)                  :: wm1,wm2
    wm1 = pi/beta ; wm2=3d0*pi/beta
    do ispin=1,Nspin
       do iorb=1,Norb
          simp(iorb,ispin) = dimag(impSmats(ispin,ispin,iorb,iorb,1)) - &
               wm1*(dimag(impSmats(ispin,ispin,iorb,iorb,2))-dimag(impSmats(ispin,ispin,iorb,iorb,1)))/(wm2-wm1)
          zimp(iorb,ispin)   = 1.d0/( 1.d0 + abs( dimag(impSmats(ispin,ispin,iorb,iorb,1))/wm1 ))
       enddo
    enddo
  end subroutine get_szr



  !+-------------------------------------------------------------------+
  !PURPOSE  : write legend, i.e. info about columns 
  !+-------------------------------------------------------------------+
  subroutine write_legend()
    integer :: unit,iorb,jorb,ispin
    unit = free_unit()
    open(unit,file="observables_info.ed")
    write(unit,"(A1,90(A10,6X))")"#",&
         (reg(txtfy(iorb))//"dens_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(Norb+iorb))//"docc_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(2*Norb+iorb))//"nup_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(3*Norb+iorb))//"ndw_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(4*Norb+iorb))//"mag_"//reg(txtfy(iorb)),iorb=1,Norb),&
         reg(txtfy(5*Norb+1))//"s2",&
         reg(txtfy(5*Norb+2))//"egs",&
         ((reg(txtfy(5*Norb+2+(iorb-1)*Norb+jorb))//"sz2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
         ((reg(txtfy((5+Norb)*Norb+2+(iorb-1)*Norb+jorb))//"n2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
         ((reg(txtfy((5+2*Norb)*Norb+2+(ispin-1)*Nspin+iorb))//"z_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin),&
         ((reg(txtfy((6+2*Norb)*Norb+2+Nspin+(ispin-1)*Nspin+iorb))//"sig_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)

    close(unit)
    !
    unit = free_unit()
    open(unit,file="parameters_info.ed")
    write(unit,"(A1,90(A14,1X))")"#","1xmu","2beta",&
         (reg(txtfy(2+iorb))//"U_"//reg(txtfy(iorb)),iorb=1,Norb),&
         reg(txtfy(2+Norb+1))//"U'",reg(txtfy(2+Norb+2))//"Jh"
    close(unit)
    !
    iolegend=.false.
  end subroutine write_legend

  subroutine write_energy_info()
    integer :: unit
    unit = free_unit()
    open(unit,file="energy_info.ed")
    write(unit,"(A1,90(A14,1X))")"#",&
         reg(txtfy(1))//"<Hi>",&
         reg(txtfy(2))//"<V>=<Hi-Ehf>",&
         reg(txtfy(3))//"<Eloc>",&
         reg(txtfy(4))//"<Ehf>",&
         reg(txtfy(5))//"<Dst>",&
         reg(txtfy(6))//"<Dnd>",&
         reg(txtfy(7))//"<Dse>",&
         reg(txtfy(8))//"<Dph>"
    close(unit)
  end subroutine write_energy_info


  !+-------------------------------------------------------------------+
  !PURPOSE  : write observables to file
  !+-------------------------------------------------------------------+
  subroutine write_observables()
    integer :: unit
    integer :: iorb,jorb,ispin
    unit = free_unit()
    open(unit,file="observables_all"//reg(ed_file_suffix)//".ed",position='append')
    write(unit,"(90(F15.9,1X))")&
         (dens(iorb),iorb=1,Norb),&
         (docc(iorb),iorb=1,Norb),&
         (dens_up(iorb),iorb=1,Norb),&
         (dens_dw(iorb),iorb=1,Norb),&
         (magz(iorb),iorb=1,Norb),&
         s2tot,egs,&
         ((sz2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
         ((n2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
         ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
         ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
    close(unit)    
    !
    unit = free_unit()
    open(unit,file="parameters_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90F15.9)")xmu,beta,(uloc(iorb),iorb=1,Norb),Ust,Jh,Jx,Jp
    close(unit)
    !
    unit = free_unit()
    open(unit,file="observables_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90(F15.9,1X))")&
         (dens(iorb),iorb=1,Norb),&
         (docc(iorb),iorb=1,Norb),&
         (dens_up(iorb),iorb=1,Norb),&
         (dens_dw(iorb),iorb=1,Norb),&
         (magz(iorb),iorb=1,Norb),&
         s2tot,egs,&
         ((sz2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
         ((n2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
         ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
         ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
    close(unit)         
  end subroutine write_observables

  subroutine write_energy()
    integer :: unit
    unit = free_unit()
    open(unit,file="energy_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90F15.9)")ed_Epot,ed_Epot-ed_Ehartree,ed_Eknot,ed_Ehartree,ed_Dust,ed_Dund,ed_Dse,ed_Dph
    close(unit)
  end subroutine write_energy



end MODULE ED_OBSERVABLES
