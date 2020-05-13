MODULE ED_OBSERVABLES
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  USE SF_ARRAYS, only: arange
  USE SF_LINALG
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_SETUP
  USE ED_HAMILTONIAN
  USE ED_BATH
  USE ED_AUX_FUNX

  implicit none
  private
  !
  public :: observables_impurity
  public :: local_energy_impurity



  logical,save                       :: iolegend=.true.
  real(8),dimension(:),allocatable   :: dens,dens_up,dens_dw
  real(8),dimension(:),allocatable   :: docc
  real(8),dimension(:),allocatable   :: magz
  real(8),dimension(:,:),allocatable :: sz2,n2
  real(8),dimension(:,:),allocatable :: zimp,simp
  real(8)                            :: dens_ph
  real(8)                            :: s2tot
  real(8)                            :: Egs
  real(8)                            :: Ei
  real(8),dimension(:),allocatable   :: Prob
  real(8),dimension(:),allocatable   :: prob_ph
  real(8),dimension(:),allocatable   :: pdf_ph
  real(8)                            :: w_ph
  !
  integer                            :: iorb,jorb,iorb1,jorb1
  integer                            :: ispin,jspin
  integer                            :: isite,jsite
  integer                            :: ibath
  integer                            :: r,m,k,k1,k2,k3,k4
  integer                            :: iup,idw
  integer                            :: jup,jdw
  integer                            :: mup,mdw
  integer                            :: iph,i_el
  real(8)                            :: sgn,sgn1,sgn2,sg1,sg2,sg3,sg4
  real(8)                            :: gs_weight
  !
  real(8)                            :: peso
  real(8)                            :: norm
  !
  integer                            :: i,j,ii
  integer                            :: isector,jsector
  integer                            :: idim,idimUP,idimDW
  !

  real(8),dimension(:),pointer       :: state_cvec
  logical                            :: Jcondition
  !



contains 

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate and print out many interesting physical qties
  !+-------------------------------------------------------------------+
  subroutine observables_impurity()
    select case(ed_diag_type)
    case default
       call lanc_observables()
    case ("full")
       call full_observables()
    end select
  end subroutine observables_impurity



  subroutine local_energy_impurity()
    select case(ed_diag_type)
    case default
       call lanc_local_energy()
    case ("full")
       call full_local_energy()
    end select
  end subroutine local_energy_impurity




  !+-------------------------------------------------------------------+
  !PURPOSE  : Lanc method
  !+-------------------------------------------------------------------+
  subroutine lanc_observables()
    integer                             :: iprob,istate,Nud(2,Ns),iud(2),jud(2)
    integer,dimension(2*Ns_Ud)          :: Indices,Jndices
    integer,dimension(Ns_Ud)            :: iDimUps,iDimDws
    integer,dimension(Ns_Ud,Ns_Orb)     :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)               :: IbUp,IbDw  ![Ns]
    real(8),dimension(Norb)             :: nup,ndw,Sz,nt
    type(sector_map),dimension(2*Ns_Ud) :: HI
    !
    allocate(dens(Norb),dens_up(Norb),dens_dw(Norb))
    allocate(docc(Norb))
    allocate(magz(Norb),sz2(Norb,Norb),n2(Norb,Norb))
    allocate(simp(Norb,Nspin),zimp(Norb,Nspin))
    allocate(Prob(3**Norb))
    allocate(prob_ph(DimPh))
    allocate(pdf_ph(Lpos))
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
    Prob    = 0.d0
    prob_ph = 0.d0
    dens_ph = 0.d0
    pdf_ph  = 0.d0
    w_ph    = w0_ph
    !
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
       !
#ifdef _MPI
       if(MpiStatus)then
          state_cvec => es_return_cvector(MpiComm,state_list,istate)
       else
          state_cvec => es_return_cvector(state_list,istate)
       endif
#else
       state_cvec => es_return_cvector(state_list,istate)
#endif
       !
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       iDim  = getdim(isector)
       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       !
       if(MpiMaster)then
          call build_sector(isector,HI)
          do i = 1,iDim
             iph = (i-1)/(iDimUp*iDimDw) + 1
             i_el = mod(i-1,iDimUp*iDimDw) + 1
             !
             call state2indices(i_el,[iDimUps,iDimDws],Indices)
             do ii=1,Ns_Ud
                mup = HI(ii)%map(Indices(ii))
                mdw = HI(ii+Ns_Ud)%map(Indices(ii+Ns_ud))
                Nups(ii,:) = Bdecomp(mup,Ns_Orb) ![Norb,1+Nbath]
                Ndws(ii,:) = Bdecomp(mdw,Ns_Orb)
             enddo
             IbUp = Breorder(Nups)
             IbDw = Breorder(Ndws)
             !
             gs_weight=peso*abs(state_cvec(i))**2
             !
             !Get operators:
             do iorb=1,Norb
                nup(iorb)= ibup(iorb)
                ndw(iorb)= ibdw(iorb)
                sz(iorb) = (nup(iorb) - ndw(iorb))/2d0
                nt(iorb) =  nup(iorb) + ndw(iorb)
             enddo
             !
             !Configuration probability
             iprob=1
             do iorb=1,Norb
                iprob=iprob+nint(nt(iorb))*3**(iorb-1)
             end do
             Prob(iprob) = Prob(iprob) + gs_weight
             !
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
             prob_ph(iph) = prob_ph(iph) + gs_weight
             dens_ph = dens_ph + (iph-1)*gs_weight
          enddo
          !
          !compute the lattice probability distribution function
          if(Dimph>1)call prob_distr_ph(state_cvec)
          !
          call delete_sector(isector,HI)
       endif
       !
#ifdef _MPI
       if(MpiStatus)then
          if(associated(state_cvec))deallocate(state_cvec)
       else
          if(associated(state_cvec))nullify(state_cvec)
       endif
#else
       if(associated(state_cvec))nullify(state_cvec)
#endif
       !
    enddo
    !
    !
    !
    !IMPURITY DENSITY MATRIX
    if(allocated(imp_density_matrix)) deallocate(imp_density_matrix)
    allocate(imp_density_matrix(Nspin,Nspin,Norb,Norb));imp_density_matrix=zero
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          state_cvec => es_return_cvector(MpiComm,state_list,istate)
       else
          state_cvec => es_return_cvector(state_list,istate)
       endif
#else
       state_cvec => es_return_cvector(state_list,istate)
#endif
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       idim  = getdim(isector)
       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       !
       if(MpiMaster)then
          call build_sector(isector,HI)
          do i=1,iDim
             iph = (i-1)/(iDimUp*iDimDw) + 1
             i_el = mod(i-1,iDimUp*iDimDw) + 1
             !
             call state2indices(i_el,[iDimUps,iDimDws],Indices)
             do ii=1,Ns_Ud
                mup = HI(ii)%map(Indices(ii))
                mdw = HI(ii+Ns_Ud)%map(Indices(ii+Ns_ud))
                Nups(ii,:) = Bdecomp(mup,Ns_Orb) ![Norb,1+Nbath]
                Ndws(ii,:) = Bdecomp(mdw,Ns_Orb)
             enddo
             Nud(1,:) = Breorder(Nups)
             Nud(2,:) = Breorder(Ndws)
             !
             !Diagonal densities
             do ispin=1,Nspin
                do iorb=1,Norb
                   imp_density_matrix(ispin,ispin,iorb,iorb) = &
                        imp_density_matrix(ispin,ispin,iorb,iorb) + &
                        peso*nud(ispin,iorb)*(state_cvec(i))*state_cvec(i)
                enddo
             enddo
             !
             !Off-diagonal
             if(ed_total_ud)then
                do ispin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         !
                         if((Nud(ispin,jorb)==1).and.(Nud(ispin,iorb)==0))then
                            iud(1) = HI(1)%map(Indices(1))
                            iud(2) = HI(2)%map(Indices(2))
                            call c(jorb,iud(ispin),r,sgn1)
                            call cdg(iorb,r,k,sgn2)
                            Jndices = Indices
                            Jndices(1+(ispin-1)*Ns_Ud) = &
                                 binary_search(HI(1+(ispin-1)*Ns_Ud)%map,k)
                            call indices2state(Jndices,[iDimUps,iDimDws],j)
                            !
                            j = j + (iph-1)*iDimUp*iDimDw
                            !
                            imp_density_matrix(ispin,ispin,iorb,jorb) = &
                                 imp_density_matrix(ispin,ispin,iorb,jorb) + &
                                 peso*sgn1*state_cvec(i)*sgn2*(state_cvec(j))
                         endif
                      enddo
                   enddo
                enddo
             endif
             !
             !
          enddo
          call delete_sector(isector,HI)         
       endif
       !
#ifdef _MPI
       if(MpiStatus)then
          if(associated(state_cvec))deallocate(state_cvec)
       else
          if(associated(state_cvec))nullify(state_cvec)
       endif
#else
       if(associated(state_cvec))nullify(state_cvec)
#endif
       !
    enddo
    !
    !
    !
    !
    if(MPIMASTER)then
       call get_szr
       if(DimPh>1) w_ph = sqrt(-2.d0*w0_ph/impDmats_ph(0)) !renormalized phonon frequency
       if(iolegend)call write_legend
       call write_observables()
       write(LOGfile,"(A,10f18.12,f18.12,A)")&
            "dens"//reg(ed_file_suffix)//"=",(dens(iorb),iorb=1,Norb),sum(dens)
       write(LOGfile,"(A,10f18.12,A)")&
            "docc"//reg(ed_file_suffix)//"=",(docc(iorb),iorb=1,Norb)
       if(Nspin==2)then
          write(LOGfile,"(A,10f18.12,A)")&
               "mag "//reg(ed_file_suffix)//"=",(magz(iorb),iorb=1,Norb)
       endif
       if(DimPh>1)call write_pdf()
       !
       do iorb=1,Norb
          ed_dens_up(iorb)=dens_up(iorb)
          ed_dens_dw(iorb)=dens_dw(iorb)
          ed_dens(iorb)   =dens(iorb)
          ed_docc(iorb)   =docc(iorb)
       enddo
    endif
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,ed_dens_up)
       call Bcast_MPI(MpiComm,ed_dens_dw)
       call Bcast_MPI(MpiComm,ed_dens)
       call Bcast_MPI(MpiComm,ed_docc)
    endif
#endif
    !
    deallocate(dens,docc,dens_up,dens_dw,magz,sz2,n2,Prob)
    deallocate(simp,zimp,prob_ph,pdf_ph)
  end subroutine lanc_observables





  !+-------------------------------------------------------------------+
  !PURPOSE  : Get internal energy from the Impurity problem.
  !+-------------------------------------------------------------------+
  subroutine lanc_local_energy()
    integer                             :: istate,iud(2),jud(2)
    integer,dimension(2*Ns_Ud)          :: Indices,Jndices
    integer,dimension(Ns_Ud)            :: iDimUps,iDimDws
    integer,dimension(Ns_Ud,Ns_Orb)     :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    real(8),dimension(Ns)               :: Nup,Ndw
    type(sector_map),dimension(2*Ns_Ud) :: H
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
    !
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          state_cvec => es_return_cvector(MpiComm,state_list,istate)
       else
          state_cvec => es_return_cvector(state_list,istate)
       endif
#else
       state_cvec => es_return_cvector(state_list,istate)
#endif
       !
       iDim  = getdim(isector)
       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       !Master:
       if(MpiMaster)then
          call build_sector(isector,H)
          do i=1,iDim
             iph = (i-1)/(iDimUp*iDimDw) + 1
             i_el = mod(i-1,iDimUp*iDimDw) + 1
             !
             call state2indices(i_el,[iDimUps,iDimDws],Indices)
             do ii=1,Ns_Ud
                mup = H(ii)%map(Indices(ii))
                mdw = H(ii+Ns_Ud)%map(Indices(ii+Ns_ud))
                Nups(ii,:) = Bdecomp(mup,Ns_Orb) ![Norb,1+Nbath]
                Ndws(ii,:) = Bdecomp(mdw,Ns_Orb)
             enddo
             Nup = Breorder(Nups)
             Ndw = Breorder(Ndws)
             !
             gs_weight=peso*abs(state_cvec(i))**2
             !
             !> H_Imp: Diagonal Elements, i.e. local part
             do iorb=1,Norb
                ed_Eknot = ed_Eknot + impHloc(1,1,iorb,iorb)*Nup(iorb)*gs_weight
                ed_Eknot = ed_Eknot + impHloc(Nspin,Nspin,iorb,iorb)*Ndw(iorb)*gs_weight
             enddo
             ! !> H_imp: Off-diagonal elements, i.e. non-local part. 
             if(ed_total_ud)then
                iup = Indices(1)  ; idw = Indices(2)
                mup = H(1)%map(iup) ; mdw = H(2)%map(idw)
                do iorb=1,Norb
                   do jorb=1,Norb
                      !UP
                      Jcondition = &
                           (impHloc(1,1,iorb,jorb)/=0d0) .AND. &
                           (Nup(jorb)==1) .AND. (Nup(iorb)==0)
                      if (Jcondition) then
                         call c(jorb,mup,k1,sg1)
                         call cdg(iorb,k1,k2,sg2)
                         jup = binary_search(H(1)%map,k2)
                         j   = jup + (idw-1)*iDimUp
                         ed_Eknot = ed_Eknot + &
                              impHloc(1,1,iorb,jorb)*sg1*sg2*state_cvec(i)*(state_cvec(j))*peso
                      endif
                      !
                      !DW
                      Jcondition = &
                           (impHloc(Nspin,Nspin,iorb,jorb)/=0d0) .AND. &
                           (ndw(jorb)==1) .AND. (ndw(iorb)==0)
                      if (Jcondition) then
                         call c(jorb,mdw,k1,sg1)
                         call cdg(iorb,k1,k2,sg2)
                         jdw = binary_search(H(2)%map,k2)
                         j   = iup + (jdw-1)*iDimUp
                         ed_Eknot = ed_Eknot + &
                              impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2*state_cvec(i)*(state_cvec(j))*peso
                      endif
                   enddo
                enddo
             ! 
             !SPIN-EXCHANGE Jx
                if(Jhflag.AND.Jx/=0d0)then
                   do iorb=1,Norb
                      do jorb=1,Norb
                         Jcondition=(&
                               (iorb/=jorb).AND.&
                               (nup(jorb)==1).AND.&
                               (ndw(iorb)==1).AND.&
                               (ndw(jorb)==0).AND.&
                               (nup(iorb)==0))
                         if(Jcondition)then
                            call c(iorb,mdw,k1,sg1)  !DW
                            call cdg(jorb,k1,k2,sg2) !DW
                            jdw=binary_search(H(2)%map,k2)
                            call c(jorb,mup,k3,sg3)  !UP
                            call cdg(iorb,k3,k4,sg4) !UP
                            jup=binary_search(H(1)%map,k4)
                            j = jup + (jdw-1)*iDimUp
                            !
                            ed_Epot = ed_Epot + Jx*sg1*sg2*sg3*sg4*state_cvec(i)*state_cvec(j)*peso
                            ed_Dse = ed_Dse + sg1*sg2*sg3*sg4*state_cvec(i)*state_cvec(j)*peso
                            !
                         endif
                      enddo
                   enddo
                endif
             !
             ! PAIR-HOPPING Jp
                if(Jhflag.AND.Jp/=0d0)then
                   do iorb=1,Norb
                      do jorb=1,Norb
                         Jcondition=(&
                               (nup(jorb)==1).AND.&
                               (ndw(jorb)==1).AND.&
                               (ndw(iorb)==0).AND.&
                               (nup(iorb)==0))
                         if(Jcondition)then
                            call c(jorb,mdw,k1,sg1)       !c_jorb_dw
                            call cdg(iorb,k1,k2,sg2)      !c^+_iorb_dw
                            jdw = binary_search(H(2)%map,k2)
                            call c(jorb,mup,k3,sg3)       !c_jorb_up
                            call cdg(iorb,k3,k4,sg4)      !c^+_iorb_up
                            jup = binary_search(H(1)%map,k4)
                            j = jup + (jdw-1)*iDimUp
                            !
                            ed_Epot = ed_Epot + Jp*sg1*sg2*sg3*sg4*state_cvec(i)*state_cvec(j)*peso
                            ed_Dph = ed_Dph + sg1*sg2*sg3*sg4*state_cvec(i)*state_cvec(j)*peso
                            !
                         endif
                      enddo
                   enddo
                endif
             endif     
             !
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
          call delete_sector(isector,H)         
       endif
       !
#ifdef _MPI
       if(MpiStatus)then
          if(associated(state_cvec))deallocate(state_cvec)
       else
          if(associated(state_cvec))nullify(state_cvec)
       endif
#else
       if(associated(state_cvec))nullify(state_cvec)
#endif
       !
    enddo
    !
    !
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,ed_Epot)
       call Bcast_MPI(MpiComm,ed_Eknot)
       call Bcast_MPI(MpiComm,ed_Ehartree)
       call Bcast_MPI(MpiComm,ed_Dust)
       call Bcast_MPI(MpiComm,ed_Dund)
    endif
#endif
    !
    ed_Epot = ed_Epot + ed_Ehartree
    !
    if(ed_verbose==3)then
       write(LOGfile,"(A,10f18.12)")"<Hint>  =",ed_Epot
       write(LOGfile,"(A,10f18.12)")"<V>     =",ed_Epot-ed_Ehartree
       write(LOGfile,"(A,10f18.12)")"<E0>    =",ed_Eknot
       write(LOGfile,"(A,10f18.12)")"<Ehf>   =",ed_Ehartree    
       write(LOGfile,"(A,10f18.12)")"Dust    =",ed_Dust
       write(LOGfile,"(A,10f18.12)")"Dund    =",ed_Dund
    endif
    !
    if(MPIMASTER)then
       call write_energy_info()
       call write_energy()
    endif
    !
    !
  end subroutine lanc_local_energy





  !####################################################################
  !####################################################################
  !####################################################################
  !####################################################################
  !####################################################################






  subroutine full_observables()
    integer                             :: iprob,i,j
    integer                             :: izero,istate
    integer                             :: isector,jsector
    integer                             :: idim,jdim
    integer                             :: iorb,jorb,ispin,jspin,isite,jsite
    integer                             :: numstates
    integer                             :: r,m,k
    real(8)                             :: sgn,sgn1,sgn2
    real(8)                             :: boltzman_weight
    real(8)                             :: state_weight
    real(8)                             :: weight
    real(8)                             :: Ei
    real(8)                             :: norm,beta_
    integer                             :: Nud(2,Ns),iud(2),jud(2)
    integer,dimension(2*Ns_Ud)          :: Indices,Jndices
    integer,dimension(Ns_Ud)            :: iDimUps,iDimDws
    integer,dimension(Ns_Ud,Ns_Orb)     :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)               :: IbUp,IbDw  ![Ns]
    real(8),dimension(Norb)             :: nup,ndw,Sz,nt
    type(sector_map),dimension(2*Ns_Ud) :: HI
    real(8),dimension(:),pointer        :: evec
    !
    !
    !LOCAL OBSERVABLES:
    allocate(dens(Norb),dens_up(Norb),dens_dw(Norb))
    allocate(docc(Norb))
    allocate(magz(Norb),sz2(Norb,Norb),n2(Norb,Norb))
    allocate(simp(Norb,Nspin),zimp(Norb,Nspin))
    allocate(Prob(3**Norb))
    allocate(prob_ph(DimPh))
    allocate(pdf_ph(Lpos))
    !
    egs     = gs_energy
    dens    = 0.d0
    dens_up = 0.d0
    dens_dw = 0.d0
    docc    = 0.d0
    magz    = 0.d0
    sz2     = 0.d0
    n2      = 0.d0
    s2tot   = 0.d0
    Prob    = 0.d0
    prob_ph = 0.d0
    dens_ph = 0.d0
    pdf_ph  = 0.d0
    w_ph    = w0_ph
    !
    beta_ = beta
    if(.not.finiteT)beta_=1000d0
    !
    do isector=1,Nsectors
       iDim  = getdim(isector)
       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       call build_sector(isector,HI)
       !
       do istate=1,iDim
          Ei=espace(isector)%e(istate)
          boltzman_weight=exp(-beta_*Ei)/zeta_function
          if(boltzman_weight < cutoff)cycle
          print*, boltzman_weight
          !
          evec => espace(isector)%M(:,istate)
          !
          do i=1,iDim
             iph = (i-1)/(iDimUp*iDimDw) + 1
             i_el = mod(i-1,iDimUp*iDimDw) + 1
             !
             call state2indices(i_el,[iDimUps,iDimDws],Indices)
             do ii=1,Ns_Ud
                mup = HI(ii)%map(Indices(ii))
                mdw = HI(ii+Ns_Ud)%map(Indices(ii+Ns_ud))
                Nups(ii,:) = Bdecomp(mup,Ns_Orb) ![Norb,1+Nbath]
                Ndws(ii,:) = Bdecomp(mdw,Ns_Orb)
             enddo
             IbUp = Breorder(Nups)
             IbDw = Breorder(Ndws)
             !
             state_weight=(evec(i))*evec(i)
             weight = boltzman_weight*state_weight
             !
             !Get operators:
             do iorb=1,Norb
                nup(iorb)= dble(IbUp(iorb))
                ndw(iorb)= dble(IbDw(iorb))
                sz(iorb) = (nup(iorb) - ndw(iorb))/2.d0
                nt(iorb) =  nup(iorb) + ndw(iorb)
             enddo
             !
             !Configuration probability
             iprob=1
             do iorb=1,Norb
                iprob=iprob+nint(nt(iorb))*3**(iorb-1)
             end do
             Prob(iprob) = Prob(iprob) + weight
             !
             !Evaluate averages of observables:
             do iorb=1,Norb
                dens(iorb)     = dens(iorb)      +  nt(iorb)*weight
                dens_up(iorb)  = dens_up(iorb)   +  nup(iorb)*weight
                dens_dw(iorb)  = dens_dw(iorb)   +  ndw(iorb)*weight
                docc(iorb)     = docc(iorb)      +  nup(iorb)*ndw(iorb)*weight
                magz(iorb)     = magz(iorb)      +  (nup(iorb)-ndw(iorb))*weight
                sz2(iorb,iorb) = sz2(iorb,iorb)  +  (sz(iorb)*sz(iorb))*weight
                n2(iorb,iorb)  = n2(iorb,iorb)   +  (nt(iorb)*nt(iorb))*weight
                do jorb=iorb+1,Norb
                   sz2(iorb,jorb) = sz2(iorb,jorb)  +  (sz(iorb)*sz(jorb))*weight
                   sz2(jorb,iorb) = sz2(jorb,iorb)  +  (sz(jorb)*sz(iorb))*weight
                   n2(iorb,jorb)  = n2(iorb,jorb)   +  (nt(iorb)*nt(jorb))*weight
                   n2(jorb,iorb)  = n2(jorb,iorb)   +  (nt(jorb)*nt(iorb))*weight
                enddo
             enddo
             s2tot = s2tot  + (sum(sz))**2*weight
             prob_ph(iph) = prob_ph(iph) + weight
             dens_ph = dens_ph + (iph-1)*weight
          enddo
          !
          !compute the lattice probability distribution function
          if(Dimph>1)call prob_distr_ph(evec)
          !
       enddo
       call delete_sector(isector,HI)
       if(associated(evec))nullify(evec)
    enddo
    !
    call get_szr
    if(DimPh>1) w_ph = sqrt(-2.d0*w0_ph/impDmats_ph(0))
    if(iolegend)call write_legend
    call write_observables()
    if(DimPh>1)call write_pdf()
    !
    write(LOGfile,"(A,10f18.12,f18.12,A)")"dens"//reg(ed_file_suffix)//"=",(dens(iorb),iorb=1,Norb),sum(dens)
    write(LOGfile,"(A,10f18.12,A)")"docc"//reg(ed_file_suffix)//"=",(docc(iorb),iorb=1,Norb)
    if(Nspin==2)write(LOGfile,"(A,10f18.12,A)") "mag "//reg(ed_file_suffix)//"=",(magz(iorb),iorb=1,Norb)
    !
    do iorb=1,Norb
       ed_dens_up(iorb)=dens_up(iorb)
       ed_dens_dw(iorb)=dens_dw(iorb)
       ed_dens(iorb)   =dens(iorb)
       ed_docc(iorb)   =docc(iorb)
    enddo
    !
    deallocate(dens,docc,dens_up,dens_dw,magz,sz2,n2,Prob)
    deallocate(simp,zimp,prob_ph,pdf_ph)    
  end subroutine full_observables




  subroutine full_local_energy()
    integer                             :: i,j
    integer                             :: izero,istate
    integer                             :: isector
    integer                             :: idim
    integer                             :: iorb,jorb,ispin
    integer                             :: numstates
    integer                             :: m,k1,k2,k3,k4
    real(8)                             :: sg1,sg2,sg3,sg4
    real(8)                             :: Ei
    real(8)                             :: boltzman_weight
    real(8)                             :: state_weight
    real(8)                             :: weight
    real(8)                             :: norm,beta_
    real(8),dimension(Nspin,Norb)       :: eloc
    real(8),dimension(:),pointer        :: evec
    integer                             :: iud(2),jud(2)
    integer,dimension(2*Ns_Ud)          :: Indices,Jndices
    integer,dimension(Ns_Ud)            :: iDimUps,iDimDws
    integer,dimension(Ns_Ud,Ns_Orb)     :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    real(8),dimension(Ns)               :: Nup,Ndw
    type(sector_map),dimension(2*Ns_Ud) :: H
    logical                             :: Jcondition
    !
    !
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
          eloc(ispin,iorb)=impHloc(ispin,ispin,iorb,iorb)
       enddo
    enddo
    !
    beta_ = beta
    if(.not.finiteT)beta_=1000d0
    !
    do isector=1,Nsectors
       iDim  = getdim(isector)
       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       call build_sector(isector,H)
       !
       do istate=1,idim
          Ei=espace(isector)%e(istate)
          boltzman_weight=exp(-beta_*Ei)/zeta_function
          if(boltzman_weight < cutoff)cycle
          !
          evec => espace(isector)%M(:,istate)
          !
          do i=1,idim
             iph = (i-1)/(iDimUp*iDimDw) + 1
             i_el = mod(i-1,iDimUp*iDimDw) + 1
             !
             call state2indices(i_el,[iDimUps,iDimDws],Indices)
             do ii=1,Ns_Ud
                mup = H(ii)%map(Indices(ii))
                mdw = H(ii+Ns_Ud)%map(Indices(ii+Ns_ud))
                Nups(ii,:) = Bdecomp(mup,Ns_Orb) ![Norb,1+Nbath]
                Ndws(ii,:) = Bdecomp(mdw,Ns_Orb)
             enddo
             Nup = Breorder(Nups)
             Ndw = Breorder(Ndws)
             !
             state_weight = (evec(i))*evec(i)
             weight = boltzman_weight*state_weight
             !
             !
             !start evaluating the Tr(H_loc) to estimate potential energy
             !
             !LOCAL ENERGY
             ed_Eknot = ed_Eknot + dot_product(eloc(1,:),nup(1:Norb))*weight + dot_product(eloc(Nspin,:),ndw(1:Norb))*weight
             !> H_imp: Off-diagonal elements, i.e. non-local part. 
             if(ed_total_ud)then
                iup = Indices(1)  ; idw = Indices(2)
                mup = H(1)%map(iup) ; mdw = H(2)%map(idw)
                do iorb=1,Norb
                   do jorb=1,Norb
                      !UP
                      Jcondition = &
                           (impHloc(1,1,iorb,jorb)/=0d0) .AND. &
                           (Nup(jorb)==1) .AND. (Nup(iorb)==0)
                      if (Jcondition) then
                         call c(jorb,mup,k1,sg1)
                         call cdg(iorb,k1,k2,sg2)
                         jup = binary_search(H(1)%map,k2)
                         j   = jup + (idw-1)*iDimUp
                         ed_Eknot = ed_Eknot + &
                              impHloc(1,1,iorb,jorb)*sg1*sg2*state_cvec(i)*(state_cvec(j))*boltzman_weight
                      endif
                      !
                      !DW
                      Jcondition = &
                           (impHloc(Nspin,Nspin,iorb,jorb)/=0d0) .AND. &
                           (ndw(jorb)==1) .AND. (ndw(iorb)==0)
                      if (Jcondition) then
                         call c(jorb,mdw,k1,sg1)
                         call cdg(iorb,k1,k2,sg2)
                         jdw = binary_search(H(2)%map,k2)
                         j   = iup + (jdw-1)*iDimUp
                         ed_Eknot = ed_Eknot + &
                              impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2*state_cvec(i)*(state_cvec(j))*boltzman_weight
                      endif
                   enddo
                enddo
             endif
             !
             !
             !DENSITY-DENSITY INTERACTION: SAME ORBITAL, OPPOSITE SPINS
             !Euloc=\sum=i U_i*(n_u*n_d)_i
             !ed_Epot = ed_Epot + dot_product(uloc,nup*ndw)*weight
             do iorb=1,Norb
                ed_Epot = ed_Epot + Uloc(iorb)*nup(iorb)*ndw(iorb)*weight
             enddo
             !
             !DENSITY-DENSITY INTERACTION: DIFFERENT ORBITALS, OPPOSITE SPINS
             !Eust=\sum_ij Ust*(n_up_i*n_dn_j + n_up_j*n_dn_i)
             !    "="\sum_ij (Uloc - 2*Jh)*(n_up_i*n_dn_j + n_up_j*n_dn_i)
             if(Norb>1)then
                do iorb=1,Norb
                   do jorb=iorb+1,Norb
                      ed_Epot = ed_Epot + Ust*(nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))*weight
                      ed_Dust = ed_Dust + (nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))*weight
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
                      ed_Epot = ed_Epot + (Ust-Jh)*(nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))*weight
                      ed_Dund = ed_Dund + (nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))*weight
                   enddo
                enddo
             endif
             !
             !HARTREE-TERMS CONTRIBUTION:
             if(hfmode)then
                do iorb=1,Norb
                   ed_Ehartree=ed_Ehartree - 0.5d0*uloc(iorb)*(nup(iorb)+ndw(iorb))*weight + 0.25d0*uloc(iorb)*weight
                enddo
                if(Norb>1)then
                   do iorb=1,Norb
                      do jorb=iorb+1,Norb
                         ed_Ehartree=ed_Ehartree - 0.5d0*Ust*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))*weight + 0.25d0*Ust*weight
                         ed_Ehartree=ed_Ehartree - 0.5d0*(Ust-Jh)*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))*weight + 0.25d0*(Ust-Jh)*weight
                      enddo
                   enddo
                endif
             endif
          enddo
       enddo
       call delete_sector(isector,H)
       if(associated(evec))nullify(evec)
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
    endif
    call write_energy_info()
    call write_energy()
    !
    !
  end subroutine full_local_energy







  !####################################################################
  !####################################################################
  !####################################################################
  !####################################################################
  !####################################################################








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
         ((reg(txtfy((6+2*Norb)*Norb+2+Nspin+(ispin-1)*Nspin+iorb))//"sig_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin),&
         reg(txtfy((6+2*Norb)*Norb+3+Nspin+(Nspin-1)*Nspin+Norb))//"nph",reg(txtfy((6+2*Norb)*Norb+4+Nspin+(Nspin-1)*Nspin+Norb))//"w_ph"

    close(unit)
    !
    unit = free_unit()
    open(unit,file="parameters_info.ed")
    write(unit,"(A1,90(A14,1X))")"#","1xmu","2beta",&
         (reg(txtfy(2+iorb))//"U_"//reg(txtfy(iorb)),iorb=1,Norb),&
         reg(txtfy(2+Norb+1))//"U'",reg(txtfy(2+Norb+2))//"Jh"
    close(unit)
    !
    unit = free_unit()
    open(unit,file="Nph_probability_info.ed")
    write(unit,"(A1,90(A10,6X))")"#",&
         (reg(txtfy(i+1))//"Nph="//reg(txtfy(i)),i=0,DimPh-1)
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
         reg(txtfy(6))//"<Dnd>"
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
         ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
         dens_ph,w_ph
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
         ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
         dens_ph,w_ph
    close(unit)         
    !
    unit = free_unit()
    open(unit,file="Occupation_prob"//reg(ed_file_suffix)//".ed")
    write(unit,"(125F15.9)")Uloc(1),Prob,sum(Prob)
    close(unit)         
    !
    unit = free_unit()
    open(unit,file="Nph_probability"//reg(ed_file_suffix)//".ed")
    write(unit,"(90(F15.9,1X))") (prob_ph(i),i=1,DimPh)
    close(unit)
  end subroutine write_observables

  subroutine write_energy()
    integer :: unit
    unit = free_unit()
    open(unit,file="energy_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90F15.9)")ed_Epot,ed_Epot-ed_Ehartree,ed_Eknot,ed_Ehartree,ed_Dust,ed_Dund,ed_Dse,ed_Dph
    close(unit)
  end subroutine write_energy

  subroutine write_pdf()
    integer :: unit,i
    real(8) :: x,dx
    unit = free_unit()
    open(unit,file="lattice_prob"//reg(ed_file_suffix)//".ed")
    dx = (xmax-xmin)/dble(Lpos)
    x = xmin
    do i=1,Lpos
       write(unit,"(2F15.9)") x,pdf_ph(i)
       x = x + dx
    enddo
    close(unit)
  end subroutine write_pdf



  !+-------------------------------------------------------------------+
  !PURPOSE  : subroutines useful for the phonons
  !+-------------------------------------------------------------------+

  !Compute the local lattice probability distribution function (PDF), i.e. the local probability of displacement
  !as a function of the displacement itself
  subroutine prob_distr_ph(vec)
    implicit none
    real(8),dimension(:),pointer          :: vec
    real(8)                               :: psi(0:DimPh-1)
    real(8)                               :: x,dx
    integer                               :: i,j,i_ph,j_ph
    integer                               :: istart,jstart,iend,jend
    !
    dx = (xmax-xmin)/dble(Lpos)
    !
    do i_ph=1,DimPh
       istart = 1 + (i_ph-1)*iDimUp*iDimDw
       iend = i_ph*iDimUp*iDimDw
       !
       do j_ph=1,DimPh
          jstart = 1 + (j_ph-1)*iDimUp*iDimDw
          jend = j_ph*iDimUp*iDimDw
          !
          x = xmin
          do i=1,Lpos
             call Hermite(x,psi)
             pdf_ph(i) = pdf_ph(i) + peso*psi(i_ph-1)*psi(j_ph-1)*dot_product(vec(istart:iend),vec(jstart:jend))
             !
             x = x + dx
          enddo
       enddo
    enddo
  end subroutine prob_distr_ph

  !Compute the Hermite functions (i.e. harmonic oscillator eigenfunctions)
  !the output is a vector with the functions up to order Dimph-1 evaluated at position x
  subroutine Hermite(x,psi)
    implicit none	
    real(8),intent(in)  ::  x
    real(8),intent(out) ::  psi(0:DimPh-1)
    integer             ::  i	
    real(8)             ::  den
    !
    den=1.331335373062335d0!pigr**(0.25d0)
    !
    psi(0)=exp(-0.5d0*x*x)/den
    psi(1)=exp(-0.5d0*x*x)*sqrt(2d0)*x/den
    !
    do i=2,DimPh-1
       psi(i)=2*x*psi(i-1)/sqrt(dble(2*i))-psi(i-2)*sqrt(dble(i-1)/dble(i))
    enddo
  end subroutine Hermite


end MODULE ED_OBSERVABLES

















