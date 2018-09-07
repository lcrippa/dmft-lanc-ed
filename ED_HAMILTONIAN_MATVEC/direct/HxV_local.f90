  do i=1,Dim
     call state2indices(i,[DimUps,DimDws],Indices)
     !
     do iud=1,Ns_Ud
        mup = Hs(iud)%map(Indices(iud))
        mdw = Hs(iud+Ns_Ud)%map(Indices(iud+Ns_ud))
        !
        Nups(iud,:) = Bdecomp(mup,Ns_Orb) ![1+Nbath]*Norb
        Ndws(iud,:) = Bdecomp(mdw,Ns_Orb) ![1+Nbath]*Norb
     enddo
     !
     Nup = Breorder(Nups)
     Ndw = Breorder(Ndws)
     !
     !
     !
     ! HxV_imp: Diagonal Elements, i.e. local part
     htmp = zero
     do iorb=1,Norb
        htmp = htmp + impHloc(1,1,iorb,iorb)*Nup(iorb)
        htmp = htmp + impHloc(Nspin,Nspin,iorb,iorb)*Ndw(iorb)
        htmp = htmp - xmu*(Nup(iorb)+Ndw(iorb))
     enddo
     !
     hv(i) = hv(i) + htmp*vin(i)
     !
     !
     !
     ! HxV_int:
     ! density-density interaction: same orbital, opposite spins:
     !  = \sum_\a U_\a*(n_{\a,up}*n_{\a,dw})
     htmp = zero
     !
     do iorb=1,Norb
        htmp = htmp + Uloc(iorb)*Nup(iorb)*Ndw(iorb)
     enddo
     if(Norb>1)then
        !density-density interaction: different orbitals, opposite spins:
        ! =   U'   *     sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
        ! =  (Uloc-2*Jh)*sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              htmp = htmp + Ust*(Nup(iorb)*Ndw(jorb) + Nup(jorb)*Ndw(iorb))
           enddo
        enddo
        !density-density interaction: different orbitals, parallel spins
        ! = \sum_{i<j}    U''     *[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
        ! = \sum_{i<j} (Uloc-3*Jh)*[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              htmp = htmp + (Ust-Jh)*(Nup(iorb)*Nup(jorb) + Ndw(iorb)*Ndw(jorb))
           enddo
        enddo
     endif
     !if using the Hartree-shifted chemical potential: mu=0 for half-filling
     !sum up the contributions of hartree terms:
     if(hfmode)then
        do iorb=1,Norb
           htmp = htmp - 0.5d0*Uloc(iorb)*(Nup(iorb)+Ndw(iorb)) + 0.25d0*uloc(iorb)
        enddo
        if(Norb>1)then
           do iorb=1,Norb
              do jorb=iorb+1,Norb
                 htmp=htmp-0.5d0*Ust*(Nup(iorb)+Ndw(iorb)+Nup(jorb)+Ndw(jorb))+0.25d0*Ust
                 htmp=htmp-0.5d0*(Ust-Jh)*(Nup(iorb)+Ndw(iorb)+Nup(jorb)+Ndw(jorb))+0.25d0*(Ust-Jh)
              enddo
           enddo
        endif
     endif
     !
     hv(i) = hv(i) + htmp*vin(i)
     !
     !
     !
     ! HxV_bath:
     !
     if(bath_type/="replica") then
        !
        !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
        htmp=zero
        do iorb=1,size(dmft_bath%e,2)
           do kp=1,Nbath
              alfa = getBathStride(iorb,kp)
              htmp =htmp + dmft_bath%e(1    ,iorb,kp)*Nup(alfa) !UP
              htmp =htmp + dmft_bath%e(Nspin,iorb,kp)*Ndw(alfa) !DW
           enddo
        enddo
        !
        hv(i) = hv(i) + htmp*vin(i)
        !
     else
        !
        !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
        htmp=zero
        do kp=1,Nbath
           do iorb=1,Norb
              alfa = getBathStride(iorb,kp)
              htmp = htmp + dmft_bath%h(1    ,    1,iorb,iorb,kp)*Nup(alfa) !UP
              htmp = htmp + dmft_bath%h(Nspin,Nspin,iorb,iorb,kp)*Ndw(alfa) !DW
           enddo
        enddo
        !
        hv(i) = hv(i) + htmp*vin(i)
        !
     endif
     !
  enddo



