  do i=MpiIstart,MpiIend
     call state2indices(i,[DimUps,DimDws],Indices)
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
     !Diagonal Elements, i.e. local part
     htmp = zero
     do iorb=1,Norb
        htmp = htmp + impHloc(1,1,iorb,iorb)*Nup(iorb)
        htmp = htmp + impHloc(Nspin,Nspin,iorb,iorb)*Ndw(iorb)
        htmp = htmp - xmu*(Nup(iorb)+Ndw(iorb))
     enddo
     !
     select case(MpiStatus)
     case (.true.)
        call sp_insert_element(MpiComm,spH0d,htmp,i,i)
     case (.false.)
        call sp_insert_element(spH0d,htmp,i,i)
     end select
     !
  enddo


