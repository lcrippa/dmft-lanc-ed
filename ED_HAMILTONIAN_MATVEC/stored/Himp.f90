  ! iup = iup_index(i,DimUp)
  ! idw = idw_index(i,DimUp)
  ! !
  ! mup = Hs(1)%map(iup)
  ! mdw = Hs(2)%map(idw)
  ! !
  ! Nups(1,:) = Bdecomp(mup,Ns)
  ! Ndws(1,:) = Bdecomp(mdw,Ns)
  ! !
  ! Nup = Breorder(Nups)
  ! Ndw = Breorder(Ndws)  
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
     !Diagonal Elements, i.e. local part
     htmp = zero
     do iorb=1,Norb
        htmp = htmp + impHloc(1,1,iorb,iorb)*Nup(iorb)
        htmp = htmp + impHloc(Nspin,Nspin,iorb,iorb)*Ndw(iorb)
        htmp = htmp - xmu*( Nup(iorb)+Ndw(iorb) )
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



  !Off-diagonal elements, i.e. non-local part
  !this loop considers only the orbital off-diagonal terms
  !because iorb=jorb can not have simultaneously
  !occupation 0 and 1, as required by this if Jcondition:
  if(ed_total_ud)then
     !
     !UP
     do iup=1,DimUp
        mup        = Hs(1)%map(iup)
        Nups(1,:)  = Bdecomp(mup,Ns)
        Nup        = Breorder(Nups)
        !
        do iorb=1,Norb
           do jorb=1,Norb
              Jcondition = &
                   (impHloc(1,1,iorb,jorb)/=zero) .AND. &
                   (Nup(jorb)==1) .AND. (Nup(iorb)==0)
              if (Jcondition) then
                 call c(jorb,mup,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 jup = binary_search(Hs(1)%map,k2)
                 htmp = impHloc(1,1,iorb,jorb)*sg1*sg2
                 !
                 call sp_insert_element(spH0ups(1),htmp,iup,jup)
                 !
              endif
           enddo
        enddo
        !
     end do
     !
     !DW
     do idw=1,DimDw
        mdw        = Hs(2)%map(idw)
        Ndws(1,:)  = bdecomp(mdw,Ns)
        Ndw        = Breorder(Ndws)
        !
        do iorb=1,Norb
           do jorb=1,Norb
              Jcondition = &
                   (impHloc(Nspin,Nspin,iorb,jorb)/=zero) .AND. &
                   (Ndw(jorb)==1) .AND. (Ndw(iorb)==0)
              if (Jcondition) then
                 call c(jorb,mdw,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 jdw = binary_search(Hs(2)%map,k2)
                 htmp = impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2
                 !
                 call sp_insert_element(spH0dws(1),htmp,idw,jdw)
                 !
              endif
           enddo
        enddo
        !
     enddo
  endif
