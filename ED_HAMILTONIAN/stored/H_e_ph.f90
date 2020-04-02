!We build the electronic part of the electron-phonon interaction: Sum_iorb g_iorb n_iorb
!The phononic part will be dealt on-the-fly since it is very simple.
  do i=MpiIstart,MpiIend
     iup = iup_index(i,DimUp)
     idw = idw_index(i,DimUp)
     !
     mup = Hs(1)%map(iup)
     mdw = Hs(2)%map(idw)
     !
     nup = bdecomp(mup,Ns)
     ndw = bdecomp(mdw,Ns)
     !
     htmp = zero
     do iorb=1,Norb
        htmp = htmp + g_ph(iorb)*(nup(iorb)+ndw(iorb))
     enddo
     !
     select case(MpiStatus)
     case (.true.)
        call sp_insert_element(MpiComm,spH0_eph,htmp,i,i)
     case (.false.)
        call sp_insert_element(spH0_eph,htmp,i,i)
     end select
     !
  enddo

















