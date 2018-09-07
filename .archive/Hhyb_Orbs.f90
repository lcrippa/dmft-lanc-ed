  do iorb=1,Ns_Ud               !1:Norb
     !
     ! IMP IORB_UP <--> BATH IORB_UP
     do iup=1,DimUps(iorb)
        mup         = Hs(iorb)%map(iup)
        Nups(iorb,:) = bdecomp(mup,Ns_Orb)
        !
        do kp=1,Nbath
           alfa=1+kp
           if( (diag_hybr(1,iorb,kp)/=0d0) .AND. (Nups(iorb,1)==1) .AND. (Nups(iorb,alfa)==0) )then              
              call c(1,mup,k1,sg1)
              call cdg(alfa,k1,k2,sg2)
              jup = binary_search(Hs(iorb)%map,k2)
              htmp = diag_hybr(1,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0ups(iorb),htmp,iup,jup)
              !
           endif
           !
           if( (diag_hybr(1,iorb,kp)/=0d0) .AND. (Nups(iorb,1)==0) .AND. (Nups(iorb,alfa)==1) )then
              call c(alfa,mup,k1,sg1)
              call cdg(1,k1,k2,sg2)
              jup=binary_search(Hs(iorb)%map,k2)
              htmp = diag_hybr(1,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0ups(iorb),htmp,iup,jup)
              !
           endif
        enddo
        !
     enddo

     !IMP IORB_DW <--> BATH IORB_DW
     do idw=1,DimDws(iorb)
        mdw          = Hs(iorb+Ns_Ud)%map(idw)
        Ndws(iorb,:) = bdecomp(mdw,Ns_Orb)
        !
        do kp=1,Nbath
           alfa=1+kp
           if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. (Ndws(iorb,1)==1) .AND. (Ndws(iorb,alfa)==0) )then
              call c(1,mdw,k1,sg1)
              call cdg(alfa,k1,k2,sg2)
              jdw=binary_search(Hs(iorb+Ns_Ud)%map,k2)
              htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0dws(iorb),htmp,idw,jdw)
              !
           endif
           if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. (Ndws(iorb,1)==0) .AND. (Ndws(iorb,alfa)==1) )then
              call c(alfa,mdw,k1,sg1)
              call cdg(1,k1,k2,sg2)
              jdw=binary_search(Hs(iorb+Ns_Ud)%map,k2)
              htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0dws(iorb),htmp,idw,jdw)
              !
           endif
        enddo
     enddo
     !
     !
  enddo
