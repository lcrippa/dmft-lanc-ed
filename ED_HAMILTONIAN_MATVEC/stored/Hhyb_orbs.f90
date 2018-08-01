  do iorb=1,Norb

     ! IMP IORB_UP <--> BATH IORB_UP
     do iup=1,DimUps(iorb)
        mup  = Hs(iorb)%map(iup)
        Ibup = bdecomp(mup,Ns_Orb)
        !
        do kp=1,Nbath
           alfa=1+kp
           if( (diag_hybr(1,iorb,kp)/=0d0) .AND. (ibup(1)==1) .AND. (ibup(alfa)==0) )then              
              call c(1,mup,k1,sg1)
              call cdg(alfa,k1,k2,sg2)
              jup = binary_search(Hs(iorb)%map,k2)
              htmp = diag_hybr(1,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0ups(iorb),htmp,iup,jup)
              !
           endif
           !
           if( (diag_hybr(1,iorb,kp)/=0d0) .AND. (ibup(1)==0) .AND. (ibup(alfa)==1) )then
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
        mdw  = Hs(iorb+Norb)%map(idw)
        Ibdw = bdecomp(mdw,Ns_Orb)
        !
        do kp=1,Nbath
           alfa=1+kp
           if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. (ibdw(1)==1) .AND. (ibdw(alfa)==0) )then
              call c(1,mdw,k1,sg1)
              call cdg(alfa,k1,k2,sg2)
              jdw=binary_search(Hs(iorb+Norb)%map,k2)
              htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0dws(iorb),htmp,idw,jdw)
              !
           endif
           if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. (ibdw(1)==0) .AND. (ibdw(alfa)==1) )then
              call c(alfa,mdw,k1,sg1)
              call cdg(1,k1,k2,sg2)
              jdw=binary_search(Hs(Norb+iorb)%map,k2)
              htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0dws(iorb),htmp,idw,jdw)
              !
           endif
        enddo
     enddo
     !
  enddo
