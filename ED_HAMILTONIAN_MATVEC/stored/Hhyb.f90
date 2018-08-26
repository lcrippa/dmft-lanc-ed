  do iup=1,DimUp
     mup  = Hs(1)%map(iup)
     nup  = bdecomp(mup,Ns)
     !
     do iorb=1,Norb
        do kp=1,Nbath
           alfa=getBathStride(iorb,kp)
           !
           if( (diag_hybr(1,iorb,kp)/=0d0) &
                .AND. (nup(iorb)==1) .AND. (nup(alfa)==0) )then              
              call c(iorb,mup,k1,sg1)
              call cdg(alfa,k1,k2,sg2)
              jup = binary_search(Hs(1)%map,k2)
              htmp = diag_hybr(1,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0up,htmp,iup,jup)
              !
           endif
           !
           if( (diag_hybr(1,iorb,kp)/=0d0) &
                .AND. (nup(iorb)==0) .AND. (nup(alfa)==1) )then
              call c(alfa,mup,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              jup=binary_search(Hs(1)%map,k2)
              htmp = diag_hybr(1,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0up,htmp,iup,jup)
              !
           endif
        enddo
     enddo
     !
  enddo

  !IMP DW <--> BATH DW
  do idw=1,DimDw
     mdw  = Hs(2)%map(idw)
     ndw  = bdecomp(mdw,Ns)
     !
     do iorb=1,Norb
        do kp=1,Nbath
           alfa=getBathStride(iorb,kp)
           !
           if( (diag_hybr(Nspin,iorb,kp)/=0d0) &
                .AND. (ndw(iorb)==1) .AND. (ndw(alfa)==0) )then
              call c(iorb,mdw,k1,sg1)
              call cdg(alfa,k1,k2,sg2)
              jdw=binary_search(Hs(2)%map,k2)
              htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0dw,htmp,idw,jdw)
              !
           endif
           if( (diag_hybr(Nspin,iorb,kp)/=0d0) &
                .AND. (ndw(iorb)==0) .AND. (ndw(alfa)==1) )then
              call c(alfa,mdw,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              jdw=binary_search(Hs(2)%map,k2)
              htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0dw,htmp,idw,jdw)
              !
           endif
        enddo
     enddo
     !
  enddo
