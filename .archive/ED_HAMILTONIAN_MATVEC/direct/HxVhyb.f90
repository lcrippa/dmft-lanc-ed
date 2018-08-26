  ! IMP UP <--> BATH UP
  ! do idw=map_first_state_dw(1),map_last_state_dw(1)
  !    do iup=map_first_state_up(idw),map_last_state_up(idw)
  do iup=first_state_up,last_state_up
     mup  = Hs(1)%map(iup)
     nup  = bdecomp(mup,Ns)
     do idw=1,DimDw
        !
        !MPI Shifts
        i    = iup + (idw-1)*dimUp
        ! impi = i - ishift
        !
        do iorb=1,Norb
           do kp=1,Nbath
              alfa=getBathStride(iorb,kp)
              !
              if( (diag_hybr(1,iorb,kp)/=0d0) .AND. (nup(iorb)==1) .AND. (nup(alfa)==0) )then              
                 call c(iorb,mup,k1,sg1)
                 call cdg(alfa,k1,k2,sg2)
                 jup = binary_search(Hs(1)%map,k2)
                 j   = jup + (idw-1)*DimUp
                 !
                 htmp = diag_hybr(1,iorb,kp)*sg1*sg2
                 !
                 ! hv(impi) = hv(impi) + htmp*vin(j)
                 vout(i) = vout(i) + htmp*vin(j)
                 !
              endif
              !
              if( (diag_hybr(1,iorb,kp)/=0d0) .AND. (nup(iorb)==0) .AND. (nup(alfa)==1) )then
                 call c(alfa,mup,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 jup = binary_search(Hs(1)%map,k2)
                 j   = jup + (idw-1)*DimUp
                 !
                 htmp = diag_hybr(1,iorb,kp)*sg1*sg2
                 !
                 ! hv(impi) = hv(impi) + htmp*vin(j)
                 vout(i) = vout(i) + htmp*vin(j)
                 !
              endif
           enddo
        enddo
        !
     end do
  end do



  ! IMP DW <--> BATH DW
  ! do idw=map_first_state_dw(1),map_last_state_dw(1)
  !    do iup=map_first_state_up(idw),map_last_state_up(idw)
  do idw=first_state_dw,last_state_dw
     mdw  = Hs(2)%map(idw)
     ndw  = bdecomp(mdw,Ns)
     do iup=1,DimUp
        !
        !MPI Shifts
        i    = iup + (idw-1)*dimUp
        impi = i - ishift
        !
        do iorb=1,Norb
           do kp=1,Nbath
              alfa=getBathStride(iorb,kp)
              !
              if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. (ndw(iorb)==1) .AND. (ndw(alfa)==0) )then
                 call c(iorb,mdw,k1,sg1)
                 call cdg(alfa,k1,k2,sg2)
                 jdw = binary_search(Hs(2)%map,k2)
                 j   = iup + (jdw-1)*DimUp
                 !
                 htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
                 !
                 ! hv(impi) = hv(impi) + htmp*vin(j)
                 vout(i) = vout(i) + htmp*vin(j)
                 !
              endif
              if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. (ndw(iorb)==0) .AND. (ndw(alfa)==1) )then
                 call c(alfa,mdw,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 jdw = binary_search(Hs(2)%map,k2)
                 j   = iup + (jdw-1)*DimUp
                 !
                 htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
                 !
                 ! hv(impi) = hv(impi) + htmp*vin(j)
                 vout(i) = vout(i) + htmp*vin(j)
                 !
              endif
           enddo
        enddo
        !
     end do
  end do
