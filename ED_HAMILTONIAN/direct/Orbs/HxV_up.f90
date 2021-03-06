  do j=1,Nloc
     i_el = mod(j-1,DimUp*DimDw) + 1	!electronic index
     iph = (j-1)/(DimUp*DimDw) + 1
     !
     call state2indices(i_el,[DimUps,DimDws],Jndices)
     !
     !>H_hyb: hopping terms for a given spin (imp <--> bath)
     do iorb=1,Ns_Ud
        mup = Hs(iorb)%map(Jndices(iorb))
        Nups(iorb,:) = Bdecomp(mup,Ns_Orb)
        !
        do kp=1,Nbath
           ialfa=1+kp
           if( (diag_hybr(1,iorb,kp)/=0d0) &
                .AND. (Nups(iorb,1)==1) .AND. (Nups(iorb,ialfa)==0) )then
              call c(1,mup,k1,sg1)
              call cdg(ialfa,k1,k2,sg2)
              Indices       = Jndices
              Indices(iorb) = binary_search(Hs(iorb)%map,k2)
              call indices2state(Indices,[DimUps,DimDws],i)
              htmp  = diag_hybr(1,iorb,kp)*sg1*sg2
              !
              i = i + (iph-1)*DimUp*DimDw
              Hv(i) = Hv(i) + htmp*Vin(j)
              !
           endif
           !
           if( (diag_hybr(1,iorb,kp)/=0d0) &
                .AND. (Nups(iorb,1)==0) .AND. (Nups(iorb,ialfa)==1) )then
              call c(ialfa,mup,k1,sg1)
              call cdg(1,k1,k2,sg2)
              Indices       = Jndices
              Indices(iorb) = binary_search(Hs(iorb)%map,k2)
              call indices2state(Indices,[DimUps,DimDws],i)
              htmp  = diag_hybr(1,iorb,kp)*sg1*sg2
              !
              i = i + (iph-1)*DimUp*DimDw
              Hv(i) = Hv(i) + htmp*Vin(j)
              !
           endif
        enddo
        !
     enddo
  enddo






