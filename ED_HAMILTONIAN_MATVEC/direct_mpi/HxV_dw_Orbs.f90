  do iidw=1,MpiQdw
     do iiup=1,DimUp
        i = iidw + (iiup-1)*DimDw
        call state2indices(i,[DimDws,DimUps],Indices)
        do iorb=1,Ns_Ud
           mdw          = Hs(iorb+Ns_Ud)%map(Indices(iorb+Ns_Ud))
           Ndws(iorb,:) = bdecomp(mdw,Ns_Orb)
           !
           !HxV_Hyb: IMP DW <--> BATH DW hoppings
           do kp=1,Nbath
              alfa=1+kp
              if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. (Ndws(iorb,1)==1) .AND. (Ndws(iorb,alfa)==0) )then
                 call c(1,mdw,k1,sg1)
                 call cdg(alfa,k1,k2,sg2)                 
                 Jndices       = Indices
                 Jndices(iorb) = binary_search(Hs(iorb+Ns_Ud)%map,k2)
                 call indices2state(Jndices,[DimUps,DimDws],j)
                 htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
                 !
                 Hvt(i) = Hvt(i) + htmp*Vt(j)
              endif
              !
              if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. (Ndws(iorb,1)==0) .AND. (Ndws(iorb,alfa)==1) )then
                 call c(alfa,mdw,k1,sg1)
                 call cdg(1,k1,k2,sg2)
                 Jndices       = Indices
                 Jndices(iorb) = binary_search(Hs(iorb+Ns_Ud)%map,k2)
                 call indices2state(Jndices,[DimUps,DimDws],j)
                 htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
                 !
                 Hvt(i) = Hvt(i) + htmp*Vt(j)
              endif
           enddo
           !
        enddo
     enddo
  enddo






