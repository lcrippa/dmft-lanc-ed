  do i=1,Dim
     call state2indices(i,[DimUps,DimDws],Indices)
     !
     do iud=1,Ns_Ud
        mdw = Hs(iud+Ns_Ud)%map(Indices(iud+Ns_ud))
        Ndws(iud,:) = Bdecomp(mdw,Ns_Orb) ![1+Nbath]*Norb
     enddo
     !
     Ndw = Breorder(Ndws)
     !
     do iud=1,Ns_Ud             !== 1
        !
        !
        ! HxV_imp: off-diagonal terms
        do iorb=1,Norb
           do jorb=1,Norb
              Jcondition = &
                   (impHloc(Nspin,Nspin,iorb,jorb)/=zero) .AND. &
                   (ndw(jorb)==1) .AND. (ndw(iorb)==0)
              if (Jcondition) then
                 call c(jorb,mdw,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 Jndices = Indices
                 Jndices(iud+Ns_Ud) = binary_search(Hs(iud+Ns_Ud)%map,k2)
                 call indices2state(Jndices,[DimUps,DimDws],j)
                 htmp = impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2
                 !
                 Hv(i) = Hv(i) + htmp*Vin(j)
              endif
           enddo
        enddo
        !
        !
        !HxV_bath: off-diagonal terms
        if(bath_type=="replica") then
           do kp=1,Nbath
              do iorb=1,Norb
                 do jorb=1,Norb
                    !
                    alfa = getBathStride(iorb,kp)
                    beta = getBathStride(jorb,kp)
                    Jcondition = &
                         (dmft_bath%h(Nspin,Nspin,iorb,jorb,kp)/=zero) .AND. &
                         (ndw(beta)==1) .AND. (ndw(alfa)==0)
                    !
                    if (Jcondition)then
                       call c(beta,mdw,k1,sg1)
                       call cdg(alfa,k1,k2,sg2)
                       Jndices = Indices
                       Jndices(iud+Ns_Ud) = binary_search(Hs(iud+Ns_Ud)%map,k2)
                       call indices2state(Jndices,[DimUps,DimDws],j)
                       htmp = dmft_bath%h(Nspin,Nspin,iorb,jorb,kp)*sg1*sg2
                       !
                       Hv(i) = Hv(i) + htmp*Vin(j)
                    endif
                 enddo
              enddo
           enddo
        end if
        !
        !
        !HxV_Hyb: IMP DW <--> BATH DW hoppings
        do iorb=1,Norb
           do kp=1,Nbath
              alfa=getBathStride(iorb,kp)
              !
              if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. &
                   (ndw(iorb)==1) .AND. (ndw(alfa)==0) )then
                 call c(iorb,mdw,k1,sg1)
                 call cdg(alfa,k1,k2,sg2)
                 Jndices = Indices
                 Jndices(iud+Ns_Ud) = binary_search(Hs(iud+Ns_Ud)%map,k2)
                 call indices2state(Jndices,[DimUps,DimDws],j)
                 htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
                 !
                 Hv(i) = Hv(i) + htmp*Vin(j)
              endif
              !
              if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. &
                   (ndw(iorb)==0) .AND. (ndw(alfa)==1) )then
                 call c(alfa,mdw,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 Jndices = Indices
                 Jndices(iud+Ns_Ud) = binary_search(Hs(iud+Ns_Ud)%map,k2)
                 call indices2state(Jndices,[DimUps,DimDws],j)
                 htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
                 !
                 Hv(i) = Hv(i) + htmp*Vin(j)
              endif
           enddo
        enddo
        !
        !
     end do
  end do
