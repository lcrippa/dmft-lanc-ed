  do i=1,Dim
     call state2indices(i,[DimUps,DimDws],Indices)
     !
     do iud=1,Ns_Ud
        mup = Hs(iud)%map(Indices(iud))
        Nups(iud,:) = Bdecomp(mup,Ns_Orb) ![1+Nbath]*Norb
     enddo
     !
     Nup = Breorder(Nups)
     !
     do iud=1,Ns_Ud             !== 1
        !
        !
        ! HxV_imp: off-diagonal terms
        do iorb=1,Norb
           do jorb=1,Norb
              Jcondition = &
                   (impHloc(1,1,iorb,jorb)/=zero) .AND. &
                   (nup(jorb)==1) .AND. (nup(iorb)==0)
              if (Jcondition) then
                 call c(jorb,mup,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 Jndices = Indices
                 Jndices(iud) = binary_search(Hs(iud)%map,k2)
                 call indices2state(Jndices,[DimUps,DimDws],j)
                 htmp = impHloc(1,1,iorb,jorb)*sg1*sg2
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
                         (dmft_bath%h(1,1,iorb,jorb,kp)/=zero) .AND.&
                         (nup(beta)==1) .AND. (nup(alfa)==0)
                    !
                    if (Jcondition)then
                       call c(beta,mup,k1,sg1)
                       call cdg(alfa,k1,k2,sg2)
                       Jndices = Indices
                       Jndices(iud) = binary_search(Hs(iud)%map,k2)
                       call indices2state(Jndices,[DimUps,DimDws],j)
                       htmp = dmft_bath%h(1,1,iorb,jorb,kp)*sg1*sg2
                       !
                       Hv(i) = Hv(i) + htmp*Vin(j)
                    endif
                 enddo
              enddo
           enddo
        end if
        !
        !
        !HxV_Hyb: IMP UP <--> BATH UP hoppings        
        do iorb=1,Norb
           do kp=1,Nbath
              alfa=getBathStride(iorb,kp)
              !
              if( (diag_hybr(1,iorb,kp)/=0d0) .AND. &
                   (nup(iorb)==1) .AND. (nup(alfa)==0) )then              
                 call c(iorb,mup,k1,sg1)
                 call cdg(alfa,k1,k2,sg2)
                 Jndices = Indices
                 Jndices(iud) = binary_search(Hs(iud)%map,k2)
                 call indices2state(Jndices,[DimUps,DimDws],j)
                 htmp = diag_hybr(1,iorb,kp)*sg1*sg2
                 !
                 Hv(i) = Hv(i) + htmp*Vin(j)
              endif
              !
              if( (diag_hybr(1,iorb,kp)/=0d0) .AND. &
                   (nup(iorb)==0) .AND. (nup(alfa)==1) )then
                 call c(alfa,mup,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 Jndices = Indices
                 Jndices(iud) = binary_search(Hs(iud)%map,k2)
                 call indices2state(Jndices,[DimUps,DimDws],j)
                 htmp = diag_hybr(1,iorb,kp)*sg1*sg2
                 !
                 Hv(i) = Hv(i) + htmp*Vin(j)
              endif
           enddo
        enddo
        !
        !
     end do
  enddo






