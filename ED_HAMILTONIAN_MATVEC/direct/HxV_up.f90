  DWstates_loop:  do idw=1,MpiQdw
     !
     !
     do iup=1,DimUp
        mup  = Hs(1)%map(iup)
        nup  = bdecomp(mup,Ns)
        !
        i    = iup + (idw-1)*dimUp
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
                 jup = binary_search(Hs(1)%map,k2)
                 j   = jup + (idw-1)*DimUp
                 htmp = impHloc(1,1,iorb,jorb)*sg1*sg2
                 !
                 Hv(i) = Hv(i) + htmp*vin(j)
                 !
              endif
           enddo
        enddo
        !
        !
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
                       jup = binary_search(Hs(1)%map,k2)
                       j   = jup + (idw-1)*DimUp
                       !
                       htmp = dmft_bath%h(1,1,iorb,jorb,kp)*sg1*sg2
                       !
                       hv(i) = hv(i) + htmp*vin(j)
                       !
                    endif
                 enddo
              enddo
           enddo
        end if
        !
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
                 jup = binary_search(Hs(1)%map,k2)
                 j   = jup + (idw-1)*DimUp
                 !
                 htmp = diag_hybr(1,iorb,kp)*sg1*sg2
                 !
                 ! hv(impi) = hv(impi) + htmp*vin(j)
                 hv(i) = hv(i) + htmp*vin(j)
                 !
              endif
              !
              if( (diag_hybr(1,iorb,kp)/=0d0) .AND. &
                   (nup(iorb)==0) .AND. (nup(alfa)==1) )then
                 call c(alfa,mup,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 jup = binary_search(Hs(1)%map,k2)
                 j   = jup + (idw-1)*DimUp
                 !
                 htmp = diag_hybr(1,iorb,kp)*sg1*sg2
                 !
                 hv(i) = hv(i) + htmp*vin(j)
                 !
              endif
           enddo
        enddo
        !
        !
     end do
     !
     !
     !
  end do DWSTATES_LOOP






