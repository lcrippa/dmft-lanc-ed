  do iup=1,MpiQup
     do idw=1,DimDw
        mdw  = Hs(2)%map(idw)
        ndw  = bdecomp(mdw,Ns)
        i    = idw + (iup-1)*DimDw
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
                 jdw = binary_search(Hs(2)%map,k2)
                 j   = jdw + (jup-1)*DimDw
                 htmp = impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2
                 !
                 Hvt(i) = Hvt(i) + htmp*Vt(j)
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
                       jdw = binary_search(Hs(2)%map,k2)
                       j   = iup + (jdw-1)*DimUp
                       htmp = dmft_bath%h(Nspin,Nspin,iorb,jorb,kp)*sg1*sg2
                       !
                       Hvt(i) = Hvt(i) + htmp*Vt(j)
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
                 jdw = binary_search(Hs(2)%map,k2)
                 j   = jdw + (jup-1)*DimDw
                 htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
                 !
                 Hvt(i) = Hvt(i) + htmp*Vt(j)
              endif
              !
              if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. &
                   (ndw(iorb)==0) .AND. (ndw(alfa)==1) )then
                 call c(alfa,mdw,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 jdw = binary_search(Hs(2)%map,k2)
                 j   = jdw + (jup-1)*DimDw
                 htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
                 !
                 Hvt(i) = Hvt(i) + htmp*Vt(j)
              endif
           enddo
        enddo
        !
        !
     end do
  end do
