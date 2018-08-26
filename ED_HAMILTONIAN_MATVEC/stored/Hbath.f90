  do i=MpiIstart,MpiIend
     iup = iup_index(i,DimUp)
     idw = idw_index(i,DimUp)
     !
     mup = Hs(1)%map(iup)
     mdw = Hs(2)%map(idw)
     !
     nup = bdecomp(mup,Ns)
     ndw = bdecomp(mdw,Ns)
     !
     if(bath_type/="replica") then
        !
        !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
        htmp=zero
        do iorb=1,size(dmft_bath%e,2)
           do kp=1,Nbath
              alfa = getBathStride(iorb,kp)
              htmp =htmp + dmft_bath%e(1    ,iorb,kp)*nup(alfa) !UP
              htmp =htmp + dmft_bath%e(Nspin,iorb,kp)*ndw(alfa) !DW
           enddo
        enddo
        !
        select case(MpiStatus)
        case (.true.)
           call sp_insert_element(MpiComm,spH0d,htmp,i,i)
        case (.false.)
           call sp_insert_element(spH0d,htmp,i,i)
        end select
        !
     else
        !
        !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
        htmp=zero
        do kp=1,Nbath
           do iorb=1,Norb
              alfa = getBathStride(iorb,kp)
              htmp = htmp + dmft_bath%h(1    ,    1,iorb,iorb,kp)*nup(alfa) !UP
              htmp = htmp + dmft_bath%h(Nspin,Nspin,iorb,iorb,kp)*ndw(alfa) !DW
           enddo
        enddo
        !
        select case(MpiStatus)
        case (.true.)
           call sp_insert_element(MpiComm,spH0d,htmp,i,i)
        case (.false.)
           call sp_insert_element(spH0d,htmp,i,i)
        end select
        !
     endif

  enddo


  !off-diagonal elements
  !
  !this loop considers only the orbital off-diagonal terms
  !because iorb=jorb can not have simultaneously
  !occupation 0 and 1, as required by this if Jcondition:
  if(bath_type=="replica") then
     !
     !UP:
     do iup=1,DimUp
        mup  = Hs(1)%map(iup)
        nup  = bdecomp(mup,Ns)
        !
        do kp=1,Nbath
           do iorb=1,Norb
              do jorb=1,Norb
                 !
                 alfa = getBathStride(iorb,kp)
                 beta = getBathStride(jorb,kp)
                 Jcondition = &
                      (dmft_bath%h(1,1,iorb,jorb,kp)/=zero) &
                      .AND. (nup(beta)==1) .AND. (nup(alfa)==0)
                 !
                 if (Jcondition)then
                    call c(beta,mup,k1,sg1)
                    call cdg(alfa,k1,k2,sg2)
                    jup = binary_search(Hs(1)%map,k2)
                    htmp = dmft_bath%h(1,1,iorb,jorb,kp)*sg1*sg2
                    !
                    call sp_insert_element(spH0up,htmp,iup,jup)
                    !
                 endif
              enddo
           enddo
        enddo
        !
     enddo
     !
     !
     !DW:
     do idw=1,DimDw
        mdw  = Hs(2)%map(idw)
        ndw  = bdecomp(mdw,Ns)
        !
        do kp=1,Nbath
           do iorb=1,Norb
              do jorb=1,Norb
                 !
                 alfa = getBathStride(iorb,kp)
                 beta = getBathStride(jorb,kp)
                 Jcondition = &
                      (dmft_bath%h(Nspin,Nspin,iorb,jorb,kp)/=zero) &
                      .AND. (ndw(beta)==1) .AND. (ndw(alfa)==0)
                 !
                 if (Jcondition)then
                    call c(beta,mdw,k1,sg1)
                    call cdg(alfa,k1,k2,sg2)
                    jdw = binary_search(Hs(2)%map,k2)
                    htmp = dmft_bath%h(Nspin,Nspin,iorb,jorb,kp)*sg1*sg2
                    !
                    call sp_insert_element(spH0dw,htmp,idw,jdw)
                    !
                 endif
              enddo
           enddo
        enddo
        !
     enddo
     !
  endif
