  do i=first_state,last_state
     m = H%map(i)
     impi = i-ishift
     ib = bdecomp(m,2*Ns)
     !
     do iorb=1,Norb
        nup(iorb)=dble(ib(iorb))
        ndw(iorb)=dble(ib(iorb+Ns))
     enddo


     if(bath_type/="replica") then
        !
        !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
        htmp=zero
        do iorb=1,size(dmft_bath%e,2)
           do kp=1,Nbath
              alfa=getBathStride(iorb,kp)
              htmp =htmp + dmft_bath%e(1,iorb,kp)*ib(alfa)        !UP
              htmp =htmp + dmft_bath%e(Nspin,iorb,kp)*ib(alfa+Ns) !DW
           enddo
        enddo
        !
        hv(impi) = hv(impi) + htmp*vin(i)
        !
     else
        !
        !
        !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
        htmp=zero
        do kp=1,Nbath
           do iorb=1,Norb
              alfa = getBathStride(iorb,kp)
              htmp = htmp + dmft_bath%h(1,1,iorb,iorb,kp)*ib(alfa)              !UP
              htmp = htmp + dmft_bath%h(Nspin,Nspin,iorb,iorb,kp)*ib(alfa+Ns)   !DW
           enddo
        enddo
        !
        hv(impi) = hv(impi) + htmp*vin(i)
        !
        !off-diagonal elements
        !
        !1. same spin:
        do kp=1,Nbath
           do iorb=1,Norb
              do jorb=1,Norb
                 !this loop considers only the orbital off-diagonal terms
                 !because iorb=jorb can not have simultaneously
                 !occupation 0 and 1, as required by this if Jcondition:
                 !UP
                 alfa = getBathStride(iorb,kp)
                 beta = getBathStride(jorb,kp)
                 Jcondition = &
                      (dmft_bath%h(1,1,iorb,jorb,kp)/=zero)               .AND. &
                      (ib(beta)==1)                                       .AND. &
                      (ib(alfa)==0)
                 if (Jcondition)then
                    call c(beta,m,k1,sg1)
                    call cdg(alfa,k1,k2,sg2)
                    j = binary_search(H%map,k2)
                    htmp = dmft_bath%h(1,1,iorb,jorb,kp)*sg1*sg2
                    !
                    hv(impi) = hv(impi) + htmp*vin(j)
                    !
                 endif
                 !DW
                 alfa = getBathStride(iorb,kp) + Ns
                 beta = getBathStride(jorb,kp) + Ns
                 Jcondition = &
                      (dmft_bath%h(Nspin,Nspin,iorb,jorb,kp)/=zero)       .AND. &
                      (ib(beta)==1)                                       .AND. &
                      (ib(alfa)==0)
                 if (Jcondition)then
                    call c(beta,m,k1,sg1)
                    call cdg(alfa,k1,k2,sg2)
                    j = binary_search(H%map,k2)
                    htmp = dmft_bath%h(Nspin,Nspin,iorb,jorb,kp)*sg1*sg2
                    !
                    hv(impi) = hv(impi) + htmp*vin(j)
                    !
                 endif
              enddo
           enddo
        enddo
        !
     endif



  enddo
