  do idw=first_state_dw,last_state_dw
     mdw  = Hs(2)%map(idw)
     ndw  = bdecomp(mdw,Ns)
     !
     do iup=first_state_up(idw),last_state_up(idw)
        mup  = Hs(1)%map(iup)
        nup  = bdecomp(mup,Ns)
        !
        !MPI Shifts
        i    = iup + (idw-1)*dimUp
        impi = i - ishift
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
                 htmp = htmp + dmft_bath%h(1    ,    1,iorb,iorb,kp)*nup(alfa) !UP
                 htmp = htmp + dmft_bath%h(Nspin,Nspin,iorb,iorb,kp)*ndw(alfa) !DW
              enddo
           enddo
           !
           hv(impi) = hv(impi) + htmp*vin(i)
           !
        endif
     enddo
  enddo



  !off-diagonal elements
  !
  !this loop considers only the orbital off-diagonal terms
  !because iorb=jorb can not have simultaneously
  !occupation 0 and 1, as required by this if Jcondition:
  !
  if(bath_type=="replica") then
     !
     !UP:
     do idw=first_state_dw,last_state_dw
        do iup=first_state_up(idw),last_state_up(idw)
           mup  = Hs(1)%map(iup)
           nup  = bdecomp(mup,Ns)
           !
           !MPI Shifts
           i    = iup + (idw-1)*dimUp
           impi = i - ishift
           !
           do kp=1,Nbath
              do iorb=1,Norb
                 do jorb=1,Norb
                    !
                    alfa = getBathStride(iorb,kp)
                    beta = getBathStride(jorb,kp)
                    Jcondition = &
                         (dmft_bath%h(1,1,iorb,jorb,kp)/=zero) .AND. (nup(beta)==1) .AND. (nup(alfa)==0)
                    !
                    if (Jcondition)then
                       call c(beta,mup,k1,sg1)
                       call cdg(alfa,k1,k2,sg2)
                       jup = binary_search(Hs(1)%map,k2)
                       j   = jup + (idw-1)*DimUp
                       !
                       htmp = dmft_bath%h(1,1,iorb,jorb,kp)*sg1*sg2
                       !
                       hv(impi) = hv(impi) + htmp*vin(j)
                       !
                    endif
                 enddo
              enddo
           enddo
           !
        end do
     end do
     !
     !
     !DW
     do idw=first_state_dw,last_state_dw
        do iup=first_state_up(idw),last_state_up(idw)
           mdw  = Hs(2)%map(idw)
           ndw  = bdecomp(mdw,Ns)
           !
           !MPI Shifts
           i    = iup + (idw-1)*dimUp
           impi = i - ishift
           !
           do kp=1,Nbath
              do iorb=1,Norb
                 do jorb=1,Norb
                    !
                    alfa = getBathStride(iorb,kp)
                    beta = getBathStride(jorb,kp)
                    Jcondition = &
                         (dmft_bath%h(Nspin,Nspin,iorb,jorb,kp)/=zero) .AND. (ndw(beta)==1) .AND. (ndw(alfa)==0)
                    !
                    if (Jcondition)then
                       call c(beta,mdw,k1,sg1)
                       call cdg(alfa,k1,k2,sg2)
                       jdw = binary_search(Hs(2)%map,k2)
                       j   = iup + (jdw-1)*DimUp
                       !
                       htmp = dmft_bath%h(Nspin,Nspin,iorb,jorb,kp)*sg1*sg2
                       !
                       hv(impi) = hv(impi) + htmp*vin(j)
                       !
                    endif
                 enddo
              enddo
           enddo
           !
        end do
     end do



  end if
