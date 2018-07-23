  !off-diagonal elements
  !
  !this loop considers only the orbital off-diagonal terms
  !because iorb=jorb can not have simultaneously
  !occupation 0 and 1, as required by this if Jcondition:
  if(bath_type=="replica") then
     !
     !UP:
     !do iup=1,DimUp             !first_state_up,last_state_up
     do iup=first_state_up,last_state_up
        mup  = Hs(1)%map(iup)
        nup  = bdecomp(mup,Ns)
        impi_up = iup - IshiftUp
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
                    !
                    vecNnz_up(impi_up) = vecNnz_up(impi_up) + 1
                 endif
              enddo
           enddo
        enddo
        !
     enddo
     !
     !
     !DW:
     ! do idw=map_first_state_dw(1),map_last_state_dw(1)
     do idw=first_state_dw,last_state_dw
        mdw  = Hs(2)%map(idw)
        ndw  = bdecomp(mdw,Ns)
        impi_dw = idw - IshiftDw
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
                    !
                    vecNnz_dw(impi_dw) = vecNnz_dw(impi_dw) + 1
                 endif
              enddo
           enddo
        enddo
        !
     enddo
     !
  endif
