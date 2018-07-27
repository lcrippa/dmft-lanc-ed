  !Off-diagonal elements, i.e. non-local part
  !this loop considers only the orbital off-diagonal terms
  !because iorb=jorb can not have simultaneously
  !occupation 0 and 1, as required by this if Jcondition:
  !
  !UP
  do iup=first_state_up,last_state_up
     mup  = Hs(1)%map(iup)
     nup  = bdecomp(mup,Ns)
     impi_up = iup - IshiftUp
     !
     do iorb=1,Norb
        do jorb=1,Norb
           Jcondition = &
                (impHloc(1,1,iorb,jorb)/=zero) .AND. &
                (nup(jorb)==1) .AND. (nup(iorb)==0)
           if (Jcondition) then
              call c(jorb,mup,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              jup = binary_search(Hs(1)%map,k2)
              vecNnz_up(impi_up) = vecNnz_up(impi_up) + 1
           endif
        enddo
     enddo
     !
  end do

  !DW
  do idw=first_state_dw,last_state_dw
     mdw  = Hs(2)%map(idw)
     ndw  = bdecomp(mdw,Ns)
     impi_dw = idw - IshiftDw
     !
     do iorb=1,Norb
        do jorb=1,Norb
           Jcondition = &
                (impHloc(Nspin,Nspin,iorb,jorb)/=zero) .AND. &
                (ndw(jorb)==1) .AND. (ndw(iorb)==0)
           if (Jcondition) then
              call c(jorb,mdw,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              jdw = binary_search(Hs(2)%map,k2)
              vecNnz_dw(impi_dw) = vecNnz_dw(impi_dw) + 1
           endif
        enddo
     enddo
     !
  enddo

