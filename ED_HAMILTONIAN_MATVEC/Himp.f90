  ! do i=first_state,last_state
  !    iup = iup_index(i,DimUp)
  !    idw = idw_index(i,DimUp)
  !    !
  !    nup  = bdecomp(Hs(1)%map(iup),Ns)
  !    ndw  = bdecomp(Hs(2)%map(idw),Ns)
  !    !
  !    impi = i - ishift
  
  
  do idw=map_first_state_dw(1),map_last_state_dw(1)
     mdw  = Hs(2)%map(idw)
     ndw  = bdecomp(mdw,Ns)
     !
     do iup=map_first_state_up(idw),map_last_state_up(idw)
        mup  = Hs(1)%map(iup)
        nup  = bdecomp(mup,Ns)
        !
        !MPI Shifts
        i    = iup + (idw-1)*dimUp
        impi = i - ishift
        !
        !
        !Diagonal Elements, i.e. local part
        htmp = zero
        do iorb=1,Norb
           htmp = htmp + impHloc(1,1,iorb,iorb)*nup(iorb)
           htmp = htmp + impHloc(Nspin,Nspin,iorb,iorb)*ndw(iorb)
           htmp = htmp - xmu*(nup(iorb)+ndw(iorb))
        enddo
        !
        call sp_insert_element(spH0,htmp,impi,i)
        !
     enddo
  enddo



  !Off-diagonal elements, i.e. non-local part
  !this loop considers only the orbital off-diagonal terms
  !because iorb=jorb can not have simultaneously
  !occupation 0 and 1, as required by this if Jcondition:
  !
  !UP
  do iup=1,DimUp                !first_state_up,last_state_up
     mup  = Hs(1)%map(iup)
     nup  = bdecomp(mup,Ns)
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
              htmp = impHloc(1,1,iorb,jorb)*sg1*sg2
              !
              call sp_insert_element(spH0up,htmp,iup,jup)
              !
           endif
        enddo
     enddo
     !
  end do

  !DW
  do idw=map_first_state_dw(1),map_last_state_dw(1)
     mdw  = Hs(2)%map(idw)
     ndw  = bdecomp(mdw,Ns)
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
              htmp = impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2
              !
              call sp_insert_element(spH0dw,htmp,idw,jdw)
              !
           endif
        enddo
     enddo
     !
  enddo

