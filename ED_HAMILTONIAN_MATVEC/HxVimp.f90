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

        !Diagonal Elements, i.e. local part
        htmp = zero
        do iorb=1,Norb
           htmp = htmp + impHloc(1,1,iorb,iorb)*nup(iorb)
           htmp = htmp + impHloc(Nspin,Nspin,iorb,iorb)*ndw(iorb)
           htmp = htmp - xmu*(nup(iorb)+ndw(iorb))
        enddo
        !
        hv(impi) = hv(impi) + htmp*vin(i)
        !
     enddo
  enddo

  !Off-diagonal elements, i.e. non-local part
  !this loop considers only the orbital off-diagonal terms
  !because iorb=jorb can not have simultaneously
  !occupation 0 and 1, as required by this if Jcondition:
  !
  !UP
  do idw=first_state_dw,last_state_dw
     do iup=first_state_up(idw),last_state_up(idw)
        mup  = Hs(1)%map(iup)
        nup  = bdecomp(mup,Ns)
        !
        !MPI Shifts
        i    = iup + (idw-1)*dimUp
        impi = i - ishift
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
                 j   = jup + (idw-1)*DimUp
                 htmp = impHloc(1,1,iorb,jorb)*sg1*sg2
                 !
                 hv(impi) = hv(impi) + htmp*vin(j)
                 !
              endif
           enddo
        enddo
        !
     end do
  end do


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
        do iorb=1,Norb
           do jorb=1,Norb
              Jcondition = &
                   (impHloc(Nspin,Nspin,iorb,jorb)/=zero) .AND. &
                   (ndw(jorb)==1) .AND. (ndw(iorb)==0)
              if (Jcondition) then
                 call c(jorb,mdw,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 jdw = binary_search(Hs(1)%map,k2)
                 j   = iup + (jdw-1)*DimUp
                 htmp = impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2
                 !
                 hv(impi) = hv(impi) + htmp*vin(j)
                 !
              endif
           enddo
        enddo
        !
     end do
  end do
