  ! IMP UP <--> BATH UP
  ! do iup=1,DimUp!first_state_up,last_state_up
  do iup=first_state_up,last_state_up
     mup  = Hs(1)%map(iup)
     nup  = bdecomp(mup,Ns)
     impi_up = iup - IshiftUp
     !
     do iorb=1,Norb
        do kp=1,Nbath
           alfa=getBathStride(iorb,kp)
           !
           if( (diag_hybr(1,iorb,kp)/=0d0) .AND. (nup(iorb)==1) .AND. (nup(alfa)==0) )then              
              call c(iorb,mup,k1,sg1)
              call cdg(alfa,k1,k2,sg2)
              jup = binary_search(Hs(1)%map,k2)
              htmp = diag_hybr(1,iorb,kp)*sg1*sg2
              !
              select case (ed_sparse_format)
              case default
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0up,htmp,iup,jup)
                 case (.false.)
                    call sp_insert_element(spH0up,htmp,iup,jup)
                 end select
              case ("ELL")
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,dpH0up,htmp,iup,jup)
                 case (.false.)
                    call sp_insert_element(dpH0up,htmp,iup,jup)
                 end select
              end select
              !
           endif
           !
           if( (diag_hybr(1,iorb,kp)/=0d0) .AND. (nup(iorb)==0) .AND. (nup(alfa)==1) )then
              call c(alfa,mup,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              jup=binary_search(Hs(1)%map,k2)
              htmp = diag_hybr(1,iorb,kp)*sg1*sg2
              !
              select case (ed_sparse_format)
              case default
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0up,htmp,iup,jup)
                 case (.false.)
                    call sp_insert_element(spH0up,htmp,iup,jup)
                 end select
              case ("ELL")
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,dpH0up,htmp,iup,jup)
                 case (.false.)
                    call sp_insert_element(dpH0up,htmp,iup,jup)
                 end select
              end select
              !
           endif
        enddo
     enddo
     !
  enddo

  !IMP DW <--> BATH DW
  ! do idw=map_first_state_dw(1),map_last_state_dw(1)
  do idw=first_state_dw,last_state_dw
     mdw  = Hs(2)%map(idw)
     ndw  = bdecomp(mdw,Ns)
     impi_dw = idw - IshiftDw
     !
     do iorb=1,Norb
        do kp=1,Nbath
           alfa=getBathStride(iorb,kp)
           !
           if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. (ndw(iorb)==1) .AND. (ndw(alfa)==0) )then
              call c(iorb,mdw,k1,sg1)
              call cdg(alfa,k1,k2,sg2)
              jdw=binary_search(Hs(2)%map,k2)
              htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
              !
              select case (ed_sparse_format)
              case default
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0dw,htmp,idw,jdw)
                 case (.false.)
                    call sp_insert_element(spH0dw,htmp,idw,jdw)
                 end select
              case ("ELL")
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,dpH0dw,htmp,idw,jdw)
                 case (.false.)
                    call sp_insert_element(dpH0dw,htmp,idw,jdw)
                 end select
              end select
              !
           endif
           if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. (ndw(iorb)==0) .AND. (ndw(alfa)==1) )then
              call c(alfa,mdw,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              jdw=binary_search(Hs(2)%map,k2)
              htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
              !
              select case (ed_sparse_format)
              case default
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0dw,htmp,idw,jdw)
                 case (.false.)
                    call sp_insert_element(spH0dw,htmp,idw,jdw)
                 end select
              case ("ELL")
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,dpH0dw,htmp,idw,jdw)
                 case (.false.)
                    call sp_insert_element(dpH0dw,htmp,idw,jdw)
                 end select
              end select
              !
           endif
        enddo
     enddo
     !
  enddo
