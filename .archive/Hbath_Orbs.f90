  do i=MpiIstart,MpiIend
     call state2indices(i,[DimUps,DimDws],Indices)
     !
     do iud=1,Ns_Ud
        mup = Hs(iud)%map(Indices(iud))
        mdw = Hs(iud+Ns_Ud)%map(Indices(iud+Ns_Ud))
        !
        Nups(iud,:) = Bdecomp(mup,Ns_Orb)
        Ndws(iud,:) = Bdecomp(mdw,Ns_Orb)
     enddo
     !
     !
     if(bath_type/="replica") then
        !
        !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
        htmp=zero
        do iorb=1,size(dmft_bath%e,2)
           do kp=1,Nbath
              htmp =htmp + dmft_bath%e(1    ,iorb,kp)*Nups(iorb,1+kp) !UP
              htmp =htmp + dmft_bath%e(Nspin,iorb,kp)*Ndws(iorb,1+kp) !DW
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
              htmp = htmp + dmft_bath%h(1    ,    1,iorb,iorb,kp)*Nups(iorb,1+kp) !UP
              htmp = htmp + dmft_bath%h(Nspin,Nspin,iorb,iorb,kp)*Ndws(iorb,1+kp) !DW
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
