  do i=1,Dim
     call state2indices(i,[DimUps,DimDws],Indices)
     do iorb=1,Norb
        Mups(iorb) = Hs(iorb)%map(Indices(iorb))
        Mdws(iorb) = Hs(Norb+iorb)%map(Indices(Norb+iorb))
        !
        Nups(iorb,:) = Bdecomp(Mups(iorb),Ns_Orb)
        Ndws(iorb,:) = Bdecomp(Mdws(iorb),Ns_Orb)
        !
        Nup(iorb) = Nups(iorb,1)
        Ndw(iorb) = Ndws(iorb,1)
     enddo

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
        call sp_insert_element(spH0d,htmp,i,i)
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
        call sp_insert_element(spH0d,htmp,i,i)
        !
     endif
  enddo
