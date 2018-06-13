subroutine ed_get_all_sigma_real_1(Sreal)
  complex(8),dimension(2,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
  Sreal(1,:,:,:,:,:) = impSreal(:,:,:,:,:)
  Sreal(2,:,:,:,:,:) = impSAreal(:,:,:,:,:)
end subroutine ed_get_all_sigma_real_1

subroutine ed_get_all_sigma_real_2(Sreal)
  complex(8),dimension(2,Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Sreal
  integer  :: io,jo,iorb,jorb,ispin,jspin
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              io = iorb + (ispin-1)*Norb
              jo = jorb + (jspin-1)*Norb
              Sreal(1,io,jo,:) = impSreal(ispin,jspin,iorb,jorb,:)
              Sreal(2,io,jo,:) = impSAreal(ispin,jspin,iorb,jorb,:)
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_all_sigma_real_2

subroutine ed_get_all_sigma_real_3(Sreal,ispin,jspin,iorb,jorb)
  complex(8),dimension(2,Lreal),intent(inout) :: Sreal
  integer                          :: iorb,jorb,ispin,jspin
  Sreal(1,:) = impSreal(ispin,jspin,iorb,jorb,:)
  Sreal(2,:) = impSAreal(ispin,jspin,iorb,jorb,:)
end subroutine ed_get_all_sigma_real_3









subroutine ed_get_all_sigma_real_lattice_1(Sreal,Nsites)
  integer                                                                  :: Nsites
  complex(8),dimension(2,Nsites,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
  Sreal(1,1:Nsites,:,:,:,:,:) = Srealii(1:Nsites,:,:,:,:,:)
  Sreal(2,1:Nsites,:,:,:,:,:) = SArealii(1:Nsites,:,:,:,:,:)
end subroutine ed_get_all_sigma_real_lattice_1

subroutine ed_get_all_sigma_real_lattice_2(Sreal,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(2,Nsites,Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Sreal
  integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
  do ilat=1,Nsites
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 io = iorb + (ispin-1)*Norb
                 jo = jorb + (jspin-1)*Norb
                 Sreal(1,ilat,io,jo,:) = Srealii(ilat,ispin,jspin,iorb,jorb,:)
                 Sreal(2,ilat,io,jo,:) = SArealii(ilat,ispin,jspin,iorb,jorb,:)
              enddo
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_all_sigma_real_lattice_2

subroutine ed_get_all_sigma_real_lattice_3(Sreal,Nsites,ispin,jspin,iorb,jorb)
  integer                                            :: Nsites
  complex(8),dimension(2,Nsites,Lreal),intent(inout) :: Sreal
  integer                                            :: iorb,jorb,ispin,jspin
  Sreal(1,1:Nsites,:) = Srealii(1:Nsites,ispin,jspin,iorb,jorb,:)
  Sreal(2,1:Nsites,:) = SArealii(1:Nsites,ispin,jspin,iorb,jorb,:)
end subroutine ed_get_all_sigma_real_lattice_3







subroutine ed_get_all_sigma_real_lattice_11(Sreal,ilat)
  integer                                                         :: ilat
  complex(8),dimension(2,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
  Sreal(1,:,:,:,:,:) = Srealii(ilat,:,:,:,:,:)
  Sreal(2,:,:,:,:,:) = SArealii(ilat,:,:,:,:,:)
end subroutine ed_get_all_sigma_real_lattice_11

subroutine ed_get_all_sigma_real_lattice_21(Sreal,ilat)
  integer                                                         :: ilat
  complex(8),dimension(2,Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Sreal
  integer                                                         :: io,jo,iorb,jorb,ispin,jspin
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              io = iorb + (ispin-1)*Norb
              jo = jorb + (jspin-1)*Norb
              Sreal(1,io,jo,:) = Srealii(ilat,ispin,jspin,iorb,jorb,:)
              Sreal(2,io,jo,:) = SArealii(ilat,ispin,jspin,iorb,jorb,:)
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_all_sigma_real_lattice_21

subroutine ed_get_all_sigma_real_lattice_31(Sreal,ilat,ispin,jspin,iorb,jorb)
  integer                                     :: ilat
  complex(8),dimension(2,Lreal),intent(inout) :: Sreal
  integer                                     :: iorb,jorb,ispin,jspin
  Sreal(1,:) = Srealii(ilat,ispin,jspin,iorb,jorb,:)
  Sreal(2,:) = SArealii(ilat,ispin,jspin,iorb,jorb,:)
end subroutine ed_get_all_sigma_real_lattice_31



