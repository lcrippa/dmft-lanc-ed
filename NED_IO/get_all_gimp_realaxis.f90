subroutine ed_get_all_gimp_real_1(Greal)
  complex(8),dimension(2,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
  Greal(1,:,:,:,:,:) = impGreal(:,:,:,:,:)
  Greal(2,:,:,:,:,:) = impFreal(:,:,:,:,:)
end subroutine ed_get_all_gimp_real_1

subroutine ed_get_all_gimp_real_2(Greal)
  complex(8),dimension(2,Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Greal
  integer  :: io,jo,iorb,jorb,ispin,jspin
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              io = iorb + (ispin-1)*Norb
              jo = jorb + (jspin-1)*Norb
              Greal(1,io,jo,:) = impGreal(ispin,jspin,iorb,jorb,:)
              Greal(2,io,jo,:) = impFreal(ispin,jspin,iorb,jorb,:)
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_all_gimp_real_2

subroutine ed_get_all_gimp_real_3(Greal,ispin,jspin,iorb,jorb)
  complex(8),dimension(2,Lreal),intent(inout) :: Greal
  integer                          :: iorb,jorb,ispin,jspin
  Greal(1,:) = impGreal(ispin,jspin,iorb,jorb,:)
  Greal(2,:) = impFreal(ispin,jspin,iorb,jorb,:)
end subroutine ed_get_all_gimp_real_3









subroutine ed_get_all_gimp_real_lattice_1(Greal,Nsites)
  integer                                                                  :: Nsites
  complex(8),dimension(2,Nsites,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
  Greal(1,1:Nsites,:,:,:,:,:) = Grealii(1:Nsites,:,:,:,:,:)
  Greal(2,1:Nsites,:,:,:,:,:) = Frealii(1:Nsites,:,:,:,:,:)
end subroutine ed_get_all_gimp_real_lattice_1

subroutine ed_get_all_gimp_real_lattice_2(Greal,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(2,Nsites,Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Greal
  integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
  do ilat=1,Nsites
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 io = iorb + (ispin-1)*Norb
                 jo = jorb + (jspin-1)*Norb
                 Greal(1,ilat,io,jo,:) = Grealii(ilat,ispin,jspin,iorb,jorb,:)
                 Greal(2,ilat,io,jo,:) = Frealii(ilat,ispin,jspin,iorb,jorb,:)
              enddo
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_all_gimp_real_lattice_2

subroutine ed_get_all_gimp_real_lattice_3(Greal,Nsites,ispin,jspin,iorb,jorb)
  integer                                            :: Nsites
  complex(8),dimension(2,Nsites,Lreal),intent(inout) :: Greal
  integer                                            :: iorb,jorb,ispin,jspin
  Greal(1,1:Nsites,:) = Grealii(1:Nsites,ispin,jspin,iorb,jorb,:)
  Greal(2,1:Nsites,:) = Frealii(1:Nsites,ispin,jspin,iorb,jorb,:)
end subroutine ed_get_all_gimp_real_lattice_3







subroutine ed_get_all_gimp_real_lattice_11(Greal,ilat)
  integer                                                         :: ilat
  complex(8),dimension(2,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
  Greal(1,:,:,:,:,:) = Grealii(ilat,:,:,:,:,:)
  Greal(2,:,:,:,:,:) = Frealii(ilat,:,:,:,:,:)
end subroutine ed_get_all_gimp_real_lattice_11

subroutine ed_get_all_gimp_real_lattice_21(Greal,ilat)
  integer                                                         :: ilat
  complex(8),dimension(2,Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Greal
  integer                                                         :: io,jo,iorb,jorb,ispin,jspin
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              io = iorb + (ispin-1)*Norb
              jo = jorb + (jspin-1)*Norb
              Greal(1,io,jo,:) = Grealii(ilat,ispin,jspin,iorb,jorb,:)
              Greal(2,io,jo,:) = Frealii(ilat,ispin,jspin,iorb,jorb,:)
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_all_gimp_real_lattice_21

subroutine ed_get_all_gimp_real_lattice_31(Greal,ilat,ispin,jspin,iorb,jorb)
  integer                                     :: ilat
  complex(8),dimension(2,Lreal),intent(inout) :: Greal
  integer                                     :: iorb,jorb,ispin,jspin
  Greal(1,:) = Grealii(ilat,ispin,jspin,iorb,jorb,:)
  Greal(2,:) = Frealii(ilat,ispin,jspin,iorb,jorb,:)
end subroutine ed_get_all_gimp_real_lattice_31



