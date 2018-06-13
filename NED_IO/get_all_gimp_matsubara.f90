subroutine ed_get_all_gimp_matsubara_1(Gmats)
  complex(8),dimension(2,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
  Gmats(1,:,:,:,:,:) = impGmats(:,:,:,:,:)
  Gmats(2,:,:,:,:,:) = impFmats(:,:,:,:,:)
end subroutine ed_get_all_gimp_matsubara_1

subroutine ed_get_all_gimp_matsubara_2(Gmats)
  complex(8),dimension(2,Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Gmats
  integer  :: io,jo,iorb,jorb,ispin,jspin
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              io = iorb + (ispin-1)*Norb
              jo = jorb + (jspin-1)*Norb
              Gmats(1,io,jo,:) = impGmats(ispin,jspin,iorb,jorb,:)
              Gmats(2,io,jo,:) = impFmats(ispin,jspin,iorb,jorb,:)
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_all_gimp_matsubara_2

subroutine ed_get_all_gimp_matsubara_3(Gmats,ispin,jspin,iorb,jorb)
  complex(8),dimension(2,Lmats),intent(inout) :: Gmats
  integer                          :: iorb,jorb,ispin,jspin
  Gmats(1,:) = impGmats(ispin,jspin,iorb,jorb,:)
  Gmats(2,:) = impFmats(ispin,jspin,iorb,jorb,:)
end subroutine ed_get_all_gimp_matsubara_3









subroutine ed_get_all_gimp_matsubara_lattice_1(Gmats,Nsites)
  integer                                                                  :: Nsites
  complex(8),dimension(2,Nsites,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
  Gmats(1,1:Nsites,:,:,:,:,:) = Gmatsii(1:Nsites,:,:,:,:,:)
  Gmats(2,1:Nsites,:,:,:,:,:) = Fmatsii(1:Nsites,:,:,:,:,:)
end subroutine ed_get_all_gimp_matsubara_lattice_1

subroutine ed_get_all_gimp_matsubara_lattice_2(Gmats,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(2,Nsites,Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Gmats
  integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
  do ilat=1,Nsites
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 io = iorb + (ispin-1)*Norb
                 jo = jorb + (jspin-1)*Norb
                 Gmats(1,ilat,io,jo,:) = Gmatsii(ilat,ispin,jspin,iorb,jorb,:)
                 Gmats(2,ilat,io,jo,:) = Fmatsii(ilat,ispin,jspin,iorb,jorb,:)
              enddo
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_all_gimp_matsubara_lattice_2

subroutine ed_get_all_gimp_matsubara_lattice_3(Gmats,Nsites,ispin,jspin,iorb,jorb)
  integer                                            :: Nsites
  complex(8),dimension(2,Nsites,Lmats),intent(inout) :: Gmats
  integer                                            :: iorb,jorb,ispin,jspin
  Gmats(1,1:Nsites,:) = Gmatsii(1:Nsites,ispin,jspin,iorb,jorb,:)
  Gmats(2,1:Nsites,:) = Fmatsii(1:Nsites,ispin,jspin,iorb,jorb,:)
end subroutine ed_get_all_gimp_matsubara_lattice_3







subroutine ed_get_all_gimp_matsubara_lattice_11(Gmats,ilat)
  integer                                                         :: ilat
  complex(8),dimension(2,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
  Gmats(1,:,:,:,:,:) = Gmatsii(ilat,:,:,:,:,:)
  Gmats(2,:,:,:,:,:) = Fmatsii(ilat,:,:,:,:,:)
end subroutine ed_get_all_gimp_matsubara_lattice_11

subroutine ed_get_all_gimp_matsubara_lattice_21(Gmats,ilat)
  integer                                                         :: ilat
  complex(8),dimension(2,Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Gmats
  integer                                                         :: io,jo,iorb,jorb,ispin,jspin
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              io = iorb + (ispin-1)*Norb
              jo = jorb + (jspin-1)*Norb
              Gmats(1,io,jo,:) = Gmatsii(ilat,ispin,jspin,iorb,jorb,:)
              Gmats(2,io,jo,:) = Fmatsii(ilat,ispin,jspin,iorb,jorb,:)
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_all_gimp_matsubara_lattice_21

subroutine ed_get_all_gimp_matsubara_lattice_31(Gmats,ilat,ispin,jspin,iorb,jorb)
  integer                                     :: ilat
  complex(8),dimension(2,Lmats),intent(inout) :: Gmats
  integer                                     :: iorb,jorb,ispin,jspin
  Gmats(1,:) = Gmatsii(ilat,ispin,jspin,iorb,jorb,:)
  Gmats(2,:) = Fmatsii(ilat,ispin,jspin,iorb,jorb,:)
end subroutine ed_get_all_gimp_matsubara_lattice_31



