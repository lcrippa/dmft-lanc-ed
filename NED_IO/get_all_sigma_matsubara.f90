subroutine ed_get_all_sigma_matsubara_1(Smats)
  complex(8),dimension(2,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Smats
  Smats(1,:,:,:,:,:) = impSmats(:,:,:,:,:)
  Smats(2,:,:,:,:,:) = impSAmats(:,:,:,:,:)
end subroutine ed_get_all_sigma_matsubara_1

subroutine ed_get_all_sigma_matsubara_2(Smats)
  complex(8),dimension(2,Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Smats
  integer  :: io,jo,iorb,jorb,ispin,jspin
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              io = iorb + (ispin-1)*Norb
              jo = jorb + (jspin-1)*Norb
              Smats(1,io,jo,:) = impSmats(ispin,jspin,iorb,jorb,:)
              Smats(2,io,jo,:) = impSAmats(ispin,jspin,iorb,jorb,:)
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_all_sigma_matsubara_2

subroutine ed_get_all_sigma_matsubara_3(Smats,ispin,jspin,iorb,jorb)
  complex(8),dimension(2,Lmats),intent(inout) :: Smats
  integer                          :: iorb,jorb,ispin,jspin
  Smats(1,:) = impSmats(ispin,jspin,iorb,jorb,:)
  Smats(2,:) = impSAmats(ispin,jspin,iorb,jorb,:)
end subroutine ed_get_all_sigma_matsubara_3









subroutine ed_get_all_sigma_matsubara_lattice_1(Smats,Nsites)
  integer                                                                  :: Nsites
  complex(8),dimension(2,Nsites,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Smats
  Smats(1,1:Nsites,:,:,:,:,:) = Smatsii(1:Nsites,:,:,:,:,:)
  Smats(2,1:Nsites,:,:,:,:,:) = SAmatsii(1:Nsites,:,:,:,:,:)
end subroutine ed_get_all_sigma_matsubara_lattice_1

subroutine ed_get_all_sigma_matsubara_lattice_2(Smats,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(2,Nsites,Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Smats
  integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
  do ilat=1,Nsites
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 io = iorb + (ispin-1)*Norb
                 jo = jorb + (jspin-1)*Norb
                 Smats(1,ilat,io,jo,:) = Smatsii(ilat,ispin,jspin,iorb,jorb,:)
                 Smats(2,ilat,io,jo,:) = SAmatsii(ilat,ispin,jspin,iorb,jorb,:)
              enddo
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_all_sigma_matsubara_lattice_2

subroutine ed_get_all_sigma_matsubara_lattice_3(Smats,Nsites,ispin,jspin,iorb,jorb)
  integer                                            :: Nsites
  complex(8),dimension(2,Nsites,Lmats),intent(inout) :: Smats
  integer                                            :: iorb,jorb,ispin,jspin
  Smats(1,1:Nsites,:) = Smatsii(1:Nsites,ispin,jspin,iorb,jorb,:)
  Smats(2,1:Nsites,:) = SAmatsii(1:Nsites,ispin,jspin,iorb,jorb,:)
end subroutine ed_get_all_sigma_matsubara_lattice_3







subroutine ed_get_all_sigma_matsubara_lattice_11(Smats,ilat)
  integer                                                         :: ilat
  complex(8),dimension(2,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Smats
  Smats(1,:,:,:,:,:) = Smatsii(ilat,:,:,:,:,:)
  Smats(2,:,:,:,:,:) = SAmatsii(ilat,:,:,:,:,:)
end subroutine ed_get_all_sigma_matsubara_lattice_11

subroutine ed_get_all_sigma_matsubara_lattice_21(Smats,ilat)
  integer                                                         :: ilat
  complex(8),dimension(2,Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Smats
  integer                                                         :: io,jo,iorb,jorb,ispin,jspin
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              io = iorb + (ispin-1)*Norb
              jo = jorb + (jspin-1)*Norb
              Smats(1,io,jo,:) = Smatsii(ilat,ispin,jspin,iorb,jorb,:)
              Smats(2,io,jo,:) = SAmatsii(ilat,ispin,jspin,iorb,jorb,:)
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_all_sigma_matsubara_lattice_21

subroutine ed_get_all_sigma_matsubara_lattice_31(Smats,ilat,ispin,jspin,iorb,jorb)
  integer                                     :: ilat
  complex(8),dimension(2,Lmats),intent(inout) :: Smats
  integer                                     :: iorb,jorb,ispin,jspin
  Smats(1,:) = Smatsii(ilat,ispin,jspin,iorb,jorb,:)
  Smats(2,:) = SAmatsii(ilat,ispin,jspin,iorb,jorb,:)
end subroutine ed_get_all_sigma_matsubara_lattice_31



