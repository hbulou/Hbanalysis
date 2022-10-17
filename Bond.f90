module Bond
  implicit none
contains

  subroutine BOND_copy(src,dest)
    use global
    implicit none
    type(t_Bond)::src,dest
    dest%idx=src%idx
    dest%list_atoms(1)%idx=src%list_atoms(1)%idx
    dest%list_atoms(2)%idx=src%list_atoms(2)%idx
    dest%list_atoms(1)%idx_molecule=src%list_atoms(1)%idx_molecule
    dest%list_atoms(2)%idx_molecule=src%list_atoms(2)%idx_molecule
  end subroutine BOND_copy
    
    
end module Bond
