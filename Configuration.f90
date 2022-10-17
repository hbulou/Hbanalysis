module Configuration
  implicit none
contains
  ! --------------------------------------------------------------------------------------
  !
  !              UNIVERS_copy()
  !
  ! --------------------------------------------------------------------------------------
  function CONFIGURATION_copy(src) result(dest)
    use global
    use Molecule
    implicit none
    type(t_Configuration)::src,dest
    integer::imol
    dest%cells%n_molecules=src%cells%n_molecules
    allocate(dest%cells%molecules(dest%cells%n_molecules))
    do imol=1,dest%cells%n_molecules
       dest%cells%molecules(imol)=MOLECULE_copy(src%cells%molecules(imol))
    end do
  end function CONFIGURATION_copy

  ! --------------------------------------------------------------------------------------
  !
  !              CONFIGURATION_new()
  !
  ! --------------------------------------------------------------------------------------
  function CONFIGURATION_init() result(dest)
    use global
    use Cell
    implicit none
    type(t_Configuration)::dest
    dest%cells=Cell_init()
  end function CONFIGURATION_init


end module Configuration
