module Atom

  implicit none
contains
#define msg(x) print *,"### ",x,__FILE__," > "
#define debug(x,idx) print *,"### DEBUG ",idx," @ line ",x,__FILE__," > "
  !   function ATOM_copy(src) result(dest)
  !   function ATOM_rij(atmi,atmj) result(r)
  !   function ATOM_new() result(dest)
  !   subroutine ATOM_del(atom)
  !   subroutine ATOM_translate(atm,t)
  ! --------------------------------------------------------------------------------------
  !
  !              ATOM_copy()
  !
  ! --------------------------------------------------------------------------------------
  function ATOM_copy(src) result(dest)
    use global
    implicit none
    type(t_Atom)::src,dest
    integer::i
    !msg(__LINE__),"dest%idx=",dest%idx
    dest%idx=src%idx
    dest%save=src%save
    !msg(__LINE__),"HHHHH"
    dest%idx_molecule=src%idx_molecule
    dest%idx_cell=src%idx_cell
    dest%constraints=src%constraints
    !read(src%elt,*) dest%elt
    dest%elt=src%elt
    !print *,dest%elt
    !debug(__LINE__,src%n_bonds)
    dest%n_bonds=src%n_bonds
    allocate(dest%bonds(dest%n_bonds))
    if(src%n_bonds.gt.0) then
       dest%bonds=src%bonds
    end if

    
    !msg(__LINE__), allocated(src%angles),"src%n_angles=",src%n_angles
    dest%n_angles=src%n_angles
    !msg(__LINE__), allocated(dest%angles),"des%n_angles=",dest%n_angles
    !do i=1,src%n_angles
    !   msg(__LINE__), i,src%angles(i)
    !end do
    if(allocated(src%angles)) then
    !if(src%n_angles.gt.0) then
       allocate(dest%angles(dest%n_angles))
       !msg(__LINE__), allocated(dest%angles),"des%n_angles=",dest%n_angles
       dest%angles=src%angles
    end if


    !msg(__LINE__),"END"
    dest%Zato=src%Zato
    dest%q=src%q
    
  end function ATOM_copy
  ! --------------------------------------------------------------------------------------
  !
  !              ATOM_rij()
  !
  ! --------------------------------------------------------------------------------------
  function ATOM_rij(atmi,atmj) result(r)
    use global
    implicit none
    type(t_Atom)::atmi,atmj
    double precision::r(4)
    integer::k
    r(4)=0.0
    do k=1,3
       r(k)=atmj%q(k)-atmi%q(k)
       r(4)=r(4)+r(k)*r(k)
    end do
    r(4)=sqrt(r(4))
  end function ATOM_rij
  ! --------------------------------------------------------------------------------------
  !
  !              ATOM_new()
  !
  ! --------------------------------------------------------------------------------------
  function ATOM_new() result(dest)
    use global
    implicit none
    type(t_Atom)::dest
    dest%save=.True.
    dest%idx=-1
    dest%idx_molecule=-1
    dest%idx_cell=-1
    dest%n_bonds=0
    dest%constraints=(/1,1,1/)
    !allocate(dest%bonds(0))
    dest%n_angles=0
    !allocate(dest%angles(0))
  end function ATOM_new
  ! --------------------------------------------------------------------------------------
  !
  !              ATOM_del()
  !
  ! --------------------------------------------------------------------------------------
  subroutine ATOM_del(atom)
    use global
    implicit none
    type(t_Atom)::atom
    if(allocated(atom%bonds)) deallocate(atom%bonds)
    if(allocated(atom%angles)) deallocate(atom%angles)
  end subroutine  ATOM_del
  ! --------------------------------------------------------------------------------------
  !
  !              ATOM_translate()
  !
  ! --------------------------------------------------------------------------------------
  subroutine ATOM_translate(atm,t)
    use global
    implicit none
    type(t_Atom)::atm
    double precision::t(3)
    integer::i
    do i=1,3
       atm%q(i)=atm%q(i)+t(i)
    end do
  end subroutine ATOM_translate
end module Atom
