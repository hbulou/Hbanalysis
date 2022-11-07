module Molecule

  !use Atom
  implicit none
contains
#define msg(x) print *,"### ",x,__FILE__," > "
#define error(x) print *,"### ERROR ",x,__FILE__," > " 
  !  subroutine MOLECULE_add_atom(mol,elt,q)
  !  subroutine MOLECULE_add_to_cell(cell)
  !  subroutine MOLECULE_align(molref,center,axe,mol_align,pos_to_center,axe_to_align)
  !  function   MOLECULE_copy(src) result(dest)
  !  subroutine MOLECULE_del(mol)
  !  subroutine MOLECULE_find_molecules(cell,PBC)
  !  subroutine MOLECULE_element_list_update(mol)
  !  function   MOLECULE_Mass_Center(molecule) result(MC)
  !  function   MOLECULE_init() result(mol)
  !  subroutine MOLECULE_rotate(molecule,alpha,beta,gamma,center_tag)
  !  function   MOLECULE_sort_z(molecule) result(tab)
  !  subroutine MOLECULE_translate_mass_center_at(molecule,vec)



  ! --------------------------------------------------------------------------------------
  !
  !              MOLECULE_add_atom()
  !
  ! --------------------------------------------------------------------------------------
  subroutine MOLECULE_add_atom(mol,atm)
    use global
    use Element
    use Atom
    implicit none
    !double precision::q(3)
    !character(len=2)::elt
    type(t_Molecule)::mol,TMPmol
    integer::i
    integer,allocatable::TMPlist(:)
    logical::newelt
    type(t_Atom)::atm
    logical::equal

    if(mol%n_atoms.eq.0) then
       !deallocate(mol%atoms)
       allocate(mol%atoms(mol%n_atoms+1))
       mol%atoms(1)=ATOM_new()
       mol%atoms(1)%elt=atm%elt
       mol%atoms(1)%Zato=elt2Z(atm%elt)
       mol%atoms(1)%q=atm%q
       mol%atoms(1)%idx=atm%idx
       !--------------------------------------------------------------------
       ! màj des listes d'atomes contraints dans la molécules
       mol%atoms(1)%constraints=atm%constraints
       equal=atm%constraints(1)==0.and.atm%constraints(2)==0.and.atm%constraints(3)==0
       if(equal) then
          mol%constraints%n_XYZ=1
          deallocate(mol%constraints%XYZ)
          allocate(mol%constraints%XYZ(mol%constraints%n_XYZ))
          mol%constraints%XYZ(mol%constraints%n_XYZ)=1
       end if
       !--------------------------------------------------------------------
       ! màj des types d'éléments dans la molécule
       mol%n_elements=1
       deallocate(mol%elements)
       allocate(mol%elements(1))
       mol%elements(1)%Zato=mol%atoms(1)%Zato
       mol%elements(1)%elt=mol%atoms(1)%elt
       mol%elements(1)%n=1
       !--------------------------------------------------------------------
    else
       !allocate(TMPmol%atoms(mol%n_atoms))
       !print *,"DEBUG"
       TMPmol=MOLECULE_copy(mol)
       deallocate(mol%atoms)
       allocate(mol%atoms(mol%n_atoms+1))
       do i=1,TMPmol%n_atoms
          mol%atoms(i)=ATOM_copy(TMPmol%atoms(i))
       end do
       mol%atoms(mol%n_atoms+1)%elt=atm%elt
       mol%atoms(mol%n_atoms+1)%Zato=elt2Z(mol%atoms(mol%n_atoms+1)%elt)
       mol%atoms(mol%n_atoms+1)%q=atm%q
       mol%atoms(mol%n_atoms+1)%idx=atm%idx
       mol%atoms(mol%n_atoms+1)%n_bonds=0
       mol%atoms(mol%n_atoms+1)%n_angles=0
       !--------------------------------------------------------------------
       ! màj des listes d'atomes contraints dans la molécules
       mol%atoms(mol%n_atoms+1)%constraints=atm%constraints
       equal=atm%constraints(1)==0.and.atm%constraints(2)==0.and.atm%constraints(3)==0
       if(equal) then
          allocate(TMPlist(mol%constraints%n_XYZ))
          TMPlist=mol%constraints%XYZ
          deallocate(mol%constraints%XYZ)
          mol%constraints%n_XYZ=mol%constraints%n_XYZ+1
          allocate(mol%constraints%XYZ(mol%constraints%n_XYZ))
          
          mol%constraints%XYZ(1:mol%constraints%n_XYZ-1)=TMPlist
          mol%constraints%XYZ(mol%constraints%n_XYZ)=mol%n_atoms+1
          deallocate(TMPlist)
       end if
       !--------------------------------------------------------------------
       ! màj des types d'éléments dans la molécule
       newelt=.True.
       i=1
       do while(newelt.and.i<=mol%n_elements)
          if(mol%elements(i)%Zato.eq.mol%atoms(mol%n_atoms+1)%Zato) then
             newelt=.False.
          end if
          i=i+1
       end do
       if(newelt) then
          !print *,"known elt ",(mol%elements(i)%elt,i=1,mol%n_elements)
          !print *,"New elt ",elt
          deallocate(mol%elements)
          allocate(mol%elements(mol%n_elements+1))
          do i=1,TMPmol%n_elements
             mol%elements(i)=ELEMENT_copy(TMPmol%elements(i))
          end do
          mol%elements(mol%n_elements+1)%Zato=mol%atoms(mol%n_atoms+1)%Zato
          mol%elements(mol%n_elements+1)%elt=mol%atoms(mol%n_atoms+1)%elt
          mol%elements(mol%n_elements+1)%n=1
          mol%n_elements=mol%n_elements+1
       else
          mol%elements(i-1)%n=mol%elements(i-1)%n+1
       end if
       !--------------------------------------------------------------------
    end if
    mol%n_atoms=mol%n_atoms+1
    !msg(__LINE__)," n_atoms=",mol%n_atoms," n_XYZ=",mol%constraints%n_XYZ
    !msg(__LINE__)," atom idx=",mol%atoms(-1)%idx
  end subroutine MOLECULE_add_atom
  ! --------------------------------------------------------------------------------------
  !
  !              MOLECULE_add_to_cell()
  !
  ! --------------------------------------------------------------------------------------
  subroutine MOLECULE_add_to_cell(cell)
    use global
    implicit none
    type(t_Cell)::cell
    
  end subroutine MOLECULE_add_to_cell
  ! --------------------------------------------------------------------------------------
  !
  !              MOLECULE_align()
  !
  ! --------------------------------------------------------------------------------------
  subroutine MOLECULE_align(center,axe,mol_align,pos_to_center,axe_to_align,angle)
      use global
      use tools
      implicit none
      type(t_Molecule)::mol_align
      integer::k,iat
      double precision::center(3),axe(3),pos_to_center(3),axe_to_align(3)
      double precision::alpha,t(3)
      double precision::ex(3),ey(3),ez(3),angle
      double precision::M_PI
      M_PI=4*atan(1.d0)

      msg(__LINE__), "Center @",center
      msg(__LINE__), pos_to_center

      !
      ! then we build the axis for rotating the molecule
      !
      call TOOLS_build_rotation_axis(center,axe,pos_to_center,axe_to_align,ex,ey,ez,alpha)


      !
      ! translation de la molécule : point iH1mol->iH1sub
      !
      do k=1,3
         t(k)=center(k)-pos_to_center(k)
      end do
      call MOLECULE_translate(mol_align,t) 

      
      !
      ! rotation of the molecule so that to align the two axes
      !
      alpha=90-alpha
      do iat=1,mol_align%n_atoms
         mol_align%atoms(iat)%q=TOOLS_rotate_around_ez(center,ex,ey,ez,alpha,mol_align%atoms(iat)%q)
      end do
      !
      ! rotation aroung the common axe
      !

      do iat=1,mol_align%n_atoms
         mol_align%atoms(iat)%q=TOOLS_rotate_around_ez(center,ez,ex,ey,angle,mol_align%atoms(iat)%q)
      end do

      
      
    end subroutine MOLECULE_align

  ! --------------------------------------------------------------------------------------
  !
  !              MOLECULE_copy()
  !
  ! --------------------------------------------------------------------------------------
  function MOLECULE_copy(src) result(dest)
    use global
    use Atom
    implicit none
    type(t_Molecule)::src
    type(t_Molecule)::dest
    integer::iat,k
    
    !msg(__LINE__)
    dest%save=src%save
    dest%idx=src%idx
    dest%n_bonds=src%n_bonds
    dest%idx_cell=src%idx_cell
    dest%n_atoms=src%n_atoms
    dest%MassCenter=src%MassCenter

    allocate(dest%atoms(dest%n_atoms))
    allocate(dest%order_idx_save(dest%n_atoms))
    !msg(__LINE__), "dest%n_atoms=",dest%n_atoms,    allocated(dest%order_idx_save)

    do iat=1,dest%n_atoms
       dest%atoms(iat)=ATOM_copy(src%atoms(iat))
       !dest%order_idx_save(iat)=src%order_idx_save(iat)
       !print *,src%atoms(iat)%Zato
    end do

    dest%n_elements=src%n_elements
    allocate(dest%elements(dest%n_elements))
    do iat=1,dest%n_elements
       dest%elements(iat)%Zato=src%elements(iat)%Zato
       dest%elements(iat)%elt=src%elements(iat)%elt
       dest%elements(iat)%n=src%elements(iat)%n
    end do
    dest%constraints%n_XYZ=src%constraints%n_XYZ


    allocate(dest%constraints%XYZ(dest%constraints%n_XYZ))
    !print *,"MOLECULE",allocated(src%constraints%XYZ)
    if(dest%constraints%n_XYZ.gt.0)    dest%constraints%XYZ=src%constraints%XYZ
    !msg(__LINE__), "END"

  end function MOLECULE_copy
  ! --------------------------------------------------------------------------------------
  !
  !              MOLECULE_del()
  !
  ! --------------------------------------------------------------------------------------
  subroutine MOLECULE_del(mol)
    use global
    use Atom
    implicit none
    type(t_Molecule)::mol
    integer::iat
    do iat=1,mol%n_atoms
       call ATOM_del(mol%atoms(iat))
    end do
    if(allocated(mol%atoms)) deallocate(mol%atoms)
    if(allocated(mol%elements)) deallocate(mol%elements)
    if(allocated(mol%order_idx_save)) deallocate(mol%order_idx_save)
    if(allocated(mol%constraints%XYZ)) deallocate(mol%constraints%XYZ)
  end subroutine MOLECULE_del
  ! --------------------------------------------------------------------------------------
  !
  !              MOLECULE_element_list_update(mol)
  !
  ! --------------------------------------------------------------------------------------
  subroutine MOLECULE_element_list_update(mol)
    use global
    use Element
    implicit none
    logical::newelt
    integer::i,iat
    type(t_Molecule)::mol,TMPmol
    msg(__LINE__),"mol%n_elements=" , mol%n_elements
    !
    ! on passe en revue tous les atomes de la molecule
    !

    do iat=1,mol%n_atoms
       !
       ! si la liste des elements est nulle ...
       !
       if(mol%n_elements.eq.0) then
          !msg(__LINE__), "1"
          mol%n_elements=1
          if(allocated(mol%elements))   deallocate(mol%elements)
          allocate(mol%elements(1))
          mol%elements(1)%Zato=mol%atoms(iat)%Zato
          mol%elements(1)%elt=mol%atoms(iat)%elt
          mol%elements(1)%n=1
       else
          !msg(__LINE__), "2"
          !
          ! si la liste des elements n'est pas nulle,
          ! on regarde si l'element en cours est deja
          ! dans la liste ...
          !
          !msg(__LINE__),allocated(mol%elements),mol%n_elements
          newelt=.True.
          i=1
          do while(newelt.and.i<=mol%n_elements)
             if(mol%elements(i)%Zato.eq.mol%atoms(iat)%Zato) then
                newelt=.False.
             end if
             i=i+1
          end do

       
          !
          ! lorsqu'il s'agit d'un element qui n'est pas encore
          ! dans la liste ...
          !
          if(newelt) then
             !print *,"known elt ",(mol%elements(i)%elt,i=1,mol%n_elements)
             !print *,"New elt ",elt
             !print *,"DEBUG",iat,mol%n_atoms,newelt
             !msg(__LINE__), "newelt"
             TMPmol=MOLECULE_copy(mol)
             !msg(__LINE__), "newelt2"
             if(allocated(mol%elements))   deallocate(mol%elements)
             allocate(mol%elements(mol%n_elements+1))
             do i=1,TMPmol%n_elements
                mol%elements(i)=ELEMENT_copy(TMPmol%elements(i))
             end do
             mol%elements(mol%n_elements+1)%Zato=mol%atoms(iat)%Zato
             mol%elements(mol%n_elements+1)%elt=mol%atoms(iat)%elt
             mol%elements(mol%n_elements+1)%n=1
             mol%n_elements=mol%n_elements+1
             call MOLECULE_del(TMPmol)
          else
             mol%elements(i-1)%n=mol%elements(i-1)%n+1
          end if
       end if
    end do
    msg(__LINE__),"END of MOLECULE_element_list_update"
  end subroutine MOLECULE_element_list_update


  ! --------------------------------------------------------------------------------------
  !
  !              MOLECULE_find_molecules()
  !
  ! --------------------------------------------------------------------------------------
  subroutine MOLECULE_find_molecules(cell,PBC)
    !use Element
    use global
    implicit none
    type(t_Molecule),allocatable::new_molecule(:)
    type(t_Cell)::cell
    integer::i
    logical::PBC
    integer::Zi,Zj,iat,jat,k,idx_molecule,imol
    double precision::d,R(3)
    integer,allocatable::n_atoms(:)    

    ! initialement on a une cell avec une seule molécule
    !
    ! Dans un premier temps, il faut calculer les distances
    ! entre tous les atomes de la molécule initiale
    
    idx_molecule=0
    do iat=1,cell%molecules(1)%n_atoms-1
       Zi=cell%molecules(1)%atoms(iat)%Zato
       do jat=iat+1,cell%molecules(1)%n_atoms
          Zj=cell%molecules(1)%atoms(jat)%Zato
          d=0.0
          do k=1,3
             R(k)=cell%molecules(1)%atoms(jat)%q(k)-cell%molecules(1)%atoms(iat)%q(k)
             if(PBC) then
                R(k)=modulo(R(k),cell%L(k,k))
                if(R(k).gt.0.5*cell%L(k,k)) then
                   R(k)=R(k)-cell%L(k,k)
                end if
                if(-R(k).gt.0.5*cell%L(k,k)) then
                   R(k)=R(k)+cell%L(k,k)
                end if
             end if
             d=d+R(k)*R(k)
          end do
          d=sqrt(d)
          !
          ! Si la distance calculée est inférieure à la distance d'équilibre de la liaison
          ! covalente aux deux espèces chimiques en présence, les deux atomes appartiennent
          ! à la même molécule. -> molecule idx_molecule ;
          ! A noter que tous les atomes de la molécule initiale ont un idx_molécule=-1.
          ! Au fur et à mesure que l'on identifie des molécules dans la molécule initiale,
          ! les atomes vont avoir des idx_molécule >= 1
          !
          if(d<ELEMENT_deq(Zi,Zj))  then
             !
             ! si aucun des deux atomes n'a encore été idx_molécule ...
             !
             if ((cell%molecules(1)%atoms(jat)%idx_molecule.eq.-1).and.(cell%molecules(1)%atoms(iat)%idx_molecule.eq.-1)) then
                idx_molecule=idx_molecule+1
             else
                !
                ! si un des deux atomes a été idx_molécule, celui qui ne l'a pas encore été prend l'idx_molecule
                ! de l'atome déjà indexé
                !
                if (cell%molecules(1)%atoms(jat)%idx_molecule*cell%molecules(1)%atoms(iat)%idx_molecule<0) then
                   idx_molecule=max(cell%molecules(1)%atoms(jat)%idx_molecule,cell%molecules(1)%atoms(iat)%idx_molecule)
                else
                   !
                   ! dernier cas de figure : les deux atomes ont déjà idx_molécule mais il appartiennent
                   ! a deux molécule différentes -> ce cas de figure peut se présenter dans le cas de grosse
                   ! molécule -> pas fonctionnelle pour l'instant, A DEVELOPPER
                   !
                   print *,"# A same atom belongs to two different molecules"
                   print *,"# TO BE DEVELOPPED"
                   stop
                end if
             end if
             cell%molecules(1)%atoms(iat)%idx_molecule=idx_molecule
             cell%molecules(1)%atoms(jat)%idx_molecule=idx_molecule

             !write(*,'(A,I5,I5,A,F12.8,A,A,I5,I5)') &
             !     "# d(",iat,jat,")=",d,molecules(1)%atoms(jat)%elt,molecules(1)%atoms(iat)%elt,&
             !     molecules(1)%atoms(iat)%idx_molecule,molecules(1)%atoms(iat)%idx_molecule
             cell%n_bonds=cell%n_bonds+1
             cell%molecules(1)%atoms(iat)%n_bonds=cell%molecules(1)%atoms(iat)%n_bonds+1
             cell%molecules(1)%atoms(jat)%n_bonds=cell%molecules(1)%atoms(jat)%n_bonds+1




             
             cell%molecules(1)%n_bonds=cell%molecules(1)%n_bonds+1
          end if
          
       end do
    end do


    !
    ! ensuite on reconstruit une nouvelle cell avec le bon 
    ! nombre de molécules
    !
    
    allocate(n_atoms(idx_molecule))
    allocate(new_molecule(idx_molecule))
    
    do imol=1,idx_molecule
       new_molecule(imol)%n_atoms=0
       n_atoms(imol)=0
    end do
    
    do iat=1,cell%molecules(1)%n_atoms
       imol=cell%molecules(1)%atoms(iat)%idx_molecule
       new_molecule(imol)%n_atoms=new_molecule(imol)%n_atoms+1
    end do
    do imol=1,idx_molecule
       allocate(new_molecule(imol)%atoms(new_molecule(imol)%n_atoms))
    end do


    do iat=1,cell%molecules(1)%n_atoms

       imol=cell%molecules(1)%atoms(iat)%idx_molecule

       n_atoms(imol)=n_atoms(imol)+1
       !print *,iat,imol,n_atoms(imol)
       new_molecule(imol)%atoms(n_atoms(imol))%elt=cell%molecules(1)%atoms(iat)%elt
       new_molecule(imol)%atoms(n_atoms(imol))%Zato=cell%molecules(1)%atoms(iat)%Zato
       new_molecule(imol)%atoms(n_atoms(imol))%n_bonds=cell%molecules(1)%atoms(iat)%n_bonds
       do k=1,3
          new_molecule(imol)%atoms(n_atoms(imol))%q(k)=cell%molecules(1)%atoms(iat)%q(k)
       end do
       new_molecule(imol)%atoms(n_atoms(imol))%idx_molecule=imol

     end do
    print *,"# MOLECULE_find_molecules>"
     
    !do imol=1,idx_molecule
    !   new_mol(imol)%n_atoms=0
    !end do
     deallocate(n_atoms)
     !molecule=>new_molecule
     deallocate(cell%molecules)
     allocate(cell%molecules(idx_molecule))
    do imol=1,idx_molecule
       cell%molecules(imol)%n_atoms=new_molecule(imol)%n_atoms
       allocate(cell%molecules(imol)%atoms(cell%molecules(imol)%n_atoms))
       do iat=1,cell%molecules(imol)%n_atoms
          cell%molecules(imol)%atoms(iat)%elt=new_molecule(imol)%atoms(iat)%elt
          cell%molecules(imol)%atoms(iat)%Zato=new_molecule(imol)%atoms(iat)%Zato
          cell%molecules(imol)%atoms(iat)%n_bonds=new_molecule(imol)%atoms(iat)%n_bonds
          do k=1,3
             cell%molecules(imol)%atoms(iat)%q(k)=new_molecule(imol)%atoms(iat)%q(k)
          end do
          cell%molecules(imol)%atoms(iat)%idx_molecule=imol
          print *,"# MOLECULE_find_molecules> molecule ",imol," n_bonds= ",cell%molecules(imol)%atoms(iat)%n_bonds
       end do
    end do
    deallocate(new_molecule)

    cell%n_molecules=idx_molecule
    print *,"# MOLECULE_find_molecules> ",cell%n_molecules," molecule(s) found"
    do i=1,cell%n_molecules
       print *,"# MOLECULE_find_molecules> molecule ",i," -> ",cell%molecules(i)%n_bonds 
    end do
    print *,"# MOLECULE_find_molecules> ",cell%n_bonds," bonds(s) found"

  end subroutine MOLECULE_find_molecules
  ! --------------------------------------------------------------------------------------
  !
  !              MOLECULE_from_file()
  !
  ! --------------------------------------------------------------------------------------
  ! function MOLECULE_from_file(input) result(molecule)
  !   use Element
  !   use global
  !   use ATOM
  !   implicit none
  !   type(t_Input)::input
  !   type(t_Molecule)::molecule
  !   type(t_Atom)::atm
  !   character (len=1024)::line
  !   character (len=NCHARFIELD)::field(32)
  !   integer::nfield,n_atoms
  !   integer::iconf,iat,io
  !   logical:: file_exists

  !   msg(__LINE__), " Reading ",input%name

  !   INQUIRE(FILE=input%name, EXIST=file_exists)
  !   if(file_exists) then
  !      open(unit=1,file=input%name,form='formatted')
  !   else
  !      error(__LINE__), " Error in MOLECULE_from_file():", input%name," doesn't exist" ; stop
  !   end if
  !   iconf=1
  !   do while(.not.(iconf.eq.input%iconf))
  !      read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
  !      read(field(1),*) n_atoms
  !      read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field)
  !      do iat=1,n_atoms
  !         read(1,'(A)',iostat=io) line
  !      end do
  !      iconf=iconf+1
  !   end do

  !   read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
  !   read(field(1),*) n_atoms
    
  !   molecule=MOLECULE_init()
  !   atm=ATOM_new()
  !   read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field)
  !   do iat=1,n_atoms
  !      read(1,'(A)',iostat=io) line
  !         call line_parser(line,nfield,field)
  !         atm%elt=field(1)
  !         atm%Zato=elt2Z(atm%elt)
  !         read(field(2),*) atm%q(1)
  !         read(field(3),*) atm%q(2)
  !         read(field(4),*) atm%q(3)
  !         atm%idx_molecule=-1
  !         atm%idx=iat
  !         atm%n_bonds=0
  !         call MOLECULE_add_atom(molecule,atm)
  !   end do
  !   close(1)

  !   call MOLECULE_update(molecule)
    
  !   msg(__LINE__)," n_atoms= ",molecule%n_atoms
  ! end function MOLECULE_from_file
  ! --------------------------------------------------------------------------------------
  !
  !             MOLECULE_analysis()
  !
  ! --------------------------------------------------------------------------------------
  subroutine MOLECULE_analysis(molecule)
    use global
    implicit none
    type(t_Molecule)::molecule
    integer::iat,i
    do i=1,3
       molecule%limits(i,1)=molecule%atoms(1)%q(i)
       molecule%limits(i,2)=molecule%atoms(1)%q(i)
    end do
    do iat=2,molecule%n_atoms
       if(molecule%atoms(iat)%q(i)<molecule%limits(i,1)) molecule%limits(i,1)=molecule%atoms(iat)%q(i)
       if(molecule%atoms(iat)%q(i)>molecule%limits(i,2)) molecule%limits(i,2)=molecule%atoms(iat)%q(i)
    end do
    write(*,*) "# Limits of the molecule:"
    do i=1,3
       write(*,*) molecule%limits(i,1),molecule%limits(i,2)
    end do
  end subroutine MOLECULE_analysis
  ! --------------------------------------------------------------------------------------
  !
  !             MOLECULE_update()
  !
  ! --------------------------------------------------------------------------------------
  subroutine MOLECULE_update(molecule)
    use global
    use Atom
    implicit none
    type(t_Molecule)::molecule
    integer::iat,jat,Zi,Zj
    double precision::r(4)
    do iat=1,molecule%n_atoms-1
       Zi=molecule%atoms(iat)%Zato
       do jat=iat+1,molecule%n_atoms
          r=ATOM_rij(molecule%atoms(jat),molecule%atoms(iat))
          Zj=molecule%atoms(jat)%Zato
          if(r(4).le.ELEMENT_deq(Zi,Zj)) then
             print *,r(4),Zi,Zj
          end if
       end do
    end do
    
  end subroutine MOLECULE_update
  ! --------------------------------------------------------------------------------------
  !
  !              molecule_mass_center()
  !
  ! --------------------------------------------------------------------------------------
  function MOLECULE_Mass_Center(molecule) result(MC)
    use Element
    use global
    implicit none

    type(t_Molecule),intent(in)::molecule
    integer::iat,k
    double precision::m,mtot
    double precision::MC(3)
    do k=1,3
       MC(k)=0.0
    end do
    do iat=1,molecule%n_atoms
       m=PeriodicTable(molecule%atoms(iat)%Zato)%mass
       !print *,m
       mtot=mtot+m
       do k=1,3
          !print *,molecule%atoms(iat)%q(k)
          MC(k)=MC(k)+m*molecule%atoms(iat)%q(k)
       end do
    end do
    do k=1,3
       MC(k)=MC(k)/mtot
    end do
  end function MOLECULE_Mass_Center
  ! --------------------------------------------------------------------------------------
  !
  !              MOLECULE_init()
  !
  ! --------------------------------------------------------------------------------------
  function MOLECULE_init() result(mol)
    use global
    implicit none
    type(t_Molecule)::mol
    mol%save=.True. 
    mol%idx=-1
    mol%idx_cell=-1
    mol%n_atoms=0
    !allocate(mol%atoms(0))
    mol%n_bonds=0
    mol%n_angles=0
    mol%n_elements=0
    allocate(mol%elements(0))
    mol%constraints%n_XYZ=0
    allocate(mol%constraints%XYZ(mol%constraints%n_XYZ))
  end function MOLECULE_init
  ! --------------------------------------------------------------------------------------
  !
  !              molecule_rotate()
  !
  ! --------------------------------------------------------------------------------------
  subroutine MOLECULE_rotate(molecule,alpha,beta,gamma,center_tag)
    use global
    implicit none
    type(t_Molecule)::molecule
    double precision::alpha,beta,gamma
    double precision::pi,center(3),newq(3),vec(3)
    character(len=32)::center_tag
    integer::iat

    pi=4.D0*DATAN(1.D0)
    alpha=alpha*pi/180
    beta=beta*pi/180
    gamma=gamma*pi/180

    if(center_tag.eq.'molecule_MC') then
       center=MOLECULE_Mass_Center(molecule)
       print *,'center=',center
    end if

    do iat=1,molecule%n_atoms
       newq(1)=center(1)+(molecule%atoms(iat)%q(1)-center(1))*cos(gamma)-(molecule%atoms(iat)%q(2)-center(2))*sin(gamma)
       newq(2)=center(2)+(molecule%atoms(iat)%q(1)-center(1))*sin(gamma)+(molecule%atoms(iat)%q(2)-center(2))*cos(gamma)
       molecule%atoms(iat)%q(1)=newq(1)
       molecule%atoms(iat)%q(2)=newq(2)
       newq(1)=center(1)+(molecule%atoms(iat)%q(1)-center(1))*cos(beta)-(molecule%atoms(iat)%q(3)-center(3))*sin(beta)
       newq(3)=center(3)+(molecule%atoms(iat)%q(1)-center(1))*sin(beta)+(molecule%atoms(iat)%q(3)-center(3))*cos(beta)
       molecule%atoms(iat)%q(1)=newq(1)
       molecule%atoms(iat)%q(3)=newq(3)
       newq(2)=center(2)+(molecule%atoms(iat)%q(2)-center(2))*cos(alpha)-(molecule%atoms(iat)%q(3)-center(3))*sin(alpha)
       newq(3)=center(3)+(molecule%atoms(iat)%q(2)-center(2))*sin(alpha)+(molecule%atoms(iat)%q(3)-center(3))*cos(alpha)
       molecule%atoms(iat)%q(2)=newq(2)
       molecule%atoms(iat)%q(3)=newq(3)
    end do

  end subroutine MOLECULE_rotate


  ! --------------------------------------------------------------------------------------
  !
  !              MOLECULE_sort_z()
  !
  ! --------------------------------------------------------------------------------------
  subroutine MOLECULE_sort_z(mol)
    use global
    !use Element
    implicit none
    type(t_Molecule)::mol
    integer::jat,idx1,idx2,nop

    if(allocated(mol%order_idx_save)) then
       deallocate(mol%order_idx_save)
    end if
    allocate(mol%order_idx_save(mol%n_atoms))
    do jat=1,mol%n_atoms
       mol%order_idx_save(jat)=jat
    end do
    nop=1
    do while(nop.gt.0)
       nop=0
       do jat=1,mol%n_atoms-1
          idx1=mol%order_idx_save(jat)
          idx2=mol%order_idx_save(jat+1)
          if((idx2>idx1).and.(ABS(AINT(mol%atoms(idx1)%q(3))).gt.ABS(AINT(mol%atoms(idx2)%q(3))))) then
             mol%order_idx_save(jat)=idx2
             mol%order_idx_save(jat+1)=idx1
             nop=nop+1
          end if
       end do
    end do
  end subroutine MOLECULE_sort_z
  ! --------------------------------------------------------------------------------------
  !
  !              molecule_mass_center()
  !
  ! --------------------------------------------------------------------------------------
  subroutine MOLECULE_translate_mass_center_at(molecule,vec)
    use global
    use Element
    implicit none
    double precision::vec(3),MC(3)
    type(t_Molecule)::molecule
    integer::iat,k
    MC=MOLECULE_Mass_Center(molecule)
    do iat=1,molecule%n_atoms
       do k=1,3
          molecule%atoms(iat)%q(k)=molecule%atoms(iat)%q(k)-MC(k)+vec(k)
       end do
    end do
  end subroutine MOLECULE_translate_mass_center_at
  ! --------------------------------------------------------------------------------------
  !
  !              molecule_translate()
  !
  ! --------------------------------------------------------------------------------------
  subroutine MOLECULE_translate(molecule,vec)
    use global
    use Element
    implicit none
    double precision::vec(3)
    type(t_Molecule)::molecule
    integer::iat,k
    do iat=1,molecule%n_atoms
       do k=1,3
          molecule%atoms(iat)%q(k)=molecule%atoms(iat)%q(k)+vec(k)
       end do
    end do
  end subroutine MOLECULE_translate


  ! --------------------------------------------------------------------------------------
  !
  !              read_multixyz()
  !
  ! --------------------------------------------------------------------------------------
  function MOLECULE_read_multixyz(name) result(mol)
    use global
    use Atom
    use Element
    implicit none
    character(len=2048)::name
    integer::io
    character (len=1024)::line
    character (len=NCHARFIELD)::field(32)
    integer::nfield,nline,nconf,natoms,i,iat
    type(t_Molecule)::mol
    type(t_Atom)::atm

    
    msg(__LINE__)," Starting read_multixyz() function ..."
    msg(__LINE__)," Reading ",trim(name)
    io=0
    nline=0
    nconf=0
    open(unit=1,file=name,form='formatted')
    do while (io==0)
       read(1,'(A)',iostat=io) line
       if(io.eq.-1) then
          exit
       else
          nline=nline+1
       end if
       call line_parser(line,nfield,field)
       if(nfield.eq.1) then
          nconf=nconf+1
          read(field(1),*) natoms
       end if
    end do
    msg(__LINE__),nline," line(s) in ",trim(name)
    msg(__LINE__),nconf," configuration(s) in ",trim(name)
    msg(__LINE__),natoms," atom(s)/configuration in ",trim(name)
    if(.not.((natoms+2)*nconf.eq.nline))    error(__LINE__),(natoms+2)*nconf,nline
    rewind(1)
    io=0
    if(nconf.gt.1) then
       i=0
       do while (.not.(i.eq.((natoms+2)*(nconf-1))))
          read(1,'(A)',iostat=io) line ; i=i+1
       end do
    end if
    read(1,'(A)',iostat=io) line
    read(1,'(A)',iostat=io) line
    mol=MOLECULE_init()
    atm=ATOM_new()
    do iat=1,natoms
       read(1,'(A)',iostat=io) line
       call line_parser(line,nfield,field)
       atm%elt=field(1)
       atm%Zato=elt2Z(atm%elt)
       read(field(2),*) atm%q(1)
       read(field(3),*) atm%q(2)
       read(field(4),*) atm%q(3)
       atm%idx_molecule=-1
       atm%n_bonds=0
       call MOLECULE_add_atom(mol,atm)
    end do
       


    close(1)
    msg(__LINE__),"... end of read_multixyz()."
  end function MOLECULE_read_multixyz

  
end module Molecule
