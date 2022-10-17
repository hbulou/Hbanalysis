module Cell
  use Molecule
  implicit none
contains
#define msg(x) print *,"### ",x,__FILE__," > "
!#define msg(x) print *,"### ",x,' ###    ',__FILE__," > "
#define error(x) print *,"### ERROR ",x,__FILE__," > "
#define debug(x,idx) print *,"### DEBUG ",idx," @ line ",x,__FILE__," > "
  !  subroutine CELL_add_bond(cell,iat,imol,jat,jmol)
  !  subroutine CELL_add_bond_to_atom(atom,idx_bond)
  !  subroutine CELL_add_molecules(cell,molec,n,adding_type)
  !  function   CELL_init() result(dest)
  !  subroutine CELL_copy(src,dest)
  !  subroutine CELL_find_molecular_topology(cells,PBC)
  !  function   CELL_free_location(cell,molecule) result(free_location)
  !  subroutine CELL_rebox(cell)
  !  function   CELL_slab(ibrav,a,csura,n,unit_cell,evacc) result(cell)
  !  function   CELL_intermolecular_distance(cell,imol,iat,jmol,jat) result(d)

  function   CELL_intermolecular_distance(cell,imol,iat,jmol,jat) result(R)
    use global
    implicit none
    type(t_Cell)::cell
    integer::iat,jat,imol,jmol,k
    double precision::R(4)
    msg(__LINE__), "cell%n_molecules=", cell%n_molecules
    do k=1,cell%n_molecules
       msg(__LINE__), "cell%molecules(",k,")%n_atoms=", cell%molecules(k)%n_atoms
    end do
    R(4)=0.0
    do k=1,3
       R(k)=cell%molecules(imol)%atoms(iat)%q(k)-cell%molecules(jmol)%atoms(jat)%q(k)
       R(k)=modulo(R(k),cell%L(k,k))
       if(R(k).gt.0.5*cell%L(k,k)) then
          R(k)=R(k)-cell%L(k,k)
       end if
       if(-R(k).gt.0.5*cell%L(k,k)) then
          R(k)=R(k)+cell%L(k,k)
       end if
       R(4)=R(4)+R(k)*R(k)
    end do
    R(4)=sqrt(R(4))

  end function CELL_intermolecular_distance
  
  ! --------------------------------------------------------------------------------------
  !
  !              CELL_add_bond()
  !
  ! --------------------------------------------------------------------------------------
  subroutine CELL_add_bond(cell,iat,imol,jat,jmol)
    use global
    implicit none
    integer::iat,jat,imol,jmol,ibond
    type(t_Cell)::cell,TMPcell
    !print *,"# CELL_add_bond> ","########################### CELL_add_bond() ################################################"
    !print *,"# CELL_add_bond> iat=",iat,"imol=",imol,"-->",allocated(cell%molecules(imol)%atoms(iat)%bonds)
    cell%n_bonds=cell%n_bonds+1
    if(allocated(cell%bonds)) then
       call CELL_copy(cell,TMPcell)
       deallocate(cell%bonds)
       allocate(cell%bonds(cell%n_bonds))
       call CELL_copy(TMPcell,cell)
    else
       allocate(cell%bonds(cell%n_bonds))
    end if
    
    cell%bonds(cell%n_bonds)%idx=cell%n_bonds

    cell%bonds(cell%n_bonds)%list_atoms(1)%idx=iat
    cell%bonds(cell%n_bonds)%list_atoms(2)%idx=jat
    cell%bonds(cell%n_bonds)%list_atoms(1)%idx_molecule=imol
    cell%bonds(cell%n_bonds)%list_atoms(2)%idx_molecule=jmol

    !print *,"# CELL_add_bond> iat=",iat,"imol=",imol,"-->",allocated(cell%molecules(imol)%atoms(iat)%bonds)
    call CELL_add_bond_to_atom(cell%molecules(imol)%atoms(iat),cell%n_bonds)
    call CELL_add_bond_to_atom(cell%molecules(jmol)%atoms(jat),cell%n_bonds)
    !print *,"# CELL_add_bond> iat=",iat,"imol=",imol,"-->",allocated(cell%molecules(imol)%atoms(iat)%bonds)
    !print *,"# CELL_add_bond>  -------------------- CELL_add_bond() ----------------------------------------------"
  end subroutine CELL_add_bond
  ! --------------------------------------------------------------------------------------
  !
  !              CELL_add_bond_to_atom(
  !
  ! --------------------------------------------------------------------------------------
  subroutine CELL_add_bond_to_atom(atom,idx_bond)
      use global
      implicit none
      integer::idx_bond,ibond
      type(t_Atom)::atom
      integer,allocatable::TMPbonds(:)
      !print *,"# CELL_add_bond_to_atom> "," ---------------------------------------------------------------------------"
      atom%n_bonds=atom%n_bonds+1
      !print *,"# CELL_add_bond_to_atom> atom_idx",atom%idx," nbonds= ",atom%n_bonds,allocated(atom%bonds)      
      if(atom%n_bonds.gt.1) then
         allocate(TMPbonds(atom%n_bonds-1))
         TMPbonds=atom%bonds
         !do ibond=1,atom%n_bonds-1
         !   TMPbonds(ibond)=atom%bonds(ibond)
         !end do
         deallocate(atom%bonds)
         allocate(atom%bonds(atom%n_bonds))
         atom%bonds(1:atom%n_bonds-1)=TMPbonds
         !do ibond=1,atom%n_bonds-1
         !   atom%bonds(ibond)=TMPbonds(ibond)
         !end do
         deallocate(TMPbonds)
      else
         allocate(atom%bonds(atom%n_bonds))
      end if
      atom%bonds(atom%n_bonds)=idx_bond
      !do ibond=1,atom%n_bonds
      !   print *,"# CELL_add_bond_to_atom>", idx_bond,"atom_idx",atom%idx,"->",atom%bonds(ibond)
      !end do
      !print *,"# CELL_add_bond_to_atom> "," ---------------------------------------------------------------------------"
      !print *,"# CELL_add_bond_to_atom> atom_idx",atom%idx," nbonds= ",atom%n_bonds,allocated(atom%bonds)      
    end subroutine CELL_add_bond_to_atom
  ! --------------------------------------------------------------------------------------
  !
  !              CELL_add_molecules()
  !
  ! --------------------------------------------------------------------------------------
  subroutine CELL_add_molecules(cell,molec,n,adding_type,chk_free_location)
    use global
    use tools
    use Molecule
    implicit none
    type(t_Cell)::cell
    type(t_Cell)::TMPcell
    type(t_Molecule)::molec
    logical::free_location,chk_free_location
    integer::i,idum,iconf,k,iat
    integer::n,imol
    character(len=32)::center_tag
    character(len=32)::adding_type  ! random || as_it
    double precision::vec(3),alpha,beta,gamma,MC(3),dlim

    msg(__LINE__)

    MC=MOLECULE_Mass_Center(molec)
    idum=1  
    i=0
    select case(adding_type)
       ! ---------------------------------------------------- !
    case('random')
       ! ---------------------------------------------------- !
       msg(__LINE__)," Adding type: random"
       do while(i.lt.n)
          
          vec=(/ cell%L(1,1)*ran(idum),cell%L(2,2)*ran(idum),cell%L(3,3)*ran(idum) /)
          call MOLECULE_translate_mass_center_at(molec,vec)
          alpha=360*ran(idum)
          beta=360*ran(idum)
          gamma=360*ran(idum)
          center_tag='molecule_MC'
          call MOLECULE_rotate(molec,alpha,beta,gamma,center_tag)
          if(chk_free_location) then
             free_location=CELL_free_location(cell,molec)
          else
             free_location=.True.
          end if
          if(free_location) then
             call CELL_copy(cell,TMPcell)
             if(allocated(cell%molecules)) then
                deallocate(cell%molecules)
             end if
             cell%n_molecules=TMPcell%n_molecules+1
             allocate(cell%molecules(cell%n_molecules))
             call CELL_copy(TMPcell,cell)
             cell%n_atoms=TMPcell%n_atoms+molec%n_atoms
             cell%molecules(cell%n_molecules)=MOLECULE_copy(molec)
             call CELL_del(TMPcell)
             !deallocate(TMPcell%molecules)
             i=i+1
          end if
       end do
       ! ---------------------------------------------------- !
    case default ! 'as_it'
       ! ---------------------------------------------------- !
       msg(__LINE__),"Adding type: as_it"
       msg(__LINE__),cell%n_molecules,"molecule(s) in cell"
       msg(__LINE__),molec%n_atoms,"atoms in the molecule(s) to add"
       if(chk_free_location) then
          msg(__LINE__),"Free location checking ON"
          free_location=CELL_free_location(cell,molec)
       else
          msg(__LINE__),"Free location checking OFF"
          free_location=.True.
       end if

       msg(__LINE__), "Free location?",free_location
       if(free_location) then
          msg(__LINE__), "Copy cell->TMPcell"
          msg(__LINE__), "cell%n_atoms=", cell%n_atoms
          msg(__LINE__), "cell%n_molecules=", cell%n_molecules
          if(cell%n_molecules.gt.0) then
             call CELL_copy(cell,TMPcell)
             msg(__LINE__), "TMPcell%n_atoms", TMPcell%n_atoms
             if(allocated(cell%molecules)) then
                deallocate(cell%molecules)
             end if
             msg(__LINE__), "TMPcell%n_atoms", TMPcell%n_atoms
             !cell%n_molecules=TMPcell%n_molecules+1
             allocate(cell%molecules(cell%n_molecules+1))
             msg(__LINE__), "Copy TMPcell-->cell"
             call CELL_copy(TMPcell,cell)
             msg(__LINE__), "n_atoms in TMPCell=",TMPcell%n_atoms
             cell%n_atoms=TMPcell%n_atoms+molec%n_atoms
             cell%n_molecules=cell%n_molecules+1
             cell%molecules(cell%n_molecules)=MOLECULE_copy(molec)
             call CELL_del(TMPcell)
          else
             cell%n_molecules=1
             msg(__LINE__), "allocated(cell%molecules)=",allocated(cell%molecules)
             if(allocated(cell%molecules)) then
                deallocate(cell%molecules)
             end if
             allocate(cell%molecules(cell%n_molecules))
             msg(__LINE__),"molec%n_atoms=", molec%n_atoms
             cell%n_atoms=molec%n_atoms
             cell%molecules(cell%n_molecules)=MOLECULE_copy(molec)
          end if
          !deallocate(TMPcell%molecules)
       end if

    end select
    msg(__LINE__),"n_molecules in cell:",cell%n_molecules
    msg(__LINE__),"n_atoms in cell:",cell%n_atoms
    do i=1,cell%n_molecules
       msg(__LINE__),"cell%molecule(",i,")%n_atoms=",cell%molecules(i)%n_atoms
    end do
    msg(__LINE__),"----------------- END of CELL_add_molecule() ---------------"

                


  end subroutine CELL_add_molecules

  ! --------------------------------------------------------------------------------------
  !
  !              CELL_init()
  !
  ! --------------------------------------------------------------------------------------
  function CELL_init() result(dest)
    use global
    implicit none
    integer::i
    type(t_Cell)::dest
    dest%idx=-1
    dest%n_atoms=0
    dest%n_molecules=0
    !allocate(dest%molecules(0))
    dest%n_bonds=0
    !allocate(dest%bonds(0))
    dest%n_angles=0
    !allocate(dest%angles(0))
    dest%periodicity='NONE'
    dest%evacc=0.0
    dest%constraints%n_XYZ=0
    !do i=1,3
    !   dest%k(i)=1
    !end do
    !allocate(dest%constraints%XYZ(dest%constraints%n_XYZ))
  end function CELL_init

  ! --------------------------------------------------------------------------------------
  !
  !              CELL_del()
  !
  ! --------------------------------------------------------------------------------------
  subroutine CELL_del(dest)
    use global
    use Molecule
    use Bond
    implicit none
    type(t_Cell)::dest
    dest%idx=0
    dest%n_bonds=0
    dest%n_atoms=0
    dest%n_molecules=0
    if(allocated(dest%molecules)) then
       deallocate(dest%molecules)
    end if
    if(allocated(dest%bonds)) then
       deallocate(dest%bonds)
    end if
    dest%constraints%n_XYZ=0
    
  end subroutine CELL_del
    ! --------------------------------------------------------------------------------------
  !
  !              CELL_copy()
  !
  ! --------------------------------------------------------------------------------------
  subroutine CELL_copy(src,dest)
    use global
    use Molecule
    use Bond
    implicit none
    type(t_Cell)::src
    type(t_Cell)::dest
    integer::imol,ibond,i

    dest%idx=src%idx
    dest%n_bonds=src%n_bonds
    dest%L=src%L
    dest%n_atoms=src%n_atoms
    dest%periodicity=src%periodicity


    if(.not.(allocated(dest%molecules))) then
       dest%n_molecules=src%n_molecules
       allocate(dest%molecules(dest%n_molecules))
       allocate(dest%bonds(dest%n_bonds))
    end if
    do imol=1,src%n_molecules
       msg(__LINE__), "molecules(",imol,")%n_atoms=",src%molecules(imol)%n_atoms
       dest%molecules(imol)=MOLECULE_copy(src%molecules(imol))
    end do
    
    do ibond=1,src%n_bonds
       call BOND_copy(src%bonds(ibond),dest%bonds(ibond))
    end do


    if(allocated(src%constraints%XYZ))     then
       dest%constraints%XYZ=src%constraints%XYZ
       dest%constraints%n_XYZ=src%constraints%n_XYZ
    else
       dest%constraints%n_XYZ=0
    end if
  end subroutine CELL_copy
  ! --------------------------------------------------------------------------------------
  !
  !              CELL_find_molecules()
  !
  ! --------------------------------------------------------------------------------------
  subroutine CELL_find_molecules(cell)
    use global
    use Atom
    implicit none
    type(t_Cell)::cell
    type(t_Molecule)::mol
    integer::imol,iat,jat,idx,k,Zi,Zj,idx_molecule,n_bonds
    double precision::d,R(3)
    ! 1 - on rassemble tous les atomes des differentes molecules de cell
    !     au sein d'une seule et meme molecule
    mol=MOLECULE_init()
    mol%n_atoms=cell%n_atoms
    print *,mol%n_atoms
    
    allocate(mol%atoms(mol%n_atoms))
    idx=1
    do imol=1,cell%n_molecules
       do iat=1,cell%molecules(imol)%n_atoms
          mol%atoms(idx)=ATOM_copy(cell%molecules(imol)%atoms(iat))
          idx=idx+1
       end do
    end do

    ! 2 - ensuite on calcule les distances
    !     entre tous les atomes de la molécule initiale
    n_bonds=0
    idx_molecule=0
    do iat=1,mol%n_atoms-1
       Zi=mol%atoms(iat)%Zato
       do jat=iat+1,mol%n_atoms
          Zj=mol%atoms(jat)%Zato
          d=0.0
          do k=1,3
             R(k)=mol%atoms(jat)%q(k)-mol%atoms(iat)%q(k)
             R(k)=modulo(R(k),cell%L(k,k))
             if(R(k).gt.0.5*cell%L(k,k)) then
                R(k)=R(k)-cell%L(k,k)
             end if
             if(-R(k).gt.0.5*cell%L(k,k)) then
                R(k)=R(k)+cell%L(k,k)
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
             if ((mol%atoms(jat)%idx_molecule.eq.-1).and.(mol%atoms(iat)%idx_molecule.eq.-1)) then
                idx_molecule=idx_molecule+1
             else
                !
                ! si un des deux atomes a été idx_molécule, celui qui ne l'a pas encore été prend l'idx_molecule
                ! de l'atome déjà indexé
                !
                if (mol%atoms(jat)%idx_molecule*mol%atoms(iat)%idx_molecule<0) then
                   idx_molecule=max(mol%atoms(jat)%idx_molecule,mol%atoms(iat)%idx_molecule)
                else
                   !
                   ! dernier cas de figure : les deux atomes ont déjà idx_molécule mais il appartiennent
                   ! a deux molécule différentes -> ce cas de figure peut se présenter dans le cas de grosse
                   ! molécule -> pas fonctionnelle pour l'instant, A DEVELOPPER
                   !
                   print *,"# A same atom belongs to two different molecules"
                   print *,"### TO BE DEVELOPPED ### ",__FILE__," ### LINE",__LINE__ ; stop        
                end if
             end if
             mol%atoms(iat)%idx_molecule=idx_molecule
             mol%atoms(jat)%idx_molecule=idx_molecule

             !write(*,'(A,I5,I5,A,F12.8,A,A,I5,I5)') &
             !     "# d(",iat,jat,")=",d,molecules(1)%atoms(jat)%elt,molecules(1)%atoms(iat)%elt,&
             !     molecules(1)%atoms(iat)%idx_molecule,molecules(1)%atoms(iat)%idx_molecule
             n_bonds=n_bonds+1
             mol%atoms(iat)%n_bonds=mol%atoms(iat)%n_bonds+1
             mol%atoms(jat)%n_bonds=mol%atoms(jat)%n_bonds+1
             
             mol%n_bonds=mol%n_bonds+1
          end if
          
       end do
    end do
    

    
  end subroutine CELL_find_molecules
  ! --------------------------------------------------------------------------------------
  !
  !              CELL_find_bonds()
  !
  ! --------------------------------------------------------------------------------------
  subroutine CELL_find_molecular_topology(cells,PBC)
    use global
    implicit none
    type(t_Cell)::cells
    integer::imol,jmol,Zi,Zj,iat,jat,k,idxi,idxj,ibond
    double precision::d,R(3)
    logical::PBC

    !print *,"# CELL_find_molecular_topology> n_bonds= ",cells%n_bonds
    !print *,"# CELL_find_molecular_topology> n_molecules= ",cells%n_molecules
    do imol=1,cells%n_molecules
       !print *,"# CELL_find_molecular_topology> molecules ",imol," -> natoms=", cells%molecules(imol)%n_atoms
       do iat=1,cells%molecules(imol)%n_atoms
          Zi=cells%molecules(imol)%atoms(iat)%Zato
          idxi=cells%molecules(imol)%atoms(iat)%idx
          do jmol=1,cells%n_molecules
             do jat=1,cells%molecules(jmol)%n_atoms
                Zj=cells%molecules(jmol)%atoms(jat)%Zato
                !if(.not.((jmol.eq.imol).and.(iat.eq.jat))) then
                idxj=cells%molecules(jmol)%atoms(jat)%idx
                !print *,"# CELL_find_molecular_topology> molecules ",imol," -> natoms=", cells%molecules(imol)%n_atoms

                if(idxj.gt.idxi) then
                   !print *,"# CELL_find_molecular_topology> molecules idxs",idxi,idxj
                   d=0.0
                   do k=1,3
                      R(k)=cells%molecules(jmol)%atoms(jat)%q(k)-cells%molecules(imol)%atoms(iat)%q(k)
                      if(PBC) then
                         R(k)=modulo(R(k),cells%L(k,k))
                         if(R(k).gt.0.5*cells%L(k,k)) then
                            R(k)=R(k)-cells%L(k,k)
                         end if
                         if(-R(k).gt.0.5*cells%L(k,k)) then
                            R(k)=R(k)+cells%L(k,k)
                         end if
                      end if
                      d=d+R(k)*R(k)
                   end do
                   d=sqrt(d)
                   !print *,d,deq(Zi,Zj)
                   if(d<ELEMENT_deq(Zi,Zj))  then
                      !print *,"# CELL_find_molecular_topology> imol=",imol," iat=",iat
                      !print *,"# CELL_find_molecular_topology> ",allocated(cells%molecules(imol)%atoms(iat)%bonds)
                      call CELL_add_bond(cells,idxi,imol,idxj,jmol)
                      !print *,"# CELL_find_molecular_topology> ",allocated(cells%molecules(imol)%atoms(iat)%bonds)
                   end if
                end if
             end do
          end do
       end do
    end do
    print *,"# CELL_find_molecular_topology> n_bonds= ",cells%n_bonds
    !
    do imol=1,cells%n_molecules
       do iat=1,cells%molecules(imol)%n_atoms
          if(cells%molecules(imol)%atoms(iat)%n_bonds.gt.1) then
             if(cells%molecules(imol)%atoms(iat)%n_bonds.eq.2) then
                print *,"#CELL_find_molecular_topology > ",imol,iat,cells%molecules(imol)%atoms(iat)%n_bonds
                cells%n_angles=cells%n_angles+1
                !if(allocated(cells%angles)) then
                !else
                !   allocate(cells%angles(cells%n_angles)
                !end if
                do ibond=1,cells%molecules(imol)%atoms(iat)%n_bonds
                   print *,cells%bonds(cells%molecules(imol)%atoms(iat)%bonds(ibond))%list_atoms(1),&
                        cells%bonds(cells%molecules(imol)%atoms(iat)%bonds(ibond))%list_atoms(2)
                end do
                
                
             else
                print *,"CELL_find_molecular_topology> "," Part to be implemented"
                stop
             end if
          end if
       end do
    end do
    
  end subroutine CELL_find_molecular_topology
  ! --------------------------------------------------------------------------------------
  !
  !              CELL_free_location()
  ! Cette routine permet de vérifier que les atomes de la nouvelle molécule que l'on ajoute
  ! ne sont pas trop proches des atomes déjà présents.
  !
  ! --------------------------------------------------------------------------------------
  function CELL_free_location(cell,molecule) result(free_location)
    use global
    implicit none
    type(t_Cell)::cell
    type(t_Molecule)::molecule
    integer::iconf,jat,imol,k,iat
    logical::free_location
    double precision::R(3),d
    double precision,parameter::dlim=2.0

    free_location=.True.
    imol=0
    do while((imol.lt.cell%n_molecules).and.free_location)
       imol=imol+1
       jat=0
       do while((jat.lt.cell%molecules(imol)%n_atoms).and.free_location)
          jat=jat+1
          iat=0
          do while((iat.lt.molecule%n_atoms).and.free_location)
             iat=iat+1
             d=0.0
             do k=1,3
                R(k)=cell%molecules(imol)%atoms(jat)%q(k)-molecule%atoms(iat)%q(k)
                R(k)=modulo(R(k),cell%L(k,k))
                if(R(k).gt.0.5*cell%L(k,k)) then
                   R(k)=R(k)-cell%L(k,k)
                end if
                if(-R(k).gt.0.5*cell%L(k,k)) then
                   R(k)=R(k)+cell%L(k,k)
                end if
                d=d+R(k)*R(k)
             end do
             d=sqrt(d)
             !print *,"d=",d
             if(d.lt.dlim) then
                free_location=.False.
             end if
          end do
       end do
    end do
    msg(__LINE__), " Free Location?",free_location
  end function CELL_free_location

  ! --------------------------------------------------------------------------------------
  !
  !              CELL_rebox()
  !
  ! --------------------------------------------------------------------------------------
  subroutine CELL_rebox(cell)
    use global
    implicit none
    type(t_Cell)::cell
    integer::imol,iat,j
    double precision::mini(3),maxi(3)
    do j=1,3
       mini(j)=cell%molecules(1)%atoms(1)%q(j)
       maxi(j)=cell%molecules(1)%atoms(1)%q(j)
    end do
    do imol=1,cell%n_molecules
       do iat=1,cell%molecules(imol)%n_atoms
          !print *,(cell%molecules(imol)%atoms(iat)%q(j),j=1,3)
          do j=1,3
             if(cell%molecules(imol)%atoms(iat)%q(j).lt.mini(j)) then
                mini(j)=cell%molecules(imol)%atoms(iat)%q(j)
             end if
             if(cell%molecules(imol)%atoms(iat)%q(j).gt.maxi(j)) then
                maxi(j)=cell%molecules(imol)%atoms(iat)%q(j)
             end if
          end do
       end do
    end do
    print *,(mini(j),j=1,3)
    print *,(maxi(j),j=1,3)
    do j=1,3
       cell%L(j,j)=maxi(j)-mini(j)
    end do

    do imol=1,cell%n_molecules
       do iat=1,cell%molecules(imol)%n_atoms
          !print *,(cell%molecules(imol)%atoms(iat)%q(j),j=1,3)
          do j=1,3
             cell%molecules(imol)%atoms(iat)%q(j)=cell%molecules(imol)%atoms(iat)%q(j)-mini(j)
          end do
       end do
    end do

    write(*,*) "    &CELL"
    write(*,*) "      A [angstrom] ",(cell%L(1,j),j=1,3)
    write(*,*) "      B [angstrom]",(cell%L(2,j),j=1,3)
    write(*,*) "      C [angstrom]",(cell%L(3,j),j=1,3)
    write(*,*) "      PERIODIC ",cell%periodicity
    write(*,*) "    &END CELL"


    
  end subroutine CELL_rebox
  ! --------------------------------------------------------------------------------------
  !
  !              CELL_set_unit_cell()
  !
  ! --------------------------------------------------------------------------------------
  subroutine CELL_set_lattice(cell,n,a,csura,evacc,ibrav)
    use global
    implicit none
    type(t_Cell)::cell
    integer::n(3),i,j,ibrav
    double precision::evacc,a,csura
    ! ibrav=4  Hexagonal and Trigonal P        celldm(3)=c/a
    ! L1 = (a/2)(-1,1,1), L2 = (a/2)(1,-1,1),  L3 = (a/2)(1,1,-1)
    cell%L(1,1)= a ;     cell%L(1,2)= 0.0;              cell%L(1,3)=0.0 ;
    cell%L(2,1)=-a*.5 ;  cell%L(2,2)= a*sqrt(3.0)*0.5;  cell%L(2,3)=0.0 ;
    cell%L(3,1)= 0.0 ;   cell%L(3,2)= 0.0;              cell%L(3,3)=3*csura*a ;
    do j=1,3
       do i=1,3
          cell%L(j,i)=cell%L(j,i)*n(j)
       end do
    end do
    cell%L(3,3)=cell%evacc+cell%L(3,3)
  end subroutine CELL_set_lattice
  ! --------------------------------------------------------------------------------------
  !
  !              CELL_set_unit_cell()
  !
  ! --------------------------------------------------------------------------------------
  function CELL_set_unit_cell(ibrav,a,csura,elt) result(cell)
    use global
    implicit none
    integer::ibrav
    character(len=2)::elt
    double precision::a,csura
    type(t_Cell)::cell
    cell%n_molecules=1
    allocate(cell%molecules(cell%n_molecules))
    cell%molecules(1)%n_atoms=3
    allocate(cell%molecules(1)%atoms(cell%molecules(1)%n_atoms))
    cell%molecules(1)%atoms(1)%elt=elt 
    cell%molecules(1)%atoms(1)%q=(/0.0,0.0,0.0/)
    cell%molecules(1)%atoms(2)%elt=elt
    cell%molecules(1)%atoms(2)%q(1)=0.5*a
    cell%molecules(1)%atoms(2)%q(2)=0.5*sqrt(3.0)*a/3.0
    cell%molecules(1)%atoms(2)%q(3)=csura*a
    cell%molecules(1)%atoms(3)%elt=elt
    cell%molecules(1)%atoms(3)%q(1)=0.5*a
    cell%molecules(1)%atoms(3)%q(2)=-0.5*sqrt(3.0)*a/3.0
    cell%molecules(1)%atoms(3)%q(3)=2*csura*a
    
  end function CELL_set_unit_cell
  ! --------------------------------------------------------------------------------------
  !
  !              CELL_slab()
  !
  ! --------------------------------------------------------------------------------------
  function CELL_slab(ibrav,a,csura,n,unit_cell,evacc) result(cell)
    use global
    use Atom
    implicit none
    type(t_Cell)::cell,unit_cell
    integer::ibrav,n(3),n0(3),i,j,k,l,m
    double precision::a,csura,q(3),evacc,evacc0
    character(len=2)::elt
    type(t_Atom)::atm


    cell=CELL_init()
    cell%n_molecules=1
    allocate(cell%molecules(cell%n_molecules))
    cell%molecules(1)%n_atoms=0
    cell%molecules(1)=MOLECULE_init()
    cell%evacc=evacc

    evacc0=0.0
    n0=(/1,1,1/)
    call CELL_set_lattice(cell,n0,a,csura,evacc0,ibrav)

    do k=0,n(3)-1
       do j=0,n(2)-1
          do i=0,n(1)-1
             do l=1,unit_cell%molecules(1)%n_atoms
                atm=ATOM_new()
                atm%elt=unit_cell%molecules(1)%atoms(l)%elt
                do m=1,3
                   atm%q(m)=unit_cell%molecules(1)%atoms(l)%q(m)+i*cell%L(1,m)+j*cell%L(2,m)+k*cell%L(3,m)
                end do
                if(atm%q(3).lt.1.0) then
                   atm%constraints=(/0,0,0/)
                end if
                call MOLECULE_add_atom(cell%molecules(1),atm)
                print *,"# CELL_slab> n_XYZ=",cell%molecules(1)%constraints%n_XYZ
                !if(k.eq.0) then
                !   unit_cell%molecules(1)%atoms(l)%constraints=(/0,0,1/)
                !end if
                cell%n_atoms=cell%n_atoms+1
                call ATOM_del(atm)
             end do
          end do
       end do
    end do

    call CELL_set_lattice(cell,n,a,csura,evacc,ibrav)
    
    





    
    print *,"# n_atom=",cell%molecules(1)%n_atoms
    do i=1,cell%molecules(1)%n_atoms
       print *,cell%molecules(1)%atoms(i)%elt,(cell%molecules(1)%atoms(i)%q(j),j=1,3)
    end do
    
  end function CELL_slab
  ! --------------------------------------------------------------------------------------
  !
  !              CELL_update()
  !
  ! --------------------------------------------------------------------------------------
  subroutine CELL_update(cell)
    use global
    implicit none
    type(t_Cell)::cell
    integer::imol,iat,idx_cur
    idx_cur=1
    do imol=1,cell%n_molecules
       write(*,*) "# Molecule ",imol,": n_atoms=",cell%molecules(imol)%n_atoms,&
            ": n_elements=",cell%molecules(imol)%n_elements
       if(.not.(allocated(cell%elements))) then
          cell%n_elements=cell%molecules(imol)%n_elements
          allocate(cell%elements(cell%n_elements))
          cell%elements=cell%molecules(imol)%elements
       else
          print *,"# CELL_update: To be implemented"
          stop
       end if

       do iat=1,cell%molecules(imol)%n_atoms
          cell%molecules(imol)%atoms%idx=idx_cur
          cell%molecules(imol)%atoms%idx_molecule=imol
          idx_cur=idx_cur+1
       end do

       
    end do
  end subroutine CELL_update


end module Cell
