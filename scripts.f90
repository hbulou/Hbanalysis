module scripts
  implicit none
contains
#define msg(x) print *,"### ",x,__FILE__," > "
!#define msg(x) print *,"### ",x,' ###    ',__FILE__," > "
#define error(x) print *,"### ERROR ",x,__FILE__," > "
#define debug(x,idx) print *,"### DEBUG ",idx," @ line ",x,__FILE__," > "

  !  subroutine interpol(param)
  !  subroutine concat(param)
  !  subroutine crystal(param)
  !  subroutine solvent()
  !  subroutine Ti()
  !  subroutine TiOH()
  !  subroutine TiOH_addOH()
  !  subroutine TiOH_moveOH()


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !               subroutine test()
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine test(param)
    use global
    use IO
    
    implicit none
    type(t_Param)::param
    type(t_CP2K_param)::cp2k_param
    msg(__LINE__), trim(param%input(1)%name)
    call IO_read_cp2k_restart_file(param%input(1)%name,cp2k_param)
    
    call IO_save(param%output(1)%name,cp2k_param)

  end subroutine test
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !               subroutine interpol(param)
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine interpol(param)
    use global
    use Atom
    use Cell
    use Molecule
    use IO
    use Univers
    implicit none
    integer::i,j,n,k
    type(t_Param)::param
    type(t_Molecule)::mol(2),newmol
    type(t_Atom)::atm
    type(t_Univers)::univ
    logical::chk_free_location
    character(len=32)::adding_type

    param%calculation%system="essai"

    if(param%n_inputs.lt.2) then
       error(__LINE__), " You must provide two topology files at least"
       stop
    end if

    univ=UNIVERS_new()
    univ%n_configurations=0
    call UNIVERS_add_configuration(univ)
    newmol=MOLECULE_init()
    !msg(__LINE__),"newmol%n_elements=",newmol%n_elements

    msg(__LINE__), "-----------------------------------------------------"
    msg(__LINE__),"Reading molecules"
    msg(__LINE__), "-----------------------------------------------------"
    do i=1,2
       msg(__LINE__), " input ",i,":",&
            trim(param%input(i)%name)," (",trim(param%input(i)%type)," )"
       mol(i)=IO_read_molecule(param%input(i))
       call MOLECULE_element_list_update(mol(i))
       msg(__LINE__), " n_elements in ",i," = ",mol(i)%n_elements
       do j=1,mol(i)%n_elements
          msg(__LINE__), mol(i)%elements(j)%elt, mol(i)%elements(j)%n
       end do
    end do
    
    do i=1,mol(1)%n_atoms
       atm=ATOM_copy(mol(1)%atoms(i))
       atm%q=0.5*(mol(1)%atoms(i)%q+mol(2)%atoms(i)%q)
       call MOLECULE_add_atom(newmol,atm)
       call ATOM_del(atm)
    end do
    msg(__LINE__), "-----------------------------------------------------"
    msg(__LINE__), "newmol n_atoms=",newmol%n_atoms
    msg(__LINE__), "-----------------------------------------------------"
    
    
    debug(__LINE__,"1")
    msg(__LINE__), " n_atoms=",univ%configurations(1)%cells%n_atoms
    call IO_read_cell_from_cp2k_restart_file(param,param%input(3),univ%configurations(1)%cells)
    call MOLECULE_del(univ%configurations(1)%cells%molecules(1))
    univ%configurations(1)%cells%n_molecules=0
    !deallocate(univ%configurations(1)%cells%molecules)

    
    do i=1,3
       msg(__LINE__),(univ%configurations(1)%cells%L(i,k),k=1,3)
    end do

    
    n=1
    adding_type='as_it'
    chk_free_location=.FALSE.
    call CELL_add_molecules(univ%configurations(1)%cells,newmol,n,adding_type,chk_free_location)
    univ%configurations(1)%cells%n_atoms=newmol%n_atoms
    msg(__LINE__), " n_molecules=",univ%configurations(1)%cells%n_molecules
    msg(__LINE__), " n_atoms=",univ%configurations(1)%cells%n_atoms

    
    call IO_write(param,univ,n)
    call save_molecule(newmol)

  end subroutine interpol
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !               subroutine concat(param)
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine concat(param)
    use global
    use tools
    use Atom
    use Cell
    use Molecule
    use IO
    use Univers
    implicit none
    integer::i,j,n,k,imol
    type(t_Param)::param
    type(t_Molecule)::mol(2),newmol
    type(t_Atom)::atm
    type(t_Univers)::univ
    character(len=32)::center_tag
    double precision::alpha,beta,gamma,t(3),P(3),M(3),S(3),Ssqr,Psqr,Msqr,ScrossMsqr,Mp(3),Mpsqr
    double precision::ex(3),ey(3),ez(3),SdotM,ScrossM(3),PdotS,cosbeta,sinbeta,Mpx,Mpy,Mx,My,angle
    double precision::center(3),axe(3),axe_to_align(3),pos_to_center(3),pos(3)
    character(len=32)::adding_type
    integer::iO1sub,iH1sub,iO2sub,iH2sub,iO1mol,iH1mol,iO2mol,iH2mol,isub,imolecule,iat,iCmol,iO3mol,iPmol
    logical::chk_free_location



    if(param%n_inputs.lt.2) then
       error(__LINE__), " You must provide two topology files at least"
       stop
    end if

    univ=UNIVERS_new()
    univ%n_configurations=0
    call UNIVERS_add_configuration(univ)

    msg(__LINE__), "-----------------------------------------------------"
    msg(__LINE__),"Reading cp2k restart file"
    msg(__LINE__), "-----------------------------------------------------"
    !debug(__LINE__,"1")
    call IO_read_cell_from_cp2k_restart_file(param,param%input(3),univ%configurations(1)%cells)
    call MOLECULE_del(univ%configurations(1)%cells%molecules(1))
    univ%configurations(1)%cells%n_molecules=0
    univ%configurations(1)%cells%n_atoms=0
    !deallocate(univ%configurations(1)%cells%molecules)



    do i=1,3
       msg(__LINE__),(univ%configurations(1)%cells%L(i,k),k=1,3)
    end do



    msg(__LINE__), "-----------------------------------------------------"
    msg(__LINE__),"Reading molecules"
    msg(__LINE__), "-----------------------------------------------------"
    n=1
    adding_type='as_it'
    chk_free_location=.FALSE.
    do i=1,2
       msg(__LINE__), " ---------------------------------------------------------------------"
       msg(__LINE__), " input ",i,":",&
            trim(param%input(i)%name)," (",trim(param%input(i)%type)," )"
       mol(i)=IO_read_molecule(param%input(i))
       call MOLECULE_element_list_update(mol(i))
       msg(__LINE__), " n_elements in ",i," = ",mol(i)%n_elements
       do j=1,mol(i)%n_elements
          msg(__LINE__), mol(i)%elements(j)%elt, mol(i)%elements(j)%n
          !if(i.eq.1) then
          !   newmol%n_atoms=mol(i)%n_atoms
          !else
          !   if(.not.(mol(i)%elements(j)%elt.eq."Ti")) newmol%n_atoms=newmol%n_atoms+mol(i)%elements(j)%n
          !end if
       end do
       msg(__LINE__),"Going to CELL_add_molecules() subroutine ..." 
       call CELL_add_molecules(univ%configurations(1)%cells,mol(i),n,adding_type,chk_free_location) 
       msg(__LINE__), "univ%configurations(1)%cells%molecules(",i,")%n_atoms=", &
            univ%configurations(1)%cells%molecules(i)%n_atoms,mol(i)%n_atoms
    end do


    msg(__LINE__), " ----------------- Computing intermolular distance ---------------------"
    msg(__LINE__), "univ%configurations(1)%cells%n_molecules=",&
         univ%configurations(1)%cells%n_molecules

    !
    ! donnée en entrée
    !
    ! isub=1
    ! iO1sub=301   ! centre de rotation
    ! iO2sub=302   ! axe de rotation
    ! iH1sub=303   ! centre de rotation
    ! iH2sub=304   ! axe de rotation

    ! imolecule=2
    ! iO1mol=3    ! point à ramener au centre
    ! iO2mol=4     ! axe de rotation à aligner
    ! iH1mol=9    ! point à ramener au centre
    ! iH2mol=10     ! axe de rotation à aligner

    ! angle=0.0

    isub=1
    iO1sub=301   ! centre de rotation
    iO2sub=302   ! axe de rotation
    iH1sub=303   ! centre de rotation
    iH2sub=304   ! axe de rotation

    imolecule=2
    iO1mol=3    ! point à ramener au centre
    iO2mol=2     ! axe de rotation à aligner
    iH1mol=9    ! point à ramener au centre
    iH2mol=28     ! axe de rotation à aligner
    iCmol=5
    iO3mol=4
    iPmol=1

    angle=0.0

    !
    ! then we build the axis for rotating the molecule
    !
    center        = univ%configurations(1)%cells%molecules(isub)%atoms(iO1sub)%q
    axe           = univ%configurations(1)%cells%molecules(isub)%atoms(iO2sub)%q
    pos_to_center = univ%configurations(1)%cells%molecules(imolecule)%atoms(iO1mol)%q
    axe_to_align  = univ%configurations(1)%cells%molecules(imolecule)%atoms(iO2mol)%q
    call TOOLS_build_rotation_axis(center,axe,pos_to_center,axe_to_align,ex,ey,ez,alpha)
    msg(__LINE__), "ex=",ex
    msg(__LINE__), "ey=",ey
    msg(__LINE__), "ez=",ez

    !
    ! rotation of the molecule so that to align the two axes
    !
    do iat=1,univ%configurations(1)%cells%molecules(imolecule)%n_atoms
       msg(__LINE__),"------------------------------------------------------------------------"
       msg(__LINE__), univ%configurations(1)%cells%molecules(imolecule)%atoms(iat)%elt
       univ%configurations(1)%cells%molecules(imolecule)%atoms(iat)%q=&
            TOOLS_rotate_around_ez(pos_to_center,ex,ey,ez,alpha,univ%configurations(1)%cells%molecules(imolecule)%atoms(iat)%q)
    end do

    alpha=-25.00
    center        = 0.5*(univ%configurations(1)%cells%molecules(imolecule)%atoms(iO1mol)%q+&
         univ%configurations(1)%cells%molecules(imolecule)%atoms(iO2mol)%q)
    do iat=1,univ%configurations(1)%cells%molecules(imolecule)%n_atoms
       univ%configurations(1)%cells%molecules(imolecule)%atoms(iat)%q=&
            TOOLS_rotate_around_ez(center,ez,ex,ey,alpha,univ%configurations(1)%cells%molecules(imolecule)%atoms(iat)%q)
       
    end do

    
    center        = univ%configurations(1)%cells%molecules(imolecule)%atoms(iPmol)%q
    axe           = univ%configurations(1)%cells%molecules(imolecule)%atoms(iCmol)%q
    axe_to_align  = univ%configurations(1)%cells%molecules(imolecule)%atoms(iO3mol)%q
    call TOOLS_build_rotation_axis(center,axe,center,axe_to_align,ex,ey,ez,alpha)
    msg(__LINE__), "ex=",ex
    msg(__LINE__), "ey=",ey
    msg(__LINE__), "ez=",ez

    angle=180.0
    do iat=1,univ%configurations(1)%cells%molecules(imolecule)%n_atoms
       if((.not.(iat.eq.iPmol)).and.(.not.(iat.eq.iCmol)).and.&
            (.not.(iat.eq.iO1mol)).and.(.not.(iat.eq.iO2mol))&
            .and.(.not.(iat.eq.iO3mol))) then
       univ%configurations(1)%cells%molecules(imolecule)%atoms(iat)%q=&
            TOOLS_rotate_around_ez(center,ez,ex,ey,angle,univ%configurations(1)%cells%molecules(imolecule)%atoms(iat)%q)
    end if
    end do

    center        = univ%configurations(1)%cells%molecules(isub)%atoms(iO1sub)%q
    pos_to_center = univ%configurations(1)%cells%molecules(imolecule)%atoms(iO1mol)%q
    t=center-pos_to_center
    call MOLECULE_translate(univ%configurations(1)%cells%molecules(imolecule),t) 


    
    
    msg(__LINE__), univ%configurations(1)%cells%molecules(imolecule)%atoms(iH1mol)%q
    ! call MOLECULE_align(molref,center_rot,axe,mol_align,to_center,to_align,angle)
    !         double precision :: center(3)
    !         double precision :: axe(3)
    !         type(t_Molecule) :: mol_align
    !         double precision :: to_center(3)
    !         double precision :: to_align(3)
    !center=univ%configurations(1)%cells%molecules(isub)%atoms(iO1sub)%q
    !axe=univ%configurations(1)%cells%molecules(isub)%atoms(iO2sub)%q
    !pos_to_center=univ%configurations(1)%cells%molecules(imolecule)%atoms(iO1mol)%q
    !axe_to_align=univ%configurations(1)%cells%molecules(imolecule)%atoms(iO2mol)%q
    
    !call MOLECULE_align(center,axe,&
    !     univ%configurations(1)%cells%molecules(imolecule),&
    !     pos_to_center,axe_to_align,angle)
  


    ! imolecule=2
    !center=univ%configurations(1)%cells%molecules(imolecule)%atoms(1)%q
    !axe=univ%configurations(1)%cells%molecules(imolecule)%atoms(5)%q
    !axe_to_align=univ%configurations(1)%cells%molecules(imolecule)%atoms(4)%q
    !msg(__LINE__), "center=",center
    !msg(__LINE__), "axe=",axe
    !msg(__LINE__), "direction=",axe-center

    !call TOOLS_build_rotation_axis(center,axe,center,axe_to_align,ex,ey,ez,alpha)
    !msg(__LINE__), "ex=",ex
    !msg(__LINE__), "ey=",ey
    !msg(__LINE__), "ez=",ez

    !angle=180.0
    !ex=(/1.0,0.0,0.0/)
    !ey=(/0.0,1.0,0.0/)
    !ez=(/0.0,0.0,1.0/)
    
    !do iat=1,univ%configurations(1)%cells%molecules(imolecule)%n_atoms
    !   univ%configurations(1)%cells%molecules(imolecule)%atoms(iat)%q=TOOLS_rotate_around_ez(center,ex,ey,ez,angle,&
    !        univ%configurations(1)%cells%molecules(imolecule)%atoms(iat)%q)
    !end do
    
    !center=(/0.0,0.0,0.5/)
    !call MOLECULE_translate(univ%configurations(1)%cells%molecules(imolecule),center) 


    ! axe_to_align=univ%configurations(1)%cells%molecules(imolecule)%atoms(5)%q
    
    ! call MOLECULE_translate(univ%configurations(1)%cells%molecules(imolecule),center) 
    ! call TOOLS_build_rotation_axis(center,axe,center,axe_to_align,ex,ey,ez,angle)
    ! msg(__LINE__), "ex=",ex
    ! msg(__LINE__), "ey=",ey
    ! msg(__LINE__), "ez=",ez
    
    
    !alpha=45.0
    !beta=90.0
    !gamma=0.0
    !center_tag='molecule_MC'
    !call MOLECULE_rotate(univ%configurations(1)%cells%molecules(2),alpha,beta,gamma,center_tag)
    !t=(/7.0,7.0,15.0/)
    !call molecule_translate_mass_center_at(univ%configurations(1)%cells%molecules(2),t)

    !univ%configurations(1)%cells%molecules(isub)%save=.False.
    univ%configurations(1)%cells%molecules(isub)%atoms(iH2sub)%save=.False.
    univ%configurations(1)%cells%molecules(isub)%atoms(iH1sub)%save=.False.
    univ%configurations(1)%cells%molecules(isub)%atoms(iO1sub)%save=.False.
    univ%configurations(1)%cells%molecules(isub)%atoms(iO2sub)%save=.False.
    !univ%configurations(1)%cells%molecules(isub)%atoms(iO1sub)%Zato=7
    !univ%configurations(1)%cells%molecules(isub)%atoms(iO2sub)%Zato=7
    univ%configurations(1)%cells%molecules(imolecule)%atoms(iH1mol)%save=.False.
    univ%configurations(1)%cells%molecules(imolecule)%atoms(iH2mol)%save=.False.
    n=1
    call IO_write(param,univ,n)
    
    !stop

    
    ! newmol=MOLECULE_init()
    ! !msg(__LINE__),"newmol%n_elements=",newmol%n_elements
    
    ! do i=1,2
    !    msg(__LINE__), " input ",i,":",&
    !         trim(param%input(i)%name)," (",trim(param%input(i)%type)," )"

    !    !mol(i)=MOLECULE_from_file(param%input(i))
    !    mol(i)=IO_read_molecule(param%input(i))
    !    !msg(__LINE__),"newmol%n_elements=",newmol%n_elements
    !    call MOLECULE_element_list_update(mol(i))
    !    msg(__LINE__), " n_elements in ",i," = ",mol(i)%n_elements
    !    do j=1,mol(i)%n_elements
    !       msg(__LINE__), mol(i)%elements(j)%elt, mol(i)%elements(j)%n
    !       !if(i.eq.1) then
    !       !   newmol%n_atoms=mol(i)%n_atoms
    !       !else
    !       !   if(.not.(mol(i)%elements(j)%elt.eq."Ti")) newmol%n_atoms=newmol%n_atoms+mol(i)%elements(j)%n
    !       !end if
    !    end do
    ! end do

    
 


    ! iH1sub=303
    ! iH2sub=304

    ! iH1mol=10
    ! iH2mol=9
    ! write(*,*) MOLECULE_intermolecular_distance(mol(1),iH1sub,mol(2),iH1mol)
    ! write(*,*) MOLECULE_intermolecular_distance(mol(1),iH2sub,mol(2),iH2mol)

    ! stop
    ! do i=1,mol(1)%n_atoms
    !    if(mol(1)%atoms(i)%elt.eq."Ti") then
    !       atm=ATOM_copy(mol(1)%atoms(i))
    !       atm%q=0.5*(mol(1)%atoms(i)%q+mol(2)%atoms(i)%q)
    !       call MOLECULE_add_atom(newmol,atm)
    !       call ATOM_del(atm)
    !    else
    !       do j=1,2
    !          atm=ATOM_copy(mol(j)%atoms(i))
    !          call MOLECULE_add_atom(newmol,atm)
    !          call ATOM_del(atm)
    !       end do
    !    end if
    ! end do
    ! msg(__LINE__), "newmol n_atoms=",newmol%n_atoms


    ! do imol=1,2
    !    do i=1,mol(imol)%n_atoms
    !       atm=ATOM_copy(mol(imol)%atoms(i))
    !       call MOLECULE_add_atom(newmol,atm)
    !       call ATOM_del(atm)
    !    end do
    ! end do
    ! msg(__LINE__), "newmol n_atoms=",newmol%n_atoms

    
    



    
    
    ! n=1
    ! adding_type='as_it'
    ! chk_free_location=.True.
    ! call CELL_add_molecules(univ%configurations(1)%cells,newmol,n,adding_type,chk_free_location)
    ! msg(__LINE__), " n_molecules=",univ%configurations(1)%cells%n_molecules

    
    ! call IO_write(param,univ,n)
    ! !call save_molecule(newmol)

  end subroutine concat



  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !             subroutine crystal()
  !   construit un substrat 
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine crystal(param)
    use global
    use Univers
    use Atom
    use IO
    use Cell
    implicit none
    type(t_Param)::param
    type(t_Cell)::unit_cell
    type(t_Univers)::univ
    integer::ibrav,i,j,nat(3),iconf
    double precision::a,csura,evacc,q(3)
    double precision,parameter::ZERO=0.0
    !type(t_File)::file
    type(t_Atom)::atm
    logical::ads
    character(len=2)::elt
    type(t_Molecule)::rlx
    integer::idx_conf_to_read
    iconf=1
    univ=UNIVERS_new()
    
    call UNIVERS_add_configuration(univ)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !           lowest part of the slab
    !
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ibrav=4 ! Hexagonal and Trigonal P
    elt="Pd"
    a=2.89

    csura=sqrt(6.0)/3.0
    unit_cell=CELL_set_unit_cell(ibrav,a,csura,elt)
    ads=.FALSE.
    !ads=.TRUE.
    
    nat=(/10,10,1/) 
    !nat=(/10,10,1/) 
    !nat=(/8,8,1/) 
    !nat=(/6,6,1/) 
    !nat=(/4,4,1/) 
    !nat=(/1,1,1/) 
    evacc=0.0
    !evacc=10.0/nat(3)
    !evacc=1.0/nat(3)
    param%input%name="cp2k.in"
    param%input%type="cp2k_input"
    param%title="run"

    param%calculation%run_type="WAVEFUNCTION_OPTIMIZATION"
    !param%calculation%run_type="GEO_OPT"
    param%calculation%ecutoff=400.0
    param%calculation%basis_name="BASIS_MOLOPT"   ! BASIS_SET
    param%calculation%potential_name="POTENTIAL"   
    param%calculation%xc_functional="PBE"
    param%calculation%k(1)=1
    param%calculation%k(2)=param%calculation%k(1)
    param%calculation%k(3)=param%calculation%k(1)
    if(param%calculation%k(1).eq.1) then
       param%calculation%scheme="gamma"
    else
       param%calculation%scheme="MONKHORST-PACK"
    end if
    
    PeriodicTable(22)%Basis_set="DZVP-MOLOPT-SR-GTH"
    PeriodicTable(22)%Potential_set="GTH-PBE-q12"
    PeriodicTable(46)%Basis_set="DZVP-MOLOPT-SR-GTH"
    PeriodicTable(46)%Potential_set="GTH-PBE-q10"
    PeriodicTable(46)%Potential_set="GTH-PBE-q18"

    
    univ%configurations(1)%cells=CELL_slab(ibrav,a,csura,nat,unit_cell,evacc)
    !print *,"# k1= ",(univ%configurations(iconf)%cells%k(i),i=1,3)
    print *,"# MAIN> n_XYZ=",univ%configurations(1)%cells%molecules(1)%constraints%n_XYZ
    

    
    
    do i=1,3
       print *,(univ%configurations(1)%cells%L(i,j),j=1,3)
    end do
    print *,"# celldm(1)=",univ%configurations(1)%cells%L(1,1)/a0
    print *,"# celldm(3)=",univ%configurations(1)%cells%L(3,3)/univ%configurations(1)%cells%L(1,1)
    
    
    print *,"# n_atoms in new=",univ%configurations(1)%cells%molecules(1)%n_atoms
    print *,"# PERIODICITY= ",univ%configurations(iconf)%cells%periodicity
    univ%configurations(iconf)%cells%periodicity="XYZ"
    !univ%configurations(iconf)%cells%periodicity="XY"


    
    


    iconf=1 ; call IO_write(param,univ,iconf)

    !iconf=1;     param%output%name="essai.xsf" ;   call save(param%output%name,univ,iconf)
    
    
  end subroutine crystal
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !             subroutine solvent()
  !   construit un substrat et ajoute une molecule de OH dessus (si ads=.TRUE.
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine solvent(param,idx_inputs)
    use global
    use Molecule
    use Univers
    use Cell
    use IO
    implicit none
    integer::idx_inputs
    type(t_Param)::param
    integer::n,j
    character(len=32)::adding_type
    type(t_Molecule)::h2o
    type(t_Univers)::univ
    logical::chk_free_location
    !h2o=MOLECULE_from_file(param%input(idx_inputs))

    h2o=IO_read_molecule(param%input(idx_inputs))
    univ=UNIVERS_new()
    univ%n_configurations=0
    call UNIVERS_add_configuration(univ)


    n=1
    adding_type='as_it'
    chk_free_location=.TRUE.
    call CELL_add_molecules(univ%configurations(1)%cells,h2o,n,adding_type,chk_free_location)

    call CELL_find_molecules(univ%configurations(1)%cells)

    call CELL_update(univ%configurations(1)%cells)

    
    call CELL_rebox(univ%configurations(1)%cells)

    do j=1,3
       univ%configurations(1)%cells%L(j,j)=45.0
    end do
    



    param%title="H2O"
    param%calculation%run_type="GEO_OPT"
    call IO_write(param,univ,n)
    
  end subroutine solvent


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !             subroutine Ti()
  !   construit un substrat 
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Ti()
    use global
    use Univers
    use Atom
    use IO
    use Cell
    implicit none
    type(t_Param)::param
    type(t_Cell)::Ti_unit_cell
    type(t_Univers)::univ
    integer::ibrav,i,j,nat(3),iconf
    double precision::a,csura,evacc,q(3)
    double precision,parameter::ZERO=0.0
    !type(t_File)::file
    type(t_Atom)::atm
    logical::ads
    character(len=2)::elt
    type(t_Molecule)::rlx
    integer::idx_conf_to_read
    iconf=1
    univ=UNIVERS_new()
    
    call UNIVERS_add_configuration(univ)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !           lowest part of the slab
    !
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ibrav=4 ! Hexagonal and Trigonal P
    elt="Ti"
    a=2.92

    csura=sqrt(6.0)/3.0
    Ti_unit_cell=CELL_set_unit_cell(ibrav,a,csura,elt)
    ads=.FALSE.
    !ads=.TRUE.
    
    nat=(/12,12,1/) 
    !nat=(/10,10,1/) 
    !nat=(/8,8,1/) 
    !nat=(/6,6,1/) 
    !nat=(/4,4,1/) 
    !nat=(/1,1,1/) 
    evacc=20.0
    !evacc=10.0/nat(3)
    !evacc=1.0/nat(3)
    param%input%name="cp2k.in"
    param%input%type="cp2k_input"
    param%title="run"
    param%machine="jeanzay"
    param%machine="hpc"

    param%calculation%run_type="WAVEFUNCTION_OPTIMIZATION"
    param%calculation%run_type="GEO_OPT"
    param%calculation%ecutoff=400.0
    param%calculation%basis_name="BASIS_MOLOPT"   ! BASIS_SET
    param%calculation%potential_name="POTENTIAL"   
    param%calculation%xc_functional="PBE"
    param%calculation%scheme="gamma"
    !param%scheme="MONKHORST-PACK"
    param%calculation%k(1)=1
    param%calculation%k(2)=param%calculation%k(1)
    param%calculation%k(3)=1

    PeriodicTable(22)%Basis_set="DZVP-MOLOPT-SR-GTH"
    PeriodicTable(22)%Potential_set="GTH-PBE-q12"
    ! Ti
    !  Ti SZV-MOLOPT-SR-GTH SZV-MOLOPT-SR-GTH-q12
    !  O  SZV-MOLOPT-GTH SZV-MOLOPT-GTH-q6
    !  O  SZV-MOLOPT-SR-GTH SZV-MOLOPT-SR-GTH-q6
    !  H  SZV-MOLOPT-GTH SZV-MOLOPT-GTH-q1



    ! Ti GTH-BLYP-q12 GTH-BLYP
    ! Ti GTH-BP-q12 GTH-BP
    ! Ti GTH-PADE-q12 GTH-LDA-q12 GTH-PADE GTH-LDA
    ! Ti GTH-PADE-q4 GTH-LDA-q4
    ! Ti GTH-PBE-q12 GTH-PBE
    ! Ti ALLELECTRON ALL
    
    !  Ti DZVP-MOLOPT-SR-GTH DZVP-MOLOPT-SR-GTH-q12
    !  O  DZVP-MOLOPT-GTH    DZVP-MOLOPT-GTH-q6
    !  O  DZVP-MOLOPT-SR-GTH DZVP-MOLOPT-SR-GTH-q6
    !  H  DZVP-MOLOPT-GTH    DZVP-MOLOPT-GTH-q1



    
    !  O  TZVP-MOLOPT-GTH TZVP-MOLOPT-GTH-q6
    !  O  TZV2P-MOLOPT-GTH TZV2P-MOLOPT-GTH-q6
    !  O  TZV2PX-MOLOPT-GTH TZV2PX-MOLOPT-GTH-q6
    !  H  TZVP-MOLOPT-GTH TZVP-MOLOPT-GTH-q1
    !  H  TZV2P-MOLOPT-GTH TZV2P-MOLOPT-GTH-q1
    !  H  TZV2PX-MOLOPT-GTH TZV2PX-MOLOPT-GTH-q1

    
    univ%configurations(1)%cells=CELL_slab(ibrav,a,csura,nat,Ti_unit_cell,evacc)
    !print *,"# k1= ",(univ%configurations(iconf)%cells%k(i),i=1,3)
    print *,"# MAIN> n_XYZ=",univ%configurations(1)%cells%molecules(1)%constraints%n_XYZ
    

    
    
    do i=1,3
       print *,(univ%configurations(1)%cells%L(i,j),j=1,3)
    end do
    print *,"# celldm(1)=",univ%configurations(1)%cells%L(1,1)/a0
    print *,"# celldm(3)=",univ%configurations(1)%cells%L(3,3)/univ%configurations(1)%cells%L(1,1)
    
    
    print *,"# n_atoms in new=",univ%configurations(1)%cells%molecules(1)%n_atoms
    print *,"# PERIODICITY= ",univ%configurations(iconf)%cells%periodicity
    univ%configurations(iconf)%cells%periodicity="XYZ"
    univ%configurations(iconf)%cells%periodicity="XY"


    
    


    iconf=1 ; call IO_write(param,univ,iconf)

    !iconf=1;     param%output%name="essai.xsf" ;   call save(param%output%name,univ,iconf)
    
    
  end subroutine Ti
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !             subroutine TiOH()
  !   construit un substrat et ajoute une molecule de OH dessus (si ads=.TRUE.
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine TiOH(param,idx_inputs)
    use global
    use Univers
    use Atom
    use IO
    use Cell
    implicit none
    type(t_Param)::param
    integer::idx_inputs
    type(t_Cell)::Ti_unit_cell
    type(t_Univers)::univ
    integer::ibrav,i,j,nat(3),iat1,iat2  
    double precision::a,csura,evacc,q(3)
    double precision,parameter::ZERO=0.0
    !type(t_File)::file
    type(t_Atom)::atm
    logical::ads
    character(len=2)::elt
    type(t_Molecule)::rlx
    integer::idx_conf_to_read
    
    univ=UNIVERS_new()
    call UNIVERS_add_configuration(univ)

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !           lowest part of the slab
    !
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ibrav=4 ! Hexagonal and Trigonal P
    elt="Ti"
    a=2.95
    csura=sqrt(6.0)/3.0
    Ti_unit_cell=CELL_set_unit_cell(ibrav,a,csura,elt)

    ads=.FALSE.
    !ads=.TRUE.
    
    !nat=(/12,12,1/) ; iat1=197 ; iat2=234
    nat=(/10,10,1/) ; iat1=134 ; iat2=165
    !nat=(/8,8,1/) ; iat1=83 ; iat2=108
    !nat=(/6,6,1/) ; iat1=44 ; iat2=63
    !nat=(/4,4,1/) ; iat1=17 ; iat2=30
    evacc=11.0/nat(3)
    !evacc=1.0/nat(3)
    univ%configurations(1)%cells=CELL_slab(ibrav,a,csura,nat,Ti_unit_cell,evacc)
    print *,"# MAIN> n_XYZ=",univ%configurations(1)%cells%molecules(1)%constraints%n_XYZ
    
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! lit les coordonnées à partir d'un fichier xyz
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    param%input%name="TiOH_10_10_1-pos-1.xyz"
    !   param%type='multixyz'
    idx_conf_to_read=1
    rlx=IO_read_molecule(param%input(idx_inputs))
    q=(/ZERO,ZERO,univ%configurations(1)%cells%L(3,3)/evacc/)
    do i=1,rlx%n_atoms
       call ATOM_translate(rlx%atoms(i),q)
       call MOLECULE_add_atom(univ%configurations(1)%cells%molecules(1),rlx%atoms(i))
       univ%configurations(1)%cells%n_atoms=univ%configurations(1)%cells%n_atoms+1
    end do
    univ%configurations(1)%cells%L(3,3)=univ%configurations(1)%cells%L(3,3)*(1.0+1.0/evacc)
    
    do i=1,3
       print *,(univ%configurations(1)%cells%L(i,j),j=1,3)
    end do
    print *,"celldm(1)=",univ%configurations(1)%cells%L(1,1)/a0
    print *,"celldm(3)=",univ%configurations(1)%cells%L(3,3)/univ%configurations(1)%cells%L(1,1)
    
    
    print *,"# n_atoms in new=",univ%configurations(1)%cells%molecules(1)%n_atoms
    !print *,allocated(univ%configurations(1)%cells%molecules(1)%atoms)
    !do i=1,univ%configurations(1)%cells%molecules(1)%n_atoms
    !   print *,univ%configurations(1)%cells%molecules(1)%atoms(i)%elt,univ%configurations(1)%cells%molecules(1)%atoms(i)%Zato
    !end do
    
    ! ! pour créer uner surface, on ajuste univ%configurations(iconf)%cells%L(3,j) :
    ! !     surf=1 -> pas de surface
    ! !     surf=2 -> couche de vide de la même épaisseur que la couche de substrat
    ! surf=9
    ! univ%configurations(1)%cells%L(3,3)=   surf*univ%configurations(1)%cells%L(3,3)
    ! !univ%configurations(iconf)%cells%periodicity="XYZ"
    
    
    if(ads) then
       atm=ATOM_new()
       atm%elt="O"
       atm%q(1)=univ%configurations(1)%cells%molecules(1)%atoms(iat1)%q(1)
       atm%q(2)=univ%configurations(1)%cells%molecules(1)%atoms(iat1)%q(2)
       atm%q(3)=2*univ%configurations(1)%cells%molecules(1)%atoms(iat2)%q(3)-&
            univ%configurations(1)%cells%molecules(1)%atoms(iat1)%q(3)
       call MOLECULE_add_atom(univ%configurations(1)%cells%molecules(1),atm)
       univ%configurations(1)%cells%n_atoms=univ%configurations(1)%cells%n_atoms+1
       atm%elt="H"
       atm%q(1)=univ%configurations(1)%cells%molecules(1)%atoms(iat1)%q(1)
       atm%q(2)=univ%configurations(1)%cells%molecules(1)%atoms(iat1)%q(2)
       atm%q(3)=2*univ%configurations(1)%cells%molecules(1)%atoms(iat2)%q(3)-&
            univ%configurations(1)%cells%molecules(1)%atoms(iat1)%q(3)+1.0
       call MOLECULE_add_atom(univ%configurations(1)%cells%molecules(1),atm)
       univ%configurations(1)%cells%n_atoms=univ%configurations(1)%cells%n_atoms+1
    end if
    
    
    !i=1;    param%output%name="essai.xsf" ;   call save(param%output%name,univ,i)
    
    param%output%name="essai.in"
    param%output%type="qe_input"
    call IO_write(param,univ,i)
    param%output%name="cp2k.in"
    param%output%type="cp2k_input"
    param%title="TiOH_10_10_2"
    param%machine="jeanzay"
    param%calculation%run_type="GEO_OPT"
    call IO_write(param,univ,i)
    
    print *,(univ%configurations(1)%cells%molecules(1)%atoms(134)%q(j),j=1,3)
    print *,(univ%configurations(1)%cells%molecules(1)%atoms(165)%q(j),j=1,3)
    print *,"n_elements=",univ%configurations(1)%cells%molecules(1)%n_elements
    print *,(univ%configurations(1)%cells%molecules(1)%elements(i)%elt,i=1,univ%configurations(1)%cells%molecules(1)%n_elements)
    
  end subroutine TiOH



  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !             subroutine TiOH_addOH()
  !  * Lit les positions à partir d'un fichier cp2k-restart ou xyz  
  !  * ajoute une molecule de OH
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! UNIVERS_add_configuration(univ)
  ! IO_read_cell_from_cp2k_restart_file
  !
  subroutine TiOH_addOH(param,idx_input)
    use global
    use Univers
    use Atom
    use IO
    use Cell
    implicit none
    integer::idx_input
    type(t_Param)::param
    !type(t_File)::file
    type(t_Atom)::atm
    type(t_Univers)::univ
    !integer::idx_conf_to_read
    integer::iat,i,k
    
    univ=UNIVERS_new()
    call UNIVERS_add_configuration(univ)

    !-----------------------------------------------------------------------------------------
    !
    ! lit les coordonnées à partir d'un fichier xyz
    !
    !param%name="/home/bulou/hpc-login/iteney/2022/cp2k/TiOH_10_10_1/TiOH_10_10_1-pos-1.xyz"
    !param%type='xyz'
    !idx_conf_to_read=1
    !call IO_read_molecule(param,molecule_rlx,idx_conf_to_read)
    !
    !-----------------------------------------------------------------------------------------
    !
    ! lit les coordonnées à partir d'un fichier cp2k restart
    !
    !param%name="/home/bulou/hpc-login/iteney/2022/cp2k/TiOH_10_10_1/TiOH_10_10_1-1.restart"
    !param%inputname="TiOH_10_10_1_pos2-1.restart"
    !param%type='cp2k_restart_file'
    call IO_read_cell_from_cp2k_restart_file(param,param%input(idx_input),univ%configurations(1)%cells)

    do i=1,3
       print *,(univ%configurations(1)%cells%L(i,k),k=1,3)
    end do
    print *,"n_atoms in ",trim(param%input(idx_input)%name),univ%configurations(1)%cells%n_atoms,&
         univ%configurations(1)%cells%molecules(1)%n_atoms

    !O    7.3750124043822654E+00    1.1070690013854406E+01    5.5494521947615292E+00
    !H    7.3750143701077500E+00    1.1070687642823296E+01    6.5342572500751057E+00


        
    atm=ATOM_new()
    atm%elt="O"
    atm%q(1)=7.3750124043822654E+00    
    atm%q(2)=1.1070690013854406E+01+5.06/3
    atm%q(3)=5.5494521947615292E+00
    PeriodicTable(8)%Basis_set="DZVP-MOLOPT-SR-GTH"
    PeriodicTable(8)%Potential_set="GTH-PBE-q6"

    call MOLECULE_add_atom(univ%configurations(1)%cells%molecules(1),atm)
    univ%configurations(1)%cells%n_atoms=univ%configurations(1)%cells%n_atoms+1
    
    atm%elt="H"
    atm%q(1)=7.3750143701077500E+00    
    atm%q(2)=1.1070687642823296E+01+5.06/3
    atm%q(3)=6.5342572500751057E+00
    PeriodicTable(1)%Basis_set="DZVP-MOLOPT-GTH"
    PeriodicTable(1)%Potential_set="GTH-PBE-q1"
    call MOLECULE_add_atom(univ%configurations(1)%cells%molecules(1),atm)
    univ%configurations(1)%cells%n_atoms=univ%configurations(1)%cells%n_atoms+1
    call MOLECULE_element_list_update(univ%configurations(1)%cells%molecules(1))


    !i=1
    !param%outputname="essai.xsf" ;
    !call save(param%output%name,univ,i)

    
    !param%output%name="TiOH_10_10_1_pos2_pos1.in"
    !param%output%type="cp2k_input"
    param%title="TiOH_10_10_1_pos2_pos1"
    !param%machine="jeanzay"
    !param%machine="hpc"
    param%calculation%run_type="GEO_OPT"
    i=1


    print *,param%n_outputs
    call IO_write(param,univ,i)


  end subroutine TiOH_addOH
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !             subroutine TiOH_moveOH()
  !  * Lit les positions à partir d'un fichier xyz  
  !  * translate la molecule de OH d'un vecteur t=(tx,ty,tz)
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine TiOH_moveOH(param,idx_input)
    use global
    use Univers
    use Atom
    use IO
    use Cell
    implicit none
    integer::idx_input
    type(t_Param)::param
    !type(t_File)::file
    type(t_Univers)::univ
    !integer::idx_conf_to_read
    integer::iat,i,k
    double precision::t(3)
    
    
    univ=UNIVERS_new()
    call UNIVERS_add_configuration(univ)
    

    t=(/0.0,1.70,0.0/)
    !
    ! lit les coordonnées à partir d'un fichier xyz
    !
    !param%inputname="/home/bulou/hpc-login/iteney/2022/cp2k/TiOH_10_10_1/TiOH_10_10_1-pos-1.xyz"
    !param%inputname="/home/bulou/hpc-login/iteney/2022/cp2k/TiOH_10_10_1/TiOH_10_10_1-1.restart"
    param%input(idx_input)%name="TiOH_10_10_1-1.restart"
    param%input(idx_input)%type='cp2k_restart_file'
    !param%type='xyz'
    !idx_conf_to_read=1
    !call IO_read_molecule(param,molecule_rlx,idx_conf_to_read)

    call IO_read_cell_from_cp2k_restart_file(param,param%input(1),univ%configurations(1)%cells)

    do i=1,3
       print *,(univ%configurations(1)%cells%L(i,k),k=1,3)
    end do
    
    do iat=1,univ%configurations(1)%cells%molecules(1)%n_atoms
       if(univ%configurations(1)%cells%molecules(1)%atoms(iat)%elt.eq."Ti") then
          univ%configurations(1)%cells%molecules(1)%atoms(iat)%idx_molecule=1
       else
          univ%configurations(1)%cells%molecules(1)%atoms(iat)%idx_molecule=2
          do k=1,3
             univ%configurations(1)%cells%molecules(1)%atoms(iat)%q(k)=&
                  univ%configurations(1)%cells%molecules(1)%atoms(iat)%q(k)+t(k)
          end do
       end if
       print *,univ%configurations(1)%cells%molecules(1)%atoms(iat)%elt,&
            univ%configurations(1)%cells%molecules(1)%atoms(iat)%idx_molecule,&
            (univ%configurations(1)%cells%molecules(1)%atoms(iat)%q(k),k=1,3)
    end do
    



    
    param%output%name="TiOH_10_10_1_pos2.in"
    param%output%type="cp2k_input"
    param%title="TiOH_10_10_1_pos2"
    !param%machine="jeanzay"
    !param%machine="hpc"
    param%calculation%run_type="GEO_OPT"
    i=1
    call IO_write(param,univ,i)



  end subroutine TiOH_moveOH
end module scripts
