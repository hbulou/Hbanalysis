! Force Field
! N2 : modèle N2 à trois sites de charge. Dans ce modèle, il y a deux sites d'azote à interaction de Lennard-Jones et une charge ponctuelle sans masse placée
! au milieu de l'azote. La masse de la charge ponctuelle/atome sans masse (au centre de N2) est fixée à 0,000000001.
!
! pair_style lj/cut/coul/long 12 11 
! pair_coeff  1 1 0 0
! pair_coeff  1 2 0 1.655
! pair_coeff  2 2 0.0742 3.31
!
!
! bond_style harmonic
! bond_coeff 1 420 0.547
! angle_style harmonic
! angle_coeff 1 147.7 180
! kspace_style pppm 0.0001

  ! COMPASS lawrence ->  /home/bulou/src/LAMMPS/lammps-20Sep2021/tools/msi2lmp/frc_files/compass_published.frc
  !         michelson -> /home/bulou/src/LAMMPS/lammps-29Oct20/tools/msi2lmp/frc_files/compass_published.frc
! n1n  14.00670     N          nitrogen in N2
!#quartic_bond     compass
!
!> E = K2 * (R - R0)^2  +  K3 * (R - R0)^3  +  K4 * (R - R0)^4
!
!Ver  Ref     I     J          R0         K2          K3          K4
!---- ---    ----  ----     -------    --------   ---------    --------

!  1.0   5     n1n   n1n       1.0977   1651.3730  -4069.3178   5984.9629
!#nonbond(9-6)     compass
!
!> E = eps(ij) [2(r(ij)*/r(ij))**9 - 3(r(ij)*/r(ij))**6]
!> where    r(ij) = [(r(i)**6 + r(j)**6))/2]**(1/6)
!>
!>        eps(ij) = 2 sqrt(eps(i) * eps(j)) * 
!>                   r(i)^3 * r(j)^3/[r(i)^6 + r(j)^6]
!
!@combination sixth-power
!@type r-eps
!
!Ver  Ref     I          r          eps 
!---- ---    ----    ---------   ---------
!
!  1.0   5     n1n        3.8008      0.0598
!#bond_increments     compass
!
!Ver  Ref     I     J       DeltaIJ     DeltaJI
!---- ---    ----  ----     -------     -------
!
! 1.0   5     n1n   n1n       0.0000      0.0000

! units metal
! boundary p p p
! atom_style full
! pair_style lj/class2



module Element
  implicit none
  !-----------------------------------------------
  type t_PeriodicTable
     integer::Z
     character(len=2)::elt
     double precision::mass
  end type t_PeriodicTable
  integer,parameter::Nelt=20
  type (t_PeriodicTable)::PeriodicTable(Nelt)
  double precision::deq(Nelt,Nelt)
  !Nelt=20
  !allocate(PeriodicTable(Nelt))

end module Element

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program generate_configuration
  use Element
  implicit none

  type t_Bond
     ! K(r-R0)^n
     double precision::R0
     double precision::K
     integer::n
  end type t_Bond
  type t_Angle
     ! K(r-R0)^n
     double precision::theta0
     double precision::K
     integer::n
  end type t_Angle
  type t_LJ
     double precision::r0
     double precision::eps
  end type t_LJ
  type t_COMPASS_ForceField
     type(t_Bond)::bond(3)
     type(t_Angle)::angle(3)
     type(t_LJ)::vdW
  end type t_COMPASS_ForceField




  !-----------------------------------------------
  type t_File
     character(len=2048)::name
     integer::n_configuration
     integer::n_atom
     integer::n_line
     double precision,allocatable::nrj(:)

  end type t_File
  !-----------------------------------------------
  type t_Atom
     integer::idx
     integer::idx_molecule
     double precision::q(3)
     character(len=2)::elt
     integer::Zato
  end type t_Atom
  !-----------------------------------------------
  type t_Molecule
     integer:: idx
     integer:: n_atom
     type(t_Atom),allocatable::atoms(:)
     double precision::MassCenter(3)
  end type t_Molecule
  !-----------------------------------------------
  type t_EltTable
     integer::idx
     character(len=2)::elt
  end type t_EltTable
  !-----------------------------------------------
  type t_Set
     integer,allocatable::list_atoms(:)
     double precision::MC(3)
  end type t_Set
  !-----------------------------------------------
  ! les groupes sont des regroupements d'atomes (set) soit par espèces chimiques (le groupe des "espèces chimique"),
  ! soit par proximité (le groupe des "atomes voisins"), ...
  type t_Group
     integer::nelt
     type(t_EltTable),allocatable::elt_table(:)    
     integer:: n_set                             ! nombre de set dans le group

     integer,allocatable:: n_atoms_set(:)         ! nombre d'atomes/set
     type (t_Set),allocatable::set(:)    ! liste des atomes dans les set
  end type t_Group
  !-----------------------------------------------
  type t_Configuration
     integer:: n_molecule
     integer:: n_atom
     type(t_Molecule),allocatable::molecules(:)
     integer::n_element
     type(t_EltTable),allocatable::list_elements(:)
  end type t_Configuration
  !-----------------------------------------------
  type t_Box
     integer:: idx
     integer:: n_configuration                     
     type(t_Configuration),allocatable::configurations(:)
     double precision::cell(3,3)
     integer::n_atom
  end type t_Box


  !type(t_Molecule),allocatable::evol(:)
  type(t_Box)::box
  type(t_Molecule)::Ca,CO3,H2O,molecule,N2,O2,air
  integer::iconf,i,num_conf,nconf
  character (len=1024)::path,filename
  character (len=2048)::trajfile,inputfile,nrjfile
  character(len=22)::fmt,conf
  character(len=32)::run_type
  type (t_Group)::group

  double precision::vec(3)
  integer::idx,j,n
  type(t_File)::file


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  character(len=1024) :: arg 


  integer::io
  character (len=1024)::line
  character (len=32)::field(32)
  integer::nfield

  double precision::alpha,beta,gamma
  character(len=32)::center_tag
  INTEGER :: idum
  logical::new_molecule
  integer::iH2O,nH2O
  call init()



  stop
  call read_COMPASS_parameters()
  stop
  ! ------------------------------------------
  !
  ! read the commande line
  !
  !-------------------------------------------
  !do i = 1, iargc() 
  !   call getarg(i, arg) 
  !   write (6,*) "#argument ",i,trim(arg) 
  !end do
  !call getarg(1,path)
  !print *,trim(path)

  ! ------------------------------------------
  !
  ! Définition des molécules Ca, CO3, et H2O
  ! on lit les coordonnées à partir d'un fichier
  !
  ! ------------------------------------------
  Ca%idx=1
  Ca%n_atom=1
  
  allocate(Ca%atoms(Ca%n_atom))

  CO3%idx=1
  CO3%n_atom=4
  allocate(CO3%atoms(CO3%n_atom))

  open(unit=1,file='CaCO3-rlx.xyz',form='formatted')
  read(1,'(A)',iostat=io) line
  read(1,'(A)',iostat=io) line
  do i=1,4
     read(1,'(A)',iostat=io) line
     call line_parser(line,nfield,field)
     CO3%atoms(i)%elt=field(1) ; call elt2Z(CO3%atoms(i)%elt,CO3%atoms(i)%Zato)
     read(field(2),*) CO3%atoms(i)%q(1)
     read(field(3),*) CO3%atoms(i)%q(2)
     read(field(4),*) CO3%atoms(i)%q(3)
  end do
  read(1,'(A)',iostat=io) line
  call line_parser(line,nfield,field)
  Ca%atoms(1)%elt=field(1) ; call elt2Z(Ca%atoms(1)%elt,Ca%atoms(1)%Zato)
  read(field(2),*) Ca%atoms(1)%q(1)
  read(field(3),*) Ca%atoms(1)%q(2)
  read(field(4),*) Ca%atoms(1)%q(3)
  close(1)


  ! H2O%idx=1
  ! H2O%n_atom=3
  ! allocate(H2O%atoms(H2O%n_atom))
  ! open(unit=1,file='H2O.xyz',form='formatted')
  ! read(1,'(A)',iostat=io) line
  ! read(1,'(A)',iostat=io) line
  ! do i=1,3
  !    read(1,'(A)',iostat=io) line
  !    call line_parser(line,nfield,field)
  !    H2O%atoms(i)%elt=field(1) ; call elt2Z(H2O%atoms(i)%elt,H2O%atoms(i)%Zato)
  !    read(field(2),*) H2O%atoms(i)%q(1)
  !    read(field(3),*) H2O%atoms(i)%q(2)
  !    read(field(4),*) H2O%atoms(i)%q(3)
  ! end do
  ! close(1)


  !call molecule_mass_center(Ca) ;   print *,(Ca%MassCenter(i),i=1,3)
  !call molecule_mass_center(CO3) ;   print *,(CO3%MassCenter(i),i=1,3)
  !call molecule_mass_center(H2O) ;   print *,(H2O%MassCenter(i),i=1,3)

  !vec=(/ 5.0,5.0,2.5 /)
  !call molecule_translate_mass_center_at(H2O,vec) ; print *,(H2O%MassCenter(i),i=1,3)
  !vec=(/ 2.5,2.5,10.0 /)
  !call molecule_translate_mass_center_at(Ca,vec) ; print *,(H2O%MassCenter(i),i=1,3)


  ! ------------------------------------------
  !
  ! Définition de la boite de simulation
  !
  ! ------------------------------------------
  iconf=1

  call box_init(box,nconf=1)


  ! on ajoute les molécules dans la boite
  !
  !call box_add_molecule(box,CO3,1)
  !call box_add_molecule(box,Ca,iconf)

  !inputfile='run-20H2O-gopt-pos-1.xyz'
  inputfile='run-30h2o-gopt-pos-1.xyz'

  file%name=inputfile
  call file_info(file)
  call molecule_from_file(molecule,inputfile,num_conf=200)

  inputfile="H2O.xyz"
  call molecule_from_file(H2O,inputfile,num_conf=1)

  inputfile="N2-rlx.xyz"
  call molecule_from_file(N2,inputfile,num_conf=4)
  !call idx_molecule(molecule,box)
  !call box_add_molecule(box,molecule,iconf)
  iconf=1
  n=10
  !call box_add_multiple_molecules(box,N2,iconf=1,n=8)

  inputfile="O2-rlx.xyz"
  call molecule_from_file(O2,inputfile,num_conf=4)
  !call idx_molecule(molecule,box)
  !call box_add_molecule(box,molecule,iconf)
  iconf=1
  n=10
  !call box_add_multiple_molecules(box,O2,iconf=1,n=2)

  inputfile="air.xyz"
  call molecule_from_file(air,inputfile,num_conf=54)
  call box_add_molecule(box,air,iconf)
  call idx_molecule(air,box)
  call box_get_element_list(box,iconf)
  !print *,"# n_element in box = ",box%configurations(iconf)%n_element
  !print *,(box%configurations(iconf)%list_elements(i)%elt,i=1,box%configurations(iconf)%n_element)
  !run_type='gopt'
  !call write_input_file(box,iconf,run_type)
  fmt='xsf'
  call save(box,fmt,iconf)


contains
  ! --------------------------------------------------------------------------------------
  !
  !              box_add_multiple_molecules()
  !
  ! --------------------------------------------------------------------------------------
  subroutine box_add_multiple_molecules(box,molecule,iconf,n)
    implicit none
    type(t_Box)::box
    type(t_Molecule)::molecule
    logical::new_molecule
    integer::i,n,idum,iconf
    character(len=32)::center_tag
    double precision::vec(3),alpha,beta,gamma
    idum=1  
    i=0
    do while(i.lt.n)
       !call box_add_molecule(box,H2O,iconf)
       vec=(/ box%cell(1,1)*ran(idum),box%cell(2,2)*ran(idum),box%cell(3,3)*ran(idum) /)
       call molecule_translate_mass_center_at(molecule,vec) !; print *,(H2O%MassCenter(i),i=1,3)
       alpha=360*ran(idum)
       beta=360*ran(idum)
       gamma=360*ran(idum)
       center_tag='molecule_MC'
       call molecule_rotate(molecule,alpha,beta,gamma,center_tag)
       call molecule_distance(box,molecule,iconf,new_molecule)
       print *,"# New molecule= ",new_molecule
       if(new_molecule) then
          call box_add_molecule(box,molecule,iconf)
          i=i+1
       end if
    end do
  end subroutine box_add_multiple_molecules

  ! --------------------------------------------------------------------------------------
  !
  !              box_init()
  !
  ! --------------------------------------------------------------------------------------
  subroutine box_init(box,nconf)
    implicit none
    type(t_Box)::box
    integer::nconf,iconf,i
    box%n_atom=0
    box%cell=0.0
    do i=1,3
       box%cell(i,i)=20.0
    end do
    box%n_configuration=nconf
    allocate(box%configurations(box%n_configuration))
    do i=1,box%n_configuration
       box%configurations(1)%n_molecule=0
       box%configurations(1)%n_element=0
    end do
  end subroutine box_init
  ! --------------------------------------------------------------------------------------
  !
  !              idx_molecule()
  !
  ! --------------------------------------------------------------------------------------
  subroutine idx_molecule(molecule,box)
    implicit none
    type(t_Box)::box
    character(len=32),parameter::name_sub='idx_molecule'
    type(t_Molecule)::molecule
    integer::i,j,k,nOH,Zi,Zj,idx_molecule_max
    double precision::R(3),d
    write(*,'(A,A,A,I5)') "# ",trim(name_sub)," > n_atom=",molecule%n_atom
    nOH=0
    idx_molecule_max=0
    do i=1,molecule%n_atom-1
       Zi=molecule%atoms(i)%Zato
       do j=i+1,molecule%n_atom
          Zj=molecule%atoms(j)%Zato
          d=0.0
          do k=1,3
             R(k)=molecule%atoms(j)%q(k)-molecule%atoms(i)%q(k)
             R(k)=modulo(R(k),box%cell(k,k))
             if(R(k).gt.0.5*box%cell(k,k)) then
                R(k)=R(k)-box%cell(k,k)
             end if
             if(-R(k).gt.0.5*box%cell(k,k)) then
                R(k)=R(k)+box%cell(k,k)
             end if
             d=d+R(k)*R(k)
          end do
          d=sqrt(d)
          if(d<deq(Zi,Zj))  then
             write(*,'(A,A,A,I5,I5,A,F12.8,A,A)') &
                  "# ",trim(name_sub)," > d(",i,j,")=",d,molecule%atoms(j)%elt,molecule%atoms(i)%elt
             nOH=nOH+1
             if((molecule%atoms(i)%idx_molecule.eq.-1).and.(molecule%atoms(j)%idx_molecule.eq.-1)) then
                idx_molecule_max=idx_molecule_max+1
                molecule%atoms(i)%idx_molecule=idx_molecule_max
                molecule%atoms(j)%idx_molecule=idx_molecule_max
             else if((molecule%atoms(i)%idx_molecule.lt.0).and.(molecule%atoms(j)%idx_molecule.eq.-1)) then
                molecule%atoms(j)%idx_molecule=molecule%atoms(i)%idx_molecule
             else if((molecule%atoms(j)%idx_molecule.lt.0).and.(molecule%atoms(j)%idx_molecule.eq.-1)) then
                molecule%atoms(i)%idx_molecule=molecule%atoms(j)%idx_molecule
             else
                print *,"LINE307"
                stop
             end if
          end if


       end do
    end do
    print *,nOH,idx_molecule_max
  end subroutine idx_molecule

  ! --------------------------------------------------------------------------------------
  !
  !              file_info()
  !
  ! --------------------------------------------------------------------------------------

  subroutine file_info(file)
    implicit none
    type(t_File)::file
    integer::io
    character (len=1024)::line
    character (len=32)::field(32)
    integer::nfield,i,j

    io=0
    file%n_line=0
    file%n_configuration=0
    open(unit=1,file=file%name,form='formatted')


    do while (io==0)
       read(1,'(A)',iostat=io) line;       if(io.eq.-1) exit ;        file%n_line=file%n_line+1

       call line_parser(line,nfield,field)
       read(field(1),*) file%n_atom
       read(1,'(A)',iostat=io) line;       if(io.eq.-1) exit ;        file%n_line=file%n_line+1
       file%n_configuration=file%n_configuration+1
       do i=1,file%n_atom
          read(1,'(A)',iostat=io) line;       if(io.eq.-1) exit 
       end do
    enddo


    print *,"# n_line=",file%n_line
    print *,"# n_atom=",file%n_atom
    print *,"# n_configuration=",file%n_configuration

    allocate(file%nrj(file%n_configuration))
    rewind(1)
    do i=1,file%n_configuration
       read(1,'(A)',iostat=io) line
       read(1,'(A)',iostat=io) line
       call line_parser(line,nfield,field)
       read(field(6),*) file%nrj(i)

       do j=1,file%n_atom
          read(1,'(A)',iostat=io) line
       end do
    end do
    close(1)
    open(unit=1,file="nrj.dat",form='formatted',status='unknown')
    do i=1,file%n_configuration
       write(1,*) i,file%nrj(i)
    end do
    close(1)
  end subroutine file_info

  ! --------------------------------------------------------------------------------------
  !
  !              molecule_from_file()
  !
  ! --------------------------------------------------------------------------------------
  subroutine molecule_from_file(molecule,xyzfile,num_conf)
    implicit none
    character(len=32),parameter::name_sub='molecule_from_file'
    character (len=1024)::xyzfile,s
    integer::i,io,natom,iconf,num_conf
    character (len=1024)::line
    character (len=32)::field(32)
    integer::nfield
    type(t_Molecule)::molecule
    open(unit=1,file=xyzfile,form='formatted')


    iconf=-1
    do while(.not.((iconf+1).eq.num_conf))
       !print *,"iconf=",iconf
       read(1,'(A)',iostat=io) line
       call line_parser(line,nfield,field)
       read(field(1),*) natom
       !print *,"# natom=",natom
       read(1,'(A)',iostat=io) line
       call line_parser(line,nfield,field)
       !print *,"#nfield= ",nfield
       if(nfield.eq.0) then
          rewind(1)
          exit
       else
          if(iconf.eq.-1) then
             iconf=1
          else
             iconf=iconf+1
          end if
          !print *,iconf,natom
          !read(field(3),*) iconf
          do i=1,natom
             read(1,'(A)',iostat=io) line
          end do

       end if
    end do
    print *,'-----------------------------------------------------------------------------'

    read(1,'(A)',iostat=io) line
    call line_parser(line,nfield,field)
    print *,trim(line)
    read(field(1),*) natom
    molecule%idx=1
    molecule%n_atom=natom
    allocate(molecule%atoms(molecule%n_atom))


    read(1,'(A)',iostat=io) line
    do i=1,natom
       read(1,'(A)',iostat=io) line
       call line_parser(line,nfield,field)
       molecule%atoms(i)%elt=field(1)
       call elt2Z(molecule%atoms(i)%elt,molecule%atoms(i)%Zato)
       read(field(2),*) molecule%atoms(i)%q(1)
       read(field(3),*) molecule%atoms(i)%q(2)
       read(field(4),*) molecule%atoms(i)%q(3)
       molecule%atoms(i)%idx_molecule=-1
       write(*,'(A,A,A,A)') "#",trim(name_sub),"> ",trim(line)
    end do
    print *,'-----------------------------------------------------------------------------'
    close(1)

  end subroutine molecule_from_file

  ! --------------------------------------------------------------------------------------
  !
  !              molecule_distance()
  !
  ! --------------------------------------------------------------------------------------
  subroutine molecule_distance(box,molecule,iconf,new_molecule)
    implicit none
    type(t_Box)::box
    type(t_Molecule)::molecule
    integer::iconf,jat,imol,k,iat
    logical::new_molecule
    double precision::R(3),d,dlim=2.0
    new_molecule=.True.
    imol=0
    do while((imol.lt.box%configurations(iconf)%n_molecule).and.new_molecule)
       imol=imol+1
       jat=0
       do while((jat.lt.box%configurations(iconf)%molecules(imol)%n_atom).and.new_molecule)
          jat=jat+1
          iat=0
          do while((iat.lt.molecule%n_atom).and.new_molecule)
             iat=iat+1
             d=0.0
             do k=1,3
                R(k)=box%configurations(iconf)%molecules(imol)%atoms(jat)%q(k)-molecule%atoms(iat)%q(k)
                R(k)=modulo(R(k),box%cell(k,k))
                if(R(k).gt.0.5*box%cell(k,k)) then
                   R(k)=R(k)-box%cell(k,k)
                end if
                if(-R(k).gt.0.5*box%cell(k,k)) then
                   R(k)=R(k)+box%cell(k,k)
                end if
                d=d+R(k)*R(k)
             end do
             d=sqrt(d)
             !print *,"d=",d
             if(d.lt.dlim) then
                new_molecule=.False.
             end if
          end do
       end do
    end do
    print *,new_molecule
  end subroutine molecule_distance
  ! --------------------------------------------------------------------------------------
  !
  !              box_get_element_list()
  !
  ! --------------------------------------------------------------------------------------
  subroutine box_get_element_list(box,iconf)
    implicit none
    type(t_Box)::box
    integer::iconf,imol,jat,ielt
    logical::new_elt
    type(t_EltTable),allocatable::temp(:)
    print *,"# Entering box_get_element_list() ..."
    do imol=1,box%configurations(iconf)%n_molecule
       do jat=1,box%configurations(iconf)%molecules(imol)%n_atom
          new_elt=.True.
          ielt=0
          do while(new_elt.and.(ielt.lt.box%configurations(iconf)%n_element))
             ielt=ielt+1
             if(box%configurations(iconf)%molecules(imol)%atoms(jat)%elt.eq.box%configurations(iconf)%list_elements(ielt)%elt) then
                new_elt=.False.
             end if
          end do
          if(new_elt) then
             box%configurations(iconf)%n_element=box%configurations(iconf)%n_element+1
             if(box%configurations(iconf)%n_element.eq.1) then
                allocate(box%configurations(iconf)%list_elements(box%configurations(iconf)%n_element))
                box%configurations(iconf)%list_elements(box%configurations(iconf)%n_element)%elt=&
                     box%configurations(iconf)%molecules(imol)%atoms(jat)%elt
             else
                allocate(temp(box%configurations(iconf)%n_element))
                temp(1:box%configurations(iconf)%n_element)=&
                     box%configurations(iconf)%list_elements(1:box%configurations(iconf)%n_element-1)
                temp(box%configurations(iconf)%n_element)%elt=box%configurations(iconf)%molecules(imol)%atoms(jat)%elt
                deallocate(box%configurations(iconf)%list_elements)
                allocate(box%configurations(iconf)%list_elements(box%configurations(iconf)%n_element))
                box%configurations(iconf)%list_elements=temp
                deallocate(temp)
                !call move_alloc(from=temp,to=box%configurations(iconf)%list_elements)
             end if
          end if
       end do
    end do

    print *,"# ... End of box_get_element_list()."
  end subroutine box_get_element_list
  ! --------------------------------------------------------------------------------------
  !
  !              box_add_molecule()
  !
  ! --------------------------------------------------------------------------------------
  subroutine box_add_molecule(box,molecule,iconf)
    implicit none
    integer::iconf,i,j
    type(t_Box)::box
    type(t_Molecule)::molecule
    type(t_Molecule),allocatable::temp(:)
    print *,"-----------------------------------------------------------------------------"
    print *,"# Entering box_add_molecule..."
    box%configurations(iconf)%n_atom=box%configurations(iconf)%n_atom+molecule%n_atom
    box%configurations(iconf)%n_molecule=box%configurations(iconf)%n_molecule+1

    if(box%configurations(iconf)%n_molecule.eq.1) then
       allocate(box%configurations(iconf)%molecules(1))
       box%configurations(iconf)%molecules(1)=molecule
    else
       allocate(temp(box%configurations(iconf)%n_molecule))
       temp(1:box%configurations(iconf)%n_molecule-1)=&
            box%configurations(iconf)%molecules(1:box%configurations(iconf)%n_molecule-1)
       temp(box%configurations(iconf)%n_molecule)=molecule
       deallocate(box%configurations(iconf)%molecules)
       allocate(box%configurations(iconf)%molecules(box%configurations(iconf)%n_molecule))
       box%configurations(iconf)%molecules=temp
       deallocate(temp)
       !       call move_alloc(from=temp,to=box(1)%configurations(iconf)%molecules)
    end if
    print *,"# Number of molecule(s) in box:",box%configurations(iconf)%n_molecule
    print *,"# End of box_add_molecule()"
    print *,"-----------------------------------------------------------------------------"
  end subroutine box_add_molecule
  ! --------------------------------------------------------------------------------------
  !
  !              molecule_mass_center()
  !
  ! --------------------------------------------------------------------------------------
  subroutine molecule_mass_center(molecule)
    use Element
    implicit none

    type(t_Molecule)::molecule
    integer::iat,k
    double precision::m,mtot
    do k=1,3
       molecule%MassCenter(k)=0.0
    end do
    do iat=1,molecule%n_atom
       m=PeriodicTable(molecule%atoms(iat)%Zato)%mass
       !print *,m
       mtot=mtot+m
       do k=1,3
          !print *,molecule%atoms(iat)%q(k)
          molecule%MassCenter(k)=molecule%MassCenter(k)+m*molecule%atoms(iat)%q(k)
       end do
    end do
    do k=1,3
       molecule%MassCenter(k)=molecule%MassCenter(k)/mtot
    end do
  end subroutine molecule_mass_center


  FUNCTION ran(idum)
    IMPLICIT NONE
    INTEGER, PARAMETER :: K4B=selected_int_kind(9)
    INTEGER(K4B), INTENT(INOUT) :: idum
    REAL :: ran
    ! “Minimal” random number generator of Park and Miller combined with a Marsaglia shift
    ! sequence. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
    ! values). This fully portable, scalar generator has the “traditional” (not Fortran 90) calling
    ! sequence with a random deviate as the returned function value: call with idum a negative
    ! integer to initialize; thereafter, do not alter idum except to reinitialize. The period of this
    ! generator is about 3.1 × 1018.
    INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
    REAL, SAVE :: am
    INTEGER(K4B), SAVE :: ix=-1,iy=-1,k
    if (idum <= 0 .or. iy < 0) then            ! Initialize.
       am=nearest(1.0,-1.0)/IM
       iy=ior(ieor(888889999,abs(idum)),1)
       ix=ieor(777755555,abs(idum))
       idum=abs(idum)+1                        ! Set idum positive.
    end if
    ix=ieor(ix,ishft(ix,13))                   ! Marsaglia shift sequence with period 232 − 1.
    ix=ieor(ix,ishft(ix,-17))
    ix=ieor(ix,ishft(ix,5))
    k=iy/IQ                                   ! Park-Miller sequence by Schrage’s method,
    iy=IA*(iy-k*IQ)-IR*k                      ! period 231 − 2.
    if (iy < 0) iy=iy+IM
    ran=am*ior(iand(IM,ieor(ix,iy)),1)        ! Combine the two generators with masking to
    !  ensure nonzero value.
  end FUNCTION ran

  ! --------------------------------------------------------------------------------------
  !
  !              molecule_rotate()
  !
  ! --------------------------------------------------------------------------------------
  subroutine molecule_rotate(molecule,alpha,beta,gamma,center_tag)
    implicit none
    type(t_Molecule)::molecule
    double precision::alpha,beta,gamma,pi,center(3),newq(3),vec(3)
    character(len=32)::center_tag
    integer::iat

    pi=4.D0*DATAN(1.D0)
    alpha=alpha*pi/180
    beta=beta*pi/180
    gamma=gamma*pi/180

    if(center_tag.eq.'molecule_MC') then
       call molecule_mass_center(molecule)
       center=molecule%MassCenter
       print *,'center=',center
    end if

    !vec=(/0.0,0.0,0.0/)
    !call molecule_translate_mass_center_at(molecule,vec) 
    do iat=1,molecule%n_atom
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
       ! newq(1)=molecule%atoms(iat)%q(1)*cos(gamma)-molecule%atoms(iat)%q(2)*sin(gamma)
       ! newq(2)=molecule%atoms(iat)%q(1)*sin(gamma)+molecule%atoms(iat)%q(2)*cos(gamma)
       ! molecule%atoms(iat)%q(1)=newq(1)
       ! molecule%atoms(iat)%q(2)=newq(2)
       ! newq(1)=molecule%atoms(iat)%q(1)*cos(beta)-molecule%atoms(iat)%q(3)*sin(beta)
       ! newq(3)=molecule%atoms(iat)%q(1)*sin(beta)+molecule%atoms(iat)%q(3)*cos(beta)
       ! molecule%atoms(iat)%q(1)=newq(1)
       ! molecule%atoms(iat)%q(3)=newq(3)
       ! newq(2)=molecule%atoms(iat)%q(2)*cos(alpha)-molecule%atoms(iat)%q(3)*sin(alpha)
       ! newq(3)=molecule%atoms(iat)%q(2)*sin(alpha)+molecule%atoms(iat)%q(3)*cos(alpha)
       ! molecule%atoms(iat)%q(2)=newq(2)
       ! molecule%atoms(iat)%q(3)=newq(3)
    end do
    !call molecule_translate_mass_center_at(molecule,center) 
  end subroutine molecule_rotate
  ! --------------------------------------------------------------------------------------
  !
  !              molecule_mass_center()
  !
  ! --------------------------------------------------------------------------------------
  subroutine molecule_translate_mass_center_at(molecule,vec)
    use Element
    implicit none
    double precision::vec(3)
    type(t_Molecule)::molecule
    integer::iat,k

    do iat=1,molecule%n_atom
       do k=1,3
          molecule%atoms(iat)%q(k)=molecule%atoms(iat)%q(k)-molecule%MassCenter(k)+vec(k)
       end do
    end do
    call molecule_mass_center(molecule)
  end subroutine molecule_translate_mass_center_at
  ! --------------------------------------------------------------------------------------
  !
  !              read_input_file()
  !
  ! --------------------------------------------------------------------------------------
  subroutine read_input_file(filename)
    implicit none
    character(len=1024)::filename
    print *,trim(filename)
  end subroutine read_input_file
  ! --------------------------------------------------------------------------------------
  !
  !              translate()
  !
  ! --------------------------------------------------------------------------------------
  subroutine translate(box,alpha)
    implicit none
    integer::imol,iconf,k,jat
    double precision::vec(3)
    type (t_Box)::box
    double precision::alpha




    vec(1)=alpha*(box%cell(1,1)+box%cell(2,1))
    vec(2)=alpha*(box%cell(2,1)+box%cell(2,2))
    vec(3)=0.0
    do iconf=1,box%n_configuration
       do imol=1,box%configurations(iconf)%n_molecule
          do jat=1,box%configurations(iconf)%molecules(imol)%n_atom
             do k=1,3
                box%configurations(iconf)%molecules(imol)%atoms(jat)%q(k)=&
                     modulo(box%configurations(iconf)%molecules(imol)%atoms(jat)%q(k)+vec(k),box%cell(k,k))
             end do
          end do
       end do
    end do
  end subroutine translate
  ! --------------------------------------------------------------------------------------
  !
  !              center_at()
  !
  ! --------------------------------------------------------------------------------------
  subroutine center_at(box,idx)
    implicit none
    type (t_Box)::box
    integer::imol,iconf,jat,k,idx
    double precision::vec(3)
    do k=1,3
       vec(k)=0.5*box%cell(k,k)
    end do
    do iconf=1,box%n_configuration
       do imol=1,box%configurations(iconf)%n_molecule
          do jat=1,box%configurations(iconf)%molecules(imol)%n_atom
             do k=1,3
                box%configurations(iconf)%molecules(imol)%atoms(jat)%q(k)=&
                     box%configurations(iconf)%molecules(imol)%atoms(jat)%q(k)-&
                     box%configurations(iconf)%molecules(imol)%atoms(idx)%q(k)+vec(k)
                if(box%configurations(iconf)%molecules(imol)%atoms(jat)%q(k).lt.0.0) then
                   box%configurations(iconf)%molecules(imol)%atoms(jat)%q(k)=&
                        box%configurations(iconf)%molecules(imol)%atoms(jat)%q(k)+box%cell(k,k)
                end if
                if(box%configurations(iconf)%molecules(imol)%atoms(jat)%q(k).gt.box%cell(k,k)) then
                   box%configurations(iconf)%molecules(imol)%atoms(jat)%q(k)=&
                        box%configurations(iconf)%molecules(imol)%atoms(jat)%q(k)-box%cell(k,k)
                end if

             end do
          end do
       end do
    end do


  end subroutine center_at

  ! --------------------------------------------------------------------------------------
  !
  !              transformation()
  !
  ! --------------------------------------------------------------------------------------
  subroutine transformation(box,iconf)
    implicit none
    integer::jat,iconf,imol
    double precision::vec(3)
    type (t_Box)::box
    double precision::alpha

    alpha=0.15

    vec(1)=alpha*(box%cell(1,1)+box%cell(2,1))
    vec(2)=alpha*(box%cell(1,2)+box%cell(2,2))
    vec(3)=0.0
    do imol=1,box%configurations(iconf)%n_molecule
       do jat=1,box%configurations(iconf)%molecules(imol)%n_atom
          if(box%configurations(iconf)%molecules(imol)%atoms(jat)%elt.eq.'Ca') then
             box%configurations(iconf)%molecules(imol)%atoms(jat)%q(1)=&
                  box%configurations(iconf)%molecules(imol)%atoms(jat)%q(1)+vec(1)
             box%configurations(iconf)%molecules(imol)%atoms(jat)%q(2)=&
                  box%configurations(iconf)%molecules(imol)%atoms(jat)%q(2)+vec(2)
          else
             box%configurations(iconf)%molecules(imol)%atoms(jat)%q(1)=&
                  box%configurations(iconf)%molecules(imol)%atoms(jat)%q(1)-vec(1)
             box%configurations(iconf)%molecules(imol)%atoms(jat)%q(2)=&
                  box%configurations(iconf)%molecules(imol)%atoms(jat)%q(2)-vec(2)
          end if
       end do
    end do
  end subroutine transformation


  ! --------------------------------------------------------------------------------------
  !
  !              read_COMPASS_parameters()
  !
  ! --------------------------------------------------------------------------------------
  subroutine read_COMPASS_parameters()
    implicit none
    character (len=1024)::inputfile
    integer::i,io
    character (len=1024)::line
    character (len=32)::field(32),ess
    integer::nfield
    logical::quartic_bond,found,quartic_angle,non_bond
    logical::atom_type
    type(t_COMPASS_ForceField)::N2
    character(len=4)::elt
    elt='n1n'
    !elt='p4='
    ! #quartic_bond     compass
    !
    ! > E = K2 * (R - R0)^2  +  K3 * (R - R0)^3  +  K4 * (R - R0)^4
    !
    ! !Ver  Ref     I     J          R0         K2          K3          K4
    ! !---- ---    ----  ----     -------    --------   ---------    --------
    quartic_bond=.FALSE.
    quartic_angle=.FALSE.
    non_bond=.FALSE.
    atom_type=.FALSE.

    inputfile='/home/bulou/src/LAMMPS/lammps-20Sep2021/tools/msi2lmp/frc_files/compass_published.frc' ! lawrence
    !inputfile='/home/bulou/src/LAMMPS/lammps-29Oct20/tools/msi2lmp/frc_files/compass_published.frc'   ! michelson
    open(unit=1,file=inputfile,form='formatted')
    i=0
    io=0
    do while(io==0)
       read(1,'(A)',iostat=io) line
       call line_parser(line,nfield,field)
       if(nfield>0) then
          if(field(1)(1:1).eq."#") then
             if(field(1).eq."#atom_types") then
                atom_type=.TRUE.
             else
                atom_type=.FALSE.
             end if
          end if
          if(atom_type) print *,trim(line)
       end if
    end do
    rewind(1)

    
    i=0
    io=0
    do while(io==0)
       read(1,'(A)',iostat=io) line
       call line_parser(line,nfield,field)
       if(nfield>0) then
          if(field(1)(1:1).eq."#") then
             print *,trim(line)
          end if
          found=.FALSE.
          i=0
          do while((i.lt.nfield).and.(.not.found))
             i=i+1
             if(field(i).eq.elt) found=.TRUE.
          end do
          if(found) print *,trim(line)
       end if
    end do
    rewind(1)

    io=0
    i=0
    do while(io==0)
       read(1,'(A)',iostat=io) line
       call line_parser(line,nfield,field)
       
       if(nfield>0) then
          if(field(1)(1:1).eq.'#') then          !print *,trim(line)
             if(field(1).eq.'#quartic_bond') then
                quartic_bond=.TRUE.  !        print *,field(1)
             else
                quartic_bond=.FALSE.
             end if
             if(field(1).eq.'#quartic_angle') then
                quartic_angle=.TRUE.  !        print *,field(1)
             else
                quartic_angle=.FALSE.
             end if
             if(field(1).eq.'#nonbond(9-6)') then
                non_bond=.TRUE.  !        print *,field(1)
             else
                non_bond=.FALSE.
             end if
          end if
          if(quartic_bond) then
             found=.FALSE.
             i=0
             do while((.not.found).and.(i.lt.nfield))
                i=i+1
                if(field(i).eq.elt) found=.TRUE.
             end do
             if(found) then
                print *,trim(line)
                read(field(5),*) N2%bond(1)%R0
                N2%bond(2)%R0=N2%bond(1)%R0
                N2%bond(3)%R0=N2%bond(1)%R0
                read(field(6),*) N2%bond(1)%K
                read(field(7),*) N2%bond(2)%K
                read(field(8),*) N2%bond(3)%K
                N2%bond(1)%n=2
                N2%bond(2)%n=3
                N2%bond(3)%n=4
                
             end if !found
          end if ! quartic_bond
          if(quartic_angle) then
             found=.FALSE.
             i=0
             do while((.not.found).and.(i.lt.nfield))
                i=i+1
                if(field(i).eq.elt) found=.TRUE.
             end do
             if(found) then
                print *,trim(line)
                read(field(6),*) N2%angle(1)%theta0
                N2%angle(2)%theta0=N2%angle(1)%theta0
                N2%angle(3)%theta0=N2%angle(1)%theta0
                read(field(7),*) N2%angle(1)%K
                read(field(8),*) N2%angle(2)%K
                read(field(8),*) N2%angle(3)%K
                N2%angle(1)%n=2
                N2%angle(2)%n=3
                N2%angle(3)%n=4
             end if !found
          end if ! quartic_bond
          if(non_bond) then
             found=.FALSE.
             i=0
             do while((.not.found).and.(i.lt.nfield))
                i=i+1
                if(field(i).eq.elt) found=.TRUE.
             end do
             if(found) then
                print *,trim(line)
                read(field(4),*) N2%vdW%r0
                read(field(5),*) N2%vdW%eps
             end if !found
          end if ! quartic_bond
       end if
    end do
    close(1)


    print *,"bond_style class2"
    print *,"bond_coeff 1   ",N2%bond(1)%R0,N2%bond(1)%R0,N2%bond(2)%K,N2%bond(3)%K
    print *,"pair_style lj/class2"
    print *,"pair_coeff 1 1",N2%vdW%eps,N2%vdW%r0

    
  end subroutine read_COMPASS_parameters


  ! --------------------------------------------------------------------------------------
  !
  !              read_cp2k_input_file()
  !
  ! --------------------------------------------------------------------------------------
  subroutine read_cp2k_input_file(inputfile,box)
    implicit none
    character (len=1024)::inputfile
    integer::i,io
    character (len=1024)::line
    character (len=32)::field(32)
    integer::nfield
    type(t_Box)::box

    open(unit=1,file=inputfile,form='formatted')

    i=0
    io=0
    do while(io==0)
       read(1,'(A)',iostat=io) line
       i=i+1
       call line_parser(line,nfield,field)
       if(field(1).eq.'ABC') then
          print *,nfield,' --> ',(field(i),i=1,nfield)
          read(field(3),*) box%cell(1,1) ; box%cell(1,2)=0.0 ; box%cell(1,3)=0.0
          read(field(4),*) box%cell(2,2) ; box%cell(2,1)=0.0 ; box%cell(2,3)=0.0
          read(field(5),*) box%cell(3,3) ; box%cell(3,2)=0.0 ; box%cell(3,1)=0.0
       end if
    end do
    close(1)

  end subroutine read_cp2k_input_file
  ! --------------------------------------------------------------------------------------
  !
  !              cp2k_aimd()
  !
  ! --------------------------------------------------------------------------------------
  subroutine cp2k_aimd(box,trajfile)
    implicit none
    type(t_Box)::box
    character (len=1024)::trajfile
    character (len=1024)::line
    character (len=32)::field(32)
    integer::nfield
    integer::io,i,natom,j
    double precision:: x,y,z

    open(unit=1,file=trajfile,form='formatted')
    read(1,'(A)',iostat=io) line
    call line_parser(line,nfield,field)
    read(field(1),*) natom
    print *,"# natom=",natom
    close(1)

    open(unit=1,file=trajfile,form='formatted')
    io=0
    i=0
    do while(io==0)
       read(1,*,iostat=io)
       i=i+1
    end do
    close(1)
    box%n_configuration=i/(natom+2)
    print *,"# n_configuration=",box%n_configuration
    allocate(box%configurations(box%n_configuration))
    do i=1,box%n_configuration
       box%configurations(i)%n_molecule=1
       allocate(box%configurations(i)%molecules(box%configurations(i)%n_molecule))
    end do
    close(1)
    open(unit=1,file=trajfile,form='formatted')
    i=0
    io=0
    do while(i<box%n_configuration)
       i=i+1
       allocate(box%configurations(i)%molecules(1)%atoms(natom))
       box%configurations(i)%molecules(1)%n_atom=natom
       read(1,'(A)',iostat=io) line
       read(1,'(A)',iostat=io) line
       j=0
       do while(j<natom)
          j=j+1
          read(1,'(A)',iostat=io) line
          call line_parser(line,nfield,field)
          read(field(1),*) box%configurations(i)%molecules(1)%atoms(j)%elt
          read(field(2),*) x
          box%configurations(i)%molecules(1)%atoms(j)%q(1)=x
          read(field(3),*) y
          box%configurations(i)%molecules(1)%atoms(j)%q(2)=y
          read(field(4),*) z
          box%configurations(i)%molecules(1)%atoms(j)%q(3)=z
          call elt2Z(box%configurations(i)%molecules(1)%atoms(j)%elt,box%configurations(i)%molecules(1)%atoms(j)%Zato)
       end do
    end do
    close(1)

  end subroutine cp2k_aimd


  subroutine update_group(group,box,iconf)
    implicit none
    type(t_Box)::box
    type (t_Group)::group
    integer::imol,jat,i,idx_set,iconf,k,j
    logical::found
    double precision::R(4),timestep
    integer::idx,idx_grp
    character(len=2)::elt
    integer,allocatable::temp(:)

    do imol=1,box%configurations(iconf)%n_molecule
       do jat=1,box%configurations(iconf)%molecules(imol)%n_atom
          found=.False.
          i=0
          do while((i.lt.group%nelt).and.(.not.found))
             i=i+1
             if(group%elt_table(i)%elt.eq.box%configurations(iconf)%molecules(imol)%atoms(jat)%elt) then
                found=.True.
                idx_set=group%elt_table(i)%idx
             end if
          end do
          if(found) then

             group%n_atoms_set(idx_set)=group%n_atoms_set(idx_set)+1
             allocate(temp(group%n_atoms_set(idx_set)))
             temp(1:group%n_atoms_set(idx_set)-1)=group%set(idx_set)%list_atoms
             temp(group%n_atoms_set(idx_set))=jat
             call move_alloc(from=temp,to=group%set(idx_set)%list_atoms)

          end if
       end do
    end do


    !print *,iconf,(group%n_atoms_set(i),i=1,group%n_set)
    do i=1,group%n_set

       do k=1,3
          group%set(i)%MC(k)=0.0
       end do
       do jat=1,group%n_atoms_set(i)
          idx=group%set(i)%list_atoms(jat)
          do k=1,3
             group%set(i)%MC(k)=group%set(i)%MC(k)+box%configurations(iconf)%molecules(1)%atoms(idx)%q(k)
          end do
       end do
       do k=1,3
          group%set(i)%MC(k)=group%set(i)%MC(k)/group%n_atoms_set(i)
       end do
       !print *,(group%set(i)%list_atoms(jat),jat=1,group%n_atoms_set(i))
       !print *,(group%set(i)%MC(k),k=1,3)
    end do

    R(4)=0.0
    do i=1,group%n_set-1
       do j=i+1,group%n_set
          do k=1,3
             R(k)=group%set(j)%MC(k)-group%set(i)%MC(k)
             R(k)=modulo(R(k),box%cell(k,k))
             if(R(k).gt.0.5*box%cell(k,k)) then
                R(k)=R(k)-box%cell(k,k)
             end if
             if(-R(k).gt.0.5*box%cell(k,k)) then
                R(k)=R(k)+box%cell(k,k)
             end if
             R(4)=R(4)+R(k)*R(k)
          end do
       end do
    end do
    R(4)=sqrt(R(4))
    timestep=0.5 !fs
    open(61,file='R.dat',action='write',position='append')
    write(61,*) iconf*timestep,R(4)
    close(61)



  end subroutine update_group


  ! --------------------------------------------------------------------------------------
  !
  !              QE_cpmd()
  !
  ! --------------------------------------------------------------------------------------
  ! subroutine QE_cpmd()
  !   implicit none
  !   type(t_Molecule),allocatable::evol(:)
  !   type(t_Molecule)::elt

  !   character (len=1024)::trajfile
  !   character (len=1024)::xyzfile
  !   character (len=1024)::line
  !   integer::io
  !   character (len=32)::field(32)
  !   integer::nfield
  !   integer::i,j,natom,nevol
  !   xyzfile="/home/bulou/occigen/scratch/276h2o/276h2o_sorted.xyz"

  !   open(unit=1,file=xyzfile,form='formatted')

  !   i=0
  !   io=0
  !   do while(io==0)
  !      read(1,'(A)',iostat=io) line
  !      i=i+1
  !      call line_parser(line,nfield,field)
  !      !print *,nfield,' --> ',(field(i),i=1,nfield)
  !      if(i==1) then
  !         !print *,nfield,' --> ',(field(i),i=1,nfield)
  !         read(field(1),*) natom
  !         print *,"# natom=",natom
  !         allocate(elt%atoms(natom))
  !      end if
  !      if(i>2) then
  !         read(field(1),*) elt%atoms(i-2)%elt
  !      end if
  !   end do
  !   close(1)


  !   trajfile="/home/bulou/occigen/scratch/276h2o/outdir/276H2O.pos"
  !   open(unit=1,file=trajfile,form='formatted')

  !   io=0
  !   i=0
  !   do while(io==0)
  !      read(1,*,iostat=io)
  !      i=i+1
  !   end do
  !   close(1)
  !   nevol=i/(natom+1)
  !   print *,"# nevol=",nevol
  !   allocate(evol(nevol))

  !   open(unit=1,file=trajfile,form='formatted')
  !   i=0
  !   io=0
  !   do while(i<nevol)
  !      i=i+1
  !      allocate(evol(i)%atoms(natom))
  !      evol(i)%n_atom=natom
  !      read(1,'(A)',iostat=io) line
  !      j=0
  !      do while(j<natom)
  !         j=j+1
  !         read(1,'(A)',iostat=io) line
  !         call line_parser(line,nfield,field)
  !         read(field(1),*) evol(i)%atoms(j)%q(1)
  !         read(field(2),*) evol(i)%atoms(j)%q(2)
  !         read(field(3),*) evol(i)%atoms(j)%q(3)
  !         evol(i)%atoms(j)%elt=elt%atoms(j)%elt

  !      end do
  !   end do
  !   close(1)


  ! end subroutine QE_cpmd
  ! --------------------------------------------------------------------------------------
  !
  !              line_parser()
  !
  ! --------------------------------------------------------------------------------------
  subroutine line_parser(line,nfield,field)
    implicit none
    character (len=1024)::line
    character (len=32)::field(32)
    integer::nfield,i

    integer::lline
    integer::beg,fin

    nfield=0
    lline=len(trim(line))
    beg=-1
    i=1
    do while(i<=lline)
       ! nous ne sommes pas encore dans un champ
       if((len(trim(line(i:i)))==1).and.(beg==-1)) then
          beg=i
          nfield=nfield+1
       end if

       ! on vient de sortir du champ
       if((len(trim(line(i:i)))==0).and.(beg>0)) then
          fin=i-1
          field(nfield)=line(beg:fin)
          beg=-1
       end if
       i=i+1
    end do


    if(beg>0) then
       fin=i-1
       field(nfield)=line(beg:fin)
    end if
  end subroutine line_parser


  ! --------------------------------------------------------------------------------------
  !
  !              init()
  !
  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------
  !
  !              elt2Z(elt,Zato)
  !
  ! --------------------------------------------------------------------------------------
  subroutine init()
    use Element
    implicit none

    PeriodicTable(1)%Z=1 ; PeriodicTable(1)%elt='H' ; PeriodicTable(1)%mass=1.0
    PeriodicTable(2)%Z=2 ; PeriodicTable(2)%elt='He' ; PeriodicTable(2)%mass=4.0
    PeriodicTable(3)%Z=3 ; PeriodicTable(3)%elt='Li' ; PeriodicTable(3)%mass=6.9
    PeriodicTable(4)%Z=4 ; PeriodicTable(4)%elt='Be' ; PeriodicTable(4)%mass=9.0
    PeriodicTable(5)%Z=5 ; PeriodicTable(5)%elt='B' ; PeriodicTable(5)%mass=10.8
    PeriodicTable(6)%Z=6 ; PeriodicTable(6)%elt='C' ; PeriodicTable(6)%mass=12.0
    PeriodicTable(7)%Z=7 ; PeriodicTable(7)%elt='N' ; PeriodicTable(7)%mass=14.0
    PeriodicTable(8)%Z=8 ; PeriodicTable(8)%elt='O' ; PeriodicTable(8)%mass=16.0
    PeriodicTable(9)%Z=9 ; PeriodicTable(9)%elt='F' ; PeriodicTable(9)%mass=19.0
    PeriodicTable(10)%Z=10 ; PeriodicTable(10)%elt='Ne' ; PeriodicTable(10)%mass=20.2
    PeriodicTable(11)%Z=11 ; PeriodicTable(11)%elt='Na' ; PeriodicTable(11)%mass=23.0
    PeriodicTable(12)%Z=12 ; PeriodicTable(12)%elt='Mg' ; PeriodicTable(12)%mass=24.3
    PeriodicTable(13)%Z=13 ; PeriodicTable(13)%elt='Al' ; PeriodicTable(13)%mass=27.0
    PeriodicTable(14)%Z=14 ; PeriodicTable(14)%elt='Si' ; PeriodicTable(14)%mass=28.1
    PeriodicTable(15)%Z=15 ; PeriodicTable(15)%elt='P' ; PeriodicTable(15)%mass=31.0
    PeriodicTable(16)%Z=16 ; PeriodicTable(16)%elt='S' ; PeriodicTable(16)%mass=32.1
    PeriodicTable(17)%Z=17 ; PeriodicTable(17)%elt='Cl' ; PeriodicTable(17)%mass=35.5
    PeriodicTable(18)%Z=18 ; PeriodicTable(18)%elt='Ar' ; PeriodicTable(18)%mass=39.9
    PeriodicTable(19)%Z=19 ; PeriodicTable(19)%elt='K' ; PeriodicTable(19)%mass=39.1
    PeriodicTable(20)%Z=20 ; PeriodicTable(20)%elt='Ca' ; PeriodicTable(20)%mass=40.1

    deq(1,8)=1.2 ; deq(8,1)=deq(1,8)
    deq(7,7)=1.2 ;
    deq(8,8)=1.3 ;
  end subroutine init

  ! --------------------------------------------------------------------------------------
  !
  !              elt2Z(elt,Zato)
  !
  ! --------------------------------------------------------------------------------------
  subroutine elt2Z(elt,Zato)
    use Element
    implicit none
    integer::Zato,i
    character(len=2)::elt
    logical::found

    found=.False.
    i=0
    do while  ((.not.found).and.(i<Nelt))
       i=i+1
       if(elt.eq.PeriodicTable(i)%elt) then
          found=.True.
       end if
    end do
    if(found) then
       Zato=PeriodicTable(i)%Z
       !write(*,*) elt,Zato
    else
       Zato=-1
       write(*,*) elt," not found"
    end if

  end subroutine elt2Z

  ! --------------------------------------------------------------------------------------
  !
  !              save()
  !
  ! --------------------------------------------------------------------------------------
  subroutine save(box,fmt,iconf)
    implicit none
    character(len=22)::fmt,conf
    type(t_Box)::box
    integer::imol,jat,k,iconf,natom
    write(*,*) '# Saving ',fmt
    ! ----------------------------- lammps trajectory file -------------------------------
    if(fmt.eq.'lammps trajectory file') then
       open(unit=1,file="evol.lammpstrj",form='formatted',status='unknown')
       do iconf=1,box%n_configuration
          do imol=1,box%configurations(iconf)%n_molecule
             write(1,*) "ITEM: TIMESTEP"
             write(1,*) i-1
             write(1,*) "ITEM: NUMBER OF ATOMS"
             write(1,*) box%configurations(iconf)%molecules(imol)%n_atom
             write(1,*) "ITEM: BOX BOUNDS pp pp pp"
             write(1,*) 0,box%cell(1,1)
             write(1,*) 0,box%cell(2,2)
             write(1,*) 0,box%cell(3,3)
             write(1,*) "ITEM: ATOMS id mol type x y z vx vy vz fx fy fz"
             do jat=1,box%configurations(iconf)%molecules(imol)%n_atom
                write(1,*) jat,1,box%configurations(iconf)%molecules(imol)%atoms(jat)%Zato,&
                     (box%configurations(iconf)%molecules(imol)%atoms(jat)%q(k),k=1,3),0.0,0.0,0.0,0.0,0.0,0.0
             end do
          end do
       end do
       close(1)
    end if
    ! ----------------------------- XYZ --------------------------------------------------
    if(fmt.eq.'xyz') then
       open(unit=1,file="evol.xyz",form='formatted',status='unknown')
       do imol=1,box%configurations(iconf)%n_molecule
          write(1,*) box%configurations(iconf)%molecules(imol)%n_atom
          write(1,*) 
          do jat=1,box%configurations(iconf)%molecules(imol)%n_atom
             write(1,*) box%configurations(iconf)%molecules(imol)%atoms(jat)%elt,&
                  box%configurations(iconf)%molecules(imol)%atoms(jat)%q(1),&
                  box%configurations(iconf)%molecules(imol)%atoms(jat)%q(2),&
                  box%configurations(iconf)%molecules(imol)%atoms(jat)%q(3)
          end do
       end do
       close(1)
    end if
    ! ----------------------------- XSF --------------------------------------------------
    ! see http://www.xcrysden.org/doc/XSF.html
    if(fmt.eq.'xsf') then
       open(unit=1,file="evol.xsf",form='formatted',status='unknown')
       write(1,*) "CRYSTAL"
       write(1,*) "# these are primitive lattice vectors (in Angstroms)"
       write(1,*) " PRIMVEC"
       write(1,*) "      ",(box%cell(1,k),k=1,3)
       write(1,*) "      ",(box%cell(2,k),k=1,3)
       write(1,*) "      ",(box%cell(3,k),k=1,3)
       write(1,*) "# these are convetional lattice vectors (in Angstroms)"
       write(1,*) " CONVEC"
       write(1,*) "      ",(box%cell(1,k),k=1,3)
       write(1,*) "      ",(box%cell(2,k),k=1,3)
       write(1,*) "      ",(box%cell(3,k),k=1,3)
       write(1,*) "# these are atomic coordinates in a primitive unit cell "
       write(1,*) " # (in Angstroms)"
       write(1,*) "PRIMCOORD ",1
       write(1,*) box%configurations(iconf)%n_atom," 1"
       do imol=1,box%configurations(iconf)%n_molecule
          do jat=1,box%configurations(iconf)%molecules(imol)%n_atom
             write(1,*) box%configurations(iconf)%molecules(imol)%atoms(jat)%Zato,&
                  (box%configurations(iconf)%molecules(imol)%atoms(jat)%q(k),k=1,3)
          end do
       end do
       close(1)
    end if
    ! ----------------------------- AXYZ --------------------------------------------------
    if(fmt.eq.'axsf') then
       open(unit=1,file="evol.axsf",form='formatted',status='unknown')
       write(1,*) "ANIMSTEPS ",box%n_configuration
       write(1,*) "CRYSTAL"
       write(1,*) "# these are primitive lattice vectors (in Angstroms)"
       write(1,*) " PRIMVEC"
       write(1,*) "      ",(box%cell(1,k),k=1,3)
       write(1,*) "      ",(box%cell(2,k),k=1,3)
       write(1,*) "      ",(box%cell(3,k),k=1,3)
       write(1,*) "# these are convetional lattice vectors (in Angstroms)"
       write(1,*) " CONVEC"
       write(1,*) "      ",(box%cell(1,k),k=1,3)
       write(1,*) "      ",(box%cell(2,k),k=1,3)
       write(1,*) "      ",(box%cell(3,k),k=1,3)
       write(1,*) "# these are atomic coordinates in a primitive unit cell "
       write(1,*) " # (in Angstroms)"
       do iconf=1,box%n_configuration
          do imol=1,box%configurations(iconf)%n_molecule
             write(1,*) "PRIMCOORD ",iconf
             write(1,*) box%configurations(iconf)%molecules(imol)%n_atom," 1"
             do jat=1,box%configurations(iconf)%molecules(imol)%n_atom
                write(1,*) box%configurations(iconf)%molecules(imol)%atoms(jat)%Zato,&
                     (box%configurations(iconf)%molecules(imol)%atoms(jat)%q(k),k=1,3)
             end do
          end do
       end do
       close(1)
    end if


  end subroutine save
  ! --------------------------------------------------------------------------------------
  !
  !              write_input_file()
  !
  ! --------------------------------------------------------------------------------------
  subroutine write_input_file(box,iconf,run_type)
    implicit none
    type(t_Box)::box
    integer::imol,iconf,jat,k
    integer::nsteps,max_scf,max_iter
    character(len=32)::run_type
    double precision::temperature,eps_scf
    print *,"-----------------------------------------------------------------------------"
    print *,"# Entering write_input_file..."
    temperature=600.0
    nsteps=10000
    eps_scf=1.0e-5
    max_scf=500
    max_iter=2000

    select case(run_type)
    case('gopt')
       open(unit=1,file="TMP-gopt.in",form='formatted',status='unknown')
       write(1,*) '@SET SYSTEM run-gopt'
       write(1,*) '#@set LIBDIR /home/bulou/src/cp2k/data'           ! pc-hervecal
       write(1,*) '@set LIBDIR /usr/local/cp2k/cp2k-6.1.0.g8/data'   ! hpc
       write(1,*) '@set EPS_SCF ',eps_scf
       write(1,*) '@set MAX_SCF ',max_scf
       write(1,*) '@set MAX_ITER ',max_iter
       write(1,*) ''
       write(1,*) '&GLOBAL'
       write(1,*) '  PROJECT ${SYSTEM}          # Name of calculation: run-1H2O-GOPT '
       write(1,*) '  RUN_TYPE GEO_OPT           # Perform a geometry optimization'
       write(1,*) '  PRINT_LEVEL LOW            # Do not output too much information'
       write(1,*) '&END GLOBAL'
       write(1,*) ''
       write(1,*) '&FORCE_EVAL'
       write(1,*) '  METHOD QuickStep'
       write(1,*) '  &DFT                       # Use density functional theory'
       write(1,*) '    BASIS_SET_FILE_NAME ${LIBDIR}/BASIS_MOLOPT'
       write(1,*)  '   POTENTIAL_FILE_NAME ${LIBDIR}/GTH_POTENTIALS'
       write(1,*) '    &MGRID'
       write(1,*) '      CUTOFF 300             # plane-wave cutoff for the charge density [Rydbergs]'
       write(1,*) '    &END MGRID  '
       write(1,*) '    &SCF                     '
       write(1,*) '       EPS_SCF ${EPS_SCF}'
       write(1,*) '       MAX_SCF ${MAX_SCF}'
       write(1,*) '    &END SCF'
       write(1,*) '    &XC                      # parameters for exchange-correlation functional'
       write(1,*) '      &XC_FUNCTIONAL PBE'
       write(1,*) '      &END XC_FUNCTIONAL'
       write(1,*) '      &XC_GRID               # some tricks to speed up the calculation of the xc potential'
       write(1,*) '        XC_DERIV       SPLINE2_smooth'
       write(1,*) '        XC_SMOOTH_RHO  NN10'
       write(1,*) '      &END XC_GRID'
       write(1,*) '    &END XC'
       write(1,*) '  &END DFT'
       write(1,*) '  &SUBSYS'
       write(1,*) '    &CELL                    # box containing the molecule: 15x15x15 Angstroms'
       write(1,*) '      ABC [angstrom] ',(box%cell(i,i),i=1,3)
       write(1,*) '      periodic xyz'
       write(1,*) '    &END CELL'
       write(1,*) '    &COORD                   # coordinates of the atoms in the box [Angstroms]'
       do imol=1,box%configurations(iconf)%n_molecule
          print *,"# Molecule:",box%configurations(iconf)%molecules(imol)%n_atom
          do jat=1,box%configurations(iconf)%molecules(imol)%n_atom
             write(1,*) box%configurations(iconf)%molecules(imol)%atoms(jat)%elt,&
                  (box%configurations(iconf)%molecules(imol)%atoms(jat)%q(k),k=1,3)
          end do
       end do
       write(1,*) '   &END COORD'
       do k=1,box%configurations(iconf)%n_element
          select case(box%configurations(iconf)%list_elements(k)%elt)
          case('H')
             write(1,*) '    &KIND H'
             write(1,*) '      BASIS_SET DZVP-MOLOPT-GTH'
             write(1,*) '      POTENTIAL GTH-PBE-q1'
             write(1,*) '    &END KIND'
          case('C')
             write(1,*) '    &KIND C                  # basis sets and pseudo-potentials for atomic species'
             write(1,*) '      BASIS_SET DZVP-MOLOPT-GTH'
             write(1,*) '      POTENTIAL GTH-PBE-q4'
             write(1,*) '    &END KIND'
          case('N')
             write(1,*) '    &KIND N'
             write(1,*) '      BASIS_SET DZVP-MOLOPT-GTH'
             write(1,*) '      POTENTIAL GTH-PBE-q5'
             write(1,*) '    &END KIND'
          case('O')
             write(1,*) '    &KIND O'
             write(1,*) '      BASIS_SET DZVP-MOLOPT-GTH'
             write(1,*) '      POTENTIAL GTH-PBE-q6'
             write(1,*) '    &END KIND'
          case('Ca')
             write(1,*) '    &KIND Ca'
             write(1,*) '      BASIS_SET DZVP-MOLOPT-SR-GTH'
             write(1,*) '      POTENTIAL GTH-PBE-q10'
             write(1,*) '    &END KIND'
          case default
             print *,"# UNKOWN element ",box%configurations(iconf)%list_elements(k)%elt
          end select
       end do

       write(1,*) '  &END SUBSYS'
       write(1,*) '&END FORCE_EVAL'

       write(1,*) '&MOTION'
       write(1,*) '  &GEO_OPT'
       write(1,*) '    TYPE MINIMIZATION'
       write(1,*) '    MAX_DR    1.0E-03'
       write(1,*) '    MAX_FORCE 1.0E-03'
       write(1,*) '    RMS_DR    1.0E-03'
       write(1,*) '    RMS_FORCE 1.0E-03'
       write(1,*) '    MAX_ITER ${MAX_ITER}'
       write(1,*) '    OPTIMIZER CG'
       write(1,*) '    &CG'
       write(1,*) '      MAX_STEEP_STEPS  0'
       write(1,*) '      RESTART_LIMIT 9.0E-01'
       write(1,*) '    &END CG'
       write(1,*) '  &END GEO_OPT'
       write(1,*) '  #&CONSTRAINT'
       write(1,*) '  #  &FIXED_ATOMS'
       write(1,*) '  #    COMPONENTS_TO_FIX XYZ'
       write(1,*) '  #    LIST 1'
       write(1,*) '  #  &END FIXED_ATOMS'
       write(1,*) '  #&END CONSTRAINT'
       write(1,*) '&END MOTION'


       write(1,*) ''
    case('aimd')
       ! cp2k - aimd
       open(unit=1,file="TMP-aimd.in",form='formatted',status='unknown')
       write(1,*) '@SET SYSTEM run-CaCO3-AIMD'
       write(1,*) '#@set LIBDIR /home/bulou/src/cp2k/data'
       write(1,*) '@set LIBDIR /usr/local/cp2k/cp2k-6.1.0.g8/data'   ! hpc
       write(1,*) '@set TEMPERATURE',temperature
       write(1,*) '@set NSTEPS',nsteps
       write(1,*) '@set EPS_SCF ',eps_scf
       write(1,*) '@set MAX_SCF ',max_scf

       write(1,*) ''
       write(1,*) ''
       write(1,*) '&GLOBAL'
       write(1,*) '  PROJECT ${SYSTEM}'
       write(1,*) '  RUN_TYPE MD'
       write(1,*) '  PRINT_LEVEL LOW'
       write(1,*) '&END GLOBAL'
       write(1,*) ''
       write(1,*) '&FORCE_EVAL'
       write(1,*) '  METHOD QuickStep'
       write(1,*) ''
       write(1,*) '  &DFT'
       write(1,*) '    BASIS_SET_FILE_NAME ${LIBDIR}/BASIS_MOLOPT'
       write(1,*) '    POTENTIAL_FILE_NAME ${LIBDIR}/GTH_POTENTIALS'
       write(1,*) ''
       write(1,*) '    &MGRID'
       write(1,*) '      CUTOFF 300'
       write(1,*) '    &END MGRID'
       write(1,*) '    &QS'
       write(1,*) '      EPS_DEFAULT 1.0E-10'
       write(1,*) '      EXTRAPOLATION ASPC'
       write(1,*) '      EXTRAPOLATION_ORDER 4'
       write(1,*) '    &END QS'
       write(1,*) '    &SCF'
       write(1,*) '      EPS_SCF ${EPS_SCF}'
       write(1,*) '      MAX_SCF ${MAX_SCF}'
       write(1,*) '      &OT ON'
       write(1,*) '        PRECONDITIONER FULL_SINGLE_INVERSE'
       write(1,*) '        MINIMIZER DIIS'
       write(1,*) '      &END OT'
       write(1,*) '    &END SCF'
       write(1,*) ''
       write(1,*) '    &XC'
       write(1,*) '      &XC_FUNCTIONAL PBE'
       write(1,*) '      &END XC_FUNCTIONAL'
       write(1,*) '      &XC_GRID'
       write(1,*) '        XC_DERIV        SPLINE2_smooth'
       write(1,*) '        XC_SMOOTH_RHO   NN10'
       write(1,*) '      &END XC_GRID'
       write(1,*) '    &END XC'
       write(1,*) ''
       write(1,*) '&poisson'
       write(1,*) 'periodic xyz'
       write(1,*) '&end poisson'
       write(1,*) ''
       write(1,*) '  &END DFT'
       write(1,*) ''
       write(1,*) '  &SUBSYS'
       write(1,*) ''
       write(1,*) '    &CELL'
       write(1,*) '      ABC [angstrom] ',(box%cell(i,i),i=1,3)
       write(1,*) '      PERIODIC XYZ'
       write(1,*) '    &END CELL'
       write(1,*) ''
       write(1,*) '    &COORD'
       do imol=1,box%configurations(iconf)%n_molecule
          print *,"# Molecule:",box%configurations(iconf)%molecules(imol)%n_atom
          do jat=1,box%configurations(iconf)%molecules(imol)%n_atom
             write(1,*) box%configurations(iconf)%molecules(imol)%atoms(jat)%elt,&
                  (box%configurations(iconf)%molecules(imol)%atoms(jat)%q(k),k=1,3)
          end do
       end do
       write(1,*) '    &END COORD'
       write(1,*) ''

       do k=1,box%configurations(iconf)%n_element
          select case(box%configurations(iconf)%list_elements(k)%elt)
          case('H')
             write(1,*) '    &KIND H'
             write(1,*) '      BASIS_SET DZVP-MOLOPT-GTH'
             write(1,*) '      POTENTIAL GTH-PBE-q1'
             write(1,*) '    &END KIND'
          case('C')
             write(1,*) '    &KIND C                  # basis sets and pseudo-potentials for atomic species'
             write(1,*) '      BASIS_SET DZVP-MOLOPT-GTH'
             write(1,*) '      POTENTIAL GTH-PBE-q4'
             write(1,*) '    &END KIND'
          case('N')
             write(1,*) '    &KIND N'
             write(1,*) '      BASIS_SET DZVP-MOLOPT-GTH'
             write(1,*) '      POTENTIAL GTH-PBE-q5'
             write(1,*) '    &END KIND'
          case('O')
             write(1,*) '    &KIND O'
             write(1,*) '      BASIS_SET DZVP-MOLOPT-GTH'
             write(1,*) '      POTENTIAL GTH-PBE-q6'
             write(1,*) '    &END KIND'
          case('Ca')
             write(1,*) '    &KIND Ca'
             write(1,*) '      BASIS_SET DZVP-MOLOPT-SR-GTH'
             write(1,*) '      POTENTIAL GTH-PBE-q10'
             write(1,*) '    &END KIND'
          case default
             print *,"# UNKOWN element ",box%configurations(iconf)%list_elements(k)%elt
             stop
          end select
       end do


       write(1,*) '  &END SUBSYS'
       write(1,*) ''
       write(1,*) '&END FORCE_EVAL'
       write(1,*) ''
       write(1,*) '&MOTION'
       write(1,*) '  &MD'
       write(1,*) '    ENSEMBLE         NVE'
       write(1,*) '    STEPS     ${NSTEPS}       '
       write(1,*) '    TIMESTEP [fs]    0.5'
       write(1,*) '    TEMPERATURE [K]  ${TEMPERATURE} '
       write(1,*) ''
       write(1,*) '    ANGVEL_ZERO                 T'
       write(1,*) '#    ANGVEL_TOL [fs^-1deg]       0.001'
       write(1,*) '    ANGVEL_TOL [fs^-1angstrom]  0.001'
       write(1,*) '    COMVEL_TOL [fs^-1angstrom]  0.001'
       write(1,*) '  &END MD'
       write(1,*) ''
       write(1,*) '  &PRINT'
       write(1,*) '    &TRAJECTORY  SILENT'
       write(1,*) '      FILENAME =${SYSTEM}.xyz'
       write(1,*) '      &EACH'
       write(1,*) '        MD 1'
       write(1,*) '      &END EACH'
       write(1,*) '    &END TRAJECTORY'
       write(1,*) '    &VELOCITIES  SILENT'
       write(1,*) '      FILENAME =${SYSTEM}.vel'
       write(1,*) '      &EACH'
       write(1,*) '        MD 1'
       write(1,*) '      &END EACH'
       write(1,*) '    &END VELOCITIES'
       write(1,*) '    &FORCES  SILENT'
       write(1,*) '      FILENAME =${SYSTEM}.force'
       write(1,*) '      &EACH'
       write(1,*) '        MD 1'
       write(1,*) '      &END EACH'
       write(1,*) '    &END FORCES'
       write(1,*) '  &END PRINT'
       write(1,*) ''
       write(1,*) '&END MOTION'

       write(1,*) '# &EXT_RESTART'
       write(1,*) '#   RESTART_FILE_NAME run-CaCO3-AIMD-1.restart'
       write(1,*) '# &END'
    case default
       print *,"# ERROR: UNKOWN run_type", run_type
    end select


    close(1)
    print *,"# End of write_input_file()"
    print *,"-----------------------------------------------------------------------------"

  end subroutine write_input_file
end program generate_configuration
