module global
  implicit none
  
  double precision,parameter::a0=0.529177249
  integer,parameter::NCHARFIELD=256
  ! ----------------------------------------------
  type t_Param_Calculation
     character(len=32)::run_type
     character(len=32)::system
     double precision::ecutoff
     character(len=32)::potential_name
     character(len=32)::basis_name
     character(len=4)::xc_functional
     integer::added_mos
     integer::k(3)
     character(len=16)::scheme
  end type t_Param_Calculation
  ! ----------------------------------------------
  type t_Machine
     character(len=32)::name
     character(len=256)::cp2k_libdir
  end type t_Machine
  integer,parameter::N_MACHINES=2
  type(t_Machine)::List_Machines(N_MACHINES)
  ! ----------------------------------------------
  type t_PeriodicTable
     integer::Z
     character(len=2)::elt
     double precision::mass
     character(len=32)::basis_set
     character(len=32)::potential_set
       ! ----------------------------------------------
  end type t_PeriodicTable
  integer,parameter::Nelt=46
  type (t_PeriodicTable)::PeriodicTable(Nelt)
  double precision::ELEMENT_deq(Nelt,Nelt)
  !-----------------------------------------------
  type t_File
     character(len=2048)::name
     character(len=32)::type
     integer::iconf
  end type t_File
  ! ----------------------------------------------
  type t_Param
     character(len=2048)::filename
     character(len=16)::run_type
     integer::n_inputs
     integer::n_configurations
     integer::n_atoms
     integer::n_lines
     double precision,allocatable::nrj(:)
     type(t_File),allocatable::input(:)
     integer::n_outputs
     type(t_File),allocatable::output(:)
     character(len=32)::title
     character(len=8)::machine  ! jeanzay, hpc
     type(t_Param_Calculation)::calculation
  end type t_Param
  !-----------------------------------------------
  type t_Atom
     logical::save
     integer::idx
     integer::idx_molecule
     integer::idx_cell
     double precision::q(3)
     integer::constraints(3)
     character(len=2)::elt
     integer::Zato
     integer::n_bonds
     integer,allocatable::bonds(:)
     integer::n_angles
     integer,allocatable::angles(:)
  end type t_Atom
  !-----------------------------------------------
  type l_Atom
     integer::idx
     integer::idx_molecule
  end type l_Atom
  type t_Bond
     integer::idx
     type(l_Atom)::list_atoms(2)
  end type t_Bond
  type t_Angle
     integer::idx
     type(l_Atom)::list_atoms(3)
  end type t_Angle
  !-----------------------------------------------
  type t_Element
     integer::Zato
     character(len=2)::elt
     integer::n
  end type t_Element
  !-----------------------------------------------
  type t_Constraint
     !integer::n_XY
     !integer,allocatable::XY(:)
     integer::n_XYZ
     integer,allocatable::XYZ(:)
  end type t_Constraint
  !-----------------------------------------------
  type t_Molecule
     logical::save
     integer:: idx
     integer::idx_cell
     integer:: n_atoms
     integer,allocatable::order_idx_save(:)
     type(t_Atom),allocatable::atoms(:)
     double precision::MassCenter(3)
     integer::n_bonds
     integer::n_angles
     integer::n_elements
     type(t_Element),allocatable::elements(:)
     type(t_Constraint)::constraints
     double precision::limits(3,2) ! (1,1)-xmin (1,2)-xmax
  end type t_Molecule
  !-----------------------------------------------
  type t_Cell
     integer::idx
     double precision::L(3,3)
     integer:: n_atoms
     ! molecules belonging to the cell
     integer::n_molecules
     type(t_Molecule),allocatable::molecules(:)
     ! bonds belonging to the cell
     integer::n_bonds
     type(t_Bond),allocatable::bonds(:)
     integer::n_angles
     type(t_Angle),allocatable::angles(:)
     !integer::k(3)
     character(len=4)::periodicity

     double precision::evacc ! utilisé pour un slab. 0 signifie pas de vide (default)
     !                       1 signifie un couche de vide de la même épaisseur
     !                          que la couche de slabe
     type(t_Constraint)::constraints
     integer::n_elements
     type(t_Element),allocatable::elements(:)
  end type t_Cell
  !-----------------------------------------------
  type t_Configuration
     !integer::n_cells
     type(t_Cell)::cells
  end type t_Configuration
  !-----------------------------------------------
  type t_Univers
     integer::n_configurations
     type(t_Configuration),allocatable::configurations(:)
  end type t_Univers

contains
  ! --------------------------------------------------------------------------------------
  !
  !              line_parser()
  !
  ! --------------------------------------------------------------------------------------
  subroutine line_parser(line,nfield,field)
    implicit none
    character (len=1024)::line
    character (len=NCHARFIELD)::field(32)
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
  !              get_last_field(line,field)
  !
  ! --------------------------------------------------------------------------------------
  subroutine get_last_field(line,field)
    implicit none
    integer::fin,deb
    character(len=32)::field
    character(len=NCHARFIELD)::line
    deb=len(trim(line))
    print *,"<",line,">",deb
    fin=-1
    do while(deb.gt.fin)
    !    print *,deb
        if(line(deb:deb).eq."/") then
           fin=len(trim(line))
           deb=deb+1
           print *,line(deb:fin)
           field=line(deb:fin)
        end if
        deb=deb-1
    end do
  end subroutine get_last_field

  ! --------------------------------------------------------------------------------------
  !
  !              GLOBAL_init_param(param)
  !
  ! --------------------------------------------------------------------------------------
  subroutine GLOBAL_init_param(param)
    implicit none
    type(t_Param)::param
    param%calculation%system="default"
    param%calculation%added_mos=200
  end subroutine GLOBAL_init_param

end module global
