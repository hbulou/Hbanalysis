module global
  implicit none
  
  double precision,parameter::a0=0.529177249
  integer,parameter::NCHARFIELD=256

  ! ----------------------------------------------
  type t_Param_Transformation_Strain
     logical::status
     double precision::eps(3,3)
  end type t_Param_Transformation_Strain
  ! ----------------------------------------------
  type t_Param_Transformation
     type(t_Param_Transformation_Strain)::strain
  end type t_Param_Transformation
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
     integer::n_transformations
     type(t_Param_Transformation),allocatable::transformation(:)
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



  ! ----------------------------------------------
  type t_CP2K_global
     logical::state
     character(len=32)::print_level
     character(len=32)::project_name
     character(len=32)::run_type
  end type t_CP2K_global

  
  type t_CP2K_motion_constraint_fixed_atoms
     character(len=3)::COMPONENTS_TO_FIX
     integer::nrec
     integer,allocatable::list(:)
  end type t_CP2K_motion_constraint_fixed_atoms

  type t_CP2K_motion_constraint
     integer::n_types
     type(t_CP2K_motion_constraint_fixed_atoms),allocatable::fixed_atoms(:)
  end type t_CP2K_motion_constraint


  
  type t_CP2K_motion_geo_opt
     integer::step_start_val
  end type t_CP2K_motion_geo_opt

  type t_CP2K_motion
     type(t_CP2K_motion_geo_opt)::geo_opt
     type(t_CP2K_motion_constraint)::constraint
  end type t_CP2K_motion


  type t_CP2K_ext_restart
     logical::RESTART_FILE_NAME 
     logical::RESTART_COUNTERS  
     logical::RESTART_POS  
     logical::RESTART_VEL  
     logical::RESTART_RANDOMG  
     logical::RESTART_SHELL_POS  
     logical::RESTART_CORE_POS  
     logical::RESTART_OPTIMIZE_INPUT_VARIABLES  
     logical::RESTART_SHELL_VELOCITY  
     logical::RESTART_CORE_VELOCITY  
     logical::RESTART_BAROSTAT  
     logical::RESTART_BAROSTAT_THERMOSTAT  
     logical::RESTART_SHELL_THERMOSTAT  
     logical::RESTART_THERMOSTAT  
     logical::RESTART_CELL  
     logical::RESTART_METADYNAMICS  
     logical::RESTART_WALKERS  
     logical::RESTART_BAND  
     logical::RESTART_QMMM  
     logical::RESTART_CONSTRAINT  
     logical::RESTART_BSSE  
     logical::RESTART_DIMER  
     logical::RESTART_AVERAGES  
     logical::RESTART_RTP  
     logical::RESTART_PINT_POS  
     logical::RESTART_PINT_VEL  
     logical::RESTART_PINT_NOSE  
     logical::RESTART_PINT_GLE  
     logical::RESTART_HELIUM_POS  
     logical::RESTART_HELIUM_PERMUTATION  
     logical::RESTART_HELIUM_FORCE  
     logical::RESTART_HELIUM_RNG  
  end type t_CP2K_ext_restart

  
  type t_CP2K_force_eval_dft_scf_mixing
     logical::state
     character(len=32)::method
     double precision::alpha
     double precision::beta
     integer::nbuffer
  end type t_CP2K_force_eval_dft_scf_mixing
  type t_CP2K_force_eval_dft_scf_smear
     logical::state
     character(len=32)::method
     double precision::electronic_temp
  end type t_CP2K_force_eval_dft_scf_smear
   type t_CP2K_force_eval_dft_scf_diag
     logical::state
     character(len=32)::algo
  end type t_CP2K_force_eval_dft_scf_diag
  type t_CP2K_force_eval_dft_scf
     integer::max_scf
     double precision::eps_scf
     character(len=32)::scf_guess
     integer::added_mos
     type(t_CP2K_force_eval_dft_scf_diag)::diag
     type(t_CP2K_force_eval_dft_scf_smear)::smear
     type(t_CP2K_force_eval_dft_scf_mixing)::mixing
  end type t_CP2K_force_eval_dft_scf
  type t_CP2K_force_eval_dft_qs
     double precision::eps_default
  end type t_CP2K_force_eval_dft_qs
  type t_CP2K_force_eval_dft_mgrid
     integer::ngrids
     double precision::cutoff
     double precision::rel_cutoff
  end type t_CP2K_force_eval_dft_mgrid
  type t_CP2K_force_eval_dft_poisson
     character(len=32)::poisson_solver
     character(len=32)::periodic
  end type t_CP2K_force_eval_dft_poisson
  type t_CP2K_force_eval_dft_rtp
     character(len=32)::initial_wfn
  end type t_CP2K_force_eval_dft_rtp


  
  type t_CP2K_force_eval_dft_xc_func
     character(len=32)::racc
     logical::PBE
  end type t_CP2K_force_eval_dft_xc_func

  type t_CP2K_force_eval_dft_xc
     double precision::density_cutoff
     double precision::gradient_cutoff
     double precision::tau_cutoff
     type(t_CP2K_force_eval_dft_xc_func)::xcfunc
  end type t_CP2K_force_eval_dft_xc


  type t_CP2K_force_eval_dft
     logical::state
     character(len=1024)::basis_set_file_name
     character(len=1024)::potential_file_name
     type(t_CP2K_force_eval_dft_scf)::scf
     type(t_CP2K_force_eval_dft_qs)::qs
     type(t_CP2K_force_eval_dft_mgrid)::mgrid
     type(t_CP2K_force_eval_dft_poisson)::poisson
     type(t_CP2K_force_eval_dft_rtp)::rtp
     type(t_CP2K_force_eval_dft_xc)::xc
  end type t_CP2K_force_eval_dft



  type t_CP2K_force_eval_subsys_cell
     logical::state
     double precision::A(3),B(3),C(3)
     character(len=3)::periodic
     integer::multiple_unit_cell(3)
  end type t_CP2K_force_eval_subsys_cell

  type t_CP2K_force_eval_subsys_coord
     logical::state
     type(t_Molecule)::mol
  end type t_CP2K_force_eval_subsys_coord


  type t_CP2K_force_eval_subsys_topology
     logical::state
     integer::number_of_atoms,multiple_unit_cell(3)
  end type t_CP2K_force_eval_subsys_topology

  type t_GTH_Potential_non_local_projector
     ! # r        : Radius of the non-local part for angular momentum quantum number l
     ! #            defined by the Gaussian function exponents alpha_prj_ppnl
     double precision::r
     ! # nprj_ppnl: Number of the non-local projectors for the angular momentum
     ! #            quantum number l
     integer::nprj_ppnl
     ! # hprj_ppnl: Coefficients of the non-local projector functions
     double precision,allocatable::hprj_ppnl(:)
  end type t_GTH_Potential_non_local_projector
  type t_GTH_Potential
     character(len=32)::name
     ! # n_elec   : Number of electrons for each angular momentum quantum number
     ! #            (electronic configuration -> s p d ...)
     integer::n_type_elec
     integer,allocatable::n_elec(:)
     ! # r_loc    : Radius for the local part defined by the Gaussian function
     ! #            exponent alpha_erf
     double precision::r_loc
     ! # nexp_ppl : Number of the local pseudopotential functions
     integer::nexp_ppl
     ! # cexp_ppl : Coefficients of the local pseudopotential functions
     double precision,allocatable::cexp_ppl(:)
     ! # nprj     : Number of the non-local projectors => nprj = SIZE(nprj_ppnl(:))
     integer::nprj
     type(t_GTH_Potential_non_local_projector),allocatable::nlprj(:)
  end type t_GTH_Potential

  ! C GTH-PBE-q4 GTH-PBE   --> [He] 2s2 2p2
  !    2    2
  !     0.33847124    2    -8.80367398     1.33921085
  !    2
  !     0.30257575    1     9.62248665
  !     0.29150694    0

  ! Ga GTH-BP-q13 GTH-BP   -->  	[Ar] 4s2 3d10 4p1
  !   n_elec s   p
  !     2    1   10
  !     0.49000000    0
  
  !     3        
  !     0.39614555    3    12.22933993    -7.15254431     2.02087227  --> 3<=> (3/2)*(3+1) (n/2)*(n+1)
  !                                       12.52084398    -5.21786976      1 ->  2 ->  3 |   
  !                                                       4.14155572        ->  4 ->  5 | stockage
  !     0.57682875    2     1.65710194     0.27257669                             ->  6 |
  !                                       -0.32251709
  !     0.23837295    1   -16.19719645


  
  type t_Kind
     character(len=2)::element
     character(len=32)::basis_set
     type(t_GTH_Potential)::potential
  end type t_Kind


  type t_CP2K_force_eval_subsys
     logical::state
     type(t_CP2K_force_eval_subsys_cell)::cell
     type(t_CP2K_force_eval_subsys_coord)::coord
     type(t_CP2K_force_eval_subsys_topology)::topology
     integer::n_kinds
     type(t_Kind),allocatable::kinds(:)
  end type t_CP2K_force_eval_subsys

  type t_CP2K_force_eval_print_forces
     character(len=4)::state
  end type t_CP2K_force_eval_print_forces

  type t_CP2K_force_eval_print
     logical::state
     type(t_CP2K_force_eval_print_forces)::forces
  end type t_CP2K_force_eval_print


  type t_CP2K_force_eval
     character(len=32)::method
     type(t_CP2K_force_eval_dft)::dft
     type(t_CP2K_force_eval_subsys)::subsys
     type(t_CP2K_force_eval_print)::print
  end type t_CP2K_force_eval


  type t_CP2K_param
     type(t_CP2K_global)::global
     type(t_CP2K_motion)::motion
     type(t_CP2K_force_eval)::force_eval
     type(t_CP2K_ext_restart)::ext_restart
  end type t_CP2K_param

  ! ----------------------------------------------
  ! ----------------------------------------------
  ! ----------------------------------------------


  
contains
  ! --------------------------------------------------------------------------------------
  !
  !              function to_upper(strIn) result(strOut)
  !
  ! --------------------------------------------------------------------------------------
  function to_upper(strIn) result(strOut)
    ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
    ! Original author: Clive Page
    
    implicit none
    
    character(len=*), intent(in) :: strIn
    character(len=len(strIn)) :: strOut
    integer :: i,j
    
    do i = 1, len(strIn)
        j = iachar(strIn(i:i))
        if (j>= iachar("a") .and. j<=iachar("z") ) then
           strOut(i:i) = achar(iachar(strIn(i:i))-32)
        else
           strOut(i:i) = strIn(i:i)
        end if
     end do
     
   end function to_upper
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
