module IO
  use global
  implicit none
contains
#define msg(x) print *,"### ",x,__FILE__," > "
!#define msg(x) print *,"### ",x,' ###    ',__FILE__," > "
#define error(x) print *,"### ERROR ",x,__FILE__," > "
#define debug(x,idx) print *,"### DEBUG ",idx," @ line ",x,__FILE__," > "

  !  subroutine IO_read_parametrization_file(param)
  !  subroutine IO_file_info(file)
  !  subroutine lammps(cells)
  !  subroutine lammpsold(filename,cell)
  !  subroutine read_cp2k_input_file(file,cell)
  !  subroutine IO_read_basis_set()
  !  subroutine IO_read_file(file,univers)
  !  subroutine IO_read_parametrization_file(param)
  !  subroutine IO_read_molecule(file,molecule,idx_conf_to_read)
  !  subroutine IO_read_molecule_from_cp2k_restart_file(file,molecule)
  !  subroutine IO_read_molecule_from_xyz_file(file,molecule,idx_conf_to_read)
  !  subroutine IO_read_cell_from_cp2k_restart_file(file,cell)
  
  !  subroutine IO_write(file,univ,iconf)
  !  subroutine IO_write_CP2K_input_file(file,univ,iconf)
  
  !  subroutine QuantumEspresso(file,univ,iconf)
  !  subroutine read_multixyz(file,univers)
  !  subroutine save(filename,univers,idx_conf_to_save)
  !  subroutine save_molecule(molecule)
  !  subroutine xsf(filename,cell)
  !  subroutine IO_write_cell_to_xyz(filename,cell)

  ! --------------------------------------------------------------------------------------
  !
  !          subroutine IO_read_parametrization_file(param)
  !
  ! --------------------------------------------------------------------------------------
  subroutine IO_read_parametrization_file(param)
    use global
    implicit none
    type(t_Param)::param
    character(len=256)::name
    integer::io
    character (len=1024)::line
    character (len=NCHARFIELD)::field(32)
    integer::nfield,i,j,k
    logical::bool_input,bool_output,bool_machine

    param%n_inputs=0
    param%n_outputs=0
    bool_input=.FALSE.
    bool_output=.FALSE.
    bool_machine=.FALSE.

    !
    ! cette partie permet de compter le nombre de blocs d'input & d'output
    !
    
    io=0
    open(unit=1,file=param%filename,form='formatted')
    do while(io.eq.0)
       read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
       if(nfield.gt.0) then

          if(field(1).eq."&run_type")         param%run_type=field(2)
          if(field(1).eq."&system")           param%calculation%system=field(2)  
          if(field(1).eq."&input")            bool_input=.TRUE.
          if(field(1).eq."&end_input")        bool_input=.FALSE.
          if(field(1).eq."&output")           bool_output=.TRUE.
          if(field(1).eq."&end_output")       bool_output=.FALSE.
          if(field(1).eq."&machine")          bool_machine=.TRUE.
          if(field(1).eq."&end_machine")      bool_machine=.FALSE.
          write(*,*) trim(line)

          if((bool_input).and.(trim(field(1)).eq."&name"))     param%n_inputs=param%n_inputs+1
          if((bool_output).and.(trim(field(1)).eq."&name"))     param%n_outputs=param%n_outputs+1
          if((bool_machine).and.(trim(field(1)).eq."&name"))     param%machine=field(2)
       end if
    end do

    allocate(param%input(param%n_inputs))
    allocate(param%output(param%n_outputs))

    !
    ! lecture des paramètres
    !
    
    rewind(1)
    j=0
    k=0
    io=0
    do while(io.eq.0)
       read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
       if(nfield.gt.0) then
          ! ----------------------------------------------------------------
          if(field(1).eq."&input")  then
             k=k+1
             bool_input=.TRUE.
          end if
          if(field(1).eq."&end_input")          bool_input=.FALSE.
          if((bool_input).and.(trim(field(1)).eq."&name"))     then
             param%input(k)%name=field(2)
          end if
          if((bool_input).and.(trim(field(1)).eq."&type"))     then
             param%input(k)%type=field(2)
          end if
          if((bool_input).and.(trim(field(1)).eq."&config"))     then
             read(field(2),*) param%input(k)%iconf
          end if
          ! ----------------------------------------------------------------
          if(field(1).eq."&output")  then
             j=j+1
             bool_output=.TRUE.
          end if
          if(field(1).eq."&end_output")          bool_output=.FALSE.
          if((bool_output).and.(trim(field(1)).eq."&name"))     then
             param%output(j)%name=field(2)
          end if
          if((bool_output).and.(trim(field(1)).eq."&type"))     then
             param%output(j)%type=field(2)
          end if


          
       end if
    end do
    close(1)

    msg(__LINE__), " n_inputs= ",param%n_inputs
    do i=1,param%n_inputs
       msg(__LINE__), " input ",i,":",trim(param%input(i)%name)," (",trim(param%input(i)%type)," )"
    end do

    msg(__LINE__), " n_outpus= ",param%n_outputs
    do i=1,param%n_outputs
       msg(__LINE__), " output ",i,":",&
            trim(param%output(i)%name)," (",trim(param%output(i)%type)," )"
    end do

  end subroutine IO_read_parametrization_file

  ! --------------------------------------------------------------------------------------
  !
  !              IO_file_info()
  !
  ! --------------------------------------------------------------------------------------
  subroutine IO_file_info(name,file)
    implicit none
    character(len=2048)::name
    type(t_Param)::file
    integer::io
    character (len=1024)::line
    character (len=NCHARFIELD)::field(32)
    integer::nfield,i,j

    io=0
    file%n_lines=0
    file%n_configurations=0
    open(unit=1,file=name,form='formatted')


    do while (io==0)
       read(1,'(A)',iostat=io) line;       if(io.eq.-1) exit  ;file%n_lines=file%n_lines+1

       call line_parser(line,nfield,field)
       read(field(1),*) file%n_atoms                          ;file%n_lines=file%n_lines+1
       read(1,'(A)',iostat=io) line;       if(io.eq.-1) exit  ;file%n_lines=file%n_lines+1
       file%n_configurations=file%n_configurations+1
       do i=1,file%n_atoms
          read(1,'(A)',iostat=io) line;       if(io.eq.-1) exit  ;file%n_lines=file%n_lines+1
       end do
    enddo

    print *,"# IO_file_info> File ",trim(name)
    print *,"# IO_file_info> number of line(s)=",file%n_lines
    print *,"# IO_file_info> number of atom(s)=",file%n_atoms
    print *,"# IO_file_info> number of configurations=",file%n_configurations

    allocate(file%nrj(file%n_configurations))
    rewind(1)
    do i=1,file%n_configurations
       read(1,'(A)',iostat=io) line
       read(1,'(A)',iostat=io) line
       call line_parser(line,nfield,field)
       if(nfield>0) then
          read(field(6),*) file%nrj(i)
       end if
       do j=1,file%n_atoms
          read(1,'(A)',iostat=io) line
       end do
    end do
    close(1)
    open(unit=1,file="nrj.dat",form='formatted',status='unknown')
    do i=1,file%n_configurations
       write(1,*) i,file%nrj(i)
    end do
    close(1)
    print *,"# IO_file_info> End of IO_file_info()"
  end subroutine IO_file_info
  ! --------------------------------------------------------------------------------------
  !
  !              lammps()
  !
  ! --------------------------------------------------------------------------------------

  subroutine lammps(output,univ)
    use global
    implicit none
    type(t_File)::output
    type(t_Univers)::univ
    integer::i,j
    integer::ibond,idxi,idxj,imol,jmol,iat
    character(len=4)::type_bond

    msg(__LINE__), " Saving lammps input file",trim(output%name)
    msg(__LINE__)," n_bonds=", univ%configurations(output%iconf)%cells%n_bonds

    open(unit=1,file=output%name,form='formatted')

    write(1,*) "  LAMMPS 'data.' description "
    write(1,*)
    write(1,'(I8.2,X,A)') univ%configurations(output%iconf)%cells%n_atoms,"atoms"
    write(1,'(I8,X,A)') univ%configurations(output%iconf)%cells%n_elements,"atom types"

    !write(1,*) univ%configurations(iconf)%cells%L(
    
    write(1,*)    
    write(1,*) 'Atoms'
    write(1,*) ''
    do imol=1,univ%configurations(output%iconf)%cells%n_molecules
       do iat=1,univ%configurations(output%iconf)%cells%molecules(imol)%n_atoms
          write(1,'(I5,1X,A,3(F12.6))') &
               univ%configurations(output%iconf)%cells%molecules(imol)%atoms(iat)%idx,&
               univ%configurations(output%iconf)%cells%molecules(imol)%atoms(iat)%elt,&
               (univ%configurations(output%iconf)%cells%molecules(imol)%atoms(iat)%q(j),j=1,3)
       end do
    end do

    write(1,*)    
    write(1,*) 'Bonds'
    write(1,*) ''
    do ibond=1,univ%configurations(output%iconf)%cells%n_bonds
       idxi=univ%configurations(output%iconf)%cells%bonds(ibond)%list_atoms(1)%idx
       idxj=univ%configurations(output%iconf)%cells%bonds(ibond)%list_atoms(2)%idx
       imol=univ%configurations(output%iconf)%cells%bonds(ibond)%list_atoms(1)%idx_molecule
       jmol=univ%configurations(output%iconf)%cells%bonds(ibond)%list_atoms(2)%idx_molecule
       type_bond=trim(univ%configurations(output%iconf)%cells%molecules(imol)%atoms(idxi)%elt)//&
            trim(univ%configurations(output%iconf)%cells%molecules(jmol)%atoms(idxj)%elt)
       write(1,*) univ%configurations(output%iconf)%cells%bonds(ibond)%idx,&
            type_bond,&
            univ%configurations(output%iconf)%cells%bonds(ibond)%list_atoms(1)%idx,&
            univ%configurations(output%iconf)%cells%bonds(ibond)%list_atoms(2)%idx
    end do
    write(1,*) ''
    !do imol=1,cells%n_molecules
    !   do iat=1,cells%molecules(imol)%n_atoms
    !      print *,"# LAMMPS> ",cells%molecules(imol)%atoms(iat)%n_bonds
    !   end do
    !end do
    close(1)
  end subroutine lammps
  ! --------------------------------------------------------------------------------------
  !
  !              lammpsold()
  !
  ! --------------------------------------------------------------------------------------
  subroutine lammpsold(filename,cell)
    implicit none
    type(t_Cell)::cell
    character(len=1024)::filename,pref
    integer::idx
    
    idx=scan(filename,".")
    pref=trim(filename(1:idx-1))//".input.lmp"
    !print *,trim(adjustl(pref))
    !str = adjustl(str)

    !
    ! lammps.input.lmp
    !
    open(unit=1,file=trim(pref),form='formatted',status='unknown')    
    write(1,*) "   # -*- lammps -*-"
    write(1,*) ''
    write(1,*) 'variable filename string H2O'
    write(1,*) ''
    write(1,*) 'variable etol equal 0.0'
    write(1,*) 'variable ftol equal 1.0e-8'
    write(1,*) 'variable maxiter equal 100000'
    write(1,*) 'variable maxeval equal 100000'
    write(1,*) ''
    write(1,*) 'variable freq_save_traj equal 1'
    write(1,*) '#variable freq_print_screen 1'
    write(1,*) ''
    write(1,*) '# ------------------------ INITIALIZATION ----------------------------'
    write(1,*) '#units real'
    write(1,*) 'units metal'
    write(1,*) 'dimension 3'
    write(1,*) 'boundary   p p p'
    write(1,*) 'atom_style    full'
    write(1,*) 'read_data ${filename}.data.lmp'
    write(1,*) ''
    write(1,*) '# ------------------------ FORCE FIELDS/Interatomic potential ------------------------------'
    write(1,*) ''
    write(1,*) 'include ${filename}.ff.lmp'
    write(1,*) ''
    write(1,*) '# Compute and print thermodynamic info (e.g. temperature, energy, pressure) on timesteps that are a multiple '
    write(1,*) '# of N and at the beginning and end of a simulation. '
    write(1,*) 'thermo 1     # sauvegarde des info thermique toutes les 1 steps'
    write(1,*) 'thermo_style custom step lx ly lz press pxx pyy pzz pe etotal bonds  # format de ce que l on veut afficher'
    write(1,*) ''
    write(1,*) 'dump 1 all  custom ${freq_save_traj} dyn.* id type element x y z fx fy fz '
    write(1,*) ''
    write(1,*) '# minimize etol ftol maxiter maxeval'
    write(1,*) 'minimize ${etol} ${ftol} ${maxiter} ${maxeval}'
    close(1)
    !
    ! lammps.ff.lmp
    !
    pref=trim(filename(1:idx-1))//".ff.lmp"
    print *,trim(pref)
    open(unit=1,file=trim(pref),form='formatted',status='unknown')    
    write(1,*) '# Berendsen et al, in "Intermolecular forces", p. 331 (1981)'
    write(1,*) '# Charges and geometry are specified in the "data." file.'
    write(1,*) ''
    write(1,*) 'mass 1 1.00794 # H'
    write(1,*) 'mass 2 15.9994 # O'
    write(1,*) ''
    write(1,*) 'pair_style lj/cut/coul/long 10.0'
    write(1,*) 'pair_modify tail yes'
    write(1,*) 'kspace_style pppm 1.0e-5'
    write(1,*) ''
    write(1,*) 'pair_coeff 1 1 0.00000 0.000'
    write(1,*) 'pair_coeff 1 2 0.00000 0.000'
    write(1,*) 'pair_coeff 2 2 0.15535 3.166'
    write(1,*) ''
    write(1,*) 'bond_style harmonic'
    write(1,*) 'bond_coeff 1 0.0 1.0'
    write(1,*) ''
    write(1,*) 'angle_style harmonic'
    write(1,*) 'angle_coeff 1 0.0 109.47'
    close(1)
    !
    ! lammps.data.lmp
    !
    pref=trim(filename(1:idx-1))//".data.lmp"
    print *,trim(pref)
    open(unit=1,file=trim(pref),form='formatted',status='unknown')

    write(1,*) "LAMMPS 'data.' description "
    write(1,*) ''
    write(1,*) '       3 atoms'
    write(1,*) '       2 bonds'
    write(1,*) '       1 angles'
    write(1,*) ''
    write(1,*) '       2 atom types'
    write(1,*) '       1 bond types'
    write(1,*) '       1 angle types'
    write(1,*) ''
    write(1,*) '    -10 10.1      xlo xhi'
    write(1,*) '    -10 10.1      ylo yhi'
    write(1,*) '    -10 10.1      zlo zhi'
    write(1,*) ''
    write(1,*) 'Atoms'
    write(1,*) ''
    write(1,*) '      1    1  2 -0.82    1.55000    1.55000    1.50000'
    write(1,*) '      2    1  1  0.41    1.55000    2.36649    2.07736'
    write(1,*) '      3    1  1  0.41    1.55000    0.73351    2.07736'
    write(1,*) ''
    write(1,*) 'Bonds'
    write(1,*) ''
    write(1,*) '      1   1      1      2'
    write(1,*) '      2   1      1      3'
    write(1,*) ''
    write(1,*) 'Angles'
    write(1,*) ''
    write(1,*) '      1   1      2      1      3'


    
    close(1)
  end subroutine lammpsold
  ! --------------------------------------------------------------------------------------
  !
  !              read_cp2k_input_file()
  !
  ! --------------------------------------------------------------------------------------
  subroutine read_cp2k_input_file(name,cell)
    implicit none
    !    type(t_Param)::file
    character(len=2048)::name
    integer::i,io
    character (len=1024)::line
    character (len=NCHARFIELD)::field(32)
    integer::nfield
    type(t_Cell)::cell

    open(unit=1,file=name,form='formatted')

    i=0
    io=0
    do while(io==0)
       read(1,'(A)',iostat=io) line
       i=i+1
       call line_parser(line,nfield,field)
       if(field(1).eq.'ABC') then
          !print *,nfield,' --> ',(field(i),i=1,nfield)
          read(field(3),*) cell%L(1,1) ; cell%L(1,2)=0.0 ; cell%L(1,3)=0.0
          read(field(4),*) cell%L(2,2) ; cell%L(2,1)=0.0 ; cell%L(2,3)=0.0
          read(field(5),*) cell%L(3,3) ; cell%L(3,2)=0.0 ; cell%L(3,1)=0.0
       end if
    end do
    close(1)

  end subroutine read_cp2k_input_file

  ! --------------------------------------------------------------------------------------
  !
  !              IO_read_cell_from_cp2k_restart_file()
  !
  ! --------------------------------------------------------------------------------------
  subroutine IO_read_cell_from_cp2k_restart_file(param,input,cell)
    use global
    use Element
    use Molecule
    use Atom
    implicit none
    type(t_File)::input
    type(t_Param)::param
    type(t_Cell)::cell
    character (len=1024)::line
    character (len=NCHARFIELD)::field(32)
    integer::nfield
    integer::iconf,iat,io,k,i,j,l
    integer::Zelt
    if(allocated(cell%molecules))    deallocate(cell%molecules)
    allocate(cell%molecules(1))
    cell%molecules(1)=MOLECULE_init()
    
    io=0
    print *,"# Reading ",trim(input%name)
    open(unit=1,file=input%name,form='formatted')
    do while((io.eq.0).and.(.not.(cell%molecules(1)%n_atoms.gt.0)))
       read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
       if(nfield.gt.0) then
          if(field(1).eq."NUMBER_OF_ATOMS") then
             read(field(2),*) cell%molecules(1)%n_atoms
          end if
          if(field(1).eq."CUTOFF") then
             read(field(2),*) param%calculation%ecutoff
          end if
          if(field(1).eq."RUN_TYPE") then
             read(field(2),*) param%calculation%run_type
          end if
          if(field(1).eq."ADDED_MOS") then
             read(field(2),*) param%calculation%added_mos
          end if
          if(field(1).eq."BASIS_SET_FILE_NAME") then
             call get_last_field(field(2),param%calculation%basis_name)
          end if
          if(field(1).eq."POTENTIAL_FILE_NAME") then
             call get_last_field(field(2),param%calculation%potential_name)
          end if
          if(field(1).eq."&XC_FUNCTIONAL") then
             read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
             l=len(trim(field(1)))
             param%calculation%xc_functional=field(1)(2:l)
             !read(field(2),*) param%calculation%xc_functional
          end if
          if(field(1).eq."&KIND") then
             Zelt=elt2Z(field(2))
             read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
             PeriodicTable(Zelt)%Basis_set=field(2)
             read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
             read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
             PeriodicTable(Zelt)%Potential_set=field(2)
             msg(__LINE__)," <",trim(PeriodicTable(Zelt)%Basis_set),">"
             msg(__LINE__)," <",trim(PeriodicTable(Zelt)%Potential_set),">"

          end if
       end if
    end do
    print *,"# Number of atoms: ",cell%molecules(1)%n_atoms
    print *,"# Functional:", param%calculation%xc_functional
    






    
    if(allocated(cell%molecules(1)%atoms)) deallocate(cell%molecules(1)%atoms)
    allocate(cell%molecules(1)%atoms(cell%molecules(1)%n_atoms))
    REWIND (1, IOSTAT=IO)
    read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
    do while(.not.(field(1).eq."&COORD"))
       read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
    end do
    do iat=1,cell%molecules(1)%n_atoms
       read(1,'(A)',iostat=io) line
       call line_parser(line,nfield,field)
       cell%molecules(1)%atoms(iat)=ATOM_new()
       cell%molecules(1)%atoms(iat)%elt=field(1)
       cell%molecules(1)%atoms(iat)%Zato=elt2Z(cell%molecules(1)%atoms(iat)%elt)
       read(field(2),*) cell%molecules(1)%atoms(iat)%q(1)
       read(field(3),*) cell%molecules(1)%atoms(iat)%q(2)
       read(field(4),*) cell%molecules(1)%atoms(iat)%q(3)
       !print *,(cell%molecules(1)%atoms(iat)%q(k),k=1,3)
       cell%molecules(1)%atoms(iat)%idx_molecule=-1
       cell%molecules(1)%atoms(iat)%idx=iat
       cell%molecules(1)%atoms(iat)%n_bonds=0
    end do



    REWIND (1, IOSTAT=IO)
    read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
    do while((.not.(field(1).eq."LIST")).and.(io.eq.0))
       read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
    end do
    if(field(1).eq."LIST") then
       print *,trim(line),nfield
       cell%molecules(1)%constraints%n_XYZ=nfield-2
       do while(.not.(field(1).eq."&END"))
          read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
          print *,trim(line),nfield
          if(nfield.eq.11) then
             cell%molecules(1)%constraints%n_XYZ=cell%molecules(1)%constraints%n_XYZ+10
          else
             if(.not.(field(1).eq."&END")) cell%molecules(1)%constraints%n_XYZ=&
                  cell%molecules(1)%constraints%n_XYZ+nfield
          end if
       end do
    end if
    print *,cell%molecules(1)%constraints%n_XYZ
    if(allocated(cell%molecules(1)%constraints%XYZ)) deallocate(cell%molecules(1)%constraints%XYZ)
    allocate(cell%molecules(1)%constraints%XYZ(cell%molecules(1)%constraints%n_XYZ))
    REWIND (1, IOSTAT=IO)
    read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
    do while((.not.(field(1).eq."LIST")).and.(io.eq.0))
       read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
    end do
    j=1
    if(field(1).eq."LIST") then
       do i=2,nfield-1
          read(field(i),*) cell%molecules(1)%constraints%XYZ(j) ; j=j+1
       end do
       do while(.not.(field(1).eq."&END"))
          read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);

          if(nfield.eq.11) then
             do i=1,10
                read(field(i),*) cell%molecules(1)%constraints%XYZ(j) ; j=j+1
             end do
          else
             if(.not.(field(1).eq."&END")) then
                do i=1,nfield
                   read(field(i),*) cell%molecules(1)%constraints%XYZ(j) ; j=j+1
                end do
             end if
          end if
       end do
    end if

    
    
    REWIND (1, IOSTAT=IO)
    read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
    do while(.not.(field(1).eq."&CELL"))
       read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
    end do
    do k=1,3
       read(1,'(A)',iostat=io) line
       call line_parser(line,nfield,field)
       do i=1,3
          read(field(i+1),*) cell%L(k,i)
       end do
    end do
    read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
    read(field(2),'(A)') cell%periodicity
    
    close(1)


    call MOLECULE_element_list_update(cell%molecules(1))

    cell%n_atoms=cell%molecules(1)%n_atoms
    cell%n_molecules=1
    
  end subroutine IO_read_cell_from_cp2k_restart_file


  ! --------------------------------------------------------------------------------------
  !
  !              IO_write()
  !
  ! --------------------------------------------------------------------------------------
  subroutine IO_write(param,univ,iconf)
    implicit none
    type(t_Param)::param
    type(t_Univers)::univ
    integer::iconf,i,j

    !print *,"# PERIODICITY= ",univ%configurations(iconf)%cells%periodicity
    !print *,"# n_outputs= ",param%n_outputs
    !do i=1,param%n_outputs
    !   print *,"# output(",i,")=",param%output(i)%name," ",trim(param%output(i)%type)
    !end do

    if(param%n_outputs.gt.0) then
       do i=1,param%n_outputs
          msg(__LINE__)," output(",i,"): ", &
               trim(param%output(i)%name)," ", trim(param%output(i)%type)
          select case (trim(param%output(i)%type))
          case('cp2k_input')
             call IO_write_CP2K_input_file(param,univ,iconf,i)
          case('lammps_input')
             call lammps(param%output(i),univ)
          case('qe_input')
             call QuantumEspresso(param%output(i),univ)
          case ('xsf')
             msg(__LINE__), " Saving xsf file",trim(param%output(i)%name)
             call xsf(param%output(i)%name,univ%configurations(iconf)%cells)
          case ('xyz')
             msg(__LINE__), " Saving xyz file",trim(param%output(i)%name)
             call IO_write_cell_to_xyz(param%output(i)%name,univ%configurations(iconf)%cells)
          case default
             error(__LINE__)," File type ",param%output(i)%type," is not (yet) implemented"
             stop
          end select
       end do
    end if
    
  end subroutine IO_write
  ! --------------------------------------------------------------------------------------
  !
  !        IO_write_CP2K_input_file()
  !
  ! --------------------------------------------------------------------------------------
  subroutine IO_write_CP2K_input_file(param,univ,iconf,ioutput)
    use Element
    use Molecule
    use Machine
    implicit none
    !type(t_Param_Calculation)::param
    type(t_Param)::param
    type(t_Univers)::univ
    integer::iconf,i,j,ioutput,idx_machine,idx

    
    msg(__LINE__)," writing CP2K input file in ",trim(param%output(ioutput)%name)    
    msg(__LINE__), "machine=",param%machine
    idx_machine=MACHINE_search_idx(param)
    msg(__LINE__), "idx_machine=",idx_machine


    
    ! structure d'un fichier d'input pour CP2K
    ! https://manual.cp2k.org/cp2k-9_1-branch/CP2K_INPUT.html
    ! GLOBAL
    !    RUN_TYPE
    !       CELL_OPT              -> Cell optimization. Both cell vectors and atomic positions are optimised.
    !       EHRENFEST_DYN         -> Ehrenfest dynamics (using real time propagation of the wavefunction)
    !       ENERGY ou WAVEFUNCTION_OPTIMIZATION               -> Computes energy
    !       FORCE_ENERGY          -> Computes energy and forces
    !       GEOMETRY_OPTIMIZATION -> Geometry Optimization
    !       LINEAR_RESPONSE       -> Linear Response
    !       MONTECARLO            -> Monte Carlo
    !       MOLECULAR_DYNAMICS    -> Molecular Dynamics
    !       NEGF                  -> Non-equilibrium Green's function method
    !       PINT                  -> Path integral
    !       SPECTRA               -> Computes absorption Spectra
    ! FORCE_EVAL
    !       DFT
    !         KPOINTS
    !       MM
    !       PROPERTIES
    !       PW_DFT
    !       QMMM
    !       SUBSYS
    ! MOTION
    
    
    
    open(unit=1,file=param%output(ioutput)%name,form='formatted')
    
    write(1,*) "@SET SYSTEM ",param%calculation%system
    write(1,*) "@set LIBDIR ",trim(List_Machines(idx_machine)%cp2k_libdir)

    !if(param%machine.eq."hpc") then
    !   write(1,*) " @set LIBDIR /home2020/home/ipcms/bulou/src/cp2k/data/"
    !end if
    !if(param%machine.eq."jeanzay") then
    !   write(1,*) " @set LIBDIR /gpfslocalsup/spack_soft/cp2k/8.2/intel-19.1.3-2l5irg7v7rchgpfxbgo66hp2rizcvo5f/share/data"
    !end if
    !write(1,*) " #@set LIBDIR /home/bulou/src/cp2k/data"
    !write(1,*) " @set LIBDIR /usr/local/cp2k/cp2k-6.1.0.g8/data"
    
    
    write(1,*) " @set ECUTOFF ",param%calculation%ecutoff
    write(1,*) " @set EPS_SCF    1.0e-10"
    write(1,*) " @set MAX_SCF    500"
    write(1,*) " @set MAX_ITER   2000"
    if(param%calculation%scheme.eq.'MONKHORST-PACK') then
       write(1,*) " @set NKX ",param%calculation%k(1)
       write(1,*) " @set NKY ",param%calculation%k(2)
       write(1,*) " @set NKZ ",param%calculation%k(3)
    end if
    write(1,*) " @set MO ",param%calculation%added_mos
    write(1,*) " @set BASIS_NAME ",param%calculation%basis_name
    write(1,*) " @set POTENTIAL_NAME ",param%calculation%potential_name
    write(1,*) " @set RUN_TYPE ",param%calculation%run_type  !   GEO_OPT CELL_OPT WAVEFUNCTION_OPTIMIZATION 
    !       write(1,*) " @set RUN_TYPE   CELL_OPT "
    !       write(1,*) " @set RUN_TYPE   WAVEFUNCTION_OPTIMIZATION "
    write(1,*) " @set ALGORITHM STANDARD "   ! DAVIDSON, FILTER_MATRIX, LANCZOS, OT, STANDARD
    ! --------------------------------------------------------------------------------------------------
    ! -------------- GLOBAL part -----------------------------------------------------------------------
    !  contains general options for the CP2K run, such as the name of the job, the type of run etc.
    !       -> for a calculation needing energy and force, we must set ''RUN_TYPE'' to ENERGY_FORCE.
    ! The keyword ''METHOD'' chooses the method for evaluating the forces on atoms to QUICKSTEP,
    !      i.e. Density Functional Theory using the Gaussian and Planewaves (GPW) method. 
    write(1,*) " &GLOBAL"
    write(1,*) "  PROJECT ${SYSTEM}"
    write(1,*) "  RUN_TYPE ${RUN_TYPE}"   ! ENERGY_FORCE ENERGY
    write(1,*) "  PRINT_LEVEL MEDIUM"
    write(1,*) "&END GLOBAL"
    
    ! --------------------------------------------------------------------------------------------------
    ! -------------- FORCE_EVAL part -------------------------------------------------------------------
    !  contains all parameters associated with the evaluation of forces on atoms, this includes the initial atomic coordinates.
    write(1,*) "&FORCE_EVAL"
    write(1,*) "  METHOD Quickstep"
    ! The subsection ''SUBSYS'' defines the simulation unit cell and the initial coordinates of atoms in the calculation
    write(1,*) "  &SUBSYS"
    ! The subsection ''KIND'' gives definitions of elements in the calculation.
    ! There must be one KIND subsection per element.
    ! In this example, for Si, we have defined the basis set to be used: DZVP-GTH-PADE (double-ζ with polarisation
    ! basis optimised for Geodecker-Teter-Hutter PADE LDA pseudopotential);
    ! and the pseudopotential: GTH-PADE-q4 (Geodecker-Teter-Hutter PADE LDA pseudopotential with 4 valence electrons).
    ! &KIND Si
    ! ELEMENT   Si
    ! BASIS_SET DZVP-GTH-PADE
    ! POTENTIAL GTH-PADE-q4
    ! &END KIND
    ! The basis set and pseudopotential names must correspond to an existing entry in the corresponding basis set
    ! and pseudopotential files defined by ''BASIS_SET_FILE_NAME'' and ''POTENTIAL_FILE_NAME'' keywords in ''DFT'' subsection,
    ! in FORCE_EVAL section. The chosen basis for Si corresponds to parameters:
    !
    ! Si DZVP-GTH-PADE
    !   2      
    !  3  0  1  4  2  2
    !        1.2032422345   0.3290350445   0.0000000000   0.0474539126   0.0000000000
    !        0.4688409786  -0.2533118323   0.0000000000  -0.2594473573   0.0000000000
    !        0.1679863234  -0.7870946277   0.0000000000  -0.5440929303   0.0000000000
    !        0.0575619526  -0.1909898479   1.0000000000  -0.3624010364   1.0000000000
    !  3  2  2  1  1
    !        0.4500000000   1.0000000000
    !
    !in file BASIS_SET; and the chosen pseudopotential corresponds to parameters:
    !
    !Si GTH-PADE-q4 GTH-LDA-q4
    !    2    2
    !     0.44000000    1    -7.33610297
    !    2
    !     0.42273813    2     5.90692831    -1.26189397
    !                                        3.25819622
    !     0.48427842    1     2.72701346
    !
    !in file GTH_POTENTIALS. 

    do i=1,univ%configurations(iconf)%cells%molecules(1)%n_elements

       write(1,*) "    &KIND ",univ%configurations(iconf)%cells%molecules(1)%elements(i)%elt
       write(1,*) "      ELEMENT   ",univ%configurations(iconf)%cells%molecules(1)%elements(i)%elt
       write(1,*) "        BASIS_SET ",PeriodicTable(elt2Z(univ%configurations(iconf)%cells%molecules(1)%elements(i)%elt))%Basis_set
       write(1,*) "        POTENTIAL ",&
            PeriodicTable(elt2Z(univ%configurations(iconf)%cells%molecules(1)%elements(i)%elt))%Potential_set
       ! if(univ%configurations(iconf)%cells%molecules(1)%elements(i)%elt.eq."H") then
       !    !          write(1,*) "      BASIS_SET DZV-GTH-PADE"
       !    write(1,*) "      POTENTIAL GTH-PADE-q1"
       ! end if
       ! if(univ%configurations(iconf)%cells%molecules(1)%elements(i)%elt.eq."O") then
       !    !write(1,*) "      BASIS_SET DZVP-GTH-PADE"
       !    write(1,*) "      POTENTIAL GTH-PADE-q6"
       ! end if
       ! if(univ%configurations(iconf)%cells%molecules(1)%elements(i)%elt.eq."Ti") then
       !    !write(1,*) "      BASIS_SET DZV-GTH-PADE"
       !    write(1,*) "      POTENTIAL GTH-PADE-q12"
       ! end if
       write(1,*) "    &END KIND"
    end do

    
    !  The subsection ''CELL'' defines the simulation unit cell used in a calculation.
    ! In this example, we define the unit cell as cubic, with lattice constant equal to 5.4306975 Angstroms.
    ! “Angstrom” is the default unit for cell vectors. ''A'', ''B'' and ''C'' are the first,
    ! second and third lattice (cell) vectors. There are many ways to define the cell, see ''CP2K'' input reference manual for more details.
    !
    !The initial atomic coordinates are specified in the ''COORD'' subsection. The default input format for atomic coordinates in CP2K is:
    !
    !<ATOM_KIND> X Y Z
    !
    !where X, Y and Z are Cartesian coordinates in Angstroms. This can be changed by configuring keyword ''SCALED'' to .TRUE.,
    !in the COORD subsection, which makes the coordinate input X Y Z to be fractional with respect to the lattice vectors.
    !One can also change the unit for the Cartesian coordinates by setting the keyword ''UNIT'' with in the subsection. <ATOM_KIND>
    !should be a label that corresponds to the definition of the elements in the ''KIND'' subsections.
    

    write(1,*) "    &CELL"
    write(1,*) "      A [angstrom] ",(univ%configurations(iconf)%cells%L(1,j),j=1,3)
    write(1,*) "      B [angstrom]",(univ%configurations(iconf)%cells%L(2,j),j=1,3)
    write(1,*) "      C [angstrom]",(univ%configurations(iconf)%cells%L(3,j),j=1,3)
    write(1,*) "      PERIODIC ",univ%configurations(iconf)%cells%periodicity
    write(1,*) "    &END CELL"

    call MOLECULE_sort_z(univ%configurations(iconf)%cells%molecules(1))
    write(1,*) "    &COORD"
    do i=1,univ%configurations(iconf)%cells%molecules(1)%n_atoms
       idx=univ%configurations(iconf)%cells%molecules(1)%order_idx_save(i)
       !write(*,*) idx
       write(1,*) univ%configurations(iconf)%cells%molecules(1)%atoms(idx)%elt,&
            (univ%configurations(iconf)%cells%molecules(1)%atoms(idx)%q(j),j=1,3)
       !               (univ%configurations(iconf)%cells%molecules(1)%atoms(i)%constraints(j),j=1,3)
    end do
    write(1,*) "    &END COORD"
    write(1,*) "  &END SUBSYS"



    !After the SUBSYS section  follows the ''DFT'' subsection, which controls all aspects of the self-consistent
    ! Kohn-Sham Density Functional Theory calculation.
    ! This subsection is only relevant if and only if the METHOD keyword in FORCE_EVAL is set to QUICKSTEP.        
    write(1,*) "  &DFT"
    write(1,*) "    BASIS_SET_FILE_NAME  ${LIBDIR}/${BASIS_NAME}"
    write(1,*) "    POTENTIAL_FILE_NAME  ${LIBDIR}/${POTENTIAL_NAME}"


    if(univ%configurations(iconf)%cells%periodicity.eq."XYZ") then
       ! KPOINT  	KPOINT {Real} {Real} {Real} {Real}	Specify kpoint coordinates and weight.
       ! SCHEME {word} SCHEME {Word} ...  Kpoint scheme to be used.
       !Available options: NONE, GAMMA, MONKHORST-PACK, MACDONALD, and GENERAL.
       !For MONKHORST-PACK and MACDONALD the number of k points in all 3 dimensions has to be supplied along with the keyword. E.g.
       !
       !MONKHORST-PACK 12 12 8.
       write(1,*) "    &KPOINTS"
       select case(param%calculation%scheme)
          case('MONKHORST-PACK')
             write(1,*) "       SCHEME MONKHORST-PACK ${NKX} ${NKY} ${NKZ}"
          case('gamma')
             write(1,*) "       SCHEME GAMMA"
          case default
             print *,"# unknown kpoint scheme: ",param%calculation%scheme
             stop
          end select
       write(1,*) "       SYMMETRY .FALSE."
       write(1,*) "    &END KPOINTS"
    end if

    ! The ''QS'' subsection contains general control parameters used by QUICKSTEP.
    ! ''EPS_DEFAULT'' sets the default value for all tolerances used within QUICKSTEP.
    ! The individual tolerances (EPS_*) can be set, and they will override the EPS_DEFAULT value. 
    write(1,*) "    &QS"
    write(1,*) "      EPS_DEFAULT ${EPS_SCF}"
    write(1,*) "    &END QS"

    ! The ''MGRID'' subsection defines how the integration grid used in QUICKSTEP calculations should be setup.
    ! QUICKSTEP uses a multi-grid method for representing Gaussian functions numerically on the grid.
    ! Narrow and sharp Gaussians are mapped onto a finer grid than wider and smoother Gaussians.
    ! In this case, we are telling the code to set up 4 levels of multi-grids, with the planewave cutoff
    ! of the finest grid set to be 300 Ry, and with the grid spacing underneath any Gaussian functions to
    ! be finer than the equivalent planewave cutoff of 60 Ry.
    ! The users should read the tutorial “Converging the CUTOFF and REL_CUTOFF” for details on how
    ! these parameters affect the grid constructed, and how to define a sufficient grid for their calculation.
    ! In this example, the grid defined has already been found to be sufficient for the energy and force calculation.
    !       &MGRID
    !  NGRIDS 4      -> 4 level of multi-grids
    !  CUTOFF 300    -> plane wave cutoff in Ry of the finest grid
    !  REL_CUTOFF 60 -> grid sapcing with the other grids
    !  &END MGRID
    write(1,*) "    &MGRID"
    write(1,*) "      NGRIDS 4"
    write(1,*) "      CUTOFF ${ECUTOFF}"
    write(1,*) "      REL_CUTOFF 60"
    write(1,*) "    &END MGRID"

    ! https://manual.cp2k.org/cp2k-9_1-branch/CP2K_INPUT/FORCE_EVAL/DFT/POISSON.html
    write(1,*) "    &POISSON"
    write(1,*) "      PERIODIC ",univ%configurations(iconf)%cells%periodicity
    if(univ%configurations(iconf)%cells%periodicity.eq."XYZ") then
       write(1,*) "      POISSON_SOLVER PERIODIC"
    else
       write(1,*) "      POISSON_SOLVER ANALYTIC"
    end if
    write(1,*) "    &END POISSON"

    ! https://manual.cp2k.org/cp2k-9_1-branch/CP2K_INPUT/FORCE_EVAL/DFT/XC.html
    ! Parameters needed for the calculation of the eXchange and Correlation potential
    ! &XC_FUNCTIONAL {Keyword}  ->	Shortcut for the most common functional combinations.
    !   B3LYP, BLYP, LDA, NONE, PADE, PBE, PBE0
    ! 
    write(1,*) "    &XC"
    write(1,*) "      &XC_FUNCTIONAL ",param%calculation%xc_functional
    write(1,*) "      &END XC_FUNCTIONAL"
    write(1,*) "    &END XC"


    ! The ''SCF'' subsection defines all the settings related to methods used to find a self-consistent
    ! solution of the Kohn-Sham DFT formalism.
    ! ''SCF_GUESS'' sets how the initial trial electron density function ρ(→r) is to be generated.
    ! In this example (ATOMIC), the initial density is to be generated using overlapping of atomic charge densities.
    ! A good starting point for the electron density in the self-consistency loop is important in obtaining a convergent
    ! result quickly. ''EPS_SCF'' sets the tolerance of the charge density residual. This overrides the EPS_DEFAULT value
    ! set in QS subsection. ''MAX_SCF'' sets the maximum number of self-consistency loops QUICKSTEP is allowed to perform
    ! for each ground-state energy calculation. 
    write(1,*) "    &SCF"
    write(1,*) "      SCF_GUESS ATOMIC"
    write(1,*) "      EPS_SCF 1.0E-6"
    write(1,*) "      MAX_SCF ${MAX_SCF}"

    ! The ''DIAGONALIZATION '' subsection tells the code to use the traditional diagonalisation method for finding
    ! the ground state Kohn-Sham energy and electron density. The subsection heading also takes an argument, and in this case
    ! is set to “ON”, which equivalent to “.TRUE.” or “T”, and indicates that the diagonalisation method is turned on. One can
    ! also omit the value of the subsection heading, which defaults to “.TRUE.”. The alternative to diagonalisation is to use the
    ! Orbital Transform (OT) method, in which case, the user should either delete the DIAGONALIZATION block or change “ON” to “OFF”
    ! (or “.FALSE.”), and add the ''OT'' subsection instead. The ''ALGORITHM'' keyword sets the algorithm to use for diagonalisation
    ! of the Kohn-Sham Hamiltonian. “STANDARD” means the standard LAPACK/SCALAPACK subroutines are to be used for diagonalisation.
    ! https://manual.cp2k.org/cp2k-9_1-branch/CP2K_INPUT/FORCE_EVAL/DFT/SCF/DIAGONALIZATION.html
    ! &DIAGONALIZATION {Logical}  -> Controls the activation of the diagonalization method
    ! ALGORITHM {Keyword}
    !    DAVIDSON       ->    Preconditioned blocked Davidson
    !    FILTER_MATRIX  ->    Filter matrix diagonalization
    !    LANCZOS        ->    Block Krylov-space approach to self-consistent diagonalisation
    !    OT             ->    Iterative diagonalization using OT method
    !    STANDARD       ->    Standard diagonalization: LAPACK methods or Jacobi
    write(1,*) "      &DIAGONALIZATION  ON"
    write(1,*) "        ALGORITHM ${ALGORITHM}"
    write(1,*) "      &END DIAGONALIZATION"

    ! The ''MIXING'' subsection contains all the parameters associated with charge mixing in a self-consistency calculation.
    ! The subsection also admits a value, which can be either .TRUE. (T) or .FALSE. (F), which switches charge mixing on or off.
    ! The default is .TRUE.. Note that this subsection only applies to the traditional diagonalisation method. The OT method uses
    ! a different approach for charge mixing, and is explained in other tutorials. The keyword ''ALPHA'' sets the mixing parameter;
    ! in this example 0.4 of the output density will be mixed with 0.6 of the input density to form the new input density in the next
    ! SCF iteration. The keyword ''METHOD'' sets the mixing method; in this case, we will use Broyden mixing. The keyword ''NBROYDEN''
    ! is an alias to the parameter NBUFFER, and it sets the number of histories to be used in the Broyden mixing algorithm. 
    write(1,*) "      &MIXING  T"
    write(1,*) "        METHOD BROYDEN_MIXING"
    write(1,*) "        ALPHA 0.4"
    write(1,*) "        BETA 1.5"
    write(1,*) "        NBROYDEN 8"
    write(1,*) "      &END MIXING"

    write(1,*) "      ADDED_MOS ${MO}" !   Extra MOs (ADDED_MOS) are required for smearing   
    write(1,*) "      &SMEAR ON"
    write(1,*) "        METHOD FERMI_DIRAC"
    write(1,*) "        ELECTRONIC_TEMPERATURE [K] 300"
    write(1,*) "      &END SMEAR"


    write(1,*) "    &END SCF"
    write(1,*) "  &END DFT"


    write(1,*) "  &PRINT"
    write(1,*) "    &FORCES ON"
    write(1,*) "    &END FORCES"
    write(1,*) "  &END PRINT"
    write(1,*) "&END FORCE_EVAL"

    if(param%calculation%run_type.eq."GEO_OPT") then
       write(1,*) "&MOTION"
       !write(1,*) "       #&FIXED_ATOMS"
       !write(1,*) "       #LIST                                1 2 3 4"
       !write(1,*) "       #COMPONENTS_TO_FIX    Z"
       !write(1,*) "       #&END FIXED_ATOMS"
       
       !write(1,*) "       #&FIXED_ATOMS"
       !write(1,*) "       #LIST                                 5 6 7 8"
       !write(1,*) "       #COMPONENTS_TO_FIX    YZ"
       !write(1,*) "       #&END FIXED_ATOMS"
       
       if(univ%configurations(iconf)%cells%molecules(1)%constraints%n_XYZ.gt.0) then
          write(1,*) "       &CONSTRAINT"
          write(1,*) "       &FIXED_ATOMS"
          write(1,*) "       LIST  ",(univ%configurations(iconf)%cells%molecules(1)%constraints%XYZ(j),&
               j=1,univ%configurations(iconf)%cells%molecules(1)%constraints%n_XYZ)
          write(1,*) "       COMPONENTS_TO_FIX    XYZ"
          write(1,*) "       &END FIXED_ATOMS"
          write(1,*) "       &END CONSTRAINT"
       end if
       
       write(1,*) "&END MOTION"
    end if

!	&FIXED_ATOMS
!	LIST 301
!	COMPONENTS_TO_FIX XY
!	&END FIXED_ATOMS



    
    write(1,*) "#&ext_restart"
    write(1,*) "#    restart_file_name ",trim(param%calculation%system),"-1.restart"
    write(1,*) "#&end ext_restart"


    close(1)
    !print *,"# iconf= ",iconf
    !print *,"# Periodicity ",univ%configurations(iconf)%cells%periodicity
    !print *,"# CP2K input file: ",trim(param%output(ioutput)%name)
    !print *,"# n_constraints_XYZ=",univ%configurations(iconf)%cells%molecules(1)%constraints%n_XYZ



    !---------------------------------------------------------------------------------------------------------------------
  end subroutine IO_write_CP2K_input_file
  ! --------------------------------------------------------------------------------------
  !
  !              IO_read_file()
  !
  ! --------------------------------------------------------------------------------------
  subroutine IO_read_file(name,type_file,univers)
    implicit none
    !type(t_Param)::file
    character(len=2048)::name
    character(len=32)::type_file
    type(t_Univers)::univers
    integer::iconf,iat


    select case (type_file)
    case('cp2k_input')
       print *,"# cp2k input file"
       do iconf=1,univers%n_configurations
          call read_cp2k_input_file(name,univers%configurations(iconf)%cells)
       end do
    case ('multixyz')
       print *,"# Multi xyz file"
       call read_multixyz(name,univers)
    case default
       print *,"# File type ",type_file," is not (yet) implemented"
       stop
    end select

  end subroutine IO_read_file
  ! --------------------------------------------------------------------------------------
  !
  !              IO_read_molecule()
  !
  ! --------------------------------------------------------------------------------------
  !  subroutine IO_read_molecule(name,type_file,molecule,idx_conf_to_read)
  function IO_read_molecule(input) result(molecule)
    use Element
    implicit none
    type(t_File)::input
    !character(len=2048)::name
    !character(len=32)::type_file
    type(t_Molecule)::molecule
    !integer::idx_conf_to_read
    
    select case(input%type)
       case('xyz')
          molecule=IO_read_molecule_from_xyz_file(input)
       case('cp2k_restart_file')
          molecule=IO_read_molecule_from_cp2k_restart_file(input)
       case default
          error(__LINE__),input%type," Unkown file format"
          error(__LINE__)," Valid type are 'xyz' or 'cp2k_restart_file' "
          stop
    end select
    msg(__LINE__),"n_atoms= ",molecule%n_atoms
  end function IO_read_molecule
  ! --------------------------------------------------------------------------------------
  !
  !              IO_read_molecule_from_cp2k_restart_file()
  !
  ! --------------------------------------------------------------------------------------
  function IO_read_molecule_from_cp2k_restart_file(input) result(mol)
    use global
    use Element
    use Molecule
    implicit none
    type(t_File)::input
    type(t_Molecule)::mol
    character (len=1024)::line
    character (len=NCHARFIELD)::field(32)
    integer::nfield
    integer::iconf,iat,io,k,i
    mol=MOLECULE_init()
    !molecule%n_atoms=-1
    io=0
    open(unit=1,file=input%name,form='formatted')
    do while((io.eq.0).and.(mol%n_atoms<0))
       read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
       if(nfield.gt.0) then
          if(field(1).eq."NUMBER_OF_ATOMS") then
             read(field(2),*) mol%n_atoms

          end if
       end if
    end do

    allocate(mol%atoms(mol%n_atoms))
    REWIND (1, IOSTAT=IO)
    read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
    do while(.not.(field(1).eq."&COORD"))
       read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
    end do
    do iat=1,mol%n_atoms
       read(1,'(A)',iostat=io) line
       call line_parser(line,nfield,field)
       mol%atoms(iat)%elt=field(1)
       mol%atoms(iat)%Zato=elt2Z(mol%atoms(iat)%elt)
       read(field(2),*) mol%atoms(iat)%q(1)
       read(field(3),*) mol%atoms(iat)%q(2)
       read(field(4),*) mol%atoms(iat)%q(3)
       print *,(mol%atoms(iat)%q(k),k=1,3)
       mol%atoms(iat)%idx_molecule=-1
       mol%atoms(iat)%idx=iat
       mol%atoms(iat)%n_bonds=0
    end do
    
    close(1)
    
    stop
  end function IO_read_molecule_from_cp2k_restart_file
  ! --------------------------------------------------------------------------------------
  !
  !              IO_read_molecule_from_xyz_file()
  !
  ! --------------------------------------------------------------------------------------
  function IO_read_molecule_from_xyz_file(input) result(mol)
    use global
    use Element
    use Molecule
    implicit none
    type(t_File)::input
    type(t_Molecule)::mol
    
    character (len=1024)::line
    character (len=NCHARFIELD)::field(32)
    integer::nfield
    integer::iconf,iat,io
    
    mol=MOLECULE_init()
    open(unit=1,file=input%name,form='formatted')
    iconf=1
    do while(.not.(iconf.eq.input%iconf))
       read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
       read(field(1),*) mol%n_atoms
       read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field)
       do iat=1,mol%n_atoms
          read(1,'(A)',iostat=io) line
       end do
       iconf=iconf+1
    end do
    
    read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
    read(field(1),*) mol%n_atoms
    allocate(mol%atoms(mol%n_atoms))
    read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field)
    do iat=1,mol%n_atoms
       read(1,'(A)',iostat=io) line
       call line_parser(line,nfield,field)
       mol%atoms(iat)%elt=field(1)
       mol%atoms(iat)%Zato=elt2Z(mol%atoms(iat)%elt)
       read(field(2),*) mol%atoms(iat)%q(1)
       read(field(3),*) mol%atoms(iat)%q(2)
       read(field(4),*) mol%atoms(iat)%q(3)
       mol%atoms(iat)%idx_molecule=-1
       mol%atoms(iat)%idx=iat
       mol%atoms(iat)%n_bonds=0
    end do
    close(1)
  end function IO_read_molecule_from_xyz_file
  ! --------------------------------------------------------------------------------------
  !
  !        QuantumEspresso()
  !
  ! --------------------------------------------------------------------------------------
  subroutine QuantumEspresso(output,univ)
    implicit none
    type(t_File)::output
    type(t_Univers)::univ
    integer::i,j

    msg(__LINE__), " Saving Quantum Espresso input file",trim(output%name)

    print *,"# writing Quantum Espress input file"
    open(unit=1,file=output%name,form='formatted')
    write(1,*) "&control"
    write(1,*) "    calculation='scf',"
    write(1,*) "    pseudo_dir='/home2020/home/ipcms/bulou/workdir/pseudo'"
    write(1,*) "    wf_collect=.true."
    write(1,*) " /"
    write(1,*) " &system"
    write(1,*) "   ibrav=4,"
    write(1,*) "   celldm(1)=",univ%configurations(output%iconf)%cells%L(1,1)/a0
    write(1,*) "   celldm(3)=",univ%configurations(output%iconf)%cells%L(3,3)/univ%configurations(output%iconf)%cells%L(1,1)
    !print *,iconf
    !write(*,*) "   celldm(3)=",univ%configurations(iconf)%cells%L(3,3)/univ%configurations(iconf)%cells%L(1,1)
    write(1,*) "   nat=",univ%configurations(output%iconf)%cells%molecules(1)%n_atoms
    !write(*,*) "   nat=",univ%configurations(iconf)%cells%n_atoms
    write(1,*) "   ntyp=1,"
    write(1,*) "   ecutwfc=50.0,"
    write(1,*) "   ecutrho=400.0"
    write(1,*) "   occupations='smearing'"
    write(1,*) "   degauss=0.1"
    write(1,*) "/"
    write(1,*) " &electrons"
    write(1,*) "   mixing_beta=0.3,"
    write(1,*) " /"
    write(1,*) "ATOMIC_SPECIES"
    write(1,*) "Ti 47.867 Ti.pz-spn-kjpaw_psl.1.0.0.UPF"
    !write(1,*) "O 15.9994 O.pz-n-kjpaw_psl.0.1.UPF"
    !write(1,*) "H 1.00794 H.pz-kjpaw_psl.0.1.UPF"
    write(1,*) ""
    write(1,*) "ATOMIC_POSITIONS angstrom"
    do i=1,univ%configurations(output%iconf)%cells%molecules(1)%n_atoms
       write(1,*) univ%configurations(output%iconf)%cells%molecules(1)%atoms(i)%elt,&
            (univ%configurations(output%iconf)%cells%molecules(1)%atoms(i)%q(j),j=1,3)
    end do
    write(1,*) ""
    write(1,*) "K_POINTS automatic"
    write(1,*) "1 1 1 0 0 0"
    close(1)
  end subroutine QuantumEspresso

  ! --------------------------------------------------------------------------------------
  !
  !              read_multixyz()
  !
  ! --------------------------------------------------------------------------------------
  subroutine read_multixyz(name,univers)
    use Element
    implicit none
    !type(t_Param)::file
    character(len=2048)::name
    type(t_Univers),target::univers
    character (len=1024)::line
    character (len=NCHARFIELD)::field(32)
    integer::nfield
    integer::iconf,iat,io
    type(t_Molecule),pointer::mol
    
    print *,"# Starting read_multixyz() subroutine ..."
    print *,"#     * Reading ",trim(name)

    open(unit=1,file=name,form='formatted')
    do iconf=1,univers%n_configurations
       univers%configurations(iconf)%cells%n_molecules=1

       mol=>univers%configurations(iconf)%cells%molecules(1)
       read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
       read(field(1),*) mol%n_atoms
       univers%configurations(iconf)%cells%n_atoms=mol%n_atoms
       read(1,'(A)',iostat=io) line ; call line_parser(line,nfield,field)
       do iat=1,mol%n_atoms
          read(1,'(A)',iostat=io) line
          call line_parser(line,nfield,field)
          mol%atoms(iat)%elt=field(1)
          mol%atoms(iat)%Zato=elt2Z(mol%atoms(iat)%elt)
          read(field(2),*) mol%atoms(iat)%q(1)
          read(field(3),*) mol%atoms(iat)%q(2)
          read(field(4),*) mol%atoms(iat)%q(3)
          mol%atoms(iat)%idx_molecule=-1
          mol%atoms(iat)%n_bonds=0
       end do
    end do

    close(1)
    print *,"# ... end of read_multixyz()."
  end subroutine read_multixyz
  ! --------------------------------------------------------------------------------------
  !
  !              save()
  !
  ! --------------------------------------------------------------------------------------
  subroutine save(filename,univers,idx_conf_to_save)
    implicit none
    type(t_Univers)::univers
    
    character(len=1024)::filename
    integer::idx
    integer,optional::idx_conf_to_save
    idx=scan(filename,".")
    select case(trim(filename(idx+1:)))
    case ('xsf') 
       print *,"# xsf file"
       print *,"# Saving configuration ",idx_conf_to_save," in ",trim(filename)
       call xsf(filename,univers%configurations(idx_conf_to_save)%cells)
    case ('xyz') 
       print *,"# xyz file"
       print *,"# Saving configuration ",idx_conf_to_save," in ",trim(filename)
       call IO_write_cell_to_xyz(filename,univers%configurations(idx_conf_to_save)%cells)
    !case ('lmp')
    !   !       call lammps(filename,univers%configurations(idx_conf_to_save)%cells)
    !   call lammps(univers%configurations(idx_conf_to_save)%cells)
    case default
       print *,"# Unkown file format"
       stop
    end select
  end subroutine save
  ! --------------------------------------------------------------------------------------
  !
  !              save_molecule()
  !
  ! --------------------------------------------------------------------------------------
  subroutine save_molecule(mol)
    use Element
    use Molecule
    implicit none
    type(t_Molecule)::mol
    integer::k,jat
    
    integer,allocatable::tab(:)

    !tab=MOLECULE_sort_z(mol)
    !do jat=1,mol%n_atoms
    !   write(*,*) tab(jat),mol%atoms(tab(jat))%q(3)
    !end do
    msg(__LINE__), " Saving into molecule.xyz"
    open(unit=1,file= "molecule.xyz",form='formatted',status='unknown')
    write(1,*) mol%n_atoms
    write(1,*)
    do jat=1,mol%n_atoms
       write(1,*) PeriodicTable(mol%atoms(jat)%Zato)%elt,&
            (mol%atoms(jat)%q(k),k=1,3)
       !,&
       !     (ABS(AINT(mol%atoms(jat)%q(k))),k=1,3)
    end do
    close(1)

    !debug(__LINE__,"1")
  end subroutine save_molecule

  ! --------------------------------------------------------------------------------------
  !
  !              xsf()
  !
  ! --------------------------------------------------------------------------------------
  subroutine xsf(filename,cell)
    implicit none
    type(t_Cell)::cell
    character(len=1024)::filename
    integer::k,imol,jat
    ! ----------------------------- XSF --------------------------------------------------
    ! see http://www.xcrysden.org/doc/XSF.html
    open(unit=1,file=trim(filename),form='formatted',status='unknown')
    write(1,*) "CRYSTAL"
    write(1,*) "# these are primitive lattice vectors (in Angstroms)"
    write(1,*) " PRIMVEC"
    write(1,*) "      ",(cell%L(1,k),k=1,3)
    write(1,*) "      ",(cell%L(2,k),k=1,3)
    write(1,*) "      ",(cell%L(3,k),k=1,3)
    write(1,*) "# these are convetional lattice vectors (in Angstroms)"
    write(1,*) " CONVEC"
    write(1,*) "      ",(cell%L(1,k),k=1,3)
    write(1,*) "      ",(cell%L(2,k),k=1,3)
    write(1,*) "      ",(cell%L(3,k),k=1,3)
    write(1,*) "# these are atomic coordinates in a primitive unit cell "
    write(1,*) " # (in Angstroms)"
    write(1,*) "PRIMCOORD ",1
    write(1,*) cell%n_atoms," 1"
    do imol=1,cell%n_molecules
       do jat=1,cell%molecules(imol)%n_atoms
          write(1,*) cell%molecules(imol)%atoms(jat)%Zato,&
               (cell%molecules(imol)%atoms(jat)%q(k),k=1,3)
       end do
    end do
    close(1)
  end subroutine xsf
  ! --------------------------------------------------------------------------------------
  !
  !              IO_write_cell_to_xyz()
  !
  ! --------------------------------------------------------------------------------------
  subroutine IO_write_cell_to_xyz(filename,cell)
    use Element
    implicit none
    type(t_Cell)::cell
    character(len=1024)::filename
    integer::k,imol,jat,n_atoms
    open(unit=1,file=trim(filename),form='formatted',status='unknown')
    n_atoms=0
    do imol=1,cell%n_molecules
       if(cell%molecules(imol)%save) then
          do  jat=1,cell%molecules(imol)%n_atoms
             if(cell%molecules(imol)%atoms(jat)%save) then
                n_atoms=n_atoms+1
             end if
          end do
       end if
    end do
    write(1,*) n_atoms
    write(1,*)
    do imol=1,cell%n_molecules
       if(cell%molecules(imol)%save) then
          msg(__LINE__) , "Saving molecule ",imol
          do jat=1,cell%molecules(imol)%n_atoms
             if(cell%molecules(imol)%atoms(jat)%save) then
                write(1,*) PeriodicTable(cell%molecules(imol)%atoms(jat)%Zato)%elt,&
                     (cell%molecules(imol)%atoms(jat)%q(k),k=1,3) !,&
                     !cell%molecules(imol)%atoms(jat)%idx
             end if
          end do
       end if
    end do
    close(1)

  end subroutine IO_write_cell_to_xyz

end module IO
