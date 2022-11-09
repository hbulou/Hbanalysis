program generate_configuration
  use Element
  use global
  use IO
  use Molecule
  use Cell
  implicit none
  character(len=1024)::filename


  type(t_File)::file
  type(t_Univers),target::univers
  integer::i,idx_conf,n,k
  type(t_Molecule),pointer::mol(:)
  type(t_Molecule)::h2o
  integer,pointer::n_atoms
  double precision::alpha,beta,gamma
  character(len=32)::adding_type  ! random || as_it
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !               INITIALISATION
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call init_element()


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !               LECTURE 
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !call system("sshfs hpc-login.u-strasbg.fr: /home/bulou/hpc-login")


  file%name="/home/bulou/hpc-login/workdir/H2O/32h2o.xyz"
  file%type='multixyz'
  call IO_file_info(file)

  idx_conf=1
  call IO_read_molecule(file,h2o,idx_conf)
  h2o%MassCenter=Mass_Center(h2o)

  

  !
  ! initialisation de l'univers
  !
  univers%n_configurations=1
  allocate(univers%configurations(univers%n_configurations))
  do i=1,univers%n_configurations
     univers%configurations(i)%cells%n_molecules=0
     univers%configurations(i)%cells%n_atoms=0
     univers%configurations(i)%cells%n_bonds=0
     allocate(univers%configurations(i)%cells%molecules(1))
     !allocate(univers%configurations(i)%cells%molecules(1)%atoms(file%n_atoms))
  end do
  file%name="/home/bulou/hpc-login/workdir/H2O/gopt-32h2o-ff.in"
  file%type='cp2k_input'
  call IO_read_file(file,univers)



  
  n=1
  adding_type='as_it'
  call CELL_add_molecules(univers%configurations(1)%cells,h2o,n,adding_type)

  !filename="essai.xyz" ;   call save(filename,univers,idx_conf)

  do i=1,univers%configurations(1)%cells%molecules(1)%n_atoms
     print *,&
          (univers%configurations(1)%cells%molecules(1)%atoms(i)%q(k),k=1,3),&
          (modulo(univers%configurations(1)%cells%molecules(1)%atoms(i)%q(k),univers%configurations(1)%cells%L(k,k)),k=1,3)
     do k=1,3
        univers%configurations(1)%cells%molecules(1)%atoms(i)%q(k)=&
             modulo(univers%configurations(1)%cells%molecules(1)%atoms(i)%q(k),univers%configurations(1)%cells%L(k,k))
             
     end do
  end do


  call MOLECULE_find_molecules(univers%configurations(univers%n_configurations)%cells,PBC=.True.)


  
  filename="essai.xsf" ;   call save(filename,univers,idx_conf)
  filename="essai.xyz" ;   call save(filename,univers,idx_conf)
  filename="essai.lmp" ;   call save(filename,univers,idx_conf)
  



  stop
  !file%name='run-30h2o-gopt-pos-1.xyz'
  file%name="/home/bulou/hpc-login/workdir/H2O/run-7h2o-gopt-pos-1.xyz"
  file%type='multixyz'
  call IO_file_info(file)
  
  univers%n_configurations=file%n_configurations
  allocate(univers%configurations(univers%n_configurations))
  do i=1,univers%n_configurations
     allocate(univers%configurations(i)%cells%molecules(1))
     allocate(univers%configurations(i)%cells%molecules(1)%atoms(file%n_atoms))
  end do
  call IO_read_file(file,univers)

  !print *,univers%configurations(univers%n_configurations)%cells%molecules(1)%atoms(1)%Zato

  file%name="/home/bulou/hpc-login/workdir/H2O/gopt-6h2o.in"
  file%type='cp2k_input'
  call IO_read_file(file,univers)

  !print *,"Zato=",univers%configurations(idx_conf)%cells%molecules(1)%atoms(1)%Zato
  print *,univers%configurations(univers%n_configurations)%cells%molecules(1)%n_atoms
  n_atoms=>univers%configurations(univers%n_configurations)%cells%molecules(1)%n_atoms
  !mol=>univers%configurations(univers%n_configurations)%cells%molecules
  !call find_molecules(mol,univers%configurations(univers%n_configurations)%cells%L)
  call MOLECULE_find_molecules(univers%configurations(univers%n_configurations)%cells,PBC=.True.)
!  print *,univers%configurations(univers%n_configurations)%cells%molecules(1)%n_atoms
!  print *,univers%configurations(univers%n_configurations)%cells%n_molecules
!  print *,univers%configurations(univers%n_configurations)%cells%n_atoms

  
  file%name="H2O.xyz"
  file%type='xyz'
  idx_conf=1

  call IO_read_molecule(file,h2o,idx_conf)
  h2o%MassCenter=Mass_Center(h2o)
!  file%type="molecule_MC"
!  alpha=0.0
!  beta=0.0
!  gamma=90.0
!  call MOLECULE_rotate(h2o,alpha,beta,gamma,file%type)
  !  call save_molecule(h2o)
  n=1
  idx_conf=univers%n_configurations
  print *,"Zato=",univers%configurations(idx_conf)%cells%molecules(1)%atoms(1)%Zato
  adding_type='random'
  call CELL_add_molecules(univers%configurations(idx_conf)%cells,h2o,n,adding_type)
  print *,"Zato=",univers%configurations(idx_conf)%cells%molecules(1)%atoms(1)%Zato

  print *,"#n_molecules=",univers%configurations(idx_conf)%cells%n_molecules
  print *,"#n_atoms=",univers%configurations(idx_conf)%cells%n_atoms


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !               SAUVEGARDES
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !filename="essai.xyz"
  !idx_conf=1
  !call save(filename,univers,idx_conf)
  filename="essai.xyz"
  call save(filename,univers,idx_conf)
  filename="essai.xsf" ;   call save(filename,univers,univers%n_configurations)
!  call system("fusermount -u /home/bulou/hpc-login")
contains     
end program generate_configuration
