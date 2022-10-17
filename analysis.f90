program analysis
  implicit none
  !-----------------------------------------------
  type t_PeriodicTable
     integer::Z
     character(len=2)::elt
  end type t_PeriodicTable
  !-----------------------------------------------
  type t_Atom
     integer::idx
     double precision::q(3)
     character(len=2)::elt
     integer::Zato
  end type t_Atom
  !-----------------------------------------------
  type t_Molecule
     integer:: idx
     integer:: n_atom
     type(t_Atom),allocatable::atoms(:)
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
     type(t_Molecule),allocatable::molecules(:)
  end type t_Configuration
  !-----------------------------------------------
  type t_Box
     integer:: idx
     integer:: n_configuration                     
     type(t_Configuration),allocatable::configurations(:)
     double precision::cell(3,3)
  end type t_Box


  !type(t_Molecule),allocatable::evol(:)
  type(t_Box),allocatable::CaCO3(:)
  integer::iconf,i
  character (len=1024)::path,filename
  character (len=2048)::trajfile,inputfile,nrjfile
  character(len=22)::fmt,conf

  type (t_Group)::group

  double precision::alpha
  integer::idx


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  character(len=1024) :: arg 
  
  do i = 1, iargc() 
     call getarg(i, arg) 
     write (6,*) "#argument ",i,trim(arg) 
  end do
  call getarg(1,path)
  print *,trim(path)

  !print *,"<---fin--->" ;  stop
  ! on alloue la boite
  allocate(CaCO3(1));  CaCO3(1)%idx=1




  ! 2 groupes d'atomes : 1-Ca et 2-CO3
  !group%nelt=3    ! trois espèces d'atomes à trier : Ca, C et O
  !allocate(group%elt_table(group%nelt)) ! 1-Ca 2-CO3 
  !group%elt_table(1)%idx=1 ; group%elt_table(1)%elt='Ca'
  !group%elt_table(2)%idx=2 ; group%elt_table(2)%elt='C'
  !group%elt_table(3)%idx=2 ; group%elt_table(3)%elt='O'


  group%nelt=2    ! trois espèces d'atomes à trier : Ca, C 
  allocate(group%elt_table(group%nelt)) ! 1-Ca 2-C
  group%elt_table(1)%idx=1 ; group%elt_table(1)%elt='Ca'
  group%elt_table(2)%idx=2 ; group%elt_table(2)%elt='C'




  group%n_set=2



  allocate(group%n_atoms_set(group%n_set))
  allocate(group%set(group%n_set))
  do i=1,group%n_set
     group%n_atoms_set(i)=0
  end do



  !path='/home/bulou/pc-hervecal/ownCloud/zim/Codes/CP2K/benchmarks/CaCO3/run_equilibre/'
  !path='/home/bulou/pc-hervecal/ownCloud/zim/Codes/CP2K/benchmarks/CaCO3/run1/'

  filename="run-CaCO3-AIMD.xyz"
  trajfile=path(1:len( trim(path) ))//filename

  print *,trim(trajfile)
  
  filename='aimd.in'
  inputfile=path(1:len( trim(path) ))//filename
  print *,trim(inputfile)

  !path='/home/bulou/pc-hervecal/ownCloud/zim/Codes/CP2K/benchmarks/CaCO3/run_equilibre/'

  !filename='run-CaCO3-AIMD-1.ener'
  !nrjfile=path(1:len( trim(path) ))//filename
  !write(*,*) trim(trajfile)
  !print *,trim(inputfile)

  ! on lit d'abord les infos sur la taille de la boite
  call read_cp2k_input_file(inputfile,CaCO3(1))
  ! puis on lit les coordonnées des atomes composant le système
  !call QE_cpmd()
  call cp2k_aimd(CaCO3(1),trajfile)


  idx=5  ! Ca
  call center_at(CaCO3(1),idx)
  !alpha=0.5;  call translate(CaCO3(1),alpha)

  ! puis on calcule la distance entre le Ca (idx=5) et le CO3 (idx=1,2,3,4)
  
  call system('rm ./R.dat')
  do iconf=1,CaCO3(1)%n_configuration
     do i=1,group%n_set
        group%n_atoms_set(i)=0
     end do
     call update_group(group,CaCO3(1),iconf)
  end do





  iconf=CaCO3(1)%n_configuration  
  !call transformation(CaCO3(1),iconf)




  !fmt='lammps trajectory file'
  !call save(CaCO3(1),fmt)
  fmt='axsf'
  !fmt='xsf'
  iconf=1
  call save(CaCO3(1),fmt,iconf)
  stop
  call write_input_file(CaCO3(1),iconf)

  !evol(nevol)%atoms(5)%q(0)=evol(nevol)%atoms(5)%q[0]+2
  !evol(nevol)%atoms(5)%q(1)=evol(nevol)%atoms(5)%q[1]+2
  !call save(evol,nevol)

contains
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
  !              elt2Z(elt,Zato)
  !
  ! --------------------------------------------------------------------------------------
  subroutine elt2Z(elt,Zato)
    implicit none
    integer::Zato,Nelt,i
    character(len=2)::elt
    type (t_PeriodicTable),allocatable::PeriodicTable(:)
    logical::found
    Nelt=20
    allocate(PeriodicTable(Nelt))
    PeriodicTable(1)%Z=1 ; PeriodicTable(1)%elt='H'
    PeriodicTable(2)%Z=2 ; PeriodicTable(2)%elt='He'
    PeriodicTable(3)%Z=3 ; PeriodicTable(3)%elt='Li'
    PeriodicTable(4)%Z=4 ; PeriodicTable(4)%elt='Be'
    PeriodicTable(5)%Z=5 ; PeriodicTable(5)%elt='B'
    PeriodicTable(6)%Z=6 ; PeriodicTable(6)%elt='C'
    PeriodicTable(7)%Z=7 ; PeriodicTable(7)%elt='N'
    PeriodicTable(8)%Z=8 ; PeriodicTable(8)%elt='O'
    PeriodicTable(9)%Z=9 ; PeriodicTable(9)%elt='F'
    PeriodicTable(10)%Z=10 ; PeriodicTable(10)%elt='Ne'
    PeriodicTable(11)%Z=11 ; PeriodicTable(11)%elt='Na'
    PeriodicTable(12)%Z=12 ; PeriodicTable(12)%elt='Mg'
    PeriodicTable(13)%Z=13 ; PeriodicTable(13)%elt='Al'
    PeriodicTable(14)%Z=14 ; PeriodicTable(14)%elt='Si'
    PeriodicTable(15)%Z=15 ; PeriodicTable(15)%elt='P'
    PeriodicTable(16)%Z=16 ; PeriodicTable(16)%elt='S'
    PeriodicTable(17)%Z=17 ; PeriodicTable(17)%elt='Cl'
    PeriodicTable(18)%Z=18 ; PeriodicTable(18)%elt='Ar'
    PeriodicTable(19)%Z=19 ; PeriodicTable(19)%elt='K'
    PeriodicTable(20)%Z=20 ; PeriodicTable(20)%elt='Ca'

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
       natom=box%configurations(iconf)%molecules(1)%n_atom
       do imol=2,box%configurations(iconf)%n_molecule
          natom=natom+box%configurations(iconf)%molecules(imol)%n_atom
       end do
       write(1,*) natom," 1"
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
  subroutine write_input_file(box,iconf)
    implicit none
    type(t_Box)::box
    integer::imol,iconf,jat,k
    integer::nsteps
    double precision::temperature
    temperature=600.0
    nsteps=1000
    ! cp2k - aimd
    open(unit=1,file="TMP-aimd.in",form='formatted',status='unknown')
    write(1,*) '@SET SYSTEM run-CaCO3-AIMD'
    write(1,*) '@set LIBDIR /home/bulou/src/cp2k/data'
    write(1,*) '@set TEMPERATURE',temperature
    write(1,*) '@set NSTEPS',nsteps
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
    write(1,*) '      EPS_SCF 1.0E-5'
    write(1,*) '      MAX_SCF 500'
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
    write(1,*) '      PERIODIC NXYZ'
    write(1,*) '    &END CELL'
    write(1,*) ''
    write(1,*) '    &COORD'
    do imol=1,box%configurations(iconf)%n_molecule
       do jat=1,box%configurations(imol)%molecules(imol)%n_atom
          write(1,*) box%configurations(iconf)%molecules(imol)%atoms(jat)%elt,&
               (box%configurations(iconf)%molecules(imol)%atoms(jat)%q(k),k=1,3)
       end do
    end do
    write(1,*) '    &END COORD'
    write(1,*) ''
    write(1,*) '    &KIND C                  # basis sets and pseudo-potentials for atomic species'
    write(1,*) '      BASIS_SET DZVP-MOLOPT-GTH'
    write(1,*) '      POTENTIAL GTH-PBE-q4'
    write(1,*) '    &END KIND'
    write(1,*) '    &KIND O'
    write(1,*) '      BASIS_SET DZVP-MOLOPT-GTH'
    write(1,*) '      POTENTIAL GTH-PBE-q6'
    write(1,*) '    &END KIND'
    write(1,*) '    &KIND Ca'
    write(1,*) '      BASIS_SET DZVP-MOLOPT-SR-GTH'
    write(1,*) '      POTENTIAL GTH-PBE-q10'
    write(1,*) '    &END KIND'
    write(1,*) '  &END SUBSYS'
    write(1,*) ''
    write(1,*) '&END FORCE_EVAL'
    write(1,*) ''
    write(1,*) '&MOTION'
    write(1,*) '  &MD'
    write(1,*) '    ENSEMBLE         NVE'
    write(1,*) '    STEPS            ',nsteps
    write(1,*) '    TIMESTEP [fs]    0.5'
    write(1,*) '    TEMPERATURE [K]  ',temperature
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

    write(1,*) ' &EXT_RESTART'
    write(1,*) '   RESTART_FILE_NAME run-CaCO3-AIMD-1.restart'
    write(1,*) ' &END'



    close(1)
  end subroutine write_input_file
end program analysis
