#define msg(x) print *,"### ",x,' ###    ',__FILE__," > "
#define error(x) print *,"### ERROR ",x,__FILE__," > "
#define debug(x,idx) print *,"### DEBUG ",idx," @ line ",x,__FILE__," > "
program new
  use global
  use Machine
  use scripts
  use Element
  use IO
  implicit none
  type(t_Param)::param
  INTEGER :: i,idx_input
  CHARACTER(len=32) :: arg


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call GLOBAL_init_param(param)

  ! lecture des parametres de la ligne de commande
  i=0
  do 
     CALL get_command_argument(i, arg)
     IF (LEN_TRIM(arg) == 0) EXIT
     if(trim(arg).eq."-i") then
        CALL get_command_argument(i+1, param%filename)
        exit
     end if
    i = i+1
  END DO
  print *,"# Input file: ",trim(param%filename)


  ! lecture du fichier d'input

  call IO_read_parametrization_file(param)

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !               INITIALISATION
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! element database
  call ELEMENT_init()
  ! momputer database
  call MACHINE_init()

  !call Ti()
  !call TiOH()
  !
  ! les subroutines correspondant aux differents cas sont dans le fichier script.f90
  !
  select case(param%run_type)
  case('test')
     call test(param)
  case('concat')
     call concat(param)
  case('crystal')
     call crystal(param)
  case('interpolation','interpol')
     call interpol(param)
  case('TiOH_addOH')
     idx_input=1
     call TiOH_addOH(param,idx_input)
  case('solvent')
     idx_input=1
     call solvent(param,idx_input)
  case('strain')
     call strain(param)
  case default
     error(__LINE__), " Keywords ",param%run_type," doesn't exist"
     stop
  end select

  !  call TiOH_moveOH()
  msg(__LINE__), " ------------------------------------------------------------------"
  msg(__LINE__)
  msg(__LINE__), "                     DONE !"
  msg(__LINE__)
  msg(__LINE__), " ------------------------------------------------------------------"
contains
end program new
