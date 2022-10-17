module Univers
  implicit none
contains
  ! subroutine UNIVERS_add_configuration(univers)
  ! subroutine UNIVERS_copy(src,dest)
  ! function UNIVERS_new() result(dest)
  ! --------------------------------------------------------------------------------------
  !
  !              UNIVERS_add_configuration(univers)
  !
  ! --------------------------------------------------------------------------------------
  subroutine UNIVERS_add_configuration(univers)
    use global
    use Configuration
    implicit none
    type(t_Configuration)::conf
    type(t_Univers)::univers,TMPunivers
    print *,"# CONFIGURATION_add_to_univers> n_configuration=",univers%n_configurations
    print *,"# CONFIGURATION_add_to_univers> n_configuration=",allocated(univers%configurations)
    if(.not.(allocated(univers%configurations))) then
    !if(univers%n_configurations.eq.0) then
       univers%n_configurations=1
       !deallocate(univers%configurations)
       allocate(univers%configurations(univers%n_configurations))
       univers%configurations(1)=CONFIGURATION_init()
       
       !univers%configurations(univers%n_configurations)%cells%n_molecules=0
    else
       call UNIVERS_copy(univers,TMPunivers)
       deallocate(univers%configurations)
       univers%n_configurations=TMPunivers%n_configurations+1
       allocate(univers%configurations(univers%n_configurations))
       call UNIVERS_copy(TMPunivers,univers)
       deallocate(TMPunivers%configurations)
    end if
    print *,"# CONFIGURATION_add_to_univers> n_configuration=",univers%n_configurations
    print *,"# CONFIGURATION_add_to_univers> n_configuration=",allocated(univers%configurations)
    print *,"-------------------------------------------------------"
  end subroutine UNIVERS_add_configuration
  ! --------------------------------------------------------------------------------------
  !
  !              UNIVERS_copy()
  !
  ! --------------------------------------------------------------------------------------
  subroutine UNIVERS_copy(src,dest)
    use global
    use Configuration
    implicit none
    type(t_Univers)::src,dest
    integer::iconf

    if(.not.(allocated(dest%configurations))) then
       print *,"# UNIVERS_copy> not allocated"
       dest%n_configurations=src%n_configurations
       allocate(dest%configurations(dest%n_configurations))
    end if
    do iconf=1,src%n_configurations
       dest%configurations(iconf)=CONFIGURATION_copy(src%configurations(iconf))
    end do
  end subroutine UNIVERS_copy
  ! --------------------------------------------------------------------------------------
  !
  !              UNIVERS_new()
  !
  ! --------------------------------------------------------------------------------------
  function UNIVERS_new() result(dest)
    use global
    use Configuration
    implicit none
    type(t_Univers)::dest
    integer::iconf

    dest%n_configurations=0
    !allocate(dest%configurations(0))
  end function UNIVERS_new

end module Univers
