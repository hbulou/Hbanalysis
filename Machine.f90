module Machine
contains
#define msg(x) print *,"### ",x,' ###    ',__FILE__," > "
#define error(x) print *,"### ERROR ",x,__FILE__," > "
#define debug(x,idx) print *,"### DEBUG ",idx," @ line ",x,__FILE__," > "

  ! --------------------------------------------------------------------------------------
  !
  !        MACHINE_init()
  !
  ! --------------------------------------------------------------------------------------
  subroutine MACHINE_init()
    use global
    implicit none
    
    List_Machines(1)%name='hpc'
    List_Machines(1)%cp2k_libdir='/home2020/home/ipcms/bulou/src/cp2k/data/'
    
    List_Machines(2)%name='JeanZay'
    List_Machines(2)%cp2k_libdir="/gpfslocalsup/spack_soft/cp2k/8.2/intel-19.1.3-2l5irg7v7rchgpfxbgo66hp2rizcvo5f/share/data"
    
  end subroutine MACHINE_init
  
  ! --------------------------------------------------------------------------------------
  !
  !        MACHINE_search_idx()
  !
  ! --------------------------------------------------------------------------------------
  function MACHINE_search_idx(param) result(idx)
    use global
    implicit none
    type(t_Param)::param
    integer::idx
    idx=1
    do while((idx.le.N_MACHINES).and.(.not.(List_Machines(idx)%name.eq.param%machine)))
       msg(__LINE__), "Machines:",List_Machines(idx)%name,param%machine
       idx=idx+1
    end do
    msg(__LINE__)," idx of ",param%machine," is ",idx
  end function MACHINE_search_idx
end module Machine
