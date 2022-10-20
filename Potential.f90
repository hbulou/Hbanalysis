module Potential

  implicit none
contains
#define msg(x) print *,"### ",x,__FILE__," > "
#define debug(x,idx) print *,"### DEBUG ",idx," @ line ",x,__FILE__," > "

  function POTENTIAL_copy(src) result(dest)
    use global
    implicit none
    integer::i
    type(t_GTH_Potential)::src,dest
    dest%name=src%name
    dest%n_type_elec=src%n_type_elec
    if(allocated(src%n_elec)) then
       allocate(dest%n_elec(dest%n_type_elec))
       dest%n_elec=src%n_elec
    end if


    dest%r_loc=src%r_loc
    dest%nexp_ppl=src%nexp_ppl
    if(allocated(src%cexp_ppl)) then
       allocate(dest%cexp_ppl(dest%nexp_ppl))
       dest%cexp_ppl=src%cexp_ppl
    end if
    dest%nprj=src%nprj
    if(allocated(src%nlprj)) then
       allocate(dest%nlprj(dest%nprj))
       do i=1,dest%nprj
          dest%nlprj(i)=POTENTIAL_copy_non_local_projector(src%nlprj(i))
       end do
    end if
  end function POTENTIAL_copy

  function POTENTIAL_copy_non_local_projector(src) result(dest)
    use global
    implicit none
    integer::k
    type(t_GTH_Potential_non_local_projector)::src,dest
    dest%r=src%r
    dest%nprj_ppnl=src%nprj_ppnl
    if(allocated(src%hprj_ppnl)) then
       k=dest%nprj_ppnl*(&
            dest%nprj_ppnl+1)/2
       allocate(dest%hprj_ppnl(k))
       dest%hprj_ppnl=src%hprj_ppnl
    end if
  end function POTENTIAL_copy_non_local_projector
  
end module Potential
