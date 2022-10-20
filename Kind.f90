module Kind
  
  implicit none
contains
#define msg(x) print *,"### ",x,__FILE__," > "
#define debug(x,idx) print *,"### DEBUG ",idx," @ line ",x,__FILE__," > "
  
  function KIND_copy(src) result(dest)
    use global
    use Potential
    implicit none
    type(t_Kind)::src,dest
    dest%element=src%element
    dest%basis_set=src%basis_set
    dest%potential=POTENTIAL_copy(src%potential)
  end function KIND_copy

  function KIND_read_potential_from_cp2k_restart_file(unit,field) result(kinds)
    use global
    implicit none
    integer::unit,io,i,j,k,l
    character (len=1024)::line
    character (len=NCHARFIELD)::field(32)
    integer::nfield
    type(t_Kind)::kinds
    
    
    kinds%element=field(2)
    do while(.not.&
         ((trim(to_upper(field(1))).eq."&END").and.&
         (trim(to_upper(field(2))).eq."KIND")))



       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       read(unit,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
       if(trim(to_upper(field(1))).eq."BASIS_SET")     then
          kinds%basis_set=field(2)
          msg(__LINE__),"FORCE EVAL > SUBSYS > KIND >  BASIS SET > ",&
               kinds%basis_set
       end if
       if(trim(to_upper(field(1))).eq."ELEMENT")     then
          kinds%element=field(2)
          msg(__LINE__),"FORCE EVAL > SUBSYS > KIND > ELEMENT > ",&
               kinds%element
       end if
       if(trim(to_upper(field(1))).eq."POTENTIAL")     then
          kinds%potential%name=field(2)
          msg(__LINE__),"FORCE EVAL > SUBSYS > KIND > POTENTIAL > ",&
               kinds%potential%name
       end if
       if(trim(to_upper(field(1))).eq."&POTENTIAL")     then
          do while(.not.&
               ((trim(to_upper(field(1))).eq."&END").and.&
               (trim(to_upper(field(2))).eq."POTENTIAL")))
             read(unit,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
             kinds%potential%n_type_elec=nfield
             allocate(kinds%potential%n_elec(nfield))
             do i=1,nfield
                read(field(i),*) kinds%potential%n_elec(i)
             end do
             msg(__LINE__),"FORCE EVAL > SUBSYS > KIND > POTENTIAL > N ELEC",&
                  kinds%potential%n_elec
             read(unit,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
             read(field(1),*) kinds%potential%r_loc
             read(field(2),*) kinds%potential%nexp_ppl
             if(kinds%potential%nexp_ppl.gt.0) then
                allocate(kinds%potential%cexp_ppl(&
                     kinds%potential%nexp_ppl))
                do i=1,kinds%potential%nexp_ppl
                   read(field(2+i),*) kinds%potential%cexp_ppl(i)
                end do
             end if
             msg(__LINE__),"FORCE EVAL > SUBSYS > KIND > POTENTIAL > r_loc > ",&
                  kinds%potential%r_loc
             msg(__LINE__),"FORCE EVAL > SUBSYS > KIND > POTENTIAL > exp_ppl > ",&
                  kinds%potential%cexp_ppl
             read(unit,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
             read(field(1),*) kinds%potential%nprj
             msg(__LINE__),"FORCE EVAL > SUBSYS > KIND > POTENTIAL > nprj > ",&
                  kinds%potential%nprj


             allocate(kinds%potential%nlprj(&
                  kinds%potential%nprj))
             do i=1,kinds%potential%nprj
                read(unit,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
                read(field(1),*) kinds%potential%nlprj(i)%r
                read(field(2),*) kinds%potential%nlprj(i)%nprj_ppnl
                msg(__LINE__),"FORCE EVAL > SUBSYS > KIND > POTENTIAL > nlprj(",i,") > r >",&
                     kinds%potential%nlprj(i)%r
                msg(__LINE__),"FORCE EVAL > SUBSYS > KIND > POTENTIAL > nlprj(",i&
                     ,") > nprj_ppnl >",&
                     kinds%potential%nlprj(i)%nprj_ppnl

                if(kinds%potential%nlprj(i)%nprj_ppnl.gt.0) then
                   k=kinds%potential%nlprj(i)%nprj_ppnl*(&
                        kinds%potential%nlprj(i)%nprj_ppnl+1)/2
                   allocate(kinds%potential%nlprj(i)%hprj_ppnl(k))
                   k=1
                   do j=1,kinds%potential%nlprj(i)%nprj_ppnl
                      read(field(2+j),*) &
                           kinds%potential%nlprj(i)%hprj_ppnl(k)
                      k=k+1
                   end do
                   if(kinds%potential%nlprj(i)%nprj_ppnl.gt.1) then
                      do l=1,kinds%potential%nlprj(i)%nprj_ppnl-1
                         read(unit,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
                         do j=1,kinds%potential%nlprj(i)%nprj_ppnl-l
                            read(field(j),*) &
                                 kinds%potential%nlprj(i)%hprj_ppnl(k)
                            k=k+1
                         end do
                      end do
                   end if
                   msg(__LINE__),"FORCE EVAL > SUBSYS > KIND > POTENTIAL > nlprj(",i&
                        ,") > nprj_ppnl >",&
                        kinds%potential%nlprj(i)%hprj_ppnl

                end if
             end do ! do i=1,nprj

             do while(.not.&
                  ((trim(to_upper(field(1))).eq."&END").and.&
                  (trim(to_upper(field(2))).eq."POTENTIAL")))

                read(unit,'(A)',iostat=io) line ; call line_parser(line,nfield,field);
                msg(__LINE__),trim(line)                                                                   
             end do



          end do  ! do while (&end potential)

       end if

       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end do
    msg(__LINE__), "End of KIND_read_potential_from_cp2k_restart_file"
  end function KIND_read_potential_from_cp2k_restart_file
  
end module Kind
