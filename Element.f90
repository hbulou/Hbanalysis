module Element

contains
#define msg(x) print *,"### ",x,' ### ',__FILE__," > "
#define error(x) print *,"### ERROR ",x,__FILE__," > "
  ! --------------------------------------------------------------------------------------
  !
  !              init()
  !
  ! --------------------------------------------------------------------------------------
  function ELEMENT_copy(src) result(dest)
    use global
    implicit none
    type(t_Element)::dest,src
    dest%Zato=src%Zato
    dest%elt=src%elt
    dest%n=src%n
  end function ELEMENT_copy
  ! --------------------------------------------------------------------------------------
  !
  !              init()
  !
  ! --------------------------------------------------------------------------------------
  subroutine ELEMENT_init()
    use global
    implicit none
    integer::i,j
    PeriodicTable(1)%Z=1 ; PeriodicTable(1)%elt='H' ; PeriodicTable(1)%mass=1.0 ; 
    PeriodicTable(2)%Z=2 ; PeriodicTable(2)%elt='He' ; PeriodicTable(2)%mass=4.0
    PeriodicTable(3)%Z=3 ; PeriodicTable(3)%elt='Li' ; PeriodicTable(3)%mass=6.9
    PeriodicTable(4)%Z=4 ; PeriodicTable(4)%elt='Be' ; PeriodicTable(4)%mass=9.0
    PeriodicTable(5)%Z=5 ; PeriodicTable(5)%elt='B' ; PeriodicTable(5)%mass=10.8
    PeriodicTable(6)%Z=6 ; PeriodicTable(6)%elt='C' ; PeriodicTable(6)%mass=12.0
    PeriodicTable(7)%Z=7 ; PeriodicTable(7)%elt='N' ; PeriodicTable(7)%mass=14.0
    PeriodicTable(8)%Z=8 ; PeriodicTable(8)%elt='O' ; PeriodicTable(8)%mass=16.0
    PeriodicTable(9)%Z=9 ; PeriodicTable(9)%elt='F' ; PeriodicTable(9)%mass=19.0
    PeriodicTable(10)%Z=10 ; PeriodicTable(10)%elt='Ne' ; PeriodicTable(10)%mass=20.2
    PeriodicTable(11)%Z=11 ; PeriodicTable(11)%elt='Na' ; PeriodicTable(11)%mass=23.0
    PeriodicTable(12)%Z=12 ; PeriodicTable(12)%elt='Mg' ; PeriodicTable(12)%mass=24.3
    PeriodicTable(13)%Z=13 ; PeriodicTable(13)%elt='Al' ; PeriodicTable(13)%mass=27.0
    PeriodicTable(14)%Z=14 ; PeriodicTable(14)%elt='Si' ; PeriodicTable(14)%mass=28.1
    PeriodicTable(15)%Z=15 ; PeriodicTable(15)%elt='P' ; PeriodicTable(15)%mass=31.0
    PeriodicTable(16)%Z=16 ; PeriodicTable(16)%elt='S' ; PeriodicTable(16)%mass=32.1
    PeriodicTable(17)%Z=17 ; PeriodicTable(17)%elt='Cl' ; PeriodicTable(17)%mass=35.5
    PeriodicTable(18)%Z=18 ; PeriodicTable(18)%elt='Ar' ; PeriodicTable(18)%mass=39.9
    PeriodicTable(19)%Z=19 ; PeriodicTable(19)%elt='K' ; PeriodicTable(19)%mass=39.1
    PeriodicTable(20)%Z=20 ; PeriodicTable(20)%elt='Ca' ; PeriodicTable(20)%mass=40.1
    PeriodicTable(21)%Z=21 ; PeriodicTable(21)%elt='Sc' ; PeriodicTable(21)%mass=44.955912
    PeriodicTable(22)%Z=22 ; PeriodicTable(22)%elt='Ti' ; PeriodicTable(22)%mass=47.867
    PeriodicTable(23)%Z=23 ; PeriodicTable(23)%elt='V' ;  PeriodicTable(23)%mass=50.942
    PeriodicTable(24)%Z=24 ; PeriodicTable(24)%elt='Cr' ; PeriodicTable(24)%mass=51.996
    PeriodicTable(25)%Z=25 ; PeriodicTable(25)%elt='Mn' ; PeriodicTable(25)%mass=54.938
    PeriodicTable(26)%Z=26 ; PeriodicTable(26)%elt='Fe' ; PeriodicTable(26)%mass=55.845
    PeriodicTable(27)%Z=27 ; PeriodicTable(27)%elt='Co' ; PeriodicTable(27)%mass=58.933
    PeriodicTable(28)%Z=28 ; PeriodicTable(28)%elt='Ni' ; PeriodicTable(28)%mass=58.693
    PeriodicTable(29)%Z=29 ; PeriodicTable(29)%elt='Cu' ; PeriodicTable(29)%mass=63.546
    PeriodicTable(30)%Z=30 ; PeriodicTable(30)%elt='Zn' ; PeriodicTable(30)%mass=65.39
    PeriodicTable(31)%Z=31 ; PeriodicTable(31)%elt='Ga' ; PeriodicTable(31)%mass=69.723
    PeriodicTable(32)%Z=32 ; PeriodicTable(32)%elt='Ge' ; PeriodicTable(32)%mass=72.64
    PeriodicTable(33)%Z=33 ; PeriodicTable(33)%elt='As' ; PeriodicTable(33)%mass=74.922
    PeriodicTable(34)%Z=34 ; PeriodicTable(34)%elt='Se' ; PeriodicTable(34)%mass=78.96
    PeriodicTable(35)%Z=35 ; PeriodicTable(35)%elt='Br' ; PeriodicTable(35)%mass=79.904

    PeriodicTable(36)%Z=36 ; PeriodicTable(36)%elt='Kr' ; PeriodicTable(36)%mass=83.798
    PeriodicTable(37)%Z=37 ; PeriodicTable(37)%elt='Rb' ; PeriodicTable(37)%mass=85.4678
    PeriodicTable(38)%Z=38 ; PeriodicTable(38)%elt='Sr' ; PeriodicTable(38)%mass=87.62
    PeriodicTable(39)%Z=39 ; PeriodicTable(39)%elt='Y' ; PeriodicTable(39)%mass= 88.90584
    PeriodicTable(40)%Z=40 ; PeriodicTable(40)%elt='Zr' ; PeriodicTable(40)%mass=91.224
    PeriodicTable(41)%Z=41 ; PeriodicTable(41)%elt='Nb' ; PeriodicTable(41)%mass=92.90637
    PeriodicTable(42)%Z=42 ; PeriodicTable(42)%elt='Mo' ; PeriodicTable(42)%mass=95.95
    PeriodicTable(43)%Z=43 ; PeriodicTable(43)%elt='Tc' ; PeriodicTable(43)%mass=98
    PeriodicTable(44)%Z=44 ; PeriodicTable(44)%elt='Ru' ; PeriodicTable(44)%mass=101.07
    PeriodicTable(45)%Z=45 ; PeriodicTable(45)%elt='Rh' ; PeriodicTable(45)%mass=102.90550
    PeriodicTable(46)%Z=46 ; PeriodicTable(46)%elt='Pd' ; PeriodicTable(46)%mass=106.42

    
    do j=1,NELT
       PeriodicTable(j)%Basis_set="NONE"
       PeriodicTable(j)%Potential_set="NONE"
       do i=1,NELT
          ELEMENT_deq(i,j)=0.0
       end do
    end do
    
    ELEMENT_deq(1,8)=1.2 ; ELEMENT_deq(8,1)=ELEMENT_deq(1,8)  ! O-H
    ELEMENT_deq(7,7)=1.2 ;                    ! N-N
    ELEMENT_deq(8,8)=1.3 ;                    ! O-O
  end subroutine ELEMENT_init

  ! --------------------------------------------------------------------------------------
  !
  !              elt2Z(elt,Zato)
  !
  ! --------------------------------------------------------------------------------------
  function elt2Z(elt) result(Zato)
    use global
    implicit none
    integer::Zato,i
    character(len=2)::elt
    logical::found

    found=.False.
    i=0
    do while((.not.found).and.(i<Nelt))
       i=i+1
       if(elt.eq.PeriodicTable(i)%elt) then
          found=.True.
       end if
    end do
    if(found) then
       Zato=PeriodicTable(i)%Z
       !msg(__LINE__), trim(elt)," found"
       !write(*,*) elt,Zato
    else
       Zato=-1
       msg(__LINE__), elt," not found"
    end if

  end function elt2Z
end module Element
