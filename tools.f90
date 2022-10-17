module tools
  implicit none
contains
#define msg(x) print *,"### ",x,__FILE__," > "
#define error(x) print *,"### ERROR ",x,__FILE__," > " 
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !                 function ran(idum)
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ran(idum)
    implicit none
    integer, parameter :: k4b=selected_int_kind(9)
    integer(k4b), intent(inout) :: idum
    real :: ran
    ! “minimal” random number generator of park and miller combined with a marsaglia shift
    ! sequence. returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
    ! values). this fully portable, scalar generator has the “traditional” (not fortran 90) calling
    ! sequence with a random deviate as the returned function value: call with idum a negative
    ! integer to initialize; thereafter, do not alter idum except to reinitialize. the period of this
    ! generator is about 3.1 × 1018.
    integer(k4b), parameter :: ia=16807,im=2147483647,iq=127773,ir=2836
    real, save :: am
    integer(k4b), save :: ix=-1,iy=-1,k
    if (idum <= 0 .or. iy < 0) then            ! initialize.
       am=nearest(1.0,-1.0)/im
       iy=ior(ieor(888889999,abs(idum)),1)
       ix=ieor(777755555,abs(idum))
       idum=abs(idum)+1                        ! set idum positive.
    end if
    ix=ieor(ix,ishft(ix,13))                   ! marsaglia shift sequence with period 232 − 1.
    ix=ieor(ix,ishft(ix,-17))
    ix=ieor(ix,ishft(ix,5))
    k=iy/iq                                   ! park-miller sequence by schrage’s method,
    iy=ia*(iy-k*iq)-ir*k                      ! period 231 − 2.
    if (iy < 0) iy=iy+im
    ran=am*ior(iand(im,ieor(ix,iy)),1)        ! combine the two generators with masking to
    !  ensure nonzero value.
  end function ran


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !                 function TOOLS_rotate_around_ez(center,ex,ey,ez,angle,q) result(Mp)
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function TOOLS_rotate_around_ez(center,ex,ey,ez,angle,q) result(Mp)
    implicit none
    double precision :: q(3),ex(3),ey(3),ez(3),center(3),M(3),Mp(3)
    double precision::x,y,z,angle,Mpx,Mpy,Mpz,Msqr,Mpsqr
    integer::k
    double precision::M_PI,anglerad
    M_PI=4*atan(1.d0)
    msg(__LINE__),"angle=",angle
    anglerad=(M_PI/180.0)*angle

    Msqr=0.0
    do k=1,3
       M(k)=q(k)-center(k)
       Msqr=Msqr+M(k)*M(k)
    end do

    x=0.0
    do k=1,3
       x=x+M(k)*ex(k)
    end do
    
    y=0.0
    do k=1,3
       y=y+M(k)*ey(k)
    end do
    
    z=0.0
    do k=1,3
       z=z+M(k)*ez(k)
    end do

    Mpx=x*cos(anglerad)-y*sin(anglerad)
    Mpy=y*cos(anglerad)+x*sin(anglerad)
    Mpz=z

    Mpsqr=0.0
    do k=1,3
       Mp(k)=Mpx*ex(k)+Mpy*ey(k)+Mpz*ez(k)
       Mpsqr=Mpsqr+Mp(k)*Mp(k)
    end do

    if(abs(Mpsqr-Msqr).gt.1.0e-12) then
       error(__LINE__), "Msqr=",Msqr,"Mpsqr=",Mpsqr,"diff=",Msqr-Mpsqr
    else
       msg(__LINE__), "Msqr=",Msqr,"Mpsqr=",Mpsqr,"diff=",Msqr-Mpsqr
    end if
    msg(__LINE__), "q=",q
    msg(__LINE__), "center=",center
    msg(__LINE__), "Mp=",Mp+center
    msg(__LINE__), "Mp-q",Mp-q+center
    msg(__LINE__), "x^2+y^2",x*x+y*y

    do k=1,3
       Mp(k)=Mp(k)+center(k)
    end do



  end function TOOLS_rotate_around_ez


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !     subroutine TOOLS_build_rotation_axis(center,axe,to_center,to_align,ex,ey,ez,alpha)
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine TOOLS_build_rotation_axis(center,axe,to_center,to_align,ex,ey,ez,alpha)
    ! INPUT
    !     double precision::center(3)
    !     double precision::axe(3)
    !     double precision::to_center(3)
    !     double precision::to_align(3)
    ! OUTPUT
    !     double precision::ex(3),ey(3),ez(3)
    !     double precision::alpha
    implicit none
    double precision::center(3),axe(3),to_center(3),to_align(3),alpha
    !double precision::alpha,beta,gamma,t(3),
    double precision:: M(3),S(3),P(3),Mp(3),Ssqr,Psqr,Msqr,ScrossMsqr,Mpsqr
    double precision::ex(3),ey(3),ez(3)
    double precision:: ScrossM(3),PcrossSdotM,SdotM,PdotS
    double precision:: cosbeta,sinbeta,Mpx,Mpy,Mx,My,Mz
    integer::k
    double precision::PcrossSsqr,PdotM,PcrossS(3)   !,M_PI,angle
    double precision::M_PI
    M_PI=4*atan(1.d0)


    Ssqr=0.0
    Msqr=0.0
    do k=1,3
       M(k)=to_align(k)-to_center(k)
       S(k)=axe(k)-center(k)
       Ssqr=Ssqr+S(k)*S(k)
       Msqr=Msqr+M(k)*M(k)
    end do

    SdotM=0.0
    do k=1,3
       SdotM=SdotM+M(k)*S(k)
    end do
    alpha=SdotM/Ssqr

    Psqr=0.0
    do k=1,3
       P(k)=M(k)-alpha*S(k)
       Psqr=Psqr+P(k)*P(k)
    end do

    PdotS=0.0
    do k=1,3
       PdotS=PdotS+P(k)*S(k)
    end do


    ScrossM(1)=S(2)*M(3)-M(2)*S(3)
    ScrossM(2)=S(3)*M(1)-M(3)*S(1)
    ScrossM(3)=S(1)*M(2)-M(1)*S(2)
    ScrossMsqr=0.0
    do k=1,3
       ScrossMsqr=ScrossMsqr+ScrossM(k)*ScrossM(k)
    end do

    cosbeta=SdotM/(sqrt(Ssqr*Msqr))
    sinbeta=sqrt(ScrossMsqr/(Ssqr*Msqr))
    Mx=SdotM/sqrt(Ssqr)
    My=sqrt(Psqr)
    Mpx=Mx*cosbeta-My*sinbeta
    Mpy=My*cosbeta+Mx*sinbeta

    Mpsqr=0.0
    do k=1,3
       Mp(k)=Mpx*P(k)/sqrt(Psqr)+Mpy*S(k)/sqrt(Ssqr)
       Mpsqr=Mpsqr+Mp(k)*Mp(k)
    end do

    PcrossS(1)=P(2)*S(3)-S(2)*P(3)
    PcrossS(2)=P(3)*S(1)-S(3)*P(1)
    PcrossS(3)=P(1)*S(2)-S(1)*P(2)
    PcrossSsqr=0.0
    do k=1,3
       PcrossSsqr=PcrossSsqr+PcrossS(k)*PcrossS(k)
    end do


    do k=1,3
       ex(k)=P(k)/sqrt(Psqr)
       ey(k)=S(k)/sqrt(Ssqr)
       ez(k)=PcrossS(k)/sqrt(PcrossSsqr)
    end do



    msg(__LINE__) ,"M=",M
    msg(__LINE__) ,"S=",S
    msg(__LINE__) ,"P=",P
    msg(__LINE__) ,"alpha=",alpha
    msg(__LINE__) ,"PdotS=",PdotS
    msg(__LINE__) ,"SdotM=",SdotM
    msg(__LINE__) ,"ScrossM=",ScrossM
    msg(__LINE__) ,"(ex,ex^2)=",ex,ex(1)*ex(1)+ex(2)*ex(2)+ex(3)*ex(3)
    msg(__LINE__) ,"(ey,ey^2)=",ey,ey(1)*ey(1)+ey(2)*ey(2)+ey(3)*ey(3)
    msg(__LINE__) ,"(ez,ez^2)=",ez,ez(1)*ez(1)+ez(2)*ez(2)+ez(3)*ez(3)
    msg(__LINE__) ,"(ex.ey)=",ex(1)*ey(1)+ex(2)*ey(2)+ex(3)*ey(3)
    msg(__LINE__) ,"(ez.ey)=",ez(1)*ey(1)+ez(2)*ey(2)+ez(3)*ey(3)
    msg(__LINE__) ,"(ex.ez)=",ex(1)*ez(1)+ex(2)*ez(2)+ex(3)*ez(3)
    msg(__LINE__) ,"(Mx,My)=",Mx,My,Mx*Mx+My*My
    msg(__LINE__) ,"(Mpx,Mpy)=",Mpx,Mpy,Mpx*Mpx+Mpy*Mpy
    msg(__LINE__) ,"cosbeta=",cosbeta
    msg(__LINE__) ,"sinbeta=",sinbeta
    msg(__LINE__) ,"cos^2(beta)+sin^2(beta)=",cosbeta*cosbeta+sinbeta*sinbeta
    msg(__LINE__) ,"Mp=",Mp
    msg(__LINE__) ,"|M|=",sqrt(MSqr)
    msg(__LINE__) ,"|Mp|=",sqrt(MpSqr)
    alpha=180.0*acos(cosbeta)/M_PI
  end subroutine TOOLS_build_rotation_axis

end module tools
