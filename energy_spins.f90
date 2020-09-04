module param
  integer nnvar   !dimension of the space
  integer dens    !points linear density in the rough search for the maximum 
  parameter (nnvar=2)!!!!changes do affect the writing format... 
                     !...in program, build_path and path_connection 
  parameter (dens=20)
  real*8 xlintol  !convergence criterion of the function variable on each linear search
  integer nitmax  !number of allowed line searchs
  real*8 ftol     !function convergence tolerance
  real*8 gtol     !gradient convergence tolerance
  parameter (xlintol=1.d-4, nitmax = 10 * nnvar, ftol = 0.d0, gtol = 1.d-7)
  real*8 mutol    !tolerance of the projection of gradient into so
  parameter (mutol=1.0d1)
  integer maxpt !maximum number of points in the mep
  parameter (maxpt=100)
  real*8 pi
  real*8 maxmin ! 1 for finding minimum (E,Grad) and -1 for maximum (-E,-Grad)
  real*8 mep(nnvar, maxpt), mepst(nnvar, maxpt) !mep and storager 
  integer nsadd !number of saddle points
  real*8, dimension(:,:), allocatable :: sadd, tsadd !storager of the saddle points and tangents
  real*8 divi !divisor of pi in the step of our "steepest descent" (step=pi/divi)
  parameter (divi=1.0d3)
  real*8, dimension(:,:), allocatable :: extrema !extrema for path connection
end module param

program auxiliar
  use param
  implicit none

  real*8 en, en2, p(nnvar), gr(nnvar), gr2(nnvar)
  real*8 energy
  external energy
  real*8 energy2
  external energy2

  pi=acos(-1.0d0)

  p(1)=pi
  p(2)=2*pi

  print *, p
  maxmin=1
  en=energy(p,gr)
  print *, p
  en2=energy2(p,gr2)
  print *, 'E=',en, en2
  print *, gr
  print *, gr2

end program auxiliar



!Energy and gradient for a finite chain of nnvar XY spins coupled 
!by exchange, dipolar and anisotropy interactions 
!Note: to transport this function to greater dimensions is non-trivial
!*************************************************************
real*8 function energy(cita,g)
  use param
  implicit none
  integer i, j, k
  real*8 cita(nnvar), g(nnvar)    
  real*8 Ee, Ea, Ed
  real*8 gee(nnvar), gea(nnvar), ged(nnvar)
  real*8 Jo(nnvar-1), KK(nnvar), kang(nnvar), Gd !model parameters
  real*8 rij

  !!!!!!!!!!! Duda: Las energias las cojo simples y el gradiente doble...?????

  !Define model parameters
  Jo=1.0d0
  Gd=0.0d0
  KK=4.0d0
  kang=0.0d0
  
  !!Calculating energy contributions
  Ee=0.0d0
  Ea=0.0d0
  Ed=0.0d0
  !exchange
  do i=1,nnvar-1
     j=i+1
     Ee=Ee-(Jo(i)*((cos(cita(i))*cos(cita(j)))+(sin(cita(i))*sin(cita(j)))))
  enddo
  !anisotropy
  do i=1,nnvar
     Ea=Ea-(KK(i)*(((cos(kang(i))*cos(cita(i)))+(sin(kang(i))*sin(cita(i))))**2))
  enddo
  !dipolar
  do i=1,nnvar-1
     do j=i+1,nnvar
        rij=real(j-i)
        Ed=Ed+(((-2*cos(cita(i))*cos(cita(j)))+(sin(cita(i))*sin(cita(j))))/rij/rij/rij)
     enddo
  enddo
  Ed=Gd*Ed

  !!Calculating gradients
  gee=0
  gea=0
  ged=0
  !exchange
  gee(1)=-Jo(1)*((-sin(cita(1))*cos(cita(2)))+(cos(cita(1))*sin(cita(2))))
  gee(nnvar)=-Jo(nnvar-1)*((-sin(cita(nnvar))*cos(cita(nnvar-1)))+(cos(cita(nnvar))*sin(cita(nnvar-1))))
  do k=2,nnvar-1
     i=k-1
     gee(k)=gee(k)-(Jo(i)*((-sin(cita(k))*cos(cita(i)))+(cos(cita(k))*sin(cita(i)))))
     i=k+1
     gee(k)=gee(k)-(Jo(k)*((-sin(cita(k))*cos(cita(i)))+(cos(cita(k))*sin(cita(i)))))
  enddo
  !anisotropy
  do k=1,nnvar
     gea(k)=-KK(k)*2*((cos(kang(k))*cos(cita(k)))+(sin(kang(k))*sin(cita(k))))*((-cos(kang(k))*sin(cita(k)))+(sin(kang(k))*cos(cita(k))))
  enddo
  !dipolar
  do k=1,nnvar
     do j=1,nnvar
        if (j.ne.k) then
           rij=real(abs(j-k))
           ged(k)=ged(k)+(((2*sin(cita(k))*cos(cita(j)))+(cos(cita(k))*sin(cita(j))))/rij/rij/rij)
        endif
     enddo
     ged(k)=Gd*ged(k)
  enddo

  !adding and inverting
  energy=maxmin*(Ee+Ea+Ed)
  g=maxmin*(gee+gea+ged) 
  return
end function energy

!Energy and gradient for the 1D xy model with anisotropy (J=1 and K=aJ) 
!**********************************************************************
real*8 function energy2(x,g)
  use param
  implicit none
  integer i
  real*8 x(nnvar), g(nnvar)
  real*8 a 
  parameter (a = 4)!!
 
  energy2 = 0.0d0
  do i=1,nnvar
     g(i)=0
  enddo
  do i=1,nnvar-1
     energy2=energy2+((cos(x(i))*cos(x(i+1)))+(sin(x(i))*sin(x(i+1))))!exchange
     energy2=energy2+(a*(cos(x(i))**2))!anisotropy
     g(i)=g(i)-((cos(x(i))*sin(x(i+1)))-(sin(x(i))*cos(x(i+1))))!exchange
     g(i+1)=g(i+1)-((sin(x(i))*cos(x(i+1)))-(cos(x(i))*sin(x(i+1))))!exchange
     g(i)=g(i)+(2*a*cos(x(i))*sin(x(i)))!anisotropy
  enddo
  energy2=energy2+(a*(cos(x(nnvar))**2))!anisotropy
  energy2=-energy2
  g(nnvar)=g(nnvar)+(2*a*cos(x(nnvar))*sin(x(nnvar)))!anisotropy
  !maxim switch
  energy2=maxmin*energy2
  g=maxmin*g
  return
end function energy2
