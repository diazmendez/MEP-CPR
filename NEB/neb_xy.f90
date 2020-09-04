!***************************************************************
!  neb.f90
!**********************************
! Este programa implementa el Nudge Elastic Band algorithm modificado 
! por Henkelman and Jonsson y encuentra la MEP entre dos minimos de energia
! de un cierto sistema.
!**********************************
! Rogelio Diaz-Mendez
! Toulouse, octubre de 2005.
!****************************************************************

!**************
program neb
  implicit none
  
  !simulaiton parameters
  integer ndim_       !dimensionality of the configurational space  
  integer nimages_    !number of images
  parameter (nimages_=8, ndim_=2)
  
  !simulation variables
  real*8   tvectors_(ndim_,nimages_+2)
  real*8   path_(ndim_,nimages_+2)
  
  integer i
  
  call set_initial_path(ndim_,nimages_,path_)
  write(*,'(2F15.8)') path_
  print *, ' '
  
  do i=1, 10 !*** En realidad esto sera un  while ()
     call set_tangent_vectors(ndim_,nimages_,path_,tvectors_)
     call vverlet(ndim_,nimages_,path_,tvectors_)
     write(*,'(2F15.8)') path_ ! Antes de la F va el numero ndim_
     print *, ' '
  end do
  
end program neb


!Calcula los versores tangentes a cada imagen para la minimizacion
!siguiendo la propuesta de J. Chem. Phys., Vol. 113, No 22 (2000).
!**********
subroutine set_tangent_vectors(nd,ni,pth,tvec)
  implicit none
  
  integer  nd
  integer  ni
  real*8   pth(nd,ni)
  real*8   tvec(nd,ni+2)
  
  real*8   energy_ 
  external energy_  
  
  integer i,k
  real*8 e
  real*8 ea
  real*8 ed  
  real*8 tmas(nd)
  real*8 tmen(nd)
  real*8 dvmax
  real*8 dvmin
  
  tvec=0
  do i=2,ni+1
     e=energy_(nd,pth(1,i))
     ea=energy_(nd,pth(1,i-1))
     ed=energy_(nd,pth(1,i+1))
     if ((ed>e).and.(e>ea)) then
        do k=1,nd
           tvec(k,i)=pth(k,i+1)-pth(k,i)
        enddo
     else
        if ((ed<e).and.(e<ea)) then
           do k=1,nd
              tvec(k,i)=pth(k,i)-pth(k,i-1)
           enddo
        else
           do k=1,nd
              tmas(k)=pth(k,i+1)-pth(k,i)
              tmen(k)=pth(k,i)-pth(k,i-1)
           enddo
           dvmax=abs(ed-e) ! esto se liga al if que le sigue...
           dvmin=abs(ea-e)
           if (dvmax<dvmin) then
              e=dvmax  ! e se redunda para ahorrar memoria (8bytes!?)
              dvmax=dvmin
              dvmin=e
           endif
           if (ed>ea) then
              do k=1,nd
                 tvec(k,i)=(dvmax*tmas(k))+(dvmin*tmen(k))
              enddo
           else
              do k=1,nd
                 tvec(k,i)=(dvmin*tmas(k))+(dvmax*tmen(k))
              enddo
           endif
        endif
     endif
     e=0  ! e vuelve a redundarse, ahora sera el modulo de la tangente
     do k=1,nd
        e=e+(tvec(k,i)**2)
     enddo
     e=sqrt(e)
     if (e.ne.0) then
        do k=1,nd
           tvec(k,i)=tvec(k,i)/e
        enddo
     endif
  enddo
  
  return
end subroutine set_tangent_vectors


! Devuelve la energia de una configuracion p. Esta funcion tiene que 
! matchear con el gradiente!!!!!
!**********
real*8 function energy_(nd_,p)  
  implicit none
  
  integer nd_
  real*8  p(nd_)
  integer i
  
  !energia para un arreglo unidimensional de espines xy interactuando
  !via exchange ferromagnetico con J=1
  energy_=0
  do i=1,nd_-1
     energy_=energy_+((cos(p(i))*cos(p(i+1)))+(sin(p(i))*sin(p(i+1))))
  enddo
  energy_=-energy_
  
  return
end function energy_


! Define la trayectoria inicial en el espacio de configuraciones
!**********
subroutine set_initial_path (nd,ni,pth)
  implicit none
  
  integer  nd
  integer  ni
  real*8   pth(nd,ni+2)
  
  integer i,j 
  real*8 dx(nd)
  
  
  !Declaracion de extremos
  pth(1:nd,1)=0                 
  pth(1:nd,ni+2)=2*3.14159
  
  !Inicializacion del path en la linea recta
  do j=1,nd
     dx(j)=(pth(j,ni+2)-pth(j,1))/(ni+1)
  enddo
  do i=2,ni+1
     do j=1,nd
        pth(j,i)=pth(j,1)+((i-1)*dx(j))+((dx(j)/6)*(-1)**(j+i))! SE ESTA CORROMPIENDO LA TRAYECTORIA
     enddo                                                     ! CON EL ULTIMO TERMINO
  enddo
  return
end subroutine set_initial_path


! Hace la dinamica de velocidades de Verlet
!*************************************************************************!
subroutine vverlet (ndim,nparts,position,tvec)
  implicit none
  
  integer ndim       
  integer nparts    
  real*8  position(ndim,nparts+2)
  real*8  tvec(ndim,nparts+2)
  
  ! simulation parameters
  integer nsteps     ! number of time steps in the simulation !this should be modified if  
  parameter(nsteps=100000)                                    !we perform a while()
  real*8 mass        ! mass of the particles
  real*8 be          ! viscossity coefficient 
  real*8 dt          ! time step
  parameter(mass=0.001,be=0.1,dt=1.0e-5)
  
  ! simulation variables
  real*8 velocity(ndim,nparts+2)
  real*8 force(ndim,nparts+2)
  real*8 accel(ndim,nparts+2)
  integer i,j,k
    
  ! set initial velocities, and accelerations
  velocity=0
  accel=0
  
  
  ! This is the main time stepping loop
  do i=1, nsteps    ! This should be changed by a while()
     call compute(nparts,ndim,position,velocity,mass,force,be,tvec)
     call update(nparts,ndim,position,velocity,force,accel,mass,dt)
  enddo
  return
end subroutine vverlet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the forces, given positions, masses,
! and velocities
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute(np,nd,pos,vel,mass,f,b,tv)
  implicit none
  
  integer np
  integer nd
  real*8  pos(nd,np+2)
  real*8  vel(nd,np+2)
  real*8  f(nd,np+2)
  real*8  mass
  real*8  b
  real*8  tv(nd,np+2)
  
  real*8 tang(nd)
  real*8 grad(nd)
  real*8 fspring(nd)
  
  ! Funciones producto interno y distancia
  real*8 dotr8
  external dotr8
  real*8 dist
  external dist
  
  integer i, k
  real*8  PI2
  parameter(PI2=3.14159265d0/2.0d0)
  
  !Clean the forces
  f=0
  
  ! HAY QUE SABER BIEN COMO DIABLOS FUNCIONA ESTO...
  ! The computation of forces and energies is fully parallel.
  !$omp  parallel do
  !$omp& default(shared)
  !$omp& private(i,j,k,rij,d)
  !$omp& reduction(+ : pot, kin)
  do i=2,np+1
     !se calcula el gradiente del potencial de energia dado en energy_()
     !Caso particular de un arreglo unidimensional de np espines xy 
     !interactuando via exchange ferromagnetico con J=1
     grad=0
     do k=1,nd
        if (k<nd) then; grad(k)=grad(k)+(cos(pos(k,i))*sin(pos(k+1,i)))-(sin(pos(k,i))*cos(pos(k+1,i))) 
        endif
        if (k>1) then; grad(k)=grad(k)+((sin(pos(k-1,i))*cos(pos(k,i)))-(cos(pos(k-1,i))*sin(pos(k,i))))
        endif
        tang(k)=tv(k,i) !updating de la tangente (ESTA LINEA NO TIENE QUE VER CON EL GRADIENTE!!!)
     enddo
     grad=-grad
     ! ...de paso, ya hemos aprovechado y YA updateamos a la tang 
     
     !Se mutila la direccion tangengial del gradiente (intrinseco de NEB)
     grad=grad-(dotr8(nd,grad,tang)*tang)
     
     !Se Calcula esta fuerza de "spring" o cuasielastica SUSTITUYENDO K POR b (ES SOLO UNA IDEA!)
     fspring=b*(dist(nd,pos(1,i+1),pos(1,i))-dist(nd,pos(1,i),pos(1,i-1)))*tang 
     
     do k=1,nd
        !se adiciona el termino debido al gradiente con la restriccion de las tangentes
        f(k,i)=f(k,i)-grad(k)
        !se adiciona el termino cuasielastico que garantiza la equidistancia
        f(k,i)=f(k,i)+fspring(k)
        !se introduce un termino viscoso para frenar en el minimo
        f(k,i)=f(k,i)-(b*vel(k,i))
     enddo
  enddo
  
  return
end subroutine compute




!Compute the norm of the difference betwen two vectors in an n-dimensional space.
real*8 function dist(n,r1,r2)
  implicit none
  
  integer n
  real*8 r1(n)
  real*8 r2(n)
  
  integer i
  
  dist = 0.0
  do i=1,n
     dist = dist + (r1(i) - r2(i))**2.
  enddo
  dist = sqrt(dist)
  
  return
end function dist


! Return the dot product between two vectors of type real*8 and length n
real*8 function dotr8(n,x,y)
  implicit none
  
  integer n
  real*8 x(n)
  real*8 y(n)
  
  integer i
  
  dotr8 = 0.0
  do i = 1,n
     dotr8 = dotr8 + x(i)*y(i)
  enddo
  
  return
end function dotr8


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Perform the time integration, using a velocity Verlet algorithm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine update(np,nd,pos,vel,f,a,mass,dt)
  implicit none
  
  integer np
  integer nd
  real*8  pos(nd,np+2)
  real*8  vel(nd,np+2)
  real*8  f(nd,np+2)
  real*8  a(nd,np+2)
  real*8  mass
  real*8  dt
  
  integer i, j
  real*8  rmass
  
  rmass = 1.0/mass
  
  ! The time integration is fully parallel
  !$omp  parallel do
  !$omp& default(shared)
  !$omp& private(i,j)
  do i = 2,np+1
     do j = 1,nd
        pos(j,i) = pos(j,i) + vel(j,i)*dt + 0.5*dt*dt*a(j,i)
        vel(j,i) = vel(j,i) + 0.5*dt*(f(j,i)*rmass + a(j,i))
        a(j,i) = f(j,i)*rmass      !** There's a better way to do this
     enddo
  enddo
  !$omp  end parallel do
  
  return
end subroutine update

