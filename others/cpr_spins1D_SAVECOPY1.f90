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

program CRP
  use param
  implicit none

  integer i, j, k, stat, re, actp, actp_c
  real*8 rr(nnvar), p_in(nnvar), p_fi(nnvar), p_mi(nnvar)
  real*8 ppi(nnvar), ppf(nnvar)
  logical done, found
  real*8 gr(nnvar)

  real*8 energy
  external energy

  pi=acos(-1.0d0)

  !File IO specification
600 format (2e20.8) !!!!!!!!!!!!!!!!!!!! this '2' must change as the dimension...
601 format (2e20.8,e20.8)!!!!!!!!!!!!!!!
  open (99, file='saddles.dat')
  write (99,*) '# List of the saddle points in the form'
  write (99,*) '# x(1) x(2)... x(n) E(x)'
  open (98, file='path.dat')
  write (98,*) '# List of the points along the MEP in the form'
  write (98,*) '# x(1) x(2)... x(n) E(x)'
  
  


  !Extrema configurations 
  !p_in(1)=0.0d0
  !p_in(2)=0.0d0
  !p_mi(1)=(pi/2)+0.1
  !p_mi(2)=(pi/2)-0.1
  !p_fi(1)=pi
  !p_fi(2)=pi

  p_in(1)=0.0d0
  p_in(2)=0.0d0
  p_fi(1)=pi
  p_fi(2)=2*pi

  !p_in(1)=0.623d0
  !p_in(2)=0.028d0 
  !p_fi(1)=-0.558d0
  !p_fi(2)=1.441d0

  !MEP initialization
  mep=-1.0d200
  mepst=-1.0d200

  mepst(:,1)=p_in(:)
  mepst(:,2)=p_fi(:)
  actp=2
  !mepst(:,1)=p_in(:)
  !mepst(:,2)=p_mi(:)
  !mepst(:,3)=p_fi(:)
  !actp=3

  mep=mepst
  actp_c=0
  ppi=p_in
  ppf=p_fi

  do i=1,maxpt
     do k=1,actp
        !write(*,'(2e20.8)') mepst(:,k)
        write(*,600) mep(:,k)
     enddo
     print *, ' '

     j=1
     found=.false.
     do while ((j<actp).and.(.not.found))
        call relax_maxim(mep(1,j),mep(1,j+1),rr,stat)
        if (stat.eq.0) then
           call act_mepst(rr,j+actp_c,actp+actp_c,1) !sobran los actp_c
           !print *, rr
           !print *, 'hi'
           actp_c=actp_c+1
           found=.true.
        elseif (j.ne.actp-1) then
           p_in(:)=mep(:,j)
           p_fi(:)=mep(:,j+2)
           rr(:)=mep(:,j+1)
           call relax_vertice(rr,p_in,p_fi,re)
           if (re.ne.-1) then
              if (re.eq.1) then
                 mepst(:,j+1)=rr
                 print *, '# Se refino el punto', j+1
              endif
              if (re.eq.0) then
                 call act_mepst(p_in,j+1+actp_c,actp+actp_c,0) !sobran los actp_c
                 actp_c=actp_c-1
              endif
              found=.true.
           endif
        endif
        !print *, j
        j=j+1
     enddo
     print *, '# Ciclo', i
     print *, '# Adiciono', actp_c, 'punto'
     print *, '#****'
     mep=mepst
     actp=actp+actp_c
     actp_c=0
     if ((j.eq.actp).and.(.not.found)) exit
  enddo
  
  if (i>=maxpt) then
     print *, '# Too large chain, something could be wrong'
     stop
  endif
  print *, '# Program finished in', i-1, 'useful cicles '
  print *, '# with an MEP of', actp, 'points '
  do k=1,actp
    !write(*,'(2e20.8)') mepst(:,k)
    write(*,600) mep(:,k)
  enddo

  maxmin=1 
  call set_sadd(actp)
  do k=1,nsadd
     write (99, 601) sadd(:,k), energy(sadd(1,k),gr)
  enddo

  call build_path()

  call path_connection(ppi,ppf)

end program CRP


!Test the good connection of the path and save minima
!****************************************************
subroutine path_connection(pin,pfi)
  use param
  implicit none
  
  real*8 pin(nnvar), pfi(nnvar)
  real*8 dp(nnvar), mdp
  real*8 ep, gr(nnvar)
  integer i
  logical conn

  real*8 energy
  external energy
  
  601 format (2e20.8,e20.8)!!!!!!!!!!!!!!! This 2 change as the dimension!

  conn=.true.

  if (nsadd>1) then
     write (99,*) '# Intermediate minima:'
  endif
   
  dp(:)=extrema(:,1)-pin(:)
  mdp=sqrt(dot_product(dp,dp))
  if (mdp>2*pi/divi) then
     print *, '# Path disconnection at the starting point!'
     conn=.false.
  endif
  dp(:)=pfi(:)-extrema(:,2*nsadd)
  mdp=sqrt(dot_product(dp,dp))
  if (mdp>2*pi/divi) then
     print *, '# Path disconnection at the ending point!'
     conn=.false.
  endif
  
  do i=1,nsadd-1
     dp(:)=extrema(:,2*i)-extrema(:,(2*(i+1))-1)
     mdp=sqrt(dot_product(dp,dp))
     if (mdp>2*pi/divi) then
        print *, '# Path disconnection betwen saddels',i,' and ',i+1 
        write (99,*) '# ...unfound minimum...'
        conn=.false.
     else
        dp(:)=extrema(:,2*i)
        ep=energy(dp,gr)
        write (99,601) dp, ep
     endif    
  enddo

  if (conn) print *, '# The path is well connected.'

  return
end subroutine path_connection



!Builds up the GDMEP and extrema
!****************************************************
subroutine build_path()
  use param
  implicit none
  
  integer i,k,n
  real*8 step, ep, p(nnvar), gp(nnvar), dp(nnvar), dpp(nnvar)
  real*8 energy
  character*8 fname
  external energy

600 format (2e20.8) !!!!!!!!!!!!!!!!!!!! this '2' must change as the dimension...
601 format (2e20.8,e20.8)!!!!!!!!!!!!!!!

  allocate (extrema(nnvar,2*nsadd))

  step=pi/divi

  do i=1,nsadd     
     !thoroug one side of the tangent
     n=2*i
     write (fname,'(I8)') n
     open(n,file=fname)
     step=pi/divi
     p(:)=sadd(:,i)+(step*tsadd(:,i))
     !do k=1,nnvar
     !   p(:)=sadd(:,i)+(step*tsadd(:,i))
     !enddo
     ep=energy(p,gp)
     write(n,601) p, ep
     write(98,601) p, ep
     dp=-step*gp/sqrt(dot_product(gp,gp))
     !gp=10
     do while (dot_product(gp,gp)>=gtol**2)
        p=p+dp
        ep=energy(p,gp)
        dpp=-step*gp/sqrt(dot_product(gp,gp))
        if (dot_product(dp,dpp)<=0) then
           p=p-dp
           exit
           !step=step/2
           !p=p-dpp
           !ep=energy(p,gp)
           !dp=-step*gp/sqrt(dot_product(gp,gp))
        else
           write(n,601) p, ep
           write(98,601) p, ep
           dp=dpp
        endif
     enddo
     extrema(:,n)=p(:)
    
     !thoroug the other side of the tangent
     n=(2*i)-1
     write (fname,'(I8)') n
     open(n,file=fname)
     step=pi/divi
     p(:)=sadd(:,i)-(step*tsadd(:,i))
     ep=energy(p,gp)
     write(n,601) p, ep
     write(98,601) p, ep
     dp=-step*gp/sqrt(dot_product(gp,gp))
     !gp=10
     do while (dot_product(gp,gp)>=gtol**2)
        p=p+dp
        ep=energy(p,gp)
        dpp=-step*gp/sqrt(dot_product(gp,gp))
        if (dot_product(dp,dpp)<=0) then
           p=p-dp
           exit
           !step=step/2
           !p=p-dpp
           !ep=energy(p,gp)
           !dp=-step*gp/sqrt(dot_product(gp,gp))
        else
           write(n,601) p, ep
           write(98,601) p, ep
           dp=dpp
        endif
     enddo
     extrema(:,n)=p(:)
  enddo
  
  return
end subroutine build_path



!Creates the array of the saddle points
!*************************************************************
subroutine set_sadd(pt)
  use param
  implicit none
    
  integer pt
  integer i
  real*8 e, ea, ed, gr(nnvar)
  real*8 p(nnvar), pa(nnvar), pd(nnvar), ta(nnvar), td(nnvar), tang(nnvar)
  real*8 tam, tdm, step
  real*8 tangtemp(nnvar,maxpt)
  real*8 energy
  external energy

  !Over all the non-extremum points
  nsadd=0
  do i=2,pt-1
     pa(:)=mep(:,i-1)
     p(:)=mep(:,i)
     pd(:)=mep(:,i+1)

     !Set and normalize tangent 
     td=pd-p
     ta=pa-p
     tdm=sqrt(dot_product(td,td))
     tam=sqrt(dot_product(ta,ta))
     td=td/tdm
     ta=ta/tam

     !Testing local maximum
     step=min(tdm,tam)/dens !step for test
     pa=p+(step*ta)
     pd=p+(step*td)
     maxmin=1  !normal function
     ea=energy(pa,gr)
     ed=energy(pd,gr)
     e=energy(p,gr)
     
     if ((ea<e).and.(e>ed)) then
        nsadd=nsadd+1
        mepst(:,nsadd)=mep(:,i)
        tang=td-ta
        tang=tang/sqrt(dot_product(tang,tang))
        tangtemp(:,nsadd)=tang(:)
     endif
  enddo
    
  allocate (sadd(nnvar,nsadd))
  allocate (tsadd(nnvar,nsadd))
  do i=1,nsadd
     sadd(:,i)=mepst(:,i)
     tsadd(:,i)=tangtemp(:,i)
  enddo
    
  return
end subroutine set_sadd


!Insert the vector r in mepst at position nn
!***************************************************************
subroutine act_mepst(r,nn,l,status)
  use param
  implicit none
  real*8 r(nnvar)
  integer nn, l
  integer i, status

  if (status.eq.0) then
    do i=nn,l-1
       mepst(:,i)=mepst(:,i+1)
    enddo
    mepst(:,l)=-1.0d200
    return 
  endif

  do i=l+1,nn+2,-1
     mepst(:,i)=mepst(:,i-1)
  enddo
  mepst(:,nn+1)=r(:)
  return
end subroutine act_mepst



!Relax the vertex if it is a non-saddle maximum
!*****************************************************************
subroutine relax_vertice(p, pa, pd, res)
  use param
  implicit none
  real*8 p(nnvar), pa(nnvar), pd(nnvar)
  integer res
  real*8 ta(nnvar), td(nnvar), tang(nnvar), tdm, tam, step
  real*8 e, ea, ed, gr(nnvar), di2, di(nnvar)
  real*8 fret
  integer niter

  real*8 energy
  external energy
  
  !Set and normalize tangent 
  td=pd-p
  ta=pa-p
  tdm=sqrt(dot_product(td,td))
  tam=sqrt(dot_product(ta,ta))
  td=td/tdm
  ta=ta/tam

  di=pd-pa
  di2=dot_product(di,di)
 
  !Testing local maximum
  step=min(tdm,tam)/dens !step for test
  pa=p+(step*ta)
  pd=p+(step*td)
  maxmin=1  !normal function
  ea=energy(pa,gr)
  ed=energy(pd,gr)
  e=energy(p,gr)

  if (((ea<e).and.(e>ed)).and.(dot_product(gr,gr)>gtol*gtol)) then
     print *, '# Refinando...' 
     tang=td-ta
     tang=tang/sqrt(dot_product(tang,tang))
     pa=p
     !Doing maximization in tang direction
     maxmin=-1 !minus the function
     call linmin(p,tang,nnvar,fret,xlintol)
     pa=pa-p 
     if (dot_product(pa,pa)>1.0d5*di2) then !condition of GOING AWAY!
        res=0                              !the vertice can not be refined, has to be removed
        return 
     endif
     !Doing the minimizations cojugated to tang
     maxmin=1 !normal function
     call conj_minim (p, nnvar, ftol, niter, fret, nitmax, xlintol, gtol, tang)  
     res=1                                 !the vertice was refined
     return
  endif

  res=-1                                   !the vertice is not a local maximum
  return
end subroutine relax_vertice





!Search a maximum if there are and returns it already cpr-relaxed 
!**************************************************************** 
subroutine relax_maxim(p_i,p_f,r,status)
  use param
  implicit none
  
  real*8 p_i(nnvar), p_f(nnvar)
  integer status
  real*8 dr(nnvar), r(nnvar), sso(nnvar)
  real*8 e, ea, ed, emax, rmax(nnvar)
  real*8 gr(nnvar)
  real*8 fret
  integer i, nspts
  integer niter
  
  real*8 energy
  external energy
  
  status=0  
    
  !Rough search for a maximum
  sso=p_f-p_i  !set direction so
  nspts=int(dens*sqrt(dot_product(sso,sso)))
  do i=1,nnvar
     dr(i)=(p_f(i)-p_i(i))/(nspts+1)
  enddo
  maxmin=1 !normal function
  r=p_i+dr
  ea=energy(p_i,gr)
  e=energy(r,gr)
  emax=-huge(emax)
  do i=1,nspts
     ed=energy(r+dr,gr)
     if ((ea<e).and.(e>ed)) then
        if (e>emax) then
           rmax=r
           emax=e
        endif
     endif
     ea=e
     e=ed
     r=r+dr
     !print *,r
  enddo
  if (emax.eq.-huge(emax)) then
     print *, '# el maximo es un extremo'
     ea=energy(p_i,gr)
     ed=energy(p_f,gr)
     if (ea<ed) then
        status=1
     else
        status=-1
     endif
     return
  endif
  sso=sso/sqrt(dot_product(sso,sso)) !normalization of so
  rmax=rmax-dr !precedent point (maybe not necessary)   

  !Doing the accurate maximization
  maxmin=-1 !minus the function for finding the maximum
  call linmin(rmax,sso,nnvar,fret,xlintol)
  !Doing the minimizations
  r=rmax
  maxmin=1 !normal function
  call conj_minim (r, nnvar, ftol, niter, fret, nitmax, xlintol, gtol, sso)  
  return
end subroutine relax_maxim



! Do the sucesive minimizations of CPR 
!**********************************************************************
subroutine conj_minim (P, N, FTOL, ITER, FRET, ITMAX, XLINTOL, gtol, so)
  ! Now modifyed for the implementation of the Conjugate
  ! Peak Refinement relation, this routine used to be:
  !
  !!! Numerical recipies' frprmn.f modified by Gus to
  !!! improve portability and handling.
  !!! Sept. 2000 
  !!!
  !!! ENERGY(P,XI) = N-variable Function to minimize (input, external) 
  !!! P  = point where ENERGY and XI are calculated
  !!! XI = gradient of FUNC 
  !!!    
  !!! FUNC = N-variable Function to minimize (input, external) OBSOLETO!
  !!! DFUNC= gradient of FUNC (input, external subroutine)     OBSOLETO!
  !!!
  !!! P = starting point (input) and solution (output)
  !!! N = number of variables
  !!! FTOL = Function convergence tolerance.
  !!! gtol = gradient convergence tolerance.
  !!! ITER = number of line searchs actually performed (output)
  !!! FRET = Function value at output P.
  !!! ITMAX = number of allowed line searchs (input)
  !!! XLINTOL = Convergence criterion of the function VARIABLE on each linear search
  !!!!!!
  !
  ! so = Direction of the maximization alrready performed (input)
  ! Hso = Estimation of $ H\cdot\vec{s_0} $
  ! ep = Auxiliar epsilon for finding hso
  ! pp, coef, soHso, gHso : are also used  

  implicit real*8 (A-H,O-Z)
  parameter (EPS=1.d-10)
  dimension P(N), G(N), H(N), XI(N), so(N)
  external energy
  
  parameter(ep=0.1) !epsilon value, it is unfortunately fixable
  dimension gr2(n), gr1(n), Hso(n), pp(n)

  !Calculation of hso
  pp=p+ep*so
  e=energy(pp,gr2)
  pp=p-ep*so
  e=energy(pp,gr1)
  Hso=(gr2-gr1)/2/ep

  !actualization of energy and gradient   
  FP = energy(P,XI)

  !Calculation of the coefficient "coef"
  gHso=0.0d0
  soHso=0.0d0
  do j=1,n
     gHso=gHso+XI(j)*Hso(j)
     soHso=soHso+so(j)*Hso(j)
  enddo
  coef=gHso/soHso

  !loop for determining the first direction
  do J=1,N
     G(J)=-XI(J)
     H(J)=G(J)+coef*so(J)
     XI(J)=H(J)
  enddo

  do ITS=1,ITMAX
     ITER=ITS
     call linmin (P, XI, N, FRET, XLINTOL)
     !        Next statement is normal Return
     if(2.d0*abs(FRET-FP).le.FTOL*(abs(FRET)+abs(FP)+EPS)) return
     FP = energy(P,XI)
     GG=0.d0
     DGG=0.d0
     ggnew = 0.d0
     do J=1,N
        ggnew = ggnew + XI(J) **2
        GG=GG+G(J)**2
        !           Fletcher-Reeves             
        !           DGG = DGG + XI(J)**2
        !           Polak-Ribiere
        DGG = DGG + (XI(J) + G(J)) * XI(J) !I'm keeping this as it is but can be wrong i have
     enddo                                 !to check the Polak-Ribiere CPR compqtibility  
     if (ggnew .lt. gtol**2) then
        !           Gradient convergence Return
        write (*,*) '# gradient return'
        FRET = FP
        return
     endif    
     if (GG .eq. 0.d0) return !!!!!!!!!!!!!!!!!!!!!!

     if (sqrt(real(nnvar))*dot_product(XI,so)/sqrt(ggnew)>mutol) then
        print *, '# Gradient*So exceeded mu tolerance. in', ITS, 'iterations'
        return
     endif

     !Calculation of the coefficient "coef"
     gHso=0.0d0
     do j=1,n
        gHso=gHso+XI(j)*Hso(j)
     enddo
     coef=gHso/soHso
     
     GAM=DGG/GG
     do J=1,N
        G(J)=-XI(J)
        H(J)=G(J)+(coef*so(J))+(GAM*H(J))
        XI(J)=H(J)
     enddo
  enddo

  write (*,*) '# FRPR maximum iterations exceeded'
  return
end subroutine conj_minim



!Numerical recipies' linmin.f follows modified to use gradient information
!******************************************************************
SUBROUTINE LINMIN (P, XI, N, FRET, XLINTOL)
  !     XLINTOL = Convergence criterion of the function VARIABLE
  !     on each linear search
  IMPLICIT real*8 (A-H,O-Z)
  !      PARAMETER (NMAX=1000,TOL=1.E-4)
  PARAMETER (NMAX=1000)
  !_no-gradient      EXTERNAL F1DIM
  EXTERNAL F1DIM, DF1DIM
  DIMENSION P(N),XI(N)
  COMMON /F1COM/ NCOM
  COMMON /F2COM/ PCOM(NMAX),XICOM(NMAX)
  !
  NCOM=N
  DO 11 J=1,N
     PCOM(J)=P(J)
     XICOM(J)=XI(J)
     !11   CONTINUE
11 enddo
  !
  !     We define the 3 starting points for bracketing such that
  !     the largest variable change is 0.1. Thought to be good 
  !     for problems involves angular variables to remain in 
  !     the closest local minimum.
  !
  gmax = 0.d0
  Do i = 1, N
     If (dabs(XICOM(i)) .gt. gmax) gmax = dabs(XICOM(i))
  EndDo
  If (gmax .eq. 0.d0) then
     write (*,*) '# Fatal in linmin: gradient equals 0 '
     stop
  EndIf
  !
  AX=0.d0
  !!      XX=1.
  !!      BX=2.
  XX=0.1d0 / gmax
  BX=0.2d0 / gmax
  !!      write (*,*) ' Entering MNBRAK '
  CALL MNBRAK(AX,XX,BX,FA,FX,FB,F1DIM)
  !!      write (*,*) ' Entering DBRENT '
  !_no-gradient      FRET=BRENT(AX,XX,BX,F1DIM,TOL,XMIN)
  FRET = DBRENT(AX, XX, BX, F1DIM, DF1DIM, XLINTOL, XMIN)
  DO 12 J=1,N
     XI(J)=XMIN*XI(J)
     P(J)=P(J)+XI(J)
     !12   CONTINUE
12 enddo
  RETURN
END SUBROUTINE LINMIN
!


!*******************************************************************
SUBROUTINE MNBRAK(AX,BX,CX,FA,FB,FC,FUNC)
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  PARAMETER (GOLD=1.618034d0, GLIMIT=100.d0, TINY=1.d-20)
  !      PARAMETER (GOLD=0.5, GLIMIT=100., TINY=1.E-20)
  FA=FUNC(AX)
  FB=FUNC(BX)
  IF(FB.GT.FA)THEN
     DUM=AX
     AX=BX
     BX=DUM
     DUM=FB
     FB=FA
     FA=DUM
  ENDIF
  CX=BX+GOLD*(BX-AX)
  FC=FUNC(CX)
1 IF(FB.GE.FC)THEN
     R=(BX-AX)*(FB-FC)
     Q=(BX-CX)*(FB-FA)
     U=BX-((BX-CX)*Q-(BX-AX)*R)/(2.d0*SIGN(MAX(ABS(Q-R),TINY),Q-R))
     ULIM=BX+GLIMIT*(CX-BX)
     IF((BX-U)*(U-CX).GT.0.d0)THEN
        FU=FUNC(U)
        IF(FU.LT.FC)THEN
           AX=BX
           FA=FB
           BX=U
           FB=FU
           GO TO 1
        ELSE IF(FU.GT.FB)THEN
           CX=U
           FC=FU
           GO TO 1
        ENDIF
        U=CX+GOLD*(CX-BX)
        FU=FUNC(U)
     ELSE IF((CX-U)*(U-ULIM).GT.0.d0)THEN
        FU=FUNC(U)
        IF(FU.LT.FC)THEN
           BX=CX
           CX=U
           U=CX+GOLD*(CX-BX)
           FB=FC
           FC=FU
           FU=FUNC(U)
        ENDIF
     ELSE IF((U-ULIM)*(ULIM-CX).GE.0.d0)THEN
        U=ULIM
        FU=FUNC(U)
     ELSE
        U=CX+GOLD*(CX-BX)
        FU=FUNC(U)
     ENDIF
     AX=BX
     BX=CX
     CX=U
     FA=FB
     FB=FC
     FC=FU
     GO TO 1
  ENDIF
  RETURN
END SUBROUTINE MNBRAK



!*******************************************************
FUNCTION DBRENT(AX,BX,CX,F,DF,TOL,XMIN)
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  PARAMETER (ITMAX=100,ZEPS=1.0d-10)
  External F, DF
  LOGICAL OK1,OK2
  A=MIN(AX,CX)
  B=MAX(AX,CX)
  V=BX
  W=V
  X=V
  E=0.d0
  FX=F(X)
  FV=FX
  FW=FX
  DX=DF(X)
  DV=DX
  DW=DX
  DO 11 ITER=1,ITMAX
     XM=0.5d0*(A+B)
     TOL1=TOL*ABS(X)+ZEPS
     TOL2=2.d0*TOL1
     IF(ABS(X-XM).LE.(TOL2-.5d0*(B-A))) GOTO 3
     IF(ABS(E).GT.TOL1) THEN
        D1=2.d0*(B-A)
        D2=D1
        IF(DW.NE.DX) D1=(W-X)*DX/(DX-DW)
        IF(DV.NE.DX) D2=(V-X)*DX/(DX-DV)
        U1=X+D1
        U2=X+D2
        OK1=((A-U1)*(U1-B).GT.0.d0).AND.(DX*D1.LE.0.d0)
        OK2=((A-U2)*(U2-B).GT.0.d0).AND.(DX*D2.LE.0.d0)
        OLDE=E
        E=D
        IF(.NOT.(OK1.OR.OK2))THEN
           GO TO 1
        ELSE IF (OK1.AND.OK2)THEN
           IF(ABS(D1).LT.ABS(D2))THEN
              D=D1
           ELSE
              D=D2
           ENDIF
        ELSE IF (OK1)THEN
           D=D1
        ELSE
           D=D2
        ENDIF
        IF(ABS(D).GT.ABS(0.5d0*OLDE))GO TO 1
        U=X+D
        IF(U-A.LT.TOL2 .OR. B-U.LT.TOL2) D=SIGN(TOL1,XM-X)
        GOTO 2
     END IF
1    IF(DX.GE.0.d0) THEN
        E=A-X
     ELSE
        E=B-X
     ENDIF
     D=0.5d0*E
2    IF(ABS(D).GE.TOL1) THEN
        U=X+D
        FU=F(U)
     ELSE
        U=X+SIGN(TOL1,D)
        FU=F(U)
        IF(FU.GT.FX)GO TO 3
     ENDIF
     DU=DF(U)
     IF(FU.LE.FX) THEN
        IF(U.GE.X) THEN
           A=X
        ELSE
           B=X
        ENDIF
        V=W
        FV=FW
        DV=DW
        W=X
        FW=FX
        DW=DX
        X=U
        FX=FU
        DX=DU
     ELSE
        IF(U.LT.X) THEN
           A=U
        ELSE
           B=U
        ENDIF
        IF(FU.LE.FW .OR. W.EQ.X) THEN
           V=W
           FV=FW
           DV=DW
           W=U
           FW=FU
           DW=DU
        ELSE IF(FU.LE.FV .OR. V.EQ.X .OR. V.EQ.W) THEN
           V=U
           FV=FU
           DV=DU
        ENDIF
     ENDIF
     !11        CONTINUE
11 enddo
  PAUSE '# DBRENT exceeded maximum iterations.'
3 XMIN=X
  DBRENT=FX
  RETURN
END FUNCTION DBRENT


!********************************************************************
FUNCTION F1DIM(X)
  IMPLICIT real*8 (A-H,O-Z)
  PARAMETER (NMAX=1000)
  COMMON /F1COM/ NCOM
  COMMON /F2COM/ PCOM(NMAX),XICOM(NMAX)
  DIMENSION XT(NMAX), gradaux(NMAX)
  DO 11 J=1,NCOM
     XT(J)=PCOM(J)+X*XICOM(J)
     !11   CONTINUE
11 enddo
  F1DIM=energy(XT,gradaux)
  RETURN
END FUNCTION F1DIM



!******************************************************************
FUNCTION DF1DIM(X)
  IMPLICIT real*8 (A-H,O-Z)
  PARAMETER (NMAX=1000)
  COMMON /F1COM/ NCOM
  COMMON /F2COM/ PCOM(NMAX),XICOM(NMAX)
  DIMENSION XT(NMAX),DF(NMAX)
  !
  DO 11 J=1,NCOM
     XT(J)=PCOM(J)+X*XICOM(J)
     !11   CONTINUE
11 enddo
  !!      CALL DFUNC(XT,DF)
  dum = energy(XT, DF)
  DF1DIM=0.d0
  DO 12 J=1,NCOM
     DF1DIM=DF1DIM+DF(J)*XICOM(J)
     !12   CONTINUE
12 enddo
  RETURN
END FUNCTION DF1DIM




!Energy and gradient for the 1D xy model with anisotropy (J=1 and K=aJ) 
!**********************************************************************
!real*8 function energy(x,g)
!  use param
!  implicit none
!  integer i
!  real*8 x(nnvar), g(nnvar)
!  real*8 a 
!  parameter (a = 4)!!
! 
!  energy = 0.0d0
!  do i=1,nnvar
!     g(i)=0
!  enddo
!  do i=1,nnvar-1
!     energy=energy+((cos(x(i))*cos(x(i+1)))+(sin(x(i))*sin(x(i+1))))!exchange
!     energy=energy+(a*(cos(x(i))**2))!anisotropy
!     g(i)=g(i)-((cos(x(i))*sin(x(i+1)))-(sin(x(i))*cos(x(i+1))))!exchange
!    g(i+1)=g(i+1)-((sin(x(i))*cos(x(i+1)))-(cos(x(i))*sin(x(i+1))))!exchange
!     g(i)=g(i)+(2*a*cos(x(i))*sin(x(i)))!anisotropy
!  enddo
!  energy=energy+(a*(cos(x(nnvar))**2))!anisotropy
!  energy=-energy
!  g(nnvar)=g(nnvar)+(2*a*cos(x(nnvar))*sin(x(nnvar)))!anisotropy
!  !maxim switch
!  energy=maxmin*energy
!  g=maxmin*g
!  return
!end function energy


!Muller-Brown potential 2D
!***************************************************************
!real*8 function energy(x,g)
!  use param
!  implicit none
!  integer i, icalle
!  real*8 x(nnvar), g(nnvar)  !nnvar has to be equal 2
!  real*8 AA(4),a(4),b(4),c(4),xo(4),yo(4)
!  real*8 temp
!  
!   
!  !3 minimus at:
!  ! x=0.623   y=0.028   -> E=-108.166
!  ! x=-0.050  y=0.466   -> E=-80.767
!  ! x=-0.558  y=1.441   -> E=-146.699
!     
!          
!  energy = 0.d0
!  do i=1,nnvar
!     g(i)=0
!  enddo
!  
!  !Muller-Brown parameters as they appear in 
!  !Quapp, J. Comp. Chem. Vol. 19, No. 9 1087-1100 (1998) 
!  AA(1)=-200.0
!  AA(2)=-100.0
!  AA(3)=-170.0
!  AA(4)=15.0
!  a(1)=-1.0
!  a(2)=-1.0
!  a(3)=-6.5
!  a(4)=0.7
!  b(1)=0.0
!  b(2)=0.0
!  b(3)=11.0
!  b(4)=0.6
!  c(1)=-10.0
!  c(2)=-10.0
!  c(3)=-6.5
!  c(4)=0.7
!  xo(1)=1.0
!  xo(2)=0.0
!  xo(3)=-0.5
!  xo(4)=-1.0
!  yo(1)=0.0
!  yo(2)=0.5
!  yo(3)=1.5
!  yo(4)=1.0
!  do i=1,4
!     temp=AA(i)*exp((a(i)*(x(1)-xo(i))**2)+(b(i)*(x(1)-xo(i))*(x(2)-yo(i)))+(c(i)*(x(2)-yo(i))**2))
!     energy=energy+temp
!     g(1)=g(1)+temp*((2*a(i)*(x(1)-xo(i)))+(b(i)*(x(2)-yo(i))))
!     g(2)=g(2)+temp*((b(i)*(x(1)-xo(i)))+(2*c(i)*(x(2)-yo(i))))
!  enddo  
!  energy=maxmin*energy
!  g=maxmin*g 
!  return
!end function energy
 

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

