      REAL FUNCTION rplam1(q,r)
      REAL*8 q
      REAL*8 r
      EXTERNAL ff1
      COMMON/z/z
      COMMON/PAWPAR/P(7)
C===  > CF proton-lambda
C---  > -------------------------> Constants --------------
      PARAMETER(pis=1.772453851) ! sqrt(pi)
      PARAMETER(hc=.1973)       ! in GeV*fm
c      PARAMETER(al=0.82)        ! suppression par. lambda
c      PARAMETER(r0=4.)          ! Emission Radius (fm)
c      PARAMETER(f0_s=2.31)      ! (-1)*(singlet scattering length)
c      PARAMETER(d0_s=3.04)      ! singlet effective radius (fm)
c      PARAMETER(f0_t=1.78)      ! (-1)*(triplet scattering length)
c      PARAMETER(d0_t=3.22)      ! triplet effective radius (fm)
c      PARAMETER(POL=0.)         ! Lambda Polarization
C---- ---------------------------------------------------
*     REAL ac(20)
*     DATA ac/4*0.,-388.,0.,388.,6*0.,-111.,111.,5*0./
C---- ---------------------------------------------------
      DATA p/1.0,10.,2.31,3.04,1.78,3.22,0./
      ak=q/2

      print *, q

      al=p(1)
c      r0=p(2)
      r0=r
      f0_s=p(3)
      d0_s=p(4)
      f0_t=p(5)
      d0_t=p(6)
      POL=p(7)

c      PRINT *, f0_s

      r0m=hc/r0

      f0m_s=(hc/f0_s+0.5*ak*ak*d0_s/hc)
      den_s=f0m_s**2+ak**2
      fre_s=f0m_s/den_s
      fim_s=ak/den_s
      fsq_s=fre_s**2+fim_s**2

      f0m_t=(hc/f0_t+0.5*ak*ak*d0_t/hc)
      den_t=f0m_t**2+ak**2
      fre_t=f0m_t/den_t
      fim_t=ak/den_t
      fsq_t=fre_t**2+fim_t**2

      z=2*ak/r0m
      zs=z*z
      ex=exp(-zs)

      call integr(ai,z)
*      ai=0.

      fsi_s=0.5*(fsq_s*r0m**2*(1-d0_s/(2*pis*r0))+4*fre_s*r0m*ai/pis)
      fsi_s=fsi_s-fim_s*r0m*(1-ex)/z
      fsi_t=0.5*(fsq_t*r0m**2*(1-d0_t/(2*pis*r0))+4*fre_t*r0m*ai/pis)
      fsi_t=fsi_t-fim_t*r0m*(1-ex)/z
*      rplam1=ai
*      return

      R_t=.25*(3+POL**2)*(1+al*fsi_t)
      R_s=.25*(1-POL**2)*(1+al*fsi_s)
      rplam1=R_t+R_s

C      print *, q,R_t,R_s,r_p_lam
      print *, rplam1
      RETURN
      END

      FUNCTION ff1(x)
      COMMON/z/z
      ff1=exp(x*x-z*z)/z
      RETURN
      END

      SUBROUTINE integr(ai,z)
C-----------------------------------------------------------------------
c      IMPLICIT REAL*8 (A-H,O-Z)
c      REAL*4 ai,z
      EXTERNAL FUN1
      zz=z
C---Setting integration accuracy
      EPS=1.0E-18
C      EPS=1.0-1
      DDD=.005
      N=100
      DE1=zz/N

C***********************************************************************
      CALL SIMPS(0.0,zz,DE1,DDD,EPS,FUN1,xx,AII,A1H,ABS)
      ai=AII
      return
      END

      FUNCTION FUN1(xx)
c      IMPLICIT REAL*8 (A-H,O-Z)
c      REAL*4 x,ff1
      x=xx
      FUN1=ff1(x)
      RETURN
      END

      SUBROUTINE SIMPS(A1,B1,H1,REPS1,AEPS1,FUNCT,X,AI,AIH,AIABS)
C SIMPS
C A1,B1 -THE LIMITS OF INTEGRATION
C H1    -AN INITIAL STEP OF INTEGRATION
C REPS1,AEPS1 - RELATIVE AND ABSOLUTE PRECISION OF INTEGRATION
C FUNCT -A NAME OF FUNCTION SUBPROGRAM FOR CALCULATION OF INTEGRAND +
C X - AN ARGUMENT OF THE INTEGRAND
C AI - THE VALUE OF INTEGRAL
C AIH- THE VALUE OF INTEGRAL WITH THE STEP OF INTEGRATION
C AIABS- THE VALUE OF INTEGRAL FOR MODULE OF THE INTEGRAND
C THIS SUBROGRAM CALCULATES THE DEFINITE INTEGRAL WITH THE RELATIVE OR
C ABSOLUTE PRECISION BY SIMPSON+S METHOD WITH THE AUTOMATICAL CHOICE
C OF THE STEP OF INTEGRATION
C IF AEPS1    IS VERY SMALL(LIKE 1.E-17),THEN CALCULATION OF INTEGRAL
C WITH REPS1,AND IF REPS1 IS VERY SMALL (LIKE 1.E-10),THEN CALCULATION
C OF INTEGRAL WITH AEPS1
C WHEN  AEPS1(REPS1(0. THEN CALCULATION WITH THE CONSTANT STEP H1
C
c      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(7),P(5)
      H=SIGN(H1,B1-A1)
      S=SIGN(1.0,H)
      A=A1
      B=B1
      AI=0.0
      AIH=0.0
      AIABS=0.0
      P(2)=4.0
      P(4)=4.0
      P(3)=2.0
      P(5)=1.0
      IF(B-A) 1,2,1
    1 REPS=ABS(REPS1)
      AEPS=ABS(AEPS1)
      DO 3 K=1,7
  3   F(K)=10.E16
      X=A
      C=0.0
      F(1)=FUNCT(X)/3
    4 X0=X
      IF((X0+4*H-B)*S) 5,5,6
    6 H=(B-X0)/4
      IF(H) 7,2,7
    7 DO 8 K=2,7
  8   F(K)=10.E16
      C=1.0
    5 DI2=F(1)
      DI3=ABS(F(1))
      DO 9 K=2,5
      X=X+H
      IF((X-B)*S) 23,24,24
   24 X=B
   23 IF(F(K)-10.E16) 10,11,10
   11 F(K)=FUNCT(X)/3
   10 DI2=DI2+P(K)*F(K)
    9 DI3=DI3+P(K)*ABS(F(K))
      DI1=(F(1)+4*F(3)+F(5))*2*H
      DI2=DI2*H
      DI3=DI3*H
      IF(REPS) 12,13,12
   13 IF(AEPS) 12,14,12
   12 EPS=ABS((AIABS+DI3)*REPS)
      IF(EPS-AEPS) 15,16,16
   15 EPS=AEPS
   16 DELTA=ABS(DI2-DI1)
      IF(DELTA-EPS) 20,21,21
   20 IF(DELTA-EPS/8) 17,14,14
   17 H=2*H
      F(1)=F(5)
      F(2)=F(6)
      F(3)=F(7)

      DO 19 K=4,7
  19  F(K)=10.E16
      GO TO 18
   14 F(1)=F(5)
      F(3)=F(6)
      F(5)=F(7)
      F(2)=10.E16
      F(4)=10.E16
      F(6)=10.E16
      F(7)=10.E16
   18 DI1=DI2+(DI2-DI1)/15
      AI=AI+DI1
      AIH=AIH+DI2
      AIABS=AIABS+DI3
      GO TO 22
   21 H=H/2
      F(7)=F(5)
      F(6)=F(4)
      F(5)=F(3)
      F(3)=F(2)
      F(2)=10.E16
      F(4)=10.E16
      X=X0
      C=0.0
      GO TO 5
   22 IF(C) 2,4,2
    2 RETURN
         END


