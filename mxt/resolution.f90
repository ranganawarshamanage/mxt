MODULE resolution_mod

REAL, PUBLIC :: a=0.,b,c,alf,bet,gam
PRIVATE
REAL(kind(1.d0)) :: as,bs,cs,cosgs,cosbs,cosas
LOGICAL :: firstcall=.TRUE.
REAL :: stol2min
REAL, ALLOCATABLE :: stol2max(:)
INTEGER :: nbin

PUBLIC :: getres, getbin, setbin, getlim, countbin

contains

SUBROUTINE getres(ih,ik,il,s2,resol)
IMPLICIT NONE
INTEGER, INTENT(IN):: ih,ik,il
REAL, INTENT(OUT)  :: s2,resol
REAL(kind(1.d0)) :: vol,sina,sinb,sing,cosa,cosb,cosg

IF (firstcall) THEN
  firstcall = .FALSE.
  if (a==0.) stop 'getres: a,b,c,alf,bet,gam unknown'
  cosa=COS(alf/57.295779513082320877d0)
  cosb=COS(bet/57.295779513082320877d0)
  cosg=COS(gam/57.295779513082320877d0)
  sina=SIN(alf/57.295779513082320877d0)
  sinb=SIN(bet/57.295779513082320877d0)
  sing=SIN(gam/57.295779513082320877d0)
  cosas=(cosb*cosg-cosa)/(sinb*sing)
  cosbs=(cosa*cosg-cosb)/(sina*sing)
  cosgs=(cosa*cosb-cosg)/(sina*sinb)
  vol=a*b*c*SQRT(1.d0-cosa**2-cosb**2-cosg**2+2.d0*cosa*cosb*cosg)
  as=b*c*sina/vol
  bs=a*c*sinb/vol
  cs=a*b*sing/vol
END IF
s2=((ih*as)**2+(ik*bs)**2+(il*cs)**2+2.d0*ih*ik*as*bs*cosgs  &
      +2.d0*ih*il*as*cs*cosbs+2.d0*ik*il*bs*cs*cosas)/4.d0
         resol = SQRT(S2)
         RESOL=1/(2.*resol)
END SUBROUTINE getres

SUBROUTINE setbin(nbin_temp,reslim1,reslim2,mode,res_arr)
  IMPLICIT NONE
  REAL, INTENT(IN) :: reslim1,reslim2
  INTEGER, INTENT(IN) :: nbin_temp,mode
  real, intent(out),dimension(:),allocatable :: res_arr
  REAL stolmax3,stolmin3,stolinc,stol3max
  INTEGER i

  nbin=nbin_temp
  if(allocated(res_arr)) deallocate(res_arr)
  allocate(res_arr(nbin))
  print*, nbin, mode, reslim1,reslim2
  if(allocated(stol2max))DEALLOCATE(stol2max)
  ALLOCATE(stol2max(nbin))
  IF (mode==1) THEN  ! equal volumes in reciprocal space
     stolmax3=(1./(2.*reslim2))**3
     stolmin3=(1./(2.*reslim1))**3
     stol2min=EXP(LOG(stolmin3)*2./3.)
  ELSE               ! linearly increasing volumes 
     stolmax3=(1./(2.*reslim2))**2
     stolmin3=(1./(2.*reslim1))**2
     stol2min=stolmin3
  END IF
  stolinc=(stolmax3-stolmin3)/nbin
  DO i=1,nbin
     stol3max=stolinc*i+stolmin3
     IF (mode==1) THEN
        stol2max(i)=EXP(LOG(stol3max)*2./3.)
     ELSE
        stol2max(i)=stol3max
     END IF
     res_arr(i) = sqrt(1./stol2max(i))/2.
     print*,'ibin, dmax=',i,sqrt(1./stol2max(i))/2.
  END DO
END SUBROUTINE setbin

SUBROUTINE getbin(s2,ibin)
! find resolution bin
IMPLICIT NONE
REAL, INTENT(IN) :: s2
INTEGER, INTENT(OUT) :: ibin
INTEGER i

DO i=1,nbin
  IF (s2 <= stol2max(i)) EXIT             
END DO
ibin=MIN(i,nbin)   ! numerically it could happen that i==nbin+1

END SUBROUTINE getbin

SUBROUTINE getlim(reslim)
! return resolution limits
IMPLICIT NONE
REAL, INTENT(OUT) :: reslim(nbin)
INTEGER i

DO i=1,nbin
  reslim(i) = sqrt(1./stol2max(i))/2         
END DO

END SUBROUTINE getlim

SUBROUTINE countbin(ispgr,shellcount)
! count unique reflections in resolution shells
USE symops_mod, ONLY:  asuput,op,ncenop
IMPLICIT NONE
INTEGER, INTENT(IN) :: ispgr
INTEGER, INTENT(OUT) :: shellcount(nbin)
INTEGER ih,ik,il,jh,jk,jl,ibin,iop,ihkl,ihmax,ikmax,ilmax,ihmin,ikmin,ilmin
INTEGER jhmax,jkmax,jlmax
REAL s2
LOGICAL(1), ALLOCATABLE :: hkl(:,:,:)

ihmin=HUGE(ihmin)
ikmin=ihmin
ilmin=ihmin
ihmax=-ihmin
ikmax=-ihmin
ilmax=-ihmin
jhmax=INT(2*SQRT(stol2max(nbin))/as)
jkmax=INT(2*SQRT(stol2max(nbin))/bs)
jlmax=INT(2*SQRT(stol2max(nbin))/cs)
!
! THIS IS NEITHER ELEGANT NOR EFFICIENT! but it should work
!
!$omp parallel do collapse(3) private(il,ik,ih,s2,jh,jk,jl,ihkl,iop) &
!$    reduction(MIN:ihmin,ikmin,ilmin) reduction(MAX:ihmax,ikmax,ilmax)
DO il=0,jlmax ! ich weiﬂ nicht, ob das bei Winkel /= 90 reicht!
DO ik=-jkmax,jkmax
DO ih=-jhmax,jhmax
  s2=((ih*as)**2+(ik*bs)**2+(il*cs)**2+2.d0*ih*ik*as*bs*cosgs  &
      +2.d0*ih*il*as*cs*cosbs+2.d0*ik*il*bs*cs*cosas)/4.d0
  IF (s2 < stol2min .OR. s2 > stol2max(nbin)) CYCLE
  CALL asuput(ih,ik,il,jh,jk,jl,ihkl,iop,op,size(op,3)/ncenop)
  ihmin=MIN(ihmin,jh)
  ihmax=MAX(ihmax,jh)
  ikmin=MIN(ikmin,jk)
  ikmax=MAX(ikmax,jk)
  ilmin=MIN(ilmin,jl)
  ilmax=MAX(ilmax,jl)
END DO
END DO
END DO
!print*,'maxhkl=',ihmax,ikmax,ilmax,ihmin,ikmin,ilmin
ALLOCATE(hkl(ihmin:ihmax,ikmin:ikmax,ilmin:ilmax))
hkl=.FALSE.
shellcount=0
!!$omp parallel do collapse(3) private(il,ik,ih,s2,jh,jk,jl,ihkl,iop,ibin)
DO il=0,jlmax ! ich weiﬂ nicht, ob das bei Winkel /= 90 reicht!
DO ik=-jkmax,jkmax
DO ih=-jhmax,jhmax
  SELECT CASE (ncenop)
  CASE(2)   
    IF (ispgr==5.OR.ispgr==20.OR.ispgr==21) THEN  ! C2, C222, C2221
      IF (MOD(ih+ik,2)/=0) CYCLE 
    ELSE ! I222, I212121,I4,I41,I422, I4122, I23, I213, I432, I4132
      IF (MOD(ih+ik+il,2)/=0) CYCLE
    END IF
  CASE(3) ! H3, H32
    IF (MOD(-ih+ik+il,3)/=0) CYCLE  
  CASE(4)   ! F222, F23, F432, F4132
    IF (MOD(ih+ik,2)/=0.OR.MOD(ih+il,2)/=0.OR.MOD(ik+il,2)/=0) CYCLE 
  END SELECT
  s2=((ih*as)**2+(ik*bs)**2+(il*cs)**2+2.d0*ih*ik*as*bs*cosgs  &
      +2.d0*ih*il*as*cs*cosbs+2.d0*ik*il*bs*cs*cosas)/4.d0
  IF (s2 < stol2min .OR. s2 > stol2max(nbin)) CYCLE
  CALL asuput(ih,ik,il,jh,jk,jl,ihkl,iop,op,size(op,3)/ncenop)
  IF (jh>ihmax.OR.jh<ihmin.OR.jk>ikmax.OR.jk<ikmin.OR.jl<ilmin.OR.jl>ilmax) THEN
    print '(a,6i4)','countbin error',jh,jk,jl,ihmax,ikmax,ilmax
  ELSE
    IF (.NOT.hkl(jh,jk,jl)) THEN
      hkl(jh,jk,jl)=.TRUE.
      CALL getbin(s2,ibin)
!racy:
      shellcount(ibin)=shellcount(ibin)+1
    END IF
  END IF
END DO
END DO
END DO
DEALLOCATE(hkl)
END SUBROUTINE countbin

END MODULE resolution_mod
