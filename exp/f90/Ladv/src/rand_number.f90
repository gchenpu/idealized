! f77 code taken from numerical recipes 

module rand_number_mod

implicit none

public :: ran, ran2, gasdev

contains

!========================= ran ======================================
FUNCTION ran(idum)
   IMPLICIT NONE
   INTEGER, PARAMETER :: K4B=selected_int_kind(9)
   INTEGER(K4B), INTENT(INOUT) :: idum
   REAL :: ran
   INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
   REAL, SAVE :: am
   INTEGER(K4B), SAVE :: ix=-1,iy=-1,k
   if (idum <= 0 .or. iy < 0) then 
      am=nearest(1.0,-1.0)/IM
      iy=ior(ieor(888889999,abs(idum)),1)
      ix=ieor(777755555,abs(idum))
      idum=abs(idum)+1 
   end if
   ix=ieor(ix,ishft(ix,13)) 
   ix=ieor(ix,ishft(ix,-17))
   ix=ieor(ix,ishft(ix,5))
   k=iy/IQ 
   iy=IA*(iy-k*IQ)-IR*k 
   if (iy < 0) iy=iy+IM
   ran=am*ior(iand(IM,ieor(ix,iy)),1) 
END FUNCTION ran

!========================= ran2 ======================================
! returns a uniform random deviate between 0.0 and 1.0
FUNCTION ran2(idum)
   INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
   REAL ran2,AM,EPS,RNMX
   PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
              IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,  &
              IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS) 
   INTEGER idum2,j,k,iv(NTAB),iy
   SAVE iv,iy,idum2
   DATA idum2/123456789/, iv/NTAB*0/, iy/0/
   if (idum.le.0) then
      idum=max(-idum,1)
      idum2=idum
      do j=NTAB+8,1,-1
         k=idum/IQ1
         idum=IA1*(idum-k*IQ1)-k*IR1
         if (idum.lt.0) idum=idum+IM1
         if (j.le.NTAB) iv(j)=idum
      enddo 
      iy=iv(1)
   endif
   k=idum/IQ1 
   idum=IA1*(idum-k*IQ1)-k*IR1 
   if(idum.lt.0) idum=idum+IM1 
   k=idum2/IQ2
   idum2=IA2*(idum2-k*IQ2)-k*IR2 
   if (idum2.lt.0) idum2=idum2+IM2
   j=1+iy/NDIV 
   iy=iv(j)-idum2 
   iv(j)=idum 
   if(iy.lt.1)iy=iy+IMM1
   ran2=min(AM*iy,RNMX)
return
END FUNCTION ran2

!========================= gasdev ======================================
FUNCTION gasdev(idum)
! Use ran
! returns a normally distributed deviate with zero mean and unit variance.
   INTEGER idum
   REAL gasdev
   INTEGER iset
   REAL fac,gset,rsq,v1,v2,ran1
   SAVE iset,gset
   DATA iset/0/   
   if (idum.lt.0) iset=0 
   if (iset.eq.0) then 
1      v1=2.*ran(idum)-1. 
       v2=2.*ran(idum)-1. 
       rsq=v1**2+v2**2 
       if(rsq.ge.1..or.rsq.eq.0.) goto 1 
       fac=sqrt(-2.*log(rsq)/rsq) 
       gset=v1*fac
       gasdev=v2*fac
       iset=1 
   else 
       gasdev=gset 
       iset=0 
   endif
   return
END FUNCTION gasdev

end module rand_number_mod
