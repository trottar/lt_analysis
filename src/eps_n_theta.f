      subroutine eps_n_theta(pid,npol,Eb,ww,qq,tt,theta_cm,eps)

c     To calculate model theta_pq in CM and epsilon. This subroutine is largely
c     based on theta_cm.f function, which in turn is based Jochen's script.

      implicit none

      character*4 pid      
      
      integer npol
      real Eb,ww,qq,tt,theta_cm,eps

      REAL s,omega,q,tmin
      REAL p1cm,p3cm,e1cm,e3cm,p1lab

      REAL m2,m3,m4
      REAL m12,m22,m32,m42

      real mp,mp2,mpi,mpi2,mn,mn2,mlamb,mlamb2,mK,mK2
      parameter (mp=.93827231)   !mp
      parameter (mp2=.88035493)  !mp^2
      parameter (mpi=.13956995)   !mpi
      parameter (mpi2=.01947977)  !mpi^2
      parameter (mn=.93956563)   !mn
      parameter (mn2=.88278357) !mn^2
      parameter (mlamb=1.115683)   !mlamb
      parameter (mlamb2=1.24474855649) !mlamb^2      
      parameter (mK=0.493677)   !mK
      parameter (mK2=0.24387)   !mK2    

      ! Check particle type and set parameters accordingly
      if (pid == "kaon") then
         m3=mK
         m32=mK2
      else if (pid == "pion") then
         m3=mpi
         m32=mpi2
      else
        print *, "*** Invalid particle type: ",pid
      endif

      if(npol.gt.0) then
         m2=mp
         m22=mp2
         m4=mlamb ! Lambda Mass
         m42=mlamb2
      else
         m2=mn
         m22=mn2
         m4=mp
         m42=mp2
      end if

      s=ww*ww
      omega=(s+qq-m22)/(2*m2)
      q=sqrt(qq+omega**2)
*     m12=qq    !error?
      m12=-qq   !mass squared of virtual photon.

      e1cm=(s+m12-m22)/(2*ww)
      e3cm=(s+m32-m42)/(2*ww)
      p1lab=q
      p1cm=p1lab*m2/ww
      p3cm=sqrt(e3cm*e3cm-m32)
      tmin=-((e1cm-e3cm)**2-(p1cm-p3cm)**2)

      if (tt.ge.tmin) then
         theta_cm=2*asin(sqrt((tt-tmin)/(4*p1cm*p3cm)))
      else
         theta_cm=-1.
         print*, 'eps_n_theta: tt=',tt,' <  tmin=',tmin
      endif
      
      eps=1.+2.*(qq+omega**2)/(4.*Eb*(Eb-omega)-qq)
      eps=1./eps

c      write(*,'(a13,7(F8.5,1x))')
c     *     'eps_n_theta: ',ww,qq,t,tmin,theta_cm,eps,omega

      end
