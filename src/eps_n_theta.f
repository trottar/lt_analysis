      subroutine eps_n_sin_theta(pid,npol,Eb,ww,qq,tt,sin_theta_cm,eps)

c     Calculates sin(theta_cm) and epsilon.
c     Unphysical theta kinematics return sin_theta_cm = -1.d0

      implicit none

      character*4 pid
      integer npol
      double precision Eb,ww,qq,tt,sin_theta_cm,eps

      double precision s,omega,q,tmin,tmax,raw,denom,tol
      double precision p1cm,p3cm,e1cm,e3cm,p1lab,p3cm_sq
      double precision sin_half_sq

      double precision m2,m3,m4
      double precision m12,m22,m32,m42

      double precision mp,mp2,mpi,mpi2,mn,mn2,mlamb,mlamb2,mK,mK2
      parameter (mp    = 0.93827231d0)
      parameter (mp2   = 0.88035493d0)
      parameter (mpi   = 0.13956995d0)
      parameter (mpi2  = 0.01947977d0)
      parameter (mn    = 0.93956563d0)
      parameter (mn2   = 0.88278357d0)
      parameter (mlamb = 1.115683d0)
      parameter (mlamb2= 1.244749d0)
      parameter (mK    = 0.493677d0)
      parameter (mK2   = 0.24387d0)
      parameter (tol   = 1.0d-10)

      sin_theta_cm = -1.d0
      eps = -1.d0

c     Check particle type and set parameters accordingly
      if (pid .eq. 'kaon') then
         m3  = mK
         m32 = mK2
      else if (pid .eq. 'pion') then
         m3  = mpi
         m32 = mpi2
      else
         write(*,*) '*** Invalid particle type: ', pid
         return
      endif

      if (npol .gt. 0) then
         m2  = mp
         m22 = mp2
         m4  = mlamb
         m42 = mlamb2
      else
         m2  = mn
         m22 = mn2
         m4  = mp
         m42 = mp2
      endif

      if (ww .le. 0.d0 .or. qq .lt. 0.d0 .or. tt .lt. 0.d0) return

      s     = ww*ww
      omega = (s + qq - m22)/(2.d0*m2)
      q     = sqrt(max(qq + omega*omega, 0.d0))
      m12   = -qq

      e1cm = (s + m12 - m22)/(2.d0*ww)
      e3cm = (s + m32 - m42)/(2.d0*ww)

      p1lab = q
      p1cm  = p1lab*m2/ww

      p3cm_sq = e3cm*e3cm - m32
      if (p1cm .le. 0.d0) return
      if (p3cm_sq .lt. -tol) return
      p3cm = sqrt(max(p3cm_sq, 0.d0))

      denom = 4.d0*p1cm*p3cm
      if (denom .le. tol) return

c     tt = -t, so tmin here is really (-t)_min
      tmin = -((e1cm - e3cm)**2 - (p1cm - p3cm)**2)
      tmax = tmin + denom

      raw = (tt - tmin)/denom

c     Reject truly unphysical points; only clip tiny numerical leakage
      if (raw .lt. -tol .or. raw .gt. 1.d0 + tol) then
         sin_theta_cm = -1.d0
      else
         sin_half_sq = min(1.d0, max(0.d0, raw))
         sin_theta_cm = 2.d0*sqrt(sin_half_sq*(1.d0 - sin_half_sq))
      endif

      eps = 1.d0 + 2.d0*(qq + omega*omega)/(4.d0*Eb*(Eb - omega) - qq)
      eps = 1.d0/eps

      return
      end