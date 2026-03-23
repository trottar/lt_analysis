      subroutine eps_n_theta(pid,npol,Eb,ww,qq,tt,sin_theta_cm,eps)

c     Calculates sin(theta_cm) and epsilon.
c     Unphysical theta kinematics return sin_theta_cm = -1.0

      implicit none

      character*4 pid
      integer npol
      real Eb,ww,qq,tt,sin_theta_cm,eps
c     Keep the interface real to match the legacy fixed-form callers.

      double precision s,omega,q,tmin,tmax,raw,denom,tol
      double precision p1cm,p3cm,e1cm,e3cm,p1lab,p3cm_sq
      double precision sin_half_sq
      double precision Eb_d,ww_d,qq_d,tt_d,eps_d

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

      sin_theta_cm = -1.0
      eps = -1.0

      Eb_d = dble(Eb)
      ww_d = dble(ww)
      qq_d = dble(qq)
      tt_d = dble(tt)

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

      if (ww_d .le. 0.d0 .or. qq_d .lt. 0.d0 .or. tt_d .lt. 0.d0)
     *     return

      s     = ww_d*ww_d
      omega = (s + qq_d - m22)/(2.d0*m2)
      q     = sqrt(max(qq_d + omega*omega, 0.d0))
      m12   = -qq_d

      e1cm = (s + m12 - m22)/(2.d0*ww_d)
      e3cm = (s + m32 - m42)/(2.d0*ww_d)

      p1lab = q
      p1cm  = p1lab*m2/ww_d

      p3cm_sq = e3cm*e3cm - m32
      if (p1cm .le. 0.d0) return
      if (p3cm_sq .lt. -tol) return
      p3cm = sqrt(max(p3cm_sq, 0.d0))

      denom = 4.d0*p1cm*p3cm
      if (denom .le. tol) return

c     tt = -t, so tmin here is really (-t)_min
      tmin = -((e1cm - e3cm)**2 - (p1cm - p3cm)**2)
      tmax = tmin + denom

      raw = (tt_d - tmin)/denom

c     Reject truly unphysical points; only clip tiny numerical leakage
      if (raw .lt. -tol .or. raw .gt. 1.d0 + tol) then
         sin_theta_cm = -1.0
      else
         sin_half_sq = min(1.d0, max(0.d0, raw))
         sin_theta_cm = sngl(2.d0*sqrt(sin_half_sq*(1.d0-sin_half_sq)))
      endif

      eps_d = 1.d0 + 2.d0*(qq_d + omega*omega)
     *     /(4.d0*Eb_d*(Eb_d - omega) - qq_d)
      eps = sngl(1.d0/eps_d)

      return
      end
