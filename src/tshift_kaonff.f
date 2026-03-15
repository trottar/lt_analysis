      program calc_tshift

      implicit none

      real*8 inp_q2, inp_w, inp_theta_cm, inp_mmshift, inp_ebeam

      write(*,*) 'Inputing Q2, W, theta_pq_cm, MM shift,'
      write(*,*) 'and beam energy:'
      read(*,*) inp_q2, inp_w, inp_theta_cm, inp_mmshift, inp_ebeam

      write(*,*) 'Q2 = ', inp_q2, ' W = ', inp_w,
     *     ' theta_cm = ', inp_theta_cm,
     *     ' MMshift(MeV) = ', inp_mmshift,
     *     ' Ebeam(MeV) = ', inp_ebeam

      call calc_tshift_kaon(inp_q2, inp_w, inp_theta_cm, inp_mmshift,
     *     inp_ebeam)

      print*, '-------------------------------------------------'
      stop
      end

*=======================================================================

      subroutine calc_tshift_kaon(q2g, wg, theta_cm_signed, mmshift,
     *     tinc)

      implicit none

      real*8 mkp,mL,mp,mpg,me,pi
      real*8 nu,nug,pgam,pgamg,pgamx,pgamy,pgamz
      real*8 einc,escat,pinc,pincx,pincz,pscat,thscat,pscatx,pscatz
      real*8 thkpcm,ekp0,ekpcm0,pkp0,pkpcm0,pkpx0,pkpy0,pkpz0
      real*8 ekp1,ekpcm1,pkp1,pkpcm1,pkpx1,pkpy1,pkpz1
      real*8 mmshift, tshift, tshiftg
      real*8 phix, thkplab0, thkplab1
      real*8 betacm, gammacm, pgamtemp, pkptemp
      real*8 t0, t1, w, ang, q2g, wg, tinc, theta_cm_signed
      real*8 thgam

c     stuff for linux_suppl.inc
      real*8 sind, cosd, tand
      real*8 asind, acosd, atand, atan2d
      external sind, cosd, tand
      external asind, acosd, atand, atan2d

      me  = 0.511d0
      mkp = 493.677d0
      mp  = 938.27d0
      mpg = mp/1000.d0
      mL  = 1115.683d0
      pi  = 3.14156d0

c     The analysis convention uses signed theta for left/right,
c     but the t-channel calculation is always forward.
c     the magnitude of theta_cm only.
      thkpcm = abs(theta_cm_signed)
      phix = 0.d0
      w = wg * 1000.d0

      nug = (wg**2 + q2g - mpg**2)/(2.d0*mpg)
      nu  = nug * 1000.d0

      pgamg = sqrt(nug**2 + q2g)
      pgam  = pgamg * 1000.d0

c     find kaon lab momentum
c     first, find speed of virtual photon+proton c.m. frame
      betacm  = pgam/(nu+mp)
      gammacm = (nu+mp)/w

      ekpcm0 = (w**2 + mkp**2 - (mL**2)) / (2.d0*w)
      ekpcm1 = (w**2 + mkp**2 - ((mL+mmshift)**2)) / (2.d0*w)
      pkpcm0 = sqrt(ekpcm0**2 - mkp**2)
      pkpcm1 = sqrt(ekpcm1**2 - mkp**2)

c     now transform to lab frame wrt q-vector
      ekp0 = gammacm * (betacm*pkpcm0*cosd(thkpcm) + ekpcm0)
      ekp1 = gammacm * (betacm*pkpcm1*cosd(thkpcm) + ekpcm1)
      pkp0 = sqrt(ekp0**2 - mkp**2)
      pkp1 = sqrt(ekp1**2 - mkp**2)
      thkplab0 = asind(sind(thkpcm)*pkpcm0/pkp0)
      thkplab1 = asind(sind(thkpcm)*pkpcm1/pkp1)

      einc = tinc + me
      pinc = sqrt(einc**2 - me**2)
      pincz = pinc
      pincx = 0.d0

      escat = einc - nu
      if (escat.lt.me) then
         write(*,*) '*** Invalid scattered electron energy'
         stop 1
      endif
      pscat = sqrt(escat**2 - me**2)

      ang = (pinc**2 + pscat**2 - pgam**2) / (2.d0*pinc*pscat)
      if (ang.lt.-1.d0 .or. ang.gt.1.d0) then
         write(*,*) '*** Invalid scattering angle in tshift_kaonff'
         stop 1
      endif
      thscat = acosd(ang)
      pscatz = pscat*cosd(thscat)
      pscatx = pscat*sind(thscat)

      pgamz = pincz - pscatz
      pgamx = pincx - pscatx
      pgamy = 0.d0

      pgamtemp = sqrt(pgamx**2 + pgamy**2 + pgamz**2)
      if (abs(pgamtemp-pgam).gt.0.2d0) write(*,150) pgam, pgamtemp
 150  format(' PGam disagreement ',2f10.3)
      thgam = atand(pgamx/pgamz)

c     now rotate to lab frame wrt e- beam
      pkpx0 = pkp0*(cosd(thkplab0)*sind(thgam) -
     *     sind(thkplab0)*cosd(thgam)*cosd(phix))
      pkpy0 = pkp0*sind(thkplab0)*sind(phix)
      pkpz0 = pkp0*(cosd(thkplab0)*cosd(thgam) +
     *     sind(thkplab0)*sind(thgam)*cosd(phix))
      pkpx1 = pkp1*(cosd(thkplab1)*sind(thgam) -
     *     sind(thkplab1)*cosd(thgam)*cosd(phix))
      pkpy1 = pkp1*sind(thkplab1)*sind(phix)
      pkpz1 = pkp1*(cosd(thkplab1)*cosd(thgam) +
     *     sind(thkplab1)*sind(thgam)*cosd(phix))

      pkptemp = sqrt(pkpx0**2 + pkpy0**2 + pkpz0**2)
      if (abs(pkptemp-pkp0).gt.0.2d0) write(*,160) pkp0, pkptemp
 160  format(' Pkp disagreement ',2f10.2)

c     equation is actually for -t
      t0 = (pgamx-pkpx0)**2 + (pgamy-pkpy0)**2 +
     *     (pgamz-pkpz0)**2 - (nu-ekp0)**2
      t1 = (pgamx-pkpx1)**2 + (pgamy-pkpy1)**2 +
     *     (pgamz-pkpz1)**2 - (nu-ekp1)**2

      tshift = t1 - t0
      tshiftg = tshift/1.d6

      write(*,170) mmshift, tshiftg
 170  format(' MM shift of ',f10.4,' MeV gives ',f12.6,' GeV^2 shift')
      write(*,171) tshiftg
 171  format(' TSHIFT_GEV2 = ',es16.8)

      return
      end

      include 'linux_suppl.inc'
