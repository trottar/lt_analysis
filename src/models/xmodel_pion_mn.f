*=======================================================================
* xmodel for pi- of FPI2
*=======================================================================      
      subroutine xmodel(pid,npol_set,Eb,q2_set,w,q2,tm,phi,
     *     eps_mod,th_mod,x_mod)

c     To calculate model cross-section, sigT+eps*sigL+ interfer._terms.

      implicit none

      character*2 prv_it
      common prv_it

      integer npol_set
      real Eb,q2_set,w,q2,tm,phi,eps_mod,th_mod,x_mod

      real targ,mp,mn
      parameter (mp=.93827231)   !mp
      parameter (mn=.93956563)   !mn
      parameter (mpipl=0.139570)
      parameter (mkpl=0.493677)
      
      real wfactor
      real thetacm

      integer i

      real sigT,sigL,sigLT,sigTT

      character*80 fn
      character*2 pol
      character*4 pid

      real par(12)
      real p,e
      real f_tm,g_W,tav,f_tav

      if(npol_set.lt.0) then
         pol='mn'
         targ=mn
      else
         pol='pl'
         targ=mp
      end if

*     Calculate model thetacm and epsilon at first.
      call eps_n_theta(pid,npol_set,Eb,w,q2,tm,thetacm,eps_mod)

*     Model fit parameters.

      write(fn,10) pid,pol,nint(q2_set*10)
 10   format(a4,'/parameters/par.',a2,'_',i2.2,'.dat')
      if (phi.lt.0.3) then
         print*, 'xmodel: fn=',fn
      endif

      open(56,file=fn)
      do while(.true.)
         read(56,*,end=9) p,e,i
         par(i)=p
         if (phi.lt.0.3) then
            write(6,101)par(i),e,i
 101        format(' xmodel: '2f11.4,i4)
         endif
c         pause
      end do

 9    close(56)

*     Model sigL, sigT, sigTT, sigLT.

      f_tm=abs(tm)/(abs(tm)+mpipl**2)**2 ! pole factor

      g_W=1./(W**2-targ**2)**2       ! W factor

* Parameterization based upon Fpi-1 pi+ IT25, 12.04.18
c      tav=(0.0735+0.028*log(q2))*q2
c      f_tav=(tm-tav)/tav
c
c      sigL=(par(1)+par(2)*log(q2))*exp((par(3)+par(4)*log(q2))*abs(tm))
c      sigT=par(5)+par(6)*log(q2)+(par(7)+par(8)*log(q2))*f_tav
c
c      sigLT=par(9)*exp(par(10)*abs(tm))*sin(thetacm)
c      sigTT=(par(11)*q2*exp(-q2)+par(12)/q2**2)*f_tm*sin(thetacm)**2

* Revised for IT22, 12.11.06
c      tav=(0.0735+0.028*log(q2_set))*q2_set
c      f_tav=(tm-tav)/tav
c
c      sigL=(par(1)+par(2)*log(q2))*exp((par(3)+par(4)*log(q2))*(abs(tm)-0.2))
c      sigT=par(5)+par(6)*log(q2)+(par(7)+par(8)*log(q2))*f_tav
c
c      sigLT=par(9)*exp(par(10)*abs(tm)+1.5/abs(tm))*sin(thetacm)
c      sigTT=(par(11)*q2*exp(-q2)+par(12)/q2**2)*f_tm*sin(thetacm)**2

* Revised for IT26, 12.11.09
      tav=(0.0735+0.028*log(q2_set))*q2_set
      f_tav=(tm-tav)/tav

      sigL=(par(1)+par(2)*log(q2))*exp((par(3)
     >     +par(4)*log(q2))*(abs(tm)-0.2))
      sigT=par(5)+par(6)*log(q2)+(par(7)+par(8)*log(q2))*f_tav

      sigLT=(par(9)*exp(par(10)*abs(tm))+par(11)/abs(tm))*sin(thetacm)
      sigTT=(par(12)*q2*exp(-q2))*f_tm*sin(thetacm)**2

** !! MODEL DEP STUDY !!
c      sigL=sigL*0.90-0.1

      x_mod=sigT+eps_mod*sigL+eps_mod*cos(2.*phi)*sigTT
     >     +sqrt(2.0*eps_mod*(1.+eps_mod))*cos(phi)*sigLT

c     Correct for W.

      wfactor=g_W
      sigL=sigL*wfactor
      sigT=sigT*wfactor
      sigTT=sigTT*wfactor
      sigLT=sigLT*wfactor
      x_mod=x_mod*wfactor

      th_mod=thetacm

      if (phi.lt.0.3) then
         write(6,102) eps_mod,tm,sigL,sigT,sigTT,sigLT,x_mod
 102     format('xmodel: eps=',f5.3,' t=',f5.3,' sigL=',f6.2,' sigT=',
     1        f6.2,' sigTT=',f5.2,' sigLT=',f5.2,' x_mod=',f5.2)
      endif

      end

*=======================================================================
