*=======================================================================
* xmodel for K^+ for KaonLT 2018-19
*=======================================================================      
      subroutine xmodel(pid,npol_set,Eb,q2_set,eps_set,w,q2,tm,phi,
     *     eps_mod,th_mod,x_mod)

c     To calculate model cross-section, sigT+eps*sigL+ interfer._terms.

      implicit none

      character*2 prv_it
      common prv_it

      integer npol_set
      real Eb,q2_set,eps_set
      real w,q2,tm,phi
      real eps_mod,th_mod,x_mod

      real targ,mp,mn,pi
      real mpipl, mkpl
      parameter (pi=3.14159)      
      parameter (mp=.93827231)   !mp
      parameter (mn=.93956563)   !mn
      parameter (mpipl=0.139570)
      parameter (mkpl=0.493677)
      
      real wfactor
      real thetacm

      integer i

      real sig,sigT,sigL,sigLT,sigTT
      real dsig,dsigT,dsigL,dsigLT,dsigTT

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
*      print*, 'param: fn=',fn

      open(56,file=fn)
      do while(.true.)
         read(56,*,end=9) p,e,i
         par(i)=p
c     Print Statements         
*         write(6,101)par(i),e,i 
* 101     format(' xmodel: '2f11.4,i4)
c         pause
      end do
      

 9    close(56)

*     Model sigL, sigT, sigTT, sigLT.

* Revised for IT26, 12.11.09
*      tav=(0.0735+0.028*log(q2_set))*q2_set
*       RLT (10/8/2023): Testing new tav parameterization
      tav=(0.1112 + 0.0066*log(q2_set))*q2_set      
      f_tav=(abs(tm)-tav)/tav
      f_tm=abs(tm)/(abs(tm)+mkpl**2)**2 ! pole factor
      
*     sigL=(par(1)+par(2)*log(q2))*exp((par(3)
*     >     +par(4)*log(q2))*(abs(tm)-0.2))
*     RLT (10/12/2023): Removed 0.2 to keep things as simple as possible for
*                       initial start parameterization
      sigL=(par(1)+par(2)*log(q2))*exp((par(3)
     >     +par(4)*log(q2))*(abs(tm)))      
      sigT=par(5)+par(6)*log(q2)+(par(7)+par(8)*log(q2))*f_tav

      sigLT=(par(9)*exp(par(10)*abs(tm))+par(11)/abs(tm))*sin(thetacm)
      sigTT=(par(12)*q2*exp(-q2))*f_tm*sin(thetacm)**2

*     RLT (9/25/2023): There are two tav parameterizations in here.
*                      I am only using the one above, for now.
*      tav=(-0.178+0.315*log(q2))*q2
      
** !! MODEL DEP STUDY !!
c      sigL=sigL*0.90-0.1

      sig=sigT+eps_mod*sigL+eps_mod*cos(2.*phi)*sigTT
     >     +sqrt(2.0*eps_mod*(1.+eps_mod))*cos(phi)*sigLT

c     Correct for W.

      g_W=1./(W**2-targ**2)**2       ! W factor
      
      wfactor=g_W
      sigL=sigL*wfactor
      sigT=sigT*wfactor
      sigTT=sigTT*wfactor
      sigLT=sigLT*wfactor

      sig=sig/2./pi/1.d+06      !dsig/dtdphicm in microbarns/MeV**2/rad

      x_mod=sig      
      
      th_mod=thetacm

***
*     RLT (9/25/2023): Temporary errors!!!!
***
      dsigL=sqrt(abs(sigL))/abs(sigL)
      dsigT=sqrt(abs(sigT))/abs(sigT)
      dsigTT=sqrt(abs(sigTT))/abs(sigTT)
      dsigLT=sqrt(abs(sigLT))/abs(sigLT)
            
      write(71,60) sigL,dsigL,sigT,dsigT,sigTT,
     *     dsigTT,sigLT,dsigLT,q2,tm
 60   format(8G15.5,2f8.5)
      
      end

*=======================================================================
