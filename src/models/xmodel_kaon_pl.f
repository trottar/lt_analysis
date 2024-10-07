*=======================================================================
* xmodel for K^+ for KaonLT 2018-19
*=======================================================================      
      subroutine xmodel(pid,npol_set,Eb,q2_set,w_set,eps_set,
     *     w,q2,tm,phi,eps_mod,th_mod,x_mod,par_fn)

c     To calculate model cross-section, sigT+eps*sigL+ interfer._terms.

      implicit none

      character*2 prv_it
      common prv_it

      integer npol_set
      real Eb,q2_set,w_set,eps_set
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

      integer i, test

      real sig,sigT,sigL,sigLT,sigTT
      real dsig,dsigT,dsigL,dsigLT,dsigTT

      character*80 par_fn
      character*2 pol
      character*4 pid

*     RLT (1/2/2024): Need to have 16 parameters (4 for L/T/LT/TT) for
*                     the xfit_in_t.py script to work. LT/TT are zeros
      real par(16)
      real p,e
      real f_tm,g_W,tav,f_tav

*     RLT (7/11/2024): Redefined functional forms of L, T, LT, TT
*                      that incorporates Q2-dep based of pi FF
      real Qdep_L, Qdep_T

      if(npol_set.lt.0) then
         pol='mn'
         targ=mn
      else
         pol='pl'
         targ=mp
      end if

      open(57, file=par_fn)
      do while (.true.)
         read(57, *, end=9) p, e, i
         par(i) = p
!     Print Statements
*         print *,"param: ", i, p, e
! You can customize the format as needed
      end do      
 9    close(57)      
      
*     Calculate model thetacm and epsilon at first.
      call eps_n_theta(pid,npol_set,Eb,w,q2,tm,
     *     thetacm,eps_mod)

*     Model sigL, sigT, sigTT, sigLT.

* Revised for IT26, 12.11.09
*      tav=(0.0735+0.028*log(q2_set))*q2_set
*       RLT (10/8/2023): Testing new tav parameterization
      tav=(0.1112 + 0.0066*log(q2_set))*q2_set      
      f_tav=(abs(tm)-tav)/tav
*     RLT (7/11/2024): Moved below for Q2dep func form
*      f_tm=abs(tm)/(abs(tm)+mkpl**2)**2 ! pole factor
      
*     sigL=(par(1)+par(2)*log(q2))*exp((par(3)
*     >     +par(4)*log(q2))*(abs(tm)-0.2))
*     RLT (10/12/2023): Removed 0.2 to keep things as simple as possible for
*                       initial start parameterization
*     RLT (2/19/2024): Adding a 0.2 term to t dependence to bring down the
*                      extreme slope at high t
*     RLT (3/09/2024): Removing +0.2 term for better parameterization of
*                      Q2=3.0, W=2.32
*      
*      sigL=(par(1)+par(2)*log(q2))*exp((par(3)
*     >     +par(4)*log(q2))*(abs(tm)))
*      sigL=(par(1)+par(2)*log(q2))*exp((par(3)
*     >     +par(4)*log(q2))*(abs(tm)+0.2))
*     RLT (4/23/2024): Marco's thesis functional forms
*     sigL=par(1)*exp(-par(2)*abs(tm))*(1.0/(1.0+par(3)*q2))
*     RLT (6/04/2024): Testing simplier exp form for L+T
**
**      sigL=(par(1)+par(2)*log(q2))*exp(par(3)*(abs(tm)))
*      sigL=(par(1)*((abs(tm)/q2)-1))*exp(par(2)*(abs(tm)))      
      
*     RLT (2/15/2024): Removing t dependence from sigT because it seems
*                        to be driving poor sep xsects results
*     RLT (2/20/2024): Added 1/Q^4 term to dampen sigT
*     RLT (2/21/2024): Using global analysis sig T model and params
*     (https://journals.aps.org/prc/pdf/10.1103/PhysRevC.85.018202)
*      
*      sigT=par(5)+par(6)*log(q2)+(par(7)+par(8)*log(q2))*f_tav
*     sigT=par(5)+par(6)*log(q2)
*     sigT=par(5)*log(q2)++par(6)/(q2**2)
*      sigT=par(5)/(1+par(6)*q2)
*     RLT (4/20/2024): Adding in t-dependence
*     sigT=(par(5)/(1+par(6)*q2))*f_tav
*     sigT=(par(5)/(1+par(6)*q2))*abs(tm)
*     RLT (4/20/2024): Exponential t-dependence
*      sigT=(par(5)/(1+par(6)*q2))*exp(par(7)*abs(tm))
*     RLT (4/23/2024): Marco's thesis functional forms
*     sigT=par(5)*exp(-par(6)*abs(tm))*(1.0/(1.0+par(7)*q2))
*     RLT (6/04/2024): Testing simplier exp form for L+T
**      
**      sigT=(par(5)*((abs(tm)/q2)-1))*exp(par(6)*(abs(tm)))
*      sigT=(par(5)+par(6)*log(q2))*exp(par(7)*(abs(tm)))            

**
      sigLT=(par(9)*exp(par(10)*abs(tm))+par(11)/abs(tm))*sin(thetacm)
*     sigLT=(par(9)+par(11)/abs(tm))*sin(thetacm)
*     RLT (4/23/2024): Marco's thesis functional forms
*      sigLT=par(9)*exp(-par(10)*abs(tm))*(1.0/(1.0+(q2**2)*par(11)))
*     RLT (1/2/2024): Need to have 16 parameters (4 for L/T/LT/TT) for
*                     the xfit_in_t.py script to work. LT/TT are zeros
*                     Therefore param 12 was also changed to 13
**      
      sigTT=(par(13)*q2*exp(-q2))*f_tm*sin(thetacm)**2
*     RLT (4/23/2024): Marco's thesis functional forms
*      sigTT=par(13)*exp(-par(14)*abs(tm))*(1.0/(1.0+(q2**2)*par(15)))

*     RLT (9/25/2023): There are two tav parameterizations in here.
*                      I am only using the one above, for now.
*      tav=(-0.178+0.315*log(q2))*q2

**      
*     RLT (7/11/2024): Redefined functional forms of L, T, LT, TT
*                      that incorporates Q2-dep based of pi FF
      f_tm=abs(tm)/(abs(tm)+mkpl**2)**2 ! pole factor
      Qdep_L=q2/(1.0+(1.77*q2)+0.12*(q2**2))
**      sigL=(par(1)*Qdep_L*f_tm)*exp(-par(2)*(abs(tm)))
      sigL=(par(1)*Qdep_L*f_tm)*exp(-par(2)*(abs(tm)))
*     sigT=(par(5)/q2)*exp(-par(6)*(q2**2))
      Qdep_T=(exp(-q2**2))/q2
**      sigT=par(5)*(par(6)+exp(-par(7)*(abs(tm))))*(Qdep_T**par(8))
      sigT=(par(5)*exp(-par(6)*(abs(tm)))+par(7)*(abs(tm)))
     >     *(Qdep_T**par(8))
***      sigLT=(par(9)/(1+q2))*sin(thetacm)
***     >     *exp(-par(10)*(abs(tm)))
*      sigLT=(par(9)*exp(par(10)*abs(tm))
*     >     +par(11)/abs(tm))*sin(thetacm)
*      sigTT=(-par(13)/(1+q2))*(sin(thetacm)**2)
*     >     *exp(-par(14)*(abs(tm)))
***      sigTT=(par(13)/(1+q2))*(sin(thetacm)**2)
***   >     *f_tm*exp(-par(14)*(q2))
****      sigTT=((-par(13)*abs(tm)+par(14))*(abs(tm)
****     >     **(q2/par(15)))-par(16)*q2)*sin(thetacm)**2
      
c     Correct for W.
      g_W=1./(W**2-targ**2)**2  ! W factor
*      g_W=1./(W**2-targ**2)**2.45 ! Q2=5.5,W=3.02, Q2=2.115, W=2.95
*      g_W=1./(W**2-targ**2)**2.65 ! W factor, Q2=4.4, W=2.74      
*     g_W=1./(W**2-targ**2)**2.25       ! W factor, Q2=3.0, W=3.14
*      g_W=1./(W**2-targ**2)**3.4       ! W factor, Q2=3.0, W=2.32

      wfactor=g_W
      sigL=sigL*wfactor
      sigT=sigT*wfactor
      sigTT=sigTT*wfactor
      sigLT=sigLT*wfactor

** !! MODEL DEP STUDY !!
c      sigL=sigL*0.90-0.1
*     RLT (2/16/2024): Moved the unsep xsect down here so that the weight
*                      factor is included in calculation
      sig=sigT+eps_mod*sigL+eps_mod*cos(2.*phi)*sigTT
     >     +sqrt(2.0*eps_mod*(1.+eps_mod))*cos(phi)*sigLT

      sig=sig/2./pi/1.d+06      !dsig/dtdphicm in microbarns/MeV**2/rad
*     sig=sig/2./pi      !dsig/dtdphicm in microbarns/GeV**2/rad      

      x_mod=sig
      
      th_mod=thetacm
      
      end

*=======================================================================
