*=======================================================================
* xmodel for K^+ for KaonLT 2018-19
*=======================================================================      
      subroutine xmodel(pid,npol_set,Eb,q2_set,w_set,eps_set,
     *     ww,qq,tt,phi,eps_mod,th_mod,x_mod,par_fn)

c     To calculate model cross-section, sig_T+eps*sig_L+ interfer._terms.

      implicit none

      character*2 prv_it
      common prv_it

      integer npol_set
      real Eb,q2_set,w_set,eps_set
      real ww,qq,tt,phi
      real eps_mod,th_mod,x_mod

      real mtar,mp,mn,pi
      real mpipl, mkpl
      parameter (pi=3.14159)
      parameter (mp=.93827231)   !mp
      parameter (mn=.93956563)   !mn
      parameter (mpipl=0.139570)
      parameter (mkpl=0.493677)
      
      real wfactor
      real theta_cm

      integer i, test

      real sig,sig_T,sig_L,sig_LT,sig_TT
      real dsig,dsig_T,dsig_L,dsig_LT,dsig_TT

      character*80 par_fn
      character*2 pol
      character*4 pid

*     RLT (1/2/2024): Need to have 16 parameters (4 for L/T/LT/TT) for
*                     the xfit_in_t.py script to work. LT/TT are zeros
      real par(16)
      real p,e
      real ft,tav,ftav

*     RLT (7/11/2024): Redefined functional forms of L, T, LT, TT
*                      that incorporates Q2-dep based of pi FF
      real Qdep_L, Qdep_T, Qdep_LT, Qdep_TT

      if(npol_set.lt.0) then
         pol='mn'
      else
         pol='pl'
      end if

*     Dynamically updated
      mtar=mp

      open(57, file=par_fn)
      do while (.true.)
         read(57, *, end=9) p, e, i
         par(i) = p
!     Print Statements
*         print *,"param: ", i, p, e
! You can customize the format as needed
      end do      
 9    close(57)      
      
*     Calculate model theta_cm and epsilon at first.
      call eps_n_theta(pid,npol_set,Eb,ww,qq,tt,
     *     theta_cm,eps_mod)

*     Model sig_L, sig_T, sig_TT, sig_LT.

*** Make sure there is a space after every equation ***
      
      tav=(0.1112 + 0.0066*log(q2_set))*q2_set
      
      ftav=(abs(tt)-tav)/tav

      ft=abs(tt)/(abs(tt)+mkpl**2)**2 ! pole factor
      
      Qdep_L=qq/(1.0+(1.77*qq)+0.12*(qq**2))
      
      sig_L=(par(1)*Qdep_L*ft)*exp(-par(2)*(abs(tt)))
      
      Qdep_T=(exp(-qq**2))/qq
      
      sig_T=(par(5)*exp(-par(6)*(abs(tt)))+par(7)*(abs(tt)))
     >     *(Qdep_T**par(8))
      
      sig_LT=(par(9)*exp(par(10)*abs(tt))+par(11)/abs(tt))*sin(theta_cm)
      
      sig_TT=((-par(13)*abs(tt)+par(14))*(abs(tt)
     >     **(qq/par(15)))-par(16)*qq)*sin(theta_cm)**2
      
c     Correct for W.
      wfactor=1./(ww**2-mtar**2)**2  ! W factor

      sig_L=sig_L*wfactor
      sig_T=sig_T*wfactor
      sig_TT=sig_TT*wfactor
      sig_LT=sig_LT*wfactor

      if (isnan(sig_L).or.(sig_L > 1.d+06)) sig_L = 0.0
      if (isnan(sig_T).or.(sig_T > 1.d+06)) sig_T = 0.0
      if (isnan(sig_LT).or.(sig_LT > 1.d+06)) sig_LT = 0.0
      if (isnan(sig_TT).or.(sig_TT > 1.d+06)) sig_TT = 0.0

      print *,"--------------"
      WRITE(*,*) 'sig_L=', sig_L
      WRITE(*,*) 'sig_T=', sig_T
      WRITE(*,*) 'sig_LT=', sig_LT
      WRITE(*,*) 'sig_TT=', sig_TT
      
      sig=sig_T+eps_mod*sig_L+eps_mod*cos(2.*phi)*sig_TT
     >     +sqrt(2.0*eps_mod*(1.+eps_mod))*cos(phi)*sig_LT

      sig=sig/2./pi/1.d+06      !dsig/dtdphicm in microbarns/MeV**2/rad
*      sig=sig/2./pi             !dsig/dtdphicm in microbarns/GeV**2/rad      

      x_mod=sig
      
      th_mod=theta_cm
      
      end

*=======================================================================
