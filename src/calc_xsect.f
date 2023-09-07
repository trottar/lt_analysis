      program calc_xsect


      implicit none

c This script computes the experimental cross section using the ratio
c DATA/MC(Yields * CS(MC)

c      character*2 prv_it
c      common prv_it

c      integer q2_bin
c      integer t_bin, phi_bin
c      common t_bin, phi_bin
      
c     Get number of the previous iteration.
      
c      if(iargc().ne.1) then
c         print*,'*** usage: calc_xsect prv_it ***'
c         stop
c      end if
c      call getarg(1,prv_it)

c     Calculate unseparated cross-sections. Now settings are for the piplus data (+)

      character*4 inp_pid
      integer inp_pol
      real inp_Q2, inp_loeps, inp_hieps      
      
      write(*,*) "Inputing particle, polarity, Q2 and both epsilons:"
      read(*,*) inp_pid, inp_pol, inp_Q2, inp_loeps, inp_hieps

      write(*,*) "PID = ",inp_pid,"POL = ",inp_pol,"Q2 = ",inp_Q2,
     *           "low_eps = ",inp_loeps,"high_eps = ",inp_hieps

      call xsect(inp_pid, inp_pol, inp_Q2, inp_loeps)
      call xsect(inp_pid, inp_pol, inp_Q2, inp_hieps)
      
      print*,  "-------------------------------------------------"
      
      stop
      end

*=======================================================================

      subroutine xsect(pid,npol_set,q2_set,eps_set)

      implicit none

      integer npol_set
      real q2_set,eps_set

      integer kset,nbin

      character*80 r_fn, kin_fn, xunsep_fn
      character*2 pol
      character*4 pid      

      integer it,ip
      real Eb,eps

      real one

      integer nt,nphi
      parameter (nt=6,nphi=16)

      real r,dr,w,dw,q2,dq2,th_only,th_cm
      real tm,tmn,tmx
      real eps_mod,th_mod,x_mod
      real x_real,dx_real

      integer ipol
      real th_pq

      real phi

      include './kin_xyz.inc'

      ipol=0
      q2=0.
      eps=0.
      tmn=0.
      tmx=0.
      kset=0
      open(55,file=trim(pid) // '/list.settings')
      do while(ipol.ne.npol_set.or.q2.ne.q2_set.or.eps.ne.eps_set)
         read(55,*) ipol,q2,eps,th_pq,tmn,tmx,nbin
c         write(6,2)ipol,q2,eps,th_pq,tmn,tmx,nbin
c 2       format(i5,5f10.5,2i5)
      end do
      close(55)
      write(6,3)tmn,tmx,kset
 3    format(' tmn, tmx: ',2f10.5,'  for kset=',i5)
      if(tmn.eq.0..or.tmx.eq.0.) 
     *     stop '*** setting is not found in list.settings'

      Eb=Ebeam(kset)

      if(npol_set.lt.0) then
         pol='mn'
      else
         pol='pl'
      end if

      write(6,4)Eb,q2,eps,pol
 4    format(' xsect: Eb=',f8.5,'   at Q2=',f7.4,
     *     '  eps=',f6.4,'  pol=',a2)

c     construct ratio data file name.

      write(r_fn,10) pid,pol,nint(q2*100),nint(eps*100)
 10   format(a4,'/averages/aver.',a2,'_',i3.3,'_',i2,'.dat')
      print*,'xsect: r_fn=',r_fn

      open(51,file=r_fn)

c     construct kinematics data file name.

      write(kin_fn,20) pid,nint(q2*100)
 20   format(a4,'/averages/avek.',i3.3,'.dat')
      print*,'xsect: kin_fn=',kin_fn

      open(52,file=kin_fn)

*     construct output file name.
      write(xunsep_fn,30) pid,pol,nint(q2_set*100),nint(eps_set*100)
 30   format(a4,'/xsects/x_unsep.',a2,'_',i3.3,'_',i2)
      print*,'xsect: xunsep_fn=',xunsep_fn
c      pause

      open(61,file=xunsep_fn,status='replace')

      do it=1,nbin

         tm=tmn+(it-0.5)*(tmx-tmn)/nbin
         read(52,*) w,dw,q2,dq2,th_cm
         write(6,32) w,dw,q2,dq2,th_cm
 32      format('xsect: ',6f10.4)

         th_cm=th_cm*3.14159D0/180.D0
         
         do ip=1,nphi

            phi=(ip-0.5)*2.*3.14159/nphi
            read(51,*) r,dr

            call xmodel(npol_set,Eb,q2_set,w,q2,tm,phi,
     *           eps_mod,th_mod,x_mod)

c angle check
            if (abs(th_mod-th_cm).gt.1.e-4) then
               write(6,*)' Angle error ',th_mod,th_cm
               stop
            endif

c ratio is data/simc - see GH logbook, p.55
            x_real=x_mod*r
            dx_real=x_mod*dr
c ratio is simc/data - see GH logbook, p.55
c            if (r.ne.0) then
c               x_real=x_mod/r
c               dx_real=x_mod*dr/r**2
c            else
c               x_real=0.
c               dx_real=0.
c            endif

            write(61,40) x_real,dx_real,x_mod,eps_mod,
     *           th_mod*180./3.14159,phi*180./3.14159,tm,w,q2
 40         format(3G15.5,f8.5,2f7.2,3f8.5)

         end do                 !phi

c        Write out kinematics for Henk.
         if(npol_set.gt.0) write(99,'(5f8.3,2x,2f6.2)')
     *   w,q2,eps_mod,th_mod*180./3.14159,tm,eps_set,q2_set

      end do                    !t

      close(51)
      close(52)
      close(61)
      print*,' '

      end

*=======================================================================

      include 'xmodel_pi_mn_fpi2.f'
      
*=======================================================================

      include 'eps_n_theta.f'
