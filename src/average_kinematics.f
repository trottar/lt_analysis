      program take_averages


      implicit none

c     This program takes averages of W and Q2 over different theta_pq
c     settings, then different polarities.
c
c     Input:  kindata/kindata.*.dat
c     Output: averages/averages.*.dat

      character*4 inp_pid
      integer inp_pol
      real inp_Q2, inp_loeps, inp_hieps
      write(*,*) "Inputing particle, polarity, Q2 and both epsilons:"
      read(*,*) inp_pid, inp_pol, inp_Q2, inp_loeps, inp_hieps

      write(*,*) "PID = ",inp_pid,"POL = ",inp_pol,"Q2 = ",inp_Q2,
     *           "low_eps = ",inp_loeps,"high_eps = ",inp_hieps
      
      call average_k(inp_pid, inp_pol, inp_Q2, inp_loeps, inp_hieps)
      print*,  "-------------------------------------------------"
      
      stop
      end

*-----------------------------------------------------------------------
      subroutine average_k(pid,pol_set,q2_set,eps_lo_set,eps_hi_set)

c     Average W and Q2 over theta_pq settings, then over low and high epsilon
c     settings, then over neg. and pos. settings,
c     and save result in averages/avek.* .

c     Fortran is annoying and can't find parameters
c     dynamically (since they must be known at compile time).
c     Therefore, I am setting is arbitrarily to allocate
c     enough space for the sets
      parameter (nbin = 10)

      real, dimension(nbin) :: aveW,errW,aveQ2,errQ2
      real, dimension(nbin) :: avett,errtt
      real, dimension(nbin,2) :: avW,erW,avQ2,erQ2
      real, dimension(nbin,2) :: avtt,ertt
      real, dimension(nbin,2,2) :: aW,eW,aQ2,eQ2
      real, dimension(nbin,2,2) :: att,ett

      real, dimension(nbin) :: thetacm_only

      real, dimension(nbin) :: eps_lo,eps_hi
      
      integer pol_set
      real q2_set

      real q2_bin
      integer t_bin, phi_bin

      integer nt,nphi
      
      real eps_set(2)

      real th_mod
      
      character*60 fn
      character*2 pol
      character*4 pid

      open (unit = 22, file =trim(pid) // "/t_bin_interval", 
     *     action='read')
      read (22,*) q2_bin, t_bin, phi_bin

      close(22)

      nt = t_bin
      nphi = phi_bin
            
      eps_set(1)=eps_lo_set
      eps_set(2)=eps_hi_set
      
      do it=1,nt

         aveW(it)=0.
         errW(it)=0.
         aveQ2(it)=0.
         errQ2(it)=0.
         avett(it)=0.
         errtt(it)=0.         

         do ip=1,1

            avW(it,ip)=0.
            erW(it,ip)=0.
            avQ2(it,ip)=0.
            erQ2(it,ip)=0.
            avtt(it,ip)=0.
            ertt(it,ip)=0.            

            do lh=1,2
               aW(it,lh,ip)=0.
               eW(it,lh,ip)=0.
               aQ2(it,lh,ip)=0.
               eQ2(it,lh,ip)=0.
               att(it,lh,ip)=0.
               ett(it,lh,ip)=0.               
            end do

         end do

      end do

c     Get low, high eps. and neg., pos. polarity data.

      do ip=1,1

         do lh=1,2

            nset=0
            open(55, file=trim(pid) // '/list.settings')
            do while(.true.)

               read(55,*,end=9) ipol,q2,eps,th_pq,tmn,tmx
               if(ipol.eq.pol_set.and.q2.eq.q2_set.and.
     &              eps.eq.eps_set(lh)) then

                  if(ipol.eq.-1) then
                     pol='mn'
                  elseif(ipol.eq.+1) then
                     pol='pl'
                  else
                     stop '*** aver: wrong pol ***'
                  endif
*                  WRITE(*,*) '------------'
*                  WRITE(*,*) 'Values read:'
*                  WRITE(*,*) '------------'
*                  WRITE(*,*) 'ipol = ', ipol
*                  WRITE(*,*) 'pol = ', pol
*                  WRITE(*,*) 'q2 = ', q2
*                  WRITE(*,*) 'eps = ', eps
*                  WRITE(*,*) 'th_pq = ', th_pq
*                  WRITE(*,*) 'tmn = ', tmn
*                  WRITE(*,*) 'tmx = ', tmx
                  write(fn,'(a4,''/kindata/kindata.'',a2,''_'',i2.2,
     *                 ''_'',i2.2,''_'',SP,i5.4,S,''.dat'')') pid, pol,
     *                 nint(q2_set*10.), nint(eps_set(lh)*100.),
     *                 nint(th_pq*1000.)
                  print*,'fn=',fn
c                 pause

                  open(66,file=fn)
                  do it=1,nt
                     read(66,*) W,dW,Q2,dQ2,tt,dtt
*                     WRITE(*,*) 'it = ', it
*                     WRITE(*,*) 'nt = ', nt
*                     WRITE(*,*) 'W = ', W
*                     WRITE(*,*) 'dW = ', dW
*                     WRITE(*,*) 'Q2 = ', Q2
*                     WRITE(*,*) 'dQ2 = ', dQ2
*                     WRITE(*,*) 'tt = ', tt
*                     WRITE(*,*) 'dtt = ', dtt                     
                     if(dW.gt.0.) then
                        aW(it,lh,ip)=aW(it,lh,ip)+W/dW**2
                        eW(it,lh,ip)=eW(it,lh,ip)+1./dW**2
                     end if
                     if(dQ2.gt.0.) then
                        aQ2(it,lh,ip)=aQ2(it,lh,ip)+Q2/dQ2**2
                        eQ2(it,lh,ip)=eQ2(it,lh,ip)+1./dQ2**2
                     end if
                     if(dtt.gt.0.) then
                        att(it,lh,ip)=att(it,lh,ip)+tt/dtt**2
                        ett(it,lh,ip)=ett(it,lh,ip)+1./dtt**2
                     end if                     
                  end do
                  close(66)

                  tmin=tmn
                  tmax=tmx
                  ntbins=nt

                  nset=nset+1

               end if           !ipol=pol_set & q2=q2_set & eps=eps_set

            end do              !while not eof.

 9          continue
            close(55)
*            WRITE(*,*) '------------'
            print*,'nset=',nset

         end do                 !lh=1,2

      end do                    !ip=1,2

c      pause

      do ip=1,1
         do lh=1,2
            do it=1,ntbins
               if (eW(it,lh,ip).gt.0.) then
                  aW(it,lh,ip)=aW(it,lh,ip)/eW(it,lh,ip)
                  eW(it,lh,ip)=1./sqrt(eW(it,lh,ip))
*****                  eW(it,lh,ip)=0.001
               end if
               if(eQ2(it,lh,ip).gt.0.) then
                  aQ2(it,lh,ip)=aQ2(it,lh,ip)/eQ2(it,lh,ip)
                  eQ2(it,lh,ip)=1./sqrt(eQ2(it,lh,ip))
*****                  eQ2(it,lh,ip)=0.001
               end if
               if(ett(it,lh,ip).gt.0.) then
                  att(it,lh,ip)=att(it,lh,ip)/ett(it,lh,ip)
                  ett(it,lh,ip)=1./sqrt(ett(it,lh,ip))
*****                  ett(it,lh,ip)=0.001
               end if               
c               write(*,'(4f8.5,2i3)') aW(it,lh,ip),eW(it,lh,ip),
c               aQ2(it,lh,ip),eQ2(it,lh,ip),it,lh
            end do
         end do
      end do
c      pause

c     Average over low and high epsilon.

      do ip=1,1
         do it=1,ntbins
            do lh=1,2
               if(eW(it,lh,ip).gt.0.) then
                  avW(it,ip)=avW(it,ip)+aW(it,lh,ip)/eW(it,lh,ip)**2
                  erW(it,ip)=erW(it,ip)+1./eW(it,lh,ip)**2
               end if
               if(eQ2(it,lh,ip).gt.0.) then
                  avQ2(it,ip)=avQ2(it,ip)+aQ2(it,lh,ip)/eQ2(it,lh,ip)**2
                  erQ2(it,ip)=erQ2(it,ip)+1./eQ2(it,lh,ip)**2
               end if
               if(ett(it,lh,ip).gt.0.) then
                  avtt(it,ip)=avtt(it,ip)+att(it,lh,ip)/ett(it,lh,ip)**2
                  ertt(it,ip)=ertt(it,ip)+1./ett(it,lh,ip)**2
               end if               
            end do
         end do
      end do

      do ip=1,1
         do it=1,ntbins
            if(erW(it,ip).gt.0.) then
               avW(it,ip)=avW(it,ip)/erW(it,ip)
               erW(it,ip)=1./sqrt(erW(it,ip))
            end if
            if(erQ2(it,ip).gt.0.) then
               avQ2(it,ip)=avQ2(it,ip)/erQ2(it,ip)
               erQ2(it,ip)=1./sqrt(erQ2(it,ip))
            end if
            if(ertt(it,ip).gt.0.) then
               avtt(it,ip)=avtt(it,ip)/ertt(it,ip)
               ertt(it,ip)=1./sqrt(ertt(it,ip))
            end if            
         end do
      end do

c     Average over neg. and pos. settings.
      do it=1,ntbins
         do ip=1,1
            if(erW(it,ip).gt.0.) then
               aveW(it)=aveW(it)+avW(it,ip)/erW(it,ip)**2
               errW(it)=errW(it)+1./erW(it,ip)**2
            end if
            if(erQ2(it,ip).gt.0.) then
               aveQ2(it)=aveQ2(it)+avQ2(it,ip)/erQ2(it,ip)**2
               errQ2(it)=errQ2(it)+1./erQ2(it,ip)**2
            end if
            if(ertt(it,ip).gt.0.) then
               avett(it)=avett(it)+avtt(it,ip)/ertt(it,ip)**2
               errtt(it)=errtt(it)+1./ertt(it,ip)**2
            end if            
         end do
      end do

      do it=1,ntbins
         aveW(it)=aveW(it)/errW(it)
         errW(it)=1./sqrt(errW(it))
         aveQ2(it)=aveQ2(it)/errQ2(it)
         errQ2(it)=1./sqrt(errQ2(it))
         avett(it)=avett(it)/errtt(it)
         errtt(it)=1./sqrt(errtt(it))         
      end do

c     Thetacm for neg. and pos. settings. It's turned out the same for
c     low and high epsilons, but different for negatives and positives.
c     So calculate for high eps., neg.-s and pos.-s.

c     Get Beam energy at first.
      Eb=0.
      open(55, file=trim(pid) // '/beam/Eb_KLT.dat')
      do while(.true.)
         read(55,*) Eb,q2,eps
         write(*,*) Eb,q2,eps
         if(q2.eq.q2_set.and.eps.eq.eps_hi_set) go to 5
      end do
 5    close(55)
c      Eb=Eb/1000.               !Mev -> Gev units.
      print*,'xsect: Eb=',Eb,'   at Q2=',q2,'  eps=',eps,'  pol=',pol

      do it=1,ntbins
         tm=tmin+(it-0.5)*(tmax-tmin)/ntbins
         call eps_n_theta(pid,pol_set,Eb,aveW(it),aveQ2(it),tm,th_mod,
     &         eps_mod)
         thetacm_only(it)=th_mod*180./3.14159   
      end do

c     Save data.

      write(fn,'(a4,''/averages/avek.'',i2.2,
     *     ''.dat'')') pid,nint(q2_set*10.)
      print*,'fn=',fn
      print*

      open(77,file=fn)
      do it=1,ntbins
         write(77,'(6f9.5,f10.2,i3)')
     *        aveW(it),errW(it),aveQ2(it),errQ2(it),
     *        avett(it), errtt(it), thetacm_only(it),it

c         write(*,'(4f8.5,i3)') aveW(it),errW(it),aveQ2(it),errQ2(it),it
      end do
      close(77)

      end

      include 'eps_n_theta.f'

