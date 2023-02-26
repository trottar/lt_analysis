      program take_averages


      implicit none

c     This program takes averages of W and Q2 over different theta_pq
c     settings, then different polarities.
c
c     Input:  kindata/kindata.*.dat
c     Output: averages/averages.*.dat

      integer inp_pol
      real inp_Q2, inp_loeps, inp_hieps
      write(*,*) "Please input your polarity, Q2 and low+high epsilon:"
      read(*,*) inp_pol, inp_Q2, inp_loeps, inp_hieps

      write(*,*) "POL = ",inp_pol,"Q2 = ",inp_Q2,
     *           "low_eps = ",inp_loeps,"high_eps = ",inp_hieps
      
      call average_k(inp_pol, inp_Q2, inp_loeps, inp_hieps)
      print*,  "-------------------------------------------------"
      
      stop
      end

*-----------------------------------------------------------------------

      subroutine average_k(npol_set,q2_set,eps_lo_set,eps_hi_set)


      implicit none

c     Average W and Q2 over theta_pq settings, then over low and high epsilon
c     settings, then over neg. and pos. settings,
c     and save result in averages/aver.* .

c      implicit none
      
      integer nu

      parameter (nu=10)

      real aveW(nu),errW(nu),aveQ2(nu),errQ2(nu)
      real avW(nu,2),erW(nu,2),avQ2(nu,2),erQ2(nu,2)
      real aW(nu,2,2),eW(nu,2,2),aQ2(nu,2,2),eQ2(nu,2,2)

      real avett(nu),    errtt(nu)    
      real avtt(nu,2),   ertt(nu,2)    
      real att(nu,2,2),  ett(nu,2,2)    


c      real thetacm_neg(nu),thetacm_pos(nu)
      real thetacm_only(nu)

      real eps_lo(nu),eps_hi(nu),um, u_min

      real eps_set(2)

      integer pol_set(1), j, npol_set

      real q2_bin, q2_set, eps_lo_set, eps_hi_set, dq2, dtt, dw, eb, eps
      real tt, wwmx, w, tmn, tmx, tmax, tmin, tm, th_mod, q2

      real one,th_pq, nset, ipol, eps_mod
      integer it, ip, lh, nbt, ntbins

      integer t_bin, phi_bin 

      real, Dimension(10) :: t_bin_boundary

      character*80:: line

      character*40 fn
      character*2 pol

      eps_set(1)=eps_lo_set
      eps_set(2)=eps_hi_set

      pol_set(1)=npol_set

      do it=1,nu

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



c c   /*--------------------------------------------------*/
c c   Read the u and phi bins 
c 

      print*, "BBBBBBBBB", q2_bin       

      open (unit = 22, file = "./t_bin_interval", action='read')
      read (22,*) q2_bin, t_bin, phi_bin


      print*, "BBBBBBBBB", q2_bin       



c      nt = t_bin
c      nphi = phi_bin 
c      read (22, '(A)') line  
c      read (22, '(A)') line  
c      print*,  line


c      read (line, end=20) 

c      print*,  trim(line)

c      read(line, *) t_bin_boundary(0), t_bin_boundary(1), t_bin_boundary(2)

c      read(line, *) (t_bin_boundary(j), j = 1, 3)

c      print*, t_bin_boundary(1)
c       read (22,*) 
c       read (22,*) q2_bin, t_bin, phi_bin
c       read (22,*) 
c 
c       print*,  t_bin, phi_bin

c      do j = 1, 3
c         print*, t_bin_boundary(j)
c      end do




      print*, q2_set

      if(q2_set.eq.1.6) then

         read (22, '(A)') line  
c         print*, line
         read(line, *) (t_bin_boundary(j), j = 1, t_bin+1)

c         print*, t_bin+1 

c         stop

c        t_bin_boundary = (/ 0.0, 0.12,  0.20, 0.40/)
c        t_bin_boundary = (/0.0, 0.10, 0.17, 0.32/)
 
      elseif(q2_set.eq.2.45) then
 
c        t_bin_boundary = (/ 0.0, 0.212, 0.33, 0.60/)
c        t_bin_boundary = (/0.0, 0.19, 0.30, 0.50/)
 
         read (22,*) 
         read (22,*) 
         read (22, '(A)') line  
         read(line, *) (t_bin_boundary(j), j = 1, t_bin+1)
         
      endif

c      stop t_bin+1

      close(22)

c      print*, t_bin , t_bin_boundary(2), t_bin_boundary(3) 
c      stop

      nbt = t_bin 

c     Get low, high eps. and neg., pos. polarity data.

      do ip=1,1

         do lh=1,2

            nset=0
            open(55,file='./list.settings')
            do while(.true.)

               read(55,*,end=9) ipol,q2,eps,th_pq,tmn,tmx,nbt
c               read(55,*,end=9) ipol,q2,eps,th_pq
               if(ipol.eq.pol_set(ip).and.q2.eq.q2_set.and.
     &              eps.eq.eps_set(lh)) then

                  if(ipol.eq.-1) then
                     pol='mn'
                  elseif(ipol.eq.+1) then
                     pol='k'
                  else
                     stop '*** aver: wrong pol ***'
                  endif


                  print*, "nbt: ", nbt

                  write(fn,'(''kindata/kindata.'',a1,''_'',i2.2,''_'',i2.2,
     *                 ''_'',SP,i5.4,S,''.dat'')')
     *                 pol,nint(q2_set*10.),nint(eps_set(lh)*100.),
     *                 nint(th_pq*1000.)
                  print*,'fn=',fn
c                 pause


c                 print*, 'aaaaaaaaaaaaaaaaaaaa '

                  open(66,file=fn)
                  read(66,*) one

                  print*, "bbbbbbbbbbbbbbbbbbbb ", one


                  do it=1,nbt
                     read(66,*) Q2,dQ2,W,dW,tt,dtt
                     print*,Q2,dQ2,W,dW,tt,dtt,it
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
                  read(66,*) one
                  close(66)

                  tmin=tmn
                  tmax=tmx
                  ntbins=nbt

                  nset=nset+1

               end if           !ipol=pol_set & q2=q2_set & eps=eps_set

            end do              !while not eof.

 9          continue
            close(55)

            print*,'nset=',nset

         end do                 !lh=1,2

      end do                    !ip=1,2




c      pause

c       open (unit = 22, file = "t_bin_interval", action='read')
c       read (22,*) q2_bin, t_bin, phi_bin
c       close(22)
c 
c 
cc c   /*--------------------------------------------------*/
cc c   Read the u and phi bins 
cc 
c
c      open (unit = 22, file = "../t_bin_interval", action='read')
c      read (22,*) q2_bin, t_bin, phi_bin
c
c        
c
c
cc      nt = t_bin
cc      nphi = phi_bin 
cc      read (22, '(A)') line  
cc      read (22, '(A)') line  
cc      print*,  line
c
c
cc      read (line, end=20) 
c
cc      print*,  trim(line)
c
cc      read(line, *) t_bin_boundary(0), t_bin_boundary(1), t_bin_boundary(2)
c
cc      read(line, *) (t_bin_boundary(j), j = 1, 3)
c
cc      print*, t_bin_boundary(1)
cc       read (22,*) 
cc       read (22,*) q2_bin, t_bin, phi_bin
cc       read (22,*) 
cc 
cc       print*,  t_bin, phi_bin
c
cc      do j = 1, 3
cc         print*, t_bin_boundary(j)
cc      end do
c
c
c
c
c      print*, q2_set
c
c      if(q2_set.eq.1.6) then
c
c         read (22, '(A)') line  
c
c         print*, line
c
c         read(line, *) (t_bin_boundary(j), j = 1, t_bin+1)
c
cc        t_bin_boundary = (/ 0.0, 0.12,  0.20, 0.40/)
cc        t_bin_boundary = (/0.0, 0.10, 0.17, 0.32/)
c 
c      elseif(q2_set.eq.2.45) then
c 
cc        t_bin_boundary = (/ 0.0, 0.212, 0.33, 0.60/)
cc        t_bin_boundary = (/0.0, 0.19, 0.30, 0.50/)
c 
c         read (22,*) 
c         read (22,*) 
c         read (22, '(A)') line  
c         read(line, *) (t_bin_boundary(j), j = 1, t_bin+1)
c
c      endif
c
c      
c      
c
c
cc      stop
c
c      close(22)
c











c      print*, t_bin, phi_bin
c      stop
      
      ntbins = t_bin

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
*****                  eQ2(it,lh,ip)=0.001
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
c      open(55,file='Beam/Eb_fpi2.dat')



      open(55,file='./beam/Eb_KLT.dat')
      do while(.true.)
         read(55,*) Eb,q2,eps
         write(*,*) Eb,q2,eps
         if(q2.eq.q2_set.and.eps.eq.eps_hi_set) go to 5
      end do
 5    close(55)

      print*,'xsect: Eb=',Eb,'   at Q2=',q2,'  eps=',eps,'  pol=',pol





c      do it=1,ntbins
c         tm=tmin+(it-0.5)*(tmax-tmin)/ntbins
c         call eps_n_theta(-1,Eb,aveW(it),aveQ2(it),tm,th_mod,eps_mod)
c         thetacm_neg(it)=th_mod*180./3.14159
c         call eps_n_theta(+1,Eb,aveW(it),aveQ2(it),tm,th_mod,eps_mod)
c         thetacm_pos(it)=th_mod*180./3.14159
c      end do

      do it=1,ntbins
         tm=tmin+(it-0.5)*(tmax-tmin)/ntbins
         um = (t_bin_boundary(it) + t_bin_boundary(it+1)) / 2


         print*, tmin, tmax, ntbins
         print*, tm, um


         print*, "~~~~~~2",q2, w, Eb, eps        

         call eps_n_theta(npol_set,Eb,aveW(it),aveQ2(it),
     *                    tm,um,u_min,th_mod,eps_mod)

         print*, "!!!!!!",th_mod


c         stop

c         call eps_n_theta(-1,Eb,aveW(it),aveQ2(it),tm,th_mod,eps_mod)
c         thetacm_neg(it)=th_mod*180./3.14159

c         call eps_n_theta(+1,Eb,aveW(it),aveQ2(it),avett(it),th_mod, 
c     *  eps_mod)







         thetacm_only(it)=th_mod*180./3.14159

      end do


c      stop




c     Save data.

      write(fn,'(''averages/avek.'',i2.2,''.dat'')')
     *     nint(q2_set*10.)
      print*,'fn=',fn
      print*

c      open(77,file=fn)
c      do it=1,ntbins
c         write(77,'(6f8.5,2f10.5,i3)')
c     *        aveW(it),errW(it),aveQ2(it),errQ2(it),
c     *        avett(it), errtt(it),
c     *        thetacm_neg(it),thetacm_pos(it),it
c



      open(77,file=fn)
      do it=1,ntbins
         write(77,'(6f9.5,f10.2,i3)')
     *        aveW(it),errW(it),aveQ2(it),errQ2(it),
     *        avett(it), errtt(it), thetacm_only(it),it



c      print*, "XXXXXX" 

c      stop


c         write(*,'(4f8.5,i3)') aveW(it),errW(it),aveQ2(it),errQ2(it),it
      end do
      close(77)

      end


      include 'eps_n_theta.f'

