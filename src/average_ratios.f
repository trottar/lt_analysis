      program average_ratios

c     Taken exp. and MC yields at each theta_pq setting, this program
c     calculates exp/simc ratios for polarity-Q2-eps setting. The idea
c     comes from Jochen.
c
c     Input: yields/yields.*.dat
c     Output: averages/aver.*.dat

      integer inp_pol
      real inp_Q2, inp_loeps, inp_hieps
      write(*,*) "Inputing particle, polarity, Q2, W and both epsilons:"
      read(*,*) inp_pid, inp_pol, inp_Q2, inp_W, inp_loeps, inp_hieps

      write(*,*) "PID = ",inp_pid,"POL = ",inp_pol,
     *     "Q2 = ",inp_Q2,"W = ",inp_W,
     *     "low_eps = ",inp_loeps,"high_eps = ",inp_hieps
      
      call average_r(inp_pid, inp_pol, inp_Q2, inp_W,
     *     inp_loeps)
      call average_r(inp_pid, inp_pol, inp_Q2, inp_W,
     *     inp_hieps)
      print*,  "-------------------------------------------------"

      
      stop
      end

*-----------------------------------------------------------------------

      subroutine average_r(pid,pol_set,q2_set,w_set,
     *     eps_set)

      implicit none

c     Aquire yields over theta_pq settings, calculate ratio, save result
c     in averages/aver.* .


c     Fortran is annoying and can't find parameters
c     dynamically (since they must be known at compile time).
c     Therefore, I am setting is arbitrarily to allocate
c     enough space for the sets
      parameter (nbin = 10)
      
      integer nt,nphi,it,ip
      integer nset,ipol
      real q2,eps,th_pq,tmn,tmx,r,e,er
      real one
      real yld
      
      integer pol_set
      real q2_set, w_set

      real q2_bin, w_bin
      integer t_bin, phi_bin

      real t_min
      
      real yrd(nphi,nt),drd(nphi,nt)
      real ymc(nphi,nt),dmc(nphi,nt)

      character*80 fn
      character*2 pol
      character*4 pid
      
      character(len=100) :: fn_t_bins
      
!     Construct the file path using a format string
      write(fn_t_bins, '(a, a, i2.2, a, i3.3, a)') trim(pid),
     *     '/t_bin_interval_Q', nint(q2_set*10), 'W', nint(w_set*100)

!     Open the file
      open (unit=22, file=fn_t_bins, action='read')
      read (22, *) q2_bin, w_bin, t_bin, phi_bin
      
      close(22)

      nt = t_bin
      nphi = phi_bin
      
      do it=1,nt
         do ip=1,nphi
            yrd(ip,it)=0.
            drd(ip,it)=0.
            ymc(ip,it)=0.
            dmc(ip,it)=0.
         end do
      end do

      nset=0
      open(55, file=trim(pid) // '/list.settings')
      do while(.true.)

         read(55,*,end=9) ipol,q2,w,eps,th_pq,tmn,tmx         
         if(ipol.eq.pol_set.and.q2.eq.q2_set.and.
     *        w.eq.w_set.and.eps.eq.eps_set) then

            if(ipol.eq.-1) then
               pol='mn'
            elseif(ipol.eq.+1) then
               pol='pl'
            else
               stop '*** aver: wrong pol ***'
            endif

            WRITE(*,*) '------------'
            WRITE(*,*) 'Values read:'
            WRITE(*,*) '------------'
            WRITE(*,*) 'ipol = ', ipol
            WRITE(*,*) 'pol = ', pol
            WRITE(*,*) 'q2 = ', q2
            WRITE(*,*) 'eps = ', eps
            WRITE(*,*) 'th_pq = ', th_pq
            WRITE(*,*) 'tmn = ', tmn
            WRITE(*,*) 'tmx = ', tmx

c     Read real data.
            write(fn,'(a4,''/yields/yield_data.'',a2,''_Q'',
     *           i2.2,''W'',i3.3,''_'',i2.2,''_'',SP,
     *           i5.4,S,''dat'')') pid, pol, 
     *           nint(q2_set*10.), nint(w_set*100.), 
     *           nint(eps_set(lh)*100.), nint(th_pq*1000.)
            print*,'fn=',fn
c            pause

            open(66,file=fn)
            read(66,*) one
            do it=1,nbin
               do ip=1,nphi
                  read(66,*) yld,er
*                  print*,yld,er
                  yrd(ip,it)=yrd(ip,it)+yld
                  drd(ip,it)=drd(ip,it)+er**2
               end do
            end do
            close(66)

c     Read real simc.
            write(fn,'(a4,''/yields/yield_simc.'',a2,''_Q'',
     *           i2.2,''W'',i3.3,''_'',i2.2,''_'',SP,
     *           i5.4,S,''dat'')') pid, pol, 
     *           nint(q2_set*10.), nint(w_set*100.), 
     *           nint(eps_set(lh)*100.), nint(th_pq*1000.)
            print*,'fn=',fn
c            pause
            
            open(66,file=fn)
            read(66,*) one
            do it=1,nbin
               do ip=1,nphi
                  read(66,*) yld,er
*                  print*,yld,er
                  ymc(ip,it)=ymc(ip,it)+yld
                  dmc(ip,it)=dmc(ip,it)+er**2
               end do
            end do
            close(66)

            do it=1,nt
               do ip=1,nphi
                  write(*,*) it,ip
                  write(*,*)'exp_yield=',yrd(ip,it),drd(ip,it)
                  write(*,*)'MC_yield=',ymc(ip,it),dmc(ip,it)
               enddo
            enddo

            nset=nset+1

         end if                 !ipol=pol_set & q2=q2_set & eps=eps_set

      end do                    !while not eof.

 9    continue
      close(55)

      print*,'nset=',nset

c      pause

      write(fn,10) pid,pol,nint(q2*10),nint(w*100),nint(eps*100)
 10   format(a4,'/averages/aver.'
     *     ,a2,'_Q',i2.2,'W',i3.3,'_',i2,'.dat')
      print*,'fn=',fn
      print*

      open(77,file=fn,status='replace')
      do it=1,nbin
         do ip=1,nphi
            r=0.
            e=0.
c ratio is data/simc - see GH logbook, p.55
            if(ymc(ip,it).ne.0.) then
               r=(yrd(ip,it))/ymc(ip,it)
*     RLT (2/11/2024): These are G.H's equations
*                      Removing to stay consistent
*                      with units of %
*     e=e+(drd(ip,it))/ymc(ip,it)**2
*     e=e+((r/ymc(ip,it))**2)*dmc(ip,it)
              e=drd(ip,it)+dmc(ip,it)
              e=sqrt(e)
              write(*,*)'     the data/simc ratio is:'
              write(*,'(2f15.5,2i3)')r,e,ip,it
              write(*,*)'--------------------------'
            end if
            write(77,'(2f15.5,2i3)') r,e,ip,it
         end do
      end do
      close(77)

      end
