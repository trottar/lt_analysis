	program iterWeight
	implicit none
	
	real pi,mtar_gev
	real mpipl, mkpl
	parameter (pi=3.14159)
	parameter (mtar_gev=0.93827231)
	parameter (mpipl=0.139570)
	parameter (mkpl=0.493677)

	real*8 q2_gev, t_gev, s, s_gev
	real*8 tav, ftav, ft
	real*8 sigl, sigt, siglt, sigtt, sig219, sig
	real*8 wfactor, wtn
	
**********************************************	
*	Read in arguments of parameters and Q2
	real*8 p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12	
	integer :: q2_set
	real*8 :: params(12)
	real*8 :: q2_sim, w_sim, t_sim, eps_sim
	real*8 :: thetacm_sim, phicm_sim
	real*8 :: sigcm_sim
	real*8 :: wt_sim
	integer :: i, argc
	character(len=30) :: arg

	! Get the number of command line arguments
	argc = COMMAND_ARGUMENT_COUNT()

	! Check if there are enough arguments
	if (argc < 21) then
	   print *, "Error: Not enough arguments provided."
	   stop
	end if

	! Get q2_set from the argument
	call GET_COMMAND_ARGUMENT(1, arg)
	read(arg, *) q2_set

	! Get q2_sim from the argument
	call GET_COMMAND_ARGUMENT(2, arg)
	read(arg, *) q2_sim

	! Get w_sim from the argument
	call GET_COMMAND_ARGUMENT(3, arg)
	read(arg, *) w_sim

	! Get t_sim from the argument
	call GET_COMMAND_ARGUMENT(4, arg)
	read(arg, *) t_sim

	! Get eps_sim from the argument
	call GET_COMMAND_ARGUMENT(5, arg)
	read(arg, *) eps_sim

	! Get thetacm_sim from the argument
	call GET_COMMAND_ARGUMENT(6, arg)
	read(arg, *) thetacm_sim

	! Get phicm_sim from the argument
	call GET_COMMAND_ARGUMENT(7, arg)
	read(arg, *) phicm_sim

	! Get sigcm_sim from the argument
	call GET_COMMAND_ARGUMENT(8, arg)
	read(arg, *) sigcm_sim

	! Get wt_sim from the argument
	call GET_COMMAND_ARGUMENT(9, arg)
	read(arg, *) wt_sim
	
	! Get params from the rest of the arguments
	do i = 10,21
	   call GET_COMMAND_ARGUMENT(i, arg)
	   read(arg, *) params(i-9)
	end do
	
	p1 = params(1)
	p2 = params(2)
	p3 = params(3)
	p4 = params(4)
	p5 = params(5)
	p6 = params(6)
	p7 = params(7)
	p8 = params(8)
	p9 = params(9)
	p10 = params(10)
	p11 = params(11)
	p12 = params(12)

c	DEBUGGING
	! Print q2_set, p1 to p12 for verification
c	print *, "q2_set =", q2_set
c	print *, "q2_sim =", q2_sim
c	print *, "w_sim =", w_sim
c	print *, "t_sim =", t_sim
c	print *, "eps_sim =", eps_sim
c	print *, "thetacm_sim =", thetacm_sim
c	print *, "phicm_sim =", phicm_sim
c	print *, "sigcm_sim =", sigcm_sim
c	print *, "wt_sim =", wt_sim
c	print *, "p1 =", p1
c	print *, "p2 =", p2
c	print *, "p3 =", p3
c	print *, "p4 =", p4
c	print *, "p5 =", p5
c	print *, "p6 =", p6
c	print *, "p7 =", p7
c	print *, "p8 =", p8
c	print *, "p9 =", p9
c	print *, "p10 =", p10
c	print *, "p11 =", p11
c	print *, "p12 =", p12	

*       ALL THIS WORKS
**********************************************
	q2_gev=q2_set/1.d6
*	t_gev=t_sim/1.d6 * This is probably not required
	t_gev=t_sim
* 	W~sqrt(s), if Mp >> E_interaction
	s = w_sim**2
	s_gev=s/1.d6
	
*	tav=(0.0735+0.028*log(q2_gev))*q2_gev
*       RLT (10/8/2023): Testing new tav parameterization
	tav=(0.1112 + 0.0066*log(q2_set))*q2_set	
	ftav=(abs(t_gev)-tav)/tav
	ft=t_gev/(abs(t_gev)+mkpl**2)**2

	sigl=(p1+p2*log(q2_gev))
     &       *exp((p3+p4*log(q2_gev))*(abs(t_gev)-0.2))
	sigt=p5+p6*log(q2_gev)
     &       +(p7+p8*log(q2_gev))*ftav

	siglt=(p9*exp(p10*abs(t_gev))
     &       +p11/abs(t_gev))*sin(thetacm_sim)
	sigtt=(p12*q2_gev*exp(-q2_gev))*ft*sin(thetacm_sim)**2

	tav=(-0.178+0.315*log(q2_gev))*q2_gev

	sig219=(sigt+eps_sim*sigl+eps_sim*cos(2.*phicm_sim)*sigtt
     &         +sqrt(2.0*eps_sim*(1.+eps_sim))
     &         *cos(phicm_sim)*siglt)/1.d0
	
	wfactor=1.D0/(s_gev-mtar_gev**2)**2
	sig=sig219*wfactor
	sig=sig/2./pi/1.d+06	!dsig/dtdphicm in microbarns/MeV**2/rad

	wtn = wt_sim*sig/sigcm_sim

        if (wtn.gt.0.0) then
          continue
        else
          wtn=0.
        endif	

	print *, "", wtn
	
	end program iterWeight
