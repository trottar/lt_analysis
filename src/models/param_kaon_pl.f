program iterWeight
	implicit none
	
	real p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12
	real, dimension(12), intent(in) :: params
	real nsigl,nsigt,nsiglt,nsigtt,tmp
	real nsig219,nsig,wtn
	real ft,tav,ftav
	
	real pi,mtar_gev,q2_gev
	real my_limit
	integer q2_set
	parameter (pi=3.14159)
	parameter (mtar_gev=0.93827231)

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
***     
*       Parameterization based upon Fpi-1 pi+ IT25, 12.04.18
*       Revised for IT21, 12.11.06
*        tav=(0.0735+0.028*log(Q2i))*Q2i
	q2_gev=float(q2_set)/100.
        tav=(0.0735+0.028*log(q2_gev))*q2_gev
        ftav=(ti-tav)/tav
	ft=ti/(ti+0.139570**2)**2

*	nsigl=(p1+p2*log(Q2i))*exp((p3+p4*log(Q2i))*ti)
	nsigl=(p1+p2*log(Q2i))*exp((p3+p4*log(Q2i))*(ti-0.2))
	nsigt=p5+p6*log(Q2i)+(p7+p8*log(Q2i))*ftav

*	nsiglt=p9*exp(p10*ti)*sin(thetacmi)
	nsiglt=(p9*exp(p10*ti)+p11/ti)*sin(thetacmi)
	nsigtt=(p12*Q2i*exp(-Q2i))*ft*sin(thetacmi)**2

	nsig219=(nsigt+epsiloni*nsigl+epsiloni*cos(2.*phicmi)*nsigtt
     1      +sqrt(2.0*epsiloni*(1.+epsiloni))*cos(phicmi)*nsiglt)/1.d0

	wfactor=1.D0/(Wcmi**2-mtar_gev**2)**2
	nsig=nsig219*wfactor

	nsig=nsig/2./pi/1.d+06   !dsig/dtdphicm in microbarns/MeV**2/rad

	wtn=Weight*nsig/dsigdt

	my_limit=0.20
        if ((wtn.lt.my_limit).and.(wtn.gt.0.0)) then
          continue
        else
          wtn=0.
        endif

!       Print the value of wtn
	print *, "wtn =", wtn

end program iterWeight
