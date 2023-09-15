	function wtn(q2_set)
	include ?

	real p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12
	real nsigl,nsigt,nsiglt,nsigtt,tmp
	real nsig219,nsig,wtn
	real ft,tav,ftav

	real pi,mtar_gev,q2_gev
	real my_limit
	integer q2_set
	parameter (pi=3.14159)
	parameter (mtar_gev=0.93827231)

	if (abs(q2_set-245).lt.1) then
	   p1=  0.25961E+02 
	   p2= -0.10000E+02 
	   p3= -0.15838E+02 
	   p4=  0.00000E+00 
	   p5=  0.46859E+02 
	   p6= -0.30000E+02 
	   p7= -0.33572E+01 
	   p8=  0.00000E+00 
	   p9=  0.10000E+04 
	   p10=-0.28000E+02 
	   p11= 0.35000E+01 
	   p12=-0.67276E+02 
	else
	   write(*,*)'wtn: q2 error ',q2_set
	endif

***
* Parameterization based upon Fpi-1 pi+ IT25, 12.04.18
* Revised for IT21, 12.11.06
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

        return
	end
