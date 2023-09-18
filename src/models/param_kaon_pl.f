	program iterWeight
	implicit none

	real nsigl,nsigt,nsiglt,nsigtt,tmp
	real nsig219,nsig,wtn
	real ft,tav,ftav
	
	real pi,mtar_gev,q2_gev
	real my_limit
	parameter (pi=3.14159)
	parameter (mtar_gev=0.93827231)

	real p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12
	
**********************************************	
*	Read in arguments of parameters and Q2
	integer :: q2_set
	real :: params(12)
	integer :: i, argc
	character(len=20) :: arg

	! Get the number of command line arguments
	argc = COMMAND_ARGUMENT_COUNT()

	! Check if there are enough arguments
	if (argc < 13) then
	   print *, "Error: Not enough arguments provided."
	   stop
	end if

	! Get q2_set from the first argument
	call GET_COMMAND_ARGUMENT(1, arg)
	read(arg, *) q2_set

	! Get params from the rest of the arguments
	do i = 2, 13
	   call GET_COMMAND_ARGUMENT(i, arg)
	   read(arg, *) params(i-1)
	end do
**********************************************
	
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

	! Print q2_set, p1 to p12 for verification
	print *, "q2_set =", q2_set
	print *, "p1 =", p1
	print *, "p2 =", p2
	print *, "p3 =", p3
	print *, "p4 =", p4
	print *, "p5 =", p5
	print *, "p6 =", p6
	print *, "p7 =", p7
	print *, "p8 =", p8
	print *, "p9 =", p9
	print *, "p10 =", p10
	print *, "p11 =", p11
	print *, "p12 =", p12	
	
	end program iterWeight
