program problem
    
    use iso_fortran_env, only: dp=>real64
    use linspace_mod
    use tdma
    
    implicit none
    !output utility
    integer, parameter :: outfile = 16
    integer, parameter :: outfile2 = 17
    
    !enter problem you want to solve here
        !problem4
    !variable declarations
    integer :: i,j,n,m,count,p,q
    real(dp), dimension(:,:),allocatable :: phi,phi_dash,gmae,gmaw,gman,gmas,Del_x,Del_y,delxe,delxw,delyn,delys,sc,sp,theta
    real(dp), dimension(:),allocatable :: z
    character(2) :: sweep
    real(dp) :: error,quee,quew,quen,ques,De,phipq,omega,Tb,U,gmac,Tmw,Nu,quec
    
    !defining the problem
    !diffusion coefficient (comment out if not constant) can be set up later as a cell property
    gmac = 10.0_dp
    !velocity of the free stream
    U = 1.0_dp
    !boundary heat flux que: e,w,n,s for the respective directions, c for common
    quec = 1000.0_dp
    quee = quec
    quew = quec
    quen = quec
    ques = quec
    !side of the square
    De = 1.0_dp
    !underrelaxation parameter
    omega = 0.5_dp
    
    !inputs
    write(*,*) "Enter number of rows"
    read(*,*) n
    write(*,*) "Enter number of columns"
    read (*,*) m
    write(*,*) "Choose the sweep direction"
    read(*,*) sweep
    write(*,*) "Choose the x-coordinate for the fixed point, must be inner"
    read(*,*) p
    write(*,*) "Choose the y-coordinate for the fixed point, must be inner"
    read(*,*) q
    write(*,*) "Choose the phi value taken by the fixed point"
    read(*,*) phipq
    
    !allocate dynamic arrays
    allocate(phi(n,m),phi_dash(n,m),theta(n,m),gmae(n,m),gmaw(n,m),gman(n,m),gmas(n,m),Del_x(n,m),Del_y(n,m),delxe(n,m),delxw(n,m),delyn(n,m),delys(n,m),sc(n,m),sp(n,m),z(m))
    
    !initialize count
    count = 0
    
    !initialize cell properties here: An example is shown here but the array can be initialized in any way
    !does not have to be a uniform mesh, but has to be a structured mesh
    !>>diffusion coefficient
    gmae = gmac
    gmaw = gmac
    gman = gmac
    gmas = gmac
    !>>cell dimensions
    Del_x = De/(1.0_dp*(m-2))
    Del_y = De/(1.0_dp*(n-2))
    delxe = De/(1.0_dp*(m-2))
    delxw = De/(1.0_dp*(m-2))
    delyn = De/(1.0_dp*(n-2))
    delys = De/(1.0_dp*(n-2))
    
    !initialize the matrix
    phi = 300.0_dp
    
    !setting up cell dimensions for the boundary
    do i = 1,m
        !south boundary
        delys(2,i) = 0.5_dp*delys(2,i)
        !north boundary
        delyn(n-1,i) = 0.5_dp*delyn(n-1,i)
	end do
    do i = 1,n
       !west boundary
       delxw(i,2) = 0.5_dp*delxw(i,2)
       !east boundary
       delxe(i,m-1) = 0.5_dp*delxe(i,m-1)
    end do
    
    !main loop start
    !source term incorporation
    do while (.true.)
        !initialize error for source term loop
        error = 0.0_dp
        !count to keep track of iterations
        count = count +1
        print*,count, 'outer'
        !backup array
        phi_dash = phi
        !declare source term here
        sc = (-4.0_dp*quec/De)
        sp = 0.0_dp
        !call newtdma2d solver to solve your problem
        call newtdma2dinmat(n,m,phi,gmae,gmaw,gman,gmas,Del_x,Del_y,delxe,delxw,delyn,delys,sc,sp,p,q,phipq,quee,quew,quen,ques,omega,sweep)
        !error calculation using eucleadian distance
        error = sqrt(sum((phi-phi_dash)**2))
        !exit condition for the main loop
        if (error<1e-10) exit
    end do
    
    !print phi to terminal
    do i = 1,n
        print*,phi(i,:)
    end do
    
    !writing phi to a file
	open(unit=outfile, file = 'phifor.dat', access='sequential',action = 'write')  
   	do i=1,n
   		write (outfile,*) (phi(i,j), j=1,m)
   	end do     
   	close(outfile)
    
    !bulk temperature
    Tb = U*sum(phi*Del_y*Del_x)/(U*sum(Del_y*Del_x))
    print*,"The bulk temperature for the given problem is", Tb
    
    !dimensionless temperature
    theta = (phi-Tb)/(quee*De/gmac)
    
    !in case you want to solve for the easiest test case, que = 0
    !theta = (phi-Tb)/Tb
    
    !midline in z dxn
    z(2:m-1) = linspace(-1.0_dp*Del_x(5,5)*(m-3)/2,Del_x(5,5)*(m-3)/2,m-2)
    z(1) = z(2)-(Del_x(5,5)/2)
    z(m) = z(m-1)+(Del_x(5,5)/2)
    
    !writing midline temperature
	open(unit=outfile2, file = 'midtemp.dat', access='sequential',action = 'write')  
   	
    if(mod(n,2)==1) then
        do i=1,m
   		    write (outfile2,*) z(i),theta((n+1)/2,i)
        end do
    end if
    if(mod(n,2)==0) then
        do i=1,m
   		    write (outfile2,*) z(i),0.5_dp*(theta((n/2),i)+theta((n/2)+1,i))
        end do
    end if
    
   	close(outfile2)
    
    !mean wall temperature
    Tmw = 0.0_dp
    do i = 1,n
        Tmw = Tmw+phi(i,1)+phi(i,m)
    end do
    do i = 1,m
        Tmw = Tmw+phi(1,i)+phi(n,i)
    end do
    Tmw = Tmw - phi(1,1) - phi(1,m) - phi(n,1) - phi(n,n)
    Tmw = Tmw/(2*(m+n)-4)
    print*, "The mean wall temperature is: ", Tmw
    
    !Nusselt number
    Nu = (quec*De/gmac)/(Tmw-Tb)
    print*, "The nusselt number is: ", Nu
    
    !deallocate the dynamic arrays
    deallocate(phi,phi_dash,theta,gmae,gmaw,gman,gmas,Del_x,Del_y,delxe,delxw,delyn,delys,sc,sp,z)
end program problem
    