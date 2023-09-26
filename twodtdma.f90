module tdma2
    use tdma
    !code for extending onedtdma to two dimensions
    subroutine twodtdma(n,m,mat,gmae,gmaw,gman,gmas,Del_x,Del_y,delxe,delxw,delyn,delys,sc,sp,sweep)
        implicit none
        !declarations
        integer, intent (in) :: n,m
        character(1), intent (in) :: sweep 
        double precision, intent (in) :: sc,sp
        double precision, dimension (1:n,1:m), intent (in) :: gmae,gmaw,gman,gmas,Del_x,Del_y,delxe,delxw,delyn,delys
        double precision, dimension (1:n,1:m), intent (inout) :: mat
        double precision, dimension (1:n,1:m) :: mat_dash
        double precision, dimension (1:m) :: lanex, ax,bx,cx,dx
        double precision, dimension (1:n) :: laney, ay,by,cy,dy
        double precision :: ae, aw, an, as, ap, delv,b, eps
        integer :: i,j
        !discretize to find a,b,c,d
        if (sweep == 'x') then
            do while (.true.)
                eps = 0.0
                !backup matrix
                mat_dash = mat
                do i = 2,m-1
                    laney = mat_dash(:,i)
                    !setting boundaries
                    ay(1) = 1.0
                    by(1) = 0.0
                    cy(1) = 0.0
                    dy(1) = mat_dash(i,1)
                    ay(n) = 1.0
                    by(n) = 0.0
                    cy(n) = 0.0
                    dy(n) = mat_dash(i,n)
                    do j = 2,n-1
                        ae = gmae(i,j)*Del_y(i,j)/delxe(i,j)
                        aw = gmaw(i,j)*Del_y(i,j)/delxw(i,j)
                        an = gman(i,j)*Del_x(i,j)/delyn(i,j)
                        as = gmas(i,j)*Del_x(i,j)/delys(i,j)
                        delv = Del_x(i,j)*Del_y(i,j)
                        b = sc*delv
                        ap = ae+aw+an+as-(sp*delv)
                        ay(j) = ap
                        by(j) = ae
                        cy(j) = aw
                        dy(j) = an*mat_dash(i,j+1)+as*mat_dash(i,j-1)+b                    
                    end do
                    call onetdma(n,ay,by,cy,dy,laney)
                    mat(:,i) = laney
                end do
                !error using eucleadian distance
		        eps = sqrt(sum((mat-mat_dash)**2))
                if (eps<1e-6) exit
            end do
        else if (sweep == 'y') then
            do while (.true.)
                eps = 0.0
                !backup matrix
                mat_dash = mat
                do i = 2,n-1
                    lanex = mat_dash(i,:)
                    !setting boundaries
                    ax(1) = 1.0
                    bx(1) = 0.0
                    cx(1) = 0.0
                    dx(1) = mat_dash(i,1)
                    ax(m) = 1.0
                    bx(m) = 0.0
                    cx(m) = 0.0
                    dx(m) = mat_dash(i,m)
                    do j = 2,m-1
                        ae = gmae(i,j)*Del_y(i,j)/delxe(i,j)
                        aw = gmaw(i,j)*Del_y(i,j)/delxw(i,j)
                        an = gman(i,j)*Del_x(i,j)/delyn(i,j)
                        as = gmas(i,j)*Del_x(i,j)/delys(i,j)
                        delv = Del_x(i,j)*Del_y(i,j)
                        b = sc*delv
                        ap = ae+aw+an+as-(sp*delv)
                        ax(j) = ap
                        bx(j) = ae
                        cx(j) = aw
                        dx(j) = an*mat_dash(i,j+1)+as*mat_dash(i,j-1)+b                    
                    end do
                    call onetdma(m,ax,bx,cx,dx,lanex)
                    mat(i,:) = lanex
                end do
                !error using eucleadian distance
		        eps = sqrt(sum((mat-mat_dash)**2))
                if (eps<1e-6) exit
            end do
        end if
    end subroutine twodtdma
end module