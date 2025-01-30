program harmonics 
    use FFT 
    integer, parameter :: n_theta=256, n_phi=256
    integer, parameter :: l_max=2
    real(kind=8), dimension(n_theta,n_phi) :: theta
    real(kind=8), dimension(n_theta,n_phi) :: phi
    complex(kind=8), dimension(0:l_max,-l_max:l_max,n_theta,n_phi) :: y_lm, y_analitc
    real(kind=8), dimension(0:l_max,-l_max:l_max,n_theta,n_phi) :: error
    character(len=20) :: filename

    integer :: i, j, l, m


    ! Mallas angulares
    do i = 1, n_theta
        do j=1, n_phi
        theta(i,j) = (dble(i-1)) * pi / dble(n_theta)
        end do
    end do
    do i = 1, n_phi
        do j=1, n_phi
            phi(i,j) = dble(j-1) * 2.0d0 * pi / dble(n_phi)
        end do
    end do

    do l=0, l_max
        do m=-l, l

    if(m.lt.0) then 
        write(filename, '(A,I1,A,I1,A)') 'Y_', l, '_-', abs(m), '.dat'
    else
        write(filename, '(A,I1,A,I1,A)') 'Y_', l, '_', m, '.dat'
    end if
    print*, filename, l, m
        
    open(1,file=trim(filename))
    
    do i=1, n_theta
        do j=1, n_theta
            y_lm(l,m,i,j) = compute_spherical_harmonic(l,m,theta(i,j),theta(i,j))
        end do
    end do

    do i=1, n_theta
        do j=1, n_theta
            write(1,*) theta(i,j), phi(i,j), dreal(y_lm(l,m,i,j)), dimag(y_lm(l,m,i,j)), abs(y_lm(l,m,i,j)) 
        end do
        write(1,*)
    end do
    write(1,*)
    write(1,*)

    close(1)

        end do
    end do


    y_analitc(0,0,:,:) = 1.0d0/dsqrt(4.0d0*pi)

    y_analitc(1,-1,:,:) = 0.5d0*dsqrt(0.5d0*3.0d0/pi)*exp(-i_unreal*phi)*dsin(theta)
    y_analitc(1,0,:,:)  = 0.5d0*dsqrt(3.0d0/pi)*dcos(theta)
    y_analitc(1,1,:,:)  = -0.5d0*dsqrt(0.5d0*3.0d0/pi)*exp(i_unreal*phi)*dsin(theta)

    y_analitc(2,-2,:,:) = 0.25d0*dsqrt(0.5d0*15.0d0/pi)*exp(-2.0d0*i_unreal*phi)*dsin(theta)**2
    y_analitc(2,-1,:,:) = 0.5d0*dsqrt(0.5d0*15.0d0/pi)*exp(-i_unreal*phi)*dsin(theta)*dcos(theta)
    y_analitc(2,0,:,:)  = 0.25d0*dsqrt(5.0d0/pi)*(3.0d0*cos(theta)**2 - 1.0d0)
    y_analitc(2,1,:,:)  = -0.25d0*dsqrt(0.5d0*15.0d0/pi)*exp(i_unreal*phi)*dsin(theta)*dcos(theta)
    y_analitc(2,2,:,:)  = 0.25d0*dsqrt(0.5d0*15.0d0/pi)*exp(2.0d0*i_unreal*phi)*dsin(theta)**2


    do l=0, l_max
        do m=-l, l
            
            error(l,m,:,:) = abs(y_lm(l,m,:,:) - y_analitc(l,m,:,:)) 

            if(m.lt.0) then 
                write(filename, '(A,I1,A,I1,A)') 'error_', l, '_-', abs(m), '.dat'
            else
                write(filename, '(A,I1,A,I1,A)') 'error_', l, '_', m, '.dat'
            end if
            print*, filename, l, m
                
            open(1,file=trim(filename))
                do i=1,n_theta
                    do j=1, n_phi
                        write(1,*) theta(i,j), phi(i,j), error(l,m,i,j)
                    end do
                    write(1,*)
                end do    
                write(1,*)
                write(1,*)
            close(1)

        end do
    end do


end program


