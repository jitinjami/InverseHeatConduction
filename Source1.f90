subroutine main(n,ch,ti)
    implicit none
    real :: le, dx, n, tt, dt, ts, k, cp, rho, q, q0, theta, s, ps, t_end, count, er, q1, t_end0, t_end1, t_end2, ch, ti, ddt
    real, dimension(n) :: it, it2, t_mat1
    real, dimension(ti) :: q2, t_mat2, tt1
    integer :: i
    
    open (10, file = 'time_intervals.txt')
    read (10,*) tt1
    ddt = tt1(2)-tt1(1)
    ts = ddt*100
    print *, ts
    if(ch==1) then
        le = 0.00218
        q = 1000000
        q0 = 50000
        k = 100
        cp = 1000
        rho = 8960
        theta = 1
        it(:) = 300
        ps = 0.0002
        s = 219
        open (12, file = 'source.txt')
        do i = 1,ti
            read (12,*) tt1(i), q2(i)
        end do
        do i = 1,ti
            if(i==1) then
                tt = tt1(i)
            else if(i>1) then
                tt = tt1(i) - tt1(i-1)
            end if
            call tdma(le, tt, q2(i), k, cp, rho, n, ts, theta, it, s, t_mat1, t_end)
            it = t_mat1
            t_mat2(i) = t_end
        end do
        print *, t_mat2
        open (28, file = 'temp_history_spherical.txt')
        write (28,*) t_mat2
        it(:) = 300
        do i = 1,ti
            t_end = t_mat2(i)
            if(i==1) then
                tt = tt1(i)
            else if(i>1) then
                tt = tt1(i) - tt1(i-1)
            end if
            call newt(q0, t_end, le, tt, k, cp, rho, n, ts, theta, it, s, ps)
            print *, s, q2(i), q0
            open(12, file = 'variable_q_1_spherical.txt')
            write(12, *) tt(i), q0
            call tdma(le, tt, q0, k, cp, rho, n, ts, theta, it, s, t_mat1, t_end)
            it = t_mat1
            q0 = 50000
        end do
    else if(ch==2) then
        le = 0.0025
        q = 1000000
        q0 = 50000
        k = 100
        cp = 1000
        rho = 8960
        theta = 1
        it(:) = 300
        ps = 0.0002
        s = 251
        open (10, file = 'source.txt')
        do i = 1,ti
            read (10,*) tt1(i), q2(i)
        end do
        do i = 1,ti
            if(i==1) then
                tt = tt1(i)
            else if(i>1) then
                tt = tt1(i) - tt1(i-1)
            end if
            call tdma(le, tt, q2(i), k, cp, rho, n, ts, theta, it, s, t_mat1, t_end)
            it = t_mat1
            t_mat2(i) = t_end
        end do
        print *, t_mat2
        open (12, file = 'temp_history_flatface.txt')
        write (12,*) t_mat2
        it(:) = 300
        do i = 1,ti
            t_end = t_mat2(i)
            if(i==1) then
                tt = tt1(i)
            else if(i>1) then
                tt = tt1(i) - tt1(i-1)
            end if
            call newt(q0, t_end, le, tt, k, cp, rho, n, ts, theta, it, s, ps)
            print *, s, q2(i), q0
            open(12, file = 'variable_q_1_flatface.txt')
            write(12, *) tt(i), q0
            call tdma(le, tt, q0, k, cp, rho, n, ts, theta, it, s, t_mat1, t_end)
            it = t_mat1
            q0 = 50000
        end do
        else if(ch==3) then
        le = 0.05
        q = 1000000
        q0 = 50000
        k = 50
        cp = 1000
        rho = 7800
        theta = 1
        it(:) = 300
        ps = 0.0002
        s = 31
        open (14, file = 'source.txt')
        do i = 1,ti
            read (14,*) tt1(i), q2(i)
        end do
        do i = 1,ti
            if(i==1) then
                tt = tt1(i)
            else if(i>1) then
                tt = tt1(i) - tt1(i-1)
            end if
            call tdma(le, tt, q2(i), k, cp, rho, n, ts, theta, it, s, t_mat1, t_end)
            it = t_mat1
            t_mat2(i) = t_end
        end do
        print *, t_mat2
        open (12, file = 'temp_history_aluminum.txt')
        write (12,*) t_mat2
        it(:) = 300
        do i = 1,ti
            t_end = t_mat2(i)
            if(i==1) then
                tt = tt1(i)
            else if(i>1) then
                tt = tt1(i) - tt1(i-1)
            end if
            call newt(q0, t_end, le, tt, k, cp, rho, n, ts, theta, it, s, ps)
            print *, s, q2(i), q0
            open(12, file = 'variable_q_1_aluminum.txt')
            write(12, *) tt(i), q0
            call tdma(le, tt, q0, k, cp, rho, n, ts, theta, it, s, t_mat1, t_end)
            it = t_mat1
            q0 = 50000
        end do
    end if
end subroutine

subroutine newt(q0, t_end, le, tt, k, cp, rho, n, ts, theta, it, s, ps)
    implicit none
    real :: q0, t_end, le, tt, k, cp, rho, n, ts, theta, s, ps, t_end0, count, er, q1
    real, dimension(n) :: t_mat1, it
    
    call tdma(le, tt, q0, k, cp, rho, n, ts, theta, it, s, t_mat1, t_end0)
    !print *,'Temperature at', s, 'node for q =', q0, 'is', t_end0
    count = 0
    do
        er = abs(t_end-t_end0)
        if (er <= 0.0001) then
            exit
        end if
        q1 = q0 - ((t_end0 - t_end)/ps)
        q0 = q1
        call tdma(le, tt, q0, k, cp, rho, n, ts, theta, it, s, t_mat1, t_end0)
        
        !print *, count ,q0, er, t_end0
        count = count + 1
    end do
    !it = t_mat2
    !print *, q0
    
end subroutine    

subroutine tdma(le, tt, q, k, cp, rho, n, ts, theta, it, s, t_mat1, t_end)
    implicit none
    
    real :: le, dx, tt, dt, q, k, cp, rho, n, ts, theta, t_end, alpha, lambda, i, j, m, s
    real, dimension(n) :: a, b, c, d, c_star, d_star, it, t_mat1
    real, dimension(n,ts) :: t_mat
    
    dx = le/(n-1)
    dt = tt/ts
    alpha = k/(cp*rho)
    lambda = (alpha*dt)/(dx**2)
    
    a(1) = 0
    b(1) = 1 + (2*lambda*theta)
    c(1) = -(2*lambda*theta)
    d(1) = (lambda*(1-theta)*2*it(1)) + ((1-(2*lambda*(1-theta)))*it(1)) + ((2*q*dt)/(rho*cp*dx))
    
    
    do i = 2,n-1
        a(i) = -(lambda*theta)
        b(i) = 1 + (2*lambda*theta)
        c(i) = -(lambda*theta)
        d(i) = (lambda*(1-theta)*it(i)*2) + ((1-(2*lambda*(1-theta)))*it(i))
    end do
        
    a(n) = -(2*lambda*theta)
    b(n) = 1 + (2*lambda*theta)
    c(n) = 0
    d(n) = (lambda*(1-theta)*2*it(n)) + ((1-(2*lambda*(1-theta)))*it(n))
    
    do j = 1,ts
        c_star(1) = c(1)/b(1)
        d_star(1) = d(1)/b(1)
                
        do i = 2,n
            m = 1/(b(i) - (a(i)*c_star(i-1)))
            c_star(i) = c(i)*m
            d_star(i) = (d(i) - (a(i)*d_star(i-1)))*m
        end do
            
        t_mat(n,j) = d_star(n)
        do i = n-1,1,-1
            t_mat(i,j) = d_star(i) - (c_star(i)*d(i+1))
        end do
            
        d(1) = (lambda*(1-theta)*2*t_mat(2,j)) + ((1-(2*lambda*(1-theta)))*t_mat(1,j)) + ((2*q*dt)/(rho*cp*dx))
        d(n) = (lambda*(1-theta)*2*t_mat(n-1,j)) + ((1-(2*lambda*(1-theta)))*t_mat(n,j))
        do i = 2,n-1
            d(i) = (lambda*(1-theta)*t_mat(i+1,j)) + ((1-(2*lambda*(1-theta)))*t_mat(i,j)) + (lambda*(1-theta)*t_mat(i-1,j))
        end do
    end do
    t_mat1(:) = t_mat(:,ts)
    t_end = t_mat(s,ts)
    
end subroutine


program iter1
    implicit none
    real :: n,ts,ch,ti
    print *, 'Type 1 for spherical, 2 for flatface, 3 for aluminum'
    read *, ch
    print *, 'Number of time intervals'
    read *, ti
    if(ch==1) then
        n = 219
        call main(n,ch,ti)
    else if(ch==2) then
        n = 251
        call main(n,ch,ti)
    else if(ch==3) then
        n = 251
        call main(n,ch,ti)
    end if
    pause
end program