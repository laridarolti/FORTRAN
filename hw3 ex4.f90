module calculation
    
        integer :: i
        real :: x, dx, dy
        real, allocatable :: a(:, :), y(:, :), adoubleprime(:, :)    
        
    contains
    
    
    subroutine derivative2(a, nx, ny, dx, dy, adoubleprime)
    integer, intent(in) :: nx, ny
    real, intent(in) :: a(nx, ny), dx, dy
    real, intent(out):: adoubleprime (nx, ny)
    integer :: i
    
    do i=1, nx
        do j=1, ny
            if (i==1 .or. j==1 .or. i==nx .or. j++ny) then
                adoubleprime (i, j) = 0
            else
                adoubleprime (i, j) = (a(i-1, j)-2*a(i, j)+a(i+1, j))/dx**2 + (a(i, j-1)-2*a(i, j)+a(i, j+1))/dy**2
            end if
        end do
    end do
 
     end subroutine derivative2
    
    end module calculation
    
    
    
    
    
    
    program hw3ex4
    
    use calculation

    implicit none
    integer:: j, nsteps, p, h
    real :: dt, deltax, deltay, time=0.
    integer :: nx, ny
    real :: L, total_time, a_ct, k=1
    real, allocatable :: T(:,:), d2Tdx2(:, :)
    
    ! reading control parameters using namelist
    namelist / inputs/ nx, ny, L, total_time, a_ct
    open (1, file='input_ex3.dat', status='old')
    read (1, inputs)
    close(1)
    write (*, inputs)
    
    allocate (T(nx, ny), d2Tdx2(nx, ny))
    call RANDOM_NUMBER(T)
    
    ! variables
    k=1
    deltax=L/(nx-1)
    deltay=L/(ny-1)
    dt=a_ct*deltax**2/k
    nsteps=total_time/dt
    
    ! calculating second derivative
    
    call derivative2(T, nx, ny, deltax, deltay, d2Tdx2)
    
    
    ! loop
    ! updating T
    
        do h= 1, nsteps
            if ( time<=total_time) then
               do i=p, nx
                   do j=1, ny
                       if(p==1 .or. j==1 .or. p==nx .or. j++ny) then
                           T (p, j) = 0
                       else
                           T(p, j)=T(p, j)+ dt*k*d2Tdx2(p, j)
                       end if
                   end do
               end do
            end if
        time=time+dt
        end do
        
 ! printing T       
        do p=1, nx
            do j=1, ny
            open(2, file='output.dat')
            write(2, *) T(p, j)
            close(2) 
            end do
    end do
    
    end program hw3ex4
    
        
    