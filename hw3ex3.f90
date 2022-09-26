module calculation
    
        integer :: i
        real :: x, dx
        real, allocatable :: a(:), y(:), d2ydx2(:)    
        
    contains
    
    
    subroutine derivative2(a, np, h, adoubleprime)
    integer, intent(in) :: np
    real, intent(in) :: a(np), h
    real, intent(out):: adoubleprime (np)
    integer :: i
    
    do i=2, np-1
        adoubleprime (i) = (a(i+1)-2*a(i)+a(i-1))/h**2
    end do
     adoubleprime (1) = 0
     adoubleprime (np) = 0
     end subroutine derivative2
    
    end module calculation
    
    
    
    
    
    
    program hw3ex3
    
    use calculation

    implicit none
    integer:: j, nsteps, N, p
    real :: dt, deltax, time=0.
    real :: L, total_time, a_ct, k=1
    real, allocatable :: T(:), d2Tdx2(:)
    
    ! reading control parameters using namelist
    namelist / inputs/ N, L, total_time, a_ct
    open (1, file='input_ex3.dat', status='old')
    read (1, inputs)
    close(1)
    write (*, inputs)
    
    allocate (T(N), d2Tdx2(N))
    call RANDOM_NUMBER(T)
    
    ! variables
    k=1
    deltax=L/(N-1)
    dt=a_ct*deltax**2/k
    nsteps=total_time/dt
    
    ! calculating second derivative
    
    call derivative2(T, N, deltax, d2Tdx2)
    
    
    ! loop
    ! updating T
    
        do j= 1, nsteps
            if ( time<=total_time) then
        T(1) = 0
        T(N) = 0
        do p=2, N-1
        T(p)=T(p)+ dt*k*d2Tdx2(p)
        end do
         end if
        time=time+dt
        end do
        
 ! printing T       
        do p=1, N
            open(2, file='output.dat')
            write(2, *) T(p)
            close(2) 
        end do
        
    
            
    
    
    end program hw3ex3
    

