    program hw2ex3

    implicit none

    integer :: n, i
    real, allocatable :: y(:), d2ydx2(:), z(:), d2zdx2(:)
    real :: x, dx
    write (*, '(a, $)') 'Input numbert of grid points:'
    read*, n
    allocate (y(n), d2ydx2(n), z(n), d2zdx2(n)) ! grid arrays
    
    dx=10.0/(n-1) ! grid spacing
    
    do i=1, n
        x=(i-1)*dx
         y(i)=sin(x)
         z(i)=x**2
    end do
    
    call derivative2 (y, n, dx, d2ydx2)
    call derivative2 (z, n, dx, d2zdx2)
    
    do i=1, n
        x=(i-1)*dx
        print*, d2ydx2(i), -sin(x), -sin(x)-d2ydx2(i), d2zdx2(i), 2, 2-d2zdx2(i)
    end do
    
    deallocate(y, d2ydx2, z, d2zdx2)
    
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
    
        
    

    end program hw2ex3

