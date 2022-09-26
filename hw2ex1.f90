program task1

    implicit none
    
    ! declaring a string of 15 characters
    character (len=15) :: f
    
    ! declaring an integer parameter equal to 5
    integer, parameter :: g=5 
    
    ! declaring a one dimensional array with indices running from -1 to +10
    ! declaring a real number e
    ! declaring a one dimensional array with 100 elements
    real :: c(-1:10), e, l(1: 100), m 
    
    ! declaring a 4 dimensional allocatable array
    integer, allocatable :: d(:, :, :, :) 
    
    ! declaring integers a and b
    ! declaring the sum ( and initialising to 0) and i as integers
    integer :: a, b, sum=0, i, j  
    
    
    ! output is nearest integer to e
    print*, 'type a real number e'
    read*, e ! reading the real number e
    print *,'the nearest integer is', nint(e) 

    
    ! printing the remainder after a is  divided by b
    print*, 'type two integers a and b'
    read*, a, b ! reading a and b
    print*, 'the remainder of a/b is', mod(a, b) 
    
    
    ! calculating the sum of even integers from 12 to 124
    do i=12, 124, 2
        sum = sum+i
    end do
    print*, 'the sum of all even numbers from 12 to 124 is', sum
    
    
   ! testing all elements of an array, if the element is positive then exit loop
    do j=1, 100
        print*, 'type element number', j, ' into the vector'
        read*, m
        l(j)=m
        if ( l(j)>0) then
         print*, 'number is positive so must exit loop'
         exit
        end if
    end do
end program task1