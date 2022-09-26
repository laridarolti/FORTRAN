module calculation
    real, allocatable :: u_fnew(:, :)
    contains
    
recursive function Vcycle_2DPoisson(u_f,rhs,h) result (resV)
    implicit none
    real resV
    real,intent(inout):: u_f(:,:)  ! arguments
    real,intent(in)   :: rhs(:,:),h
    integer         :: nx,ny,nxc,nyc, i  ! local variables
    real,allocatable:: res_c(:,:),corr_c(:,:),res_f(:,:),corr_f(:,:)
    real            :: alpha=0.7, res_rms

    nx=size(u_f,1); ny=size(u_f,2)  ! must be power of 2 plus 1
    if( nx-1/=2*((nx-1)/2) .or. ny-1/=2*((ny-1)/2) ) &
         stop 'ERROR:not a power of 2'
    nxc=1+(nx-1)/2; nyc=1+(ny-1)/2  ! coarse grid size

    if (min(nx,ny)>5) then  ! not the coarsest level
        allocate(res_f(nx,ny),corr_f(nx,ny), &
            corr_c(nxc,nyc),res_c(nxc,nyc))
        !---------- take 2 iterations on the fine grid--------------
       res_rms = iteration_2DPoisson(u_f,rhs,h,alpha) 
       res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)
       !---------- restrict the residue to the coarse grid --------
       call residue_2DPoisson(u_f,rhs,h,res_f) 
       call restrict(res_f,res_c)
       !---------- solve for the coarse grid correction -----------
       corr_c = 0.  
       res_rms = Vcycle_2DPoisson(corr_c,res_c,h*2) ! *RECURSIVE CALL*
       !---- prolongate (interpolate) the correction to the fine grid 
       call prolongate(corr_c,corr_f)
       !---------- correct the fine-grid solution -----------------
       u_f = u_f - corr_f  
       !---------- two more smoothing iterations on the fine grid---
       res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)
       res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)
       deallocate(res_f,corr_f,res_c,corr_c)
    else  
       !----- coarsest level (ny=5): iterate to get 'exact' solution
       do i = 1,100
          res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)
       end do
    end if
    resV = res_rms   ! returns the rms. residue
end function Vcycle_2DPoisson



subroutine residue_2DPoisson(u_f, rhs, h, res_f)
implicit none
real, intent(in) :: u_f(:,:), rhs(:, :)
real, intent(in):: h
real, intent(out):: res_f(:, :)
integer :: nx, ny, i, j
! WRITE YOUR CODE
nx=size(u_f, 1)
ny=size(u_f, 2)
res_f(:, :)=0
do concurrent (i=2:nx-1, j=2:ny-1)
        res_f(i, j) = (u_f(i, j+1) + u_f(i, j-1) + u_f(i+1, j) + u_f(i-1, j)-4*u_f(i, j))/h**2 - rhs(i, j)
end do
! inputs needed or calculated are u, h, f
! WRITE YOUR CODE
end subroutine residue_2DPoisson



subroutine restrict (res_f, res_c)
implicit none
real, intent(in) :: res_f(:, :)
integer :: i, j, nx, ny, nxc, nyc, i_c, j_c ! indices for restricted grid
real, allocatable, intent(out) :: res_c (:, :)
! WRITE YOUR CODE
nx=size(res_f, 1)
ny=size(res_f, 2)
! need to test if nx, ny are powers of 2 + 1?
! nx=size(u_f,1); ny=size(u_f,2)  ! must be power of 2 plus 1
   ! if( nx-1/=2*((nx-1)/2) .or. ny-1/=2*((ny-1)/2) ) &
!         stop 'ERROR:not a power of 2'
nxc=1+(nx-1)/2; nyc=1+(ny-1)/2  ! coarse grid size
allocate (res_c(nxc, nyc))
do concurrent (i=1:nx, j=1:ny)
    if (mod(i, 2)==1 .and. mod(j, 2)==1) then
        i_c = (i+1)/2
        j_c = (j+1)/2
        res_c(i_c, j_c) = res_f(i, j)
    end if
end do       
! WRITE YOUR CODE
end subroutine restrict



subroutine prolongate (res_c, res_f)
implicit none
real, intent(in) :: res_c(:, :)
real, allocatable, intent(out) :: res_f(:, :)
integer :: nx, ny, nxc, nyc, i_c, j_c, i, j ! indices for coarse, respectively fine matrix
! WRITE YOUR CODE
nxc=size(res_c, 1)
nyc=size(res_c, 2)
nx= 1+ (nxc-1)*2
ny= 1+ (nyc-1)*2
allocate(res_f(nx, ny))

! STEP 1 : copying course pts into every other point of the grid
do concurrent (i_c=1:nxc, j_c=1:nyc)
    i = i_c*2-1
    j = j_c*2-1
    res_f(i, j) = res_c(i_c, j_c)
end do

! STEP 2 : for the points located between 2 known points from STEP 1 : calculate average
do concurrent (i=1:nx, j=1:ny)
    if ( mod(i, 2)==1 .and. mod(j, 2)==0) then
    res_f(i, j) = (res_f(i, j-1) +res_f(i, j+1))/2.
    end if
    if (mod(i, 2)==0 .and. mod(j, 2)==1) then
        res_f(i, j) = (res_f(i-1, j) +res_f(i+1, j))/2.
    end if
end do

! STEP 3 : for all the remaining pts, calculate average of all 4 points surrounding it
do concurrent (i=1:nx, j=1:ny)
    if ( mod(i, 2)==0 .and. mod(j, 2)==0) then
    res_f(i, j) = (res_f(i, j-1) +res_f(i, j+1)+res_f(i-1, j) +res_f(i+1, j))/4.
    end if
end do
! WRITE YOUR CODE
end subroutine prolongate 




function iteration_2DPoisson ( u_f, rhs, h, alpha) result (res_rms)
implicit none
real, intent(in) :: rhs(:,:), u_f(:, :)
real, intent(in):: h, alpha
real :: res_rms
integer :: nx, ny, i, j
real, allocatable :: res_f(:, :)
! WRITE YOUR CODE
nx=size(u_f, 1)
ny=size(u_f, 2)
allocate (res_f(nx, ny))
res_f(:, :)=0
u_fnew(:, :)=0

! ITERATION
do concurrent(i=2:nx-1, j=2:ny-1)
        res_f(i, j) = (u_f(i, j+1) + u_f(i, j-1) + u_f(i+1, j) + u_f(i-1, j)-4*u_f(i, j))/(h**2) - rhs(i, j)
        u_fnew(i, j)=u_f(i, j)+ alpha*res_f(i, j)*(h**2)/4
end do
! CALCULATING RMS RESIDUE
res_rms=0
do i=1, nx
         do j=1, ny
         res_rms = res_rms + (res_f(i, j))**2
         end do
end do
res_rms = sqrt(res_rms/(nx*ny))
! WRITE YOUR CODE
end function iteration_2DPoisson

    
    
end module calculation

    










program LDarolti_multigrid_method
    use calculation
    implicit none
    
    
    
    ! declaring variables
    integer :: nx, ny, i, j, niter=0
    logical :: multigrid
    real :: alpha, h, res_rms, f_rms=0., r
    real, allocatable :: init(:, :), rhs(:,:), u_f(:,:)
    ! reading control parameters using namelist
    namelist / inputs / nx, ny, init, alpha, multigrid
    open (1, file='relaxIN.txt', status='old')
    read (1, inputs)
    close(1)
    write (*, inputs)
    
    !nx=2049
    !ny=2049
    !alpha=0.7
    !multigrid = .true.
    
    
    
    ! initialising variables
    h=1./(ny-1)
    allocate (rhs(nx, ny), u_f(nx, ny), init(nx, ny),  u_fnew(nx, ny))
    rhs=init
    do concurrent(i=1:nx, j=1:ny)
        if (i==nx/2+1 .and. j==ny/2+1) then
            rhs(i, j)=1/h**2
        else
            rhs(i, j)=0
        end if
    end do
    !call random_number(rhs)
    u_f(:, :)=0
    print*, "niter", niter
    
    ! calculating f_rms
    do i=1, nx
        do j=1, ny
            f_rms = f_rms + (rhs(i, j))**2
        end do
    end do
    f_rms=sqrt(f_rms/(nx*ny))
    ! iterating until error is small enough
     print*, "f_rms", f_rms
     res_rms= iteration_2DPoisson ( u_f, rhs, h, alpha)
     print*, "res_rms", res_rms
     r =res_rms/f_rms
     print*, "r", r
    do while ( r>=1e-5)
        res_rms= iteration_2DPoisson ( u_f, rhs, h, alpha)
        r =res_rms/f_rms
        niter=niter+1
        u_f = u_fnew
        
        !if (mod(niter, 1000)==0) then
            !print*, "niter", niter
            !print*, "r", r
        !end if
    end do
    print*, niter
   ! PRINTING U
    open(2, file='outputHW5.txt', status='replace')
        do j=1, ny
                write(2, '(2049(1pe16.7, " "))') u_fnew(:, j)
        end do
        close(2)
    end program LDarolti_multigrid_method
