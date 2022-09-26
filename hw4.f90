module calculation
    
        integer, private :: i, nx, ny
        real :: dx, dy
        
    contains
    
    ! subroutine for second derivative
    subroutine derivative2(T, nx, ny, dx, dy, del2T)
    integer, intent(in) :: nx, ny
    real, intent(in) :: T(nx, ny), dx, dy
    real, intent(out):: del2T (nx, ny)
    integer :: i
    
    do i=1, nx
        do j=1, ny
            if (i==1 .or. j==1 .or. i==nx .or. j==ny) then
                del2T (i, j) = 0
            else
                del2T (i, j) = (T(i-1, j)-2*T(i, j)+T(i+1, j))/dx**2 + (T(i, j-1)-2*T(i, j)+T(i, j+1))/dy**2
            end if
        end do
    end do
 
    end subroutine derivative2
    
    
    
    ! subroutine for vx and vy velocities
    subroutine x_and_y_velocities(S, nx, ny, dx, dy, vx, vy)
    integer, intent(in) :: nx, ny
    real, intent(in) :: S(nx, ny), dx, dy
    real, intent(out):: vx (nx, ny), vy (nx, ny)
    integer:: i
    
    do i=1, nx
        do j=1, ny
            if (i==1 .or. j==1 .or. i==nx .or. j==ny) then
                vx (i, j) = 0
                vy (i, j) = 0
            else
                vx (i, j) = (S(i, j+1)-S(i, j-1))/(2*dy)
                vy (i, j) = (S(i-1, j)-S(i+1, j))/(2*dx)
            end if
        end do
    end do
 
    end subroutine x_and_y_velocities
    
    
    ! subroutine for calculating vgradT
   subroutine vgradientT( nx, ny, dx, dy, vx, vy, T, vgradT)
    integer, intent(in) :: nx, ny
    real, intent(in) :: vx(nx, ny), vy(nx, ny), T(nx, ny), dx, dy
    real, intent(out):: vgradT (nx, ny)
    integer:: i
    real :: vx_dTdx(nx, ny), vy_dTdy(nx, ny)
   
    do i=1, nx
        do j=1, ny
            if (i==1 .or. j==1 .or. i==nx .or. j==ny) then
                vx_dTdx(i, j)=0
                vy_dTdy(i, j)=0
            else
            ! calculating vx_dTdx
                      if(vx(i, j)>0) then
                          vx_dTdx(i, j) = vx(i, j) * (T(i, j)-T(i-1, j))/dx
                      else if (vx(i, j)<0) then
                          vx_dTdx(i, j) = vx(i, j) * (T(i+1, j)-T(i, j))/dx
                      else
                          vx_dTdx(i, j) = 0
                      end if
                      
                      ! calculating vy_dTdy
                      if(vy(i, j)>0) then
                          vy_dTdy(i, j)=vy(i, j)*(T(i, j)-T(i, j-1))/dy
                      else if (vy(i, j)<0) then
                          vy_dTdy(i, j)=vy(i, j)*(T(i, j+1)-T(i, j))/dy
                      else
                         vy_dTdy(i, j)=0
                      end if
            end if
            ! calculating vgradT
            vgradT(i, j) = vx_dTdx(i, j) + vy_dTdy(i, j)
        end do
    end do
 
    end subroutine vgradientT
    
    end module calculation
    
    
    
    
    
    
    program hw4
    
    use calculation

    implicit none
    integer:: j, nsteps, i, nx=100, ny=100, Numb_of_iter
    real :: h, time=0.,  L, total_time, k, B, a_dif, a_adv, pi, deltax, deltay, dt_dif,dt_adv,  vx_max, vy_max, dt
    real, allocatable :: T(:,:), d2Tdx2(:, :), S(:, :), vx(:, :), vy(:, :), del2T(:, :), vgradT(:, :), vx_dTdx(:, :), vy_dTdy(:, :)
    
    ! reading control parameters using namelist
    namelist / inputs / nx, ny, B, k, a_adv, a_dif, total_time
    open (1, file='AdvDif_parameters.dat', status='old')
    read (1, inputs)
    close(1)
    write (*, inputs)
    
   allocate (T(nx, ny), d2Tdx2(nx, ny), S(nx, ny), vx(nx, ny), vy(nx, ny), del2T(nx, ny), vgradT(nx, ny), vx_dTdx(nx, ny), vy_dTdy(nx, ny))
   call RANDOM_NUMBER(T)
    
    ! variables
    !h=L/(ny-1)
    L=1
    deltax=L/(nx-1)
    deltay=L/(ny-1)
    dt_dif=a_dif*(min(deltax, deltay))**2/k
    print*, 'deltax=', deltax
    print*, 'deltay=', deltay
    print*, 'dt_dif=', dt_dif
    
    pi=DACOS(-1.D0) ! defining pi
    
    ! defining the streamfunction S
    do i=1, nx
        do j=1, ny
            S(i, j) = B*sin(pi*(i-1)/(nx-1))*sin(pi*(j-1)/(ny-1))
        end do
    end do
      
      ! calculativg vx and vy from S
      call x_and_y_velocities(S, nx, ny, deltax, deltay, vx, vy)
      vx_max=maxval(vx)
      vy_max=maxval(vy)
      
      dt_adv=a_adv*min((deltax/vx_max), (deltay/vy_max))
      dt= min(dt_dif, dt_adv)
      print*, 'dt_adv=', dt_adv
      print*, 'dt=', dt
      print*, 'vx_max', vx_max
      print*, 'vy_max', vy_max
      
      nsteps=total_time/dt
      print*,'nsteps=', nsteps
    !Numb_of_iter=0
    ! loop ! updating T
            do while(time<=total_time)
            call derivative2(T, nx, ny, deltax, deltay, del2T)
            call vgradientT( nx, ny, deltax, deltay, vx, vy, T, vgradT)
            ! calculating T
            T=T+ dt*(k*del2T - vgradT) 
            T (:, ny) = 0
            T(:, 1) = 1
            T(1, :)=T(2, :)
            T(nx, :)=T(nx-1, :)
            time=time+dt
            !print*, 'time=', time
            !Numb_of_iter=Numb_of_iter+1
            !if(mOd(Numb_of_iter, 20000)==0) then
            !print*, 'Numb_of_iter=', Numb_of_iter
            end if
            end do
    
 ! printing T       
        open(2, file='outputhw3ex4_Bis1.dat', status='replace')
        do j=1, ny
                write(2, '(100(1pe12.6, " "))') T(:, j)
        end do
        close(2)
   
    end program hw4
    
