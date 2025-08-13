!=========================
! File: reaction_module.f90
!=========================
module reaction_diffusion_mod
  implicit none
  private
  public :: initialize_fields, solve_rd_finite_diff, solve_rd_fullTime, write_fields
  public :: check_conditions

contains

  subroutine check_conditions(a, b, c, d, dt, Lx, Ly, Nx, Ny, model_flag)
    implicit none
    real(8), intent(in) :: a, b, c, d, dt
    real(8), intent(in) :: Lx, Ly
    integer, intent(in) :: Nx, Ny, model_flag
    real(8) :: apb, diff, CFL, dx, dy
    real(8) :: uo, fu, fv, gu, gv
    

    dx = Lx / dble(Nx)
    dy = Ly / dble(Ny)
    CFL = dt / (d * dx * dx)

    print *, 'CFL = ', CFL



    if(model_flag == 1)then
      apb = a + b
      diff = b - a

      fu = (b-a)/(a+b)
      fv = (a+b)**2
      gu = -2*b/(a+b)
      gv = -(a+b)**2

      print *, '------------Insatibility Conditions------------------------'
      print '(A34, L1, A, F8.4, A, F8.4)', 'Condition 1 (fu + gv < 0):  ',  &
           (fu + gv < 0), ' --> fu = ', fu, ', gv = ', gv
      
      print '(A34, L1, A, F8.4, A, F8.4)', 'Condition 2 (fu*gv - fv*gu > 0):    ',  &
           (fu*gv - fv*gu > 0), ' --> fu*gv = ', fu*gv, ', fv*gu = ', fv*gu
      
      print '(A34, L1, A, F8.4, A, F8.4)', 'Condition 3 (d*fu + gv > 0):  ',  &
           (d*fu + gv > 0), ' --> d*fu = ', d*fu, ', gv = ', gv
      
      print '(A34, L1, A, F8.4)', 'Condition 4 (Discriminant > 0):  ', &
           ((d*fu + gv)**2 - 4*d*(fu*gv - fv*gu) > 0), ' --> Discriminant = ',  &
           (d*fu + gv)**2 - 4*d*(fu*gv - fv*gu)
      print *, '-----------------------------------------------------------'

    else if(model_flag == 2)then

      uo = (a+c)/b
      fu = - b + 2*c/uo 
      fv = - c**2/uo**2
      gu =  2 * uo
      gv = -c


      print *, '------------Insatibility Conditions------------------------'
      print '(A34, L1, A, F8.4, A, F8.4)', 'Condition 1 (fu + gv < 0):  ',  &
           (fu + gv < 0), ' --> fu = ', fu, ', gv = ', gv
      
      print '(A34, L1, A, F8.4, A, F8.4)', 'Condition 2 (fu*gv - fv*gu > 0):    ',  &
           (fu*gv - fv*gu > 0), ' --> fu*gv = ', fu*gv, ', fv*gu = ', fv*gu
      
      print '(A34, L1, A, F8.4, A, F8.4)', 'Condition 3 (d*fu + gv > 0):  ',  &
           (d*fu + gv > 0), ' --> d*fu = ', d*fu, ', gv = ', gv
      
      print '(A34, L1, A, F8.4)', 'Condition 4 (Discriminant > 0):  ', &
           ((d*fu + gv)**2 - 4*d*(fu*gv - fv*gu) > 0), ' --> Discriminant = ',  &
           (d*fu + gv)**2 - 4*d*(fu*gv - fv*gu)
      print *, '-----------------------------------------------------------'


    else
      print*, "Wrong model type"
    
    end if


  end subroutine check_conditions


  subroutine initialize_fields(u, v, x, y, Nx, Ny, Lx, Ly, a, b, c,  model_flag)
    implicit none
    integer, intent(in) :: Nx, Ny, model_flag
    real(8), intent(in) :: Lx, Ly, a, b, c
    real(8), intent(out) :: u(Ny, Nx), v(Ny, Nx)
    real(8), intent(out) :: x(Nx), y(Ny)
    real(8) :: dx, dy
    integer :: ix, iy

    dx = Lx / dble(Nx)
    dy = Ly / dble(Ny)

    call random_seed()
    call random_number(u)
    call random_number(v)

    do ix = 1, Nx
      x(ix) = (ix - 1) * dx
    end do
    do iy = 1, Ny
      y(iy) = (iy - 1) * dy
    end do

!    u = a + b + 0.02d0 * (u - 0.5d0)
!    v = b / (a + b)**2 + 0.02d0 * (v - 0.5d0)

    if (model_flag == 1) then
      ! Initialization for Original Model
      u = a + b + 0.02d0 * (u - 0.5d0)
      v = b / (a + b)**2 + 0.02d0 * (v - 0.5d0)
    
    else if (model_flag == 2) then
      ! Initialization for Gierer–Meinhardt Model
      u = (a + c)/b + 0.02d0 * (u - 0.5d0)
      v = u**2 / c + 0.02d0 * (v - 0.5d0)
    
    else
      print *, "Unknown model_flag during initialization:", model_flag
      stop
    end if

!    do ix = 1, Nx
!      do iy = 1, Ny
!        u(ix,iy) = a + b + 0.01 * cos(2*x(ix)) * cos(2*y(iy))
!        v(ix,iy) = b / (a + b)**2 
!      end do
!    end do

  end subroutine initialize_fields

  subroutine solve_rd_fullTime(u, v, Nx, Ny, Lx, Ly, Nt, dt, a, b, d, r, gamm0, rate, plot_every)
    implicit none
    integer, intent(in) :: Nx, Ny, Nt, plot_every
    real(8), intent(in) :: Lx, Ly, dt, a, b, d, r, gamm0, rate
    real(8), intent(inout) :: u(Ny, Nx), v(Ny, Nx)

    real(8) :: dx, dy, gamm, f_u, g_v
    real(8), dimension(Ny, Nx) :: u_lap, v_lap
    integer :: it, ix, iy, im, ip, jm, jp
    character(len=40) :: filename1, filename2

    real(8) :: apb, diff, CFL

    dx = Lx / dble(Nx)
    dy = Ly / dble(Ny)


!    call check_conditions(a, b, d, dt, Lx, Ly, Nx, Ny)


    gamm = gamm0

    do it = 1, Nt
      write(*,'(A,F8.5,A,I0)', advance='no') achar(13)//'gamma: ', gamm, ', Step: ', it

      do iy = 1, Ny
        jm = mod(iy - 2 + Ny, Ny) + 1
        jp = mod(iy     + Ny, Ny) + 1
        do ix = 1, Nx
          im = mod(ix - 2 + Nx, Nx) + 1
          ip = mod(ix     + Nx, Nx) + 1

          u_lap(iy,ix) = (u(iy,im) + u(iy,ip) + u(jm,ix) + u(jp,ix) - &
                          4.d0*u(iy,ix)) / dx**2

          v_lap(iy,ix) = (v(iy,im) + v(iy,ip) + v(jm,ix) + v(jp,ix) - &
                          4.d0*v(iy,ix)) / dx**2
        end do
      end do

      do iy = 1, Ny
        do ix = 1, Nx
          f_u = a - u(iy,ix) + u(iy,ix)**2 * v(iy,ix)
          g_v = b - u(iy,ix)**2 * v(iy,ix)

          u(iy,ix) = u(iy,ix) + dt * (gamm*f_u - r*u(iy,ix) + &
                      exp(-2.d0*r*it*dt)*u_lap(iy,ix))

          v(iy,ix) = v(iy,ix) + dt * (gamm*g_v - r*v(iy,ix) + &
                      d*exp(-2.d0*r*it*dt)*v_lap(iy,ix))
        end do
      end do

      gamm = gamm + dt * rate

      if (mod(it, plot_every) == 0) then
        call write_fields(u, v, Nx, Ny, it)
      end if
    end do

    print *
    print *, 'Simulation complete.'
  end subroutine solve_rd_fullTime


  subroutine solve_rd_finite_diff(u, v, Nx, Ny, Lx, Ly, a, b, c, d, r, &
      growth_dist, dt, gamm, it, model_flag)
     
      implicit none
      integer, intent(in) :: Nx, Ny, it, model_flag
      real(8), intent(in) :: Lx, Ly, a, b, c, d, r, dt, gamm
      real(8), intent(inout) :: u(Ny, Nx), v(Ny, Nx)
      character(100), intent(in) :: growth_dist
      real(8) :: dx, dy
      real(8), dimension(Ny, Nx) :: u_lap, v_lap
      real(8) :: f_u, g_v, growth_rate
      integer :: ix, iy, im, ip, jm, jp
    
      dx = Lx / dble(Nx)
      dy = Ly / dble(Ny)
    
      ! Compute Laplacians with periodic boundaries
      do iy = 1, Ny
        jm = mod(iy - 2 + Ny, Ny) + 1
        jp = mod(iy     + Ny, Ny) + 1
        do ix = 1, Nx
          im = mod(ix - 2 + Nx, Nx) + 1
          ip = mod(ix     + Nx, Nx) + 1
    
          u_lap(iy,ix) = (u(iy,im) + u(iy,ip) + u(jm,ix) + u(jp,ix) & 
            - 4.d0*u(iy,ix)) / dx**2
          v_lap(iy,ix) = (v(iy,im) + v(iy,ip) + v(jm,ix) + v(jp,ix) & 
            - 4.d0*v(iy,ix)) / dx**2
        end do
      end do
    
      ! Reaction-Diffusion update
      do iy = 1, Ny
        do ix = 1, Nx
          if (model_flag == 1) then
            ! Schnakenberg model
            f_u = a - u(iy,ix) + u(iy,ix)**2 * v(iy,ix)
            g_v = b - u(iy,ix)**2 * v(iy,ix)
          else if (model_flag == 2) then
            ! Gierer–Meinhardt model
            if (v(iy,ix) > 1.d-12) then
              f_u = a - b * u(iy,ix) + (u(iy,ix)**2) / v(iy,ix)
              g_v = (u(iy,ix)**2)  - c * v(iy,ix)
            else
              f_u = 0.d0
              g_v = 0.d0
            end if
          else
            print *, 'Error: Unknown model_flag =', model_flag
            stop
          end if



          if (trim(growth_dist) == 'uniform') then
            growth_rate = r
          else
            growth_rate = r * u(iy, ix)
          end if

    
          u(iy,ix) = u(iy,ix) + &
            dt * (gamm*f_u - growth_rate*u(iy,ix) + & 
            exp(-2.d0*growth_rate*it*dt)*u_lap(iy,ix))

          v(iy,ix) = v(iy,ix) + &
            dt * (gamm*g_v - growth_rate*v(iy,ix) + & 
            d*exp(-2.d0*growth_rate*it*dt)*v_lap(iy,ix))

        end do
      end do
  end subroutine solve_rd_finite_diff


  subroutine write_fields(u, v, Nx, Ny, it)
    implicit none
    integer, intent(in) :: Nx, Ny, it
    real(8), intent(in) :: u(Ny, Nx), v(Ny, Nx)
    character(len=40) :: filename1, filename2
    integer :: ix, iy

    write(filename1, '(A,I0,A)') 'data/u_', it, '.dat'
    write(filename2, '(A,I0,A)') 'data/v_', it, '.dat'

    open(unit=10, file=filename1, status='replace')
    open(unit=11, file=filename2, status='replace')
    do iy = 1, Ny
      write(10,*) (u(iy,ix), ix=1,Nx)
      write(11,*) (v(iy,ix), ix=1,Nx)
    end do
    close(10)
    close(11)
  end subroutine write_fields

end module reaction_diffusion_mod

