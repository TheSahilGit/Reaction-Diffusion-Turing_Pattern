program main_sim
  use reaction_diffusion_mod
  implicit none

  integer, parameter :: Nx = 128, Ny = 128, Nt = 5000000, plot_every = 1000
  real(8), parameter :: pi = acos(-1.0d0)
  real(8), parameter :: Lx = 2*pi, Ly = 2*pi, dt = 2.5d-5
  real(8) :: a, b, c, d, r  
  real(8), parameter :: gamm0 = 10.0d0, gamm_rate = 0.0d0
  integer :: it
  integer, parameter :: model_flag = 2

  real(8), dimension(Ny,Nx) :: u, v
  real(8), dimension(Nx) :: x
  real(8), dimension(Ny) :: y
  character(100) :: growth_dist



  
  if(model_flag == 1)then
    print*, 'Model type: ', 'Schnakenberg model'
    a = 0.1d0; b = 0.9d0; c = 0.0d0; d = 20.0d0; r = 0.00d0

  elseif(model_flag == 2) then 
    print*, 'Model type: ', 'Gierer–Meinhardt model model'
    a = 1.12d0; b = 6.67d0; c = 4.45d0;  d = 10.0d0; r = 0.0d0

  end if


  growth_dist = 'uniform'

  print*, 'Growth rate: ', r, "Growth dist : ", trim(growth_dist)

  !growth_dist = 'uniform'


  call check_conditions(a, b, c, d, dt, Lx, Ly, Nx, Ny, model_flag)
  call initialize_fields(u, v, x, y, Nx, Ny, Lx, Ly, a, b, c, model_flag)
!  call solve_rd_fullTime(u, v, Nx, Ny, Lx, Ly, Nt, dt, a, b, d, r, gamm0, gamm_rate, plot_every)

  do it = 1, Nt
    call solve_rd(u, v, Nx, Ny, Lx, Ly, a, b, c, d, r, growth_dist, dt, gamm0, it, model_flag)

    write(*,'(A,I0,A,I0)', advance='no') achar(13)//'Step: ', it, '/', Nt

    if (mod(it, plot_every) == 0) then
      call write_fields(u, v, Nx, Ny, it)
    end if


  end do

   
end program main_sim

