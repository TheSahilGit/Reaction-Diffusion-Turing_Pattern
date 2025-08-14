program main_sim
  use reaction_diffusion_mod
  use rd_spectral_fftw

  implicit none

  integer, parameter :: Nx = 128, Ny = 128, Nt = 5000000, plot_every = 1000
  real(8), parameter :: pi = acos(-1.0d0)
  real(8), parameter :: Lx = 2*pi, Ly = 2*pi, dt = 2.5d-5
  real(8) :: a, b, c, d, r  
  real(8), parameter :: gamm0 = 10.0d0, gamm_rate = 0.0d0
  integer :: it
  integer, parameter :: model_flag = 1
  integer, parameter :: solve_method = 1

  real(8), dimension(Ny,Nx) :: u, v
  real(8), dimension(Nx) :: x
  real(8), dimension(Ny) :: y
  character(100) :: growth_dist

  type(rd_spec_plan) :: P

  



  
  if(model_flag == 1)then
    print*, 'Model type: ', 'Schnakenberg model'
    a = 0.1d0; b = 0.9d0; c = 0.0d0; d = 10.0d0; r = 0.00d0

  elseif(model_flag == 2) then 
    print*, 'Model type: ', 'Gierer–Meinhardt model'
    a = 1.12d0; b = 6.67d0; c = 4.45d0;  d = 10.0d0; r = 0.0d0

  end if


  growth_dist = 'uniform'

  print*, 'Growth rate: ', r, "Growth dist : ", trim(growth_dist)

  !growth_dist = 'uniform'


  call check_conditions(a, b, c, d, dt, Lx, Ly, Nx, Ny, model_flag)
  call initialize_fields(u, v, x, y, Nx, Ny, Lx, Ly, a, b, c, model_flag)


  if(solve_method == 1)then
    print*, "Numerical Method :  Finite Difference" 

  elseif(solve_method == 2)then
    print*, "Numerical Method :  Spectral"

    call init_rd_spec_plan(P, Nx, Ny, Lx, Ly)
    P%u = u
    P%v = v

  end if



  do it = 1, Nt

    if(solve_method == 1)then
      call solve_rd_finite_diff(u, v, Nx, Ny, Lx, Ly, a, b, c, d, r, & 
        growth_dist, dt, gamm0, it, model_flag)

    elseif(solve_method == 2)then
      call rd_step_spectral_uniform_growth(P, a, b, c, d, r, gamm0, dt, & 
        it*dt, model_flag)
      u = P%u
      v = P%v
    end if

    write(*,'(A,I0,A,I0)', advance='no') achar(13)//'Step: ', it, '/', Nt

    if (it == 1 .or. mod(it, plot_every) == 0) then
      call write_fields(u, v, Nx, Ny, it)
    end if


  end do

  if (solve_method == 2) then
    call destroy_rd_spec_plan(P)
  end if

   
end program main_sim

