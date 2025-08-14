program test_fft_derivative
  use rd_spectral_fftw
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
  integer, parameter :: dp = selected_real_kind(15,300)
  type(rd_spec_plan) :: P
  integer :: Nx, Ny, i, j, nxc
  real(dp) :: Lx, Ly, two_pi, dx, dy
  real(dp), allocatable :: u_orig(:,:), dudx_spec(:,:), dudx_fd(:,:)
  real(dp) :: err_max

  Nx = 64; Ny = 64
  Lx = 1.0_dp; Ly = 1.0_dp
  call init_rd_spec_plan(P, Nx, Ny, Lx, Ly, .false.)
  nxc = Nx/2 + 1

  allocate(u_orig(Ny,Nx))
  allocate(dudx_spec(Ny,Nx))
  allocate(dudx_fd(Ny,Nx))

  two_pi = 2.0_dp * acos(-1.0_dp)
  dx = Lx / real(Nx,dp)
  dy = Ly / real(Ny,dp)

  ! Fill u with a smooth periodic pattern
  do j = 1, Ny
     do i = 1, Nx
        P%u(j,i) = sin(two_pi*real(i-1,dp)/Nx) + 0.37_dp*cos(two_pi*real(j-1,dp)/Ny)
        u_orig(j,i) = P%u(j,i)
     end do
  end do

  ! Forward FFT
  call fftw_execute_dft_r2c(P%plan_r2c_u, P%u, P%uhat)

  ! Spectral derivative in x: multiply by i*kx
  do j = 1, Ny
     do i = 1, nxc
        P%uhat(j,i) = (0.0_dp,1.0_dp) * real(P%kx(i),dp) * P%uhat(j,i)
     end do
  end do

  ! Inverse FFT to get dudx
  call fftw_execute_dft_c2r(P%plan_c2r_u, P%uhat, dudx_spec)
  dudx_spec = dudx_spec * P%invN

  ! Finite difference derivative in x for comparison
  do j = 1, Ny
     do i = 1, Nx
        dudx_fd(j,i) = (u_orig(j,mod(i,Nx)+1) - u_orig(j,mod(i-2,Nx)+1))/(2.0_dp*dx)
     end do
  end do

  ! Maximum absolute difference
  err_max = maxval(abs(dudx_spec - dudx_fd))
  print *, "Spectral vs FD derivative max error =", err_max

  call destroy_rd_spec_plan(P)
end program test_fft_derivative

