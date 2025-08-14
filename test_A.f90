program test_fft_roundtrip
  use rd_spectral_fftw
  implicit none
  integer, parameter :: dp = selected_real_kind(15,300)
  type(rd_spec_plan) :: P
  integer :: Nx, Ny, i, j
  real(dp) :: Lx, Ly, err_max, two_pi
  real(dp), allocatable :: u_orig(:,:)

  ! grid size
  Nx = 64; Ny = 64
  Lx = 1.0_dp; Ly = 1.0_dp

  ! initialize plan (no de-aliasing for this simple check)
  call init_rd_spec_plan(P, Nx, Ny, Lx, Ly, .false.)

  allocate(u_orig(Ny,Nx))

  two_pi = 2.0_dp * acos(-1.0_dp)

  ! fill P%u with a test pattern (nontrivial)
  do j = 1, Ny
     do i = 1, Nx
        P%u(j,i) = sin(two_pi * real(i-1,dp) / real(Nx,dp)) + 0.37_dp * cos(two_pi * real(j-1,dp) / real(Ny,dp))
        u_orig(j,i) = P%u(j,i)
     end do
  end do

  ! Use the module's step routine with no reaction (gamm=0) and dt=0 so the spectral
  ! transform and inverse happen but nothing changes; this exercises the forward+inverse.
  call rd_step_spectral_uniform_growth(P, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1)

  ! compute max abs difference between original and after roundtrip
  err_max = 0.0_dp
  do j = 1, Ny
     do i = 1, Nx
        err_max = max(err_max, abs(P%u(j,i) - u_orig(j,i)))
     end do
  end do

  print*,  'FFT roundtrip max error = ', err_max
  ! expected err_max ~ 1e-12 ... 1e-14 (double precision)

  call destroy_rd_spec_plan(P)
end program test_fft_roundtrip

