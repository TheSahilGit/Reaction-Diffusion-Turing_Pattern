!============================================================
! FFTW3 pseudo-spectral reaction–diffusion solver (Fortran)
!   - Periodic BCs (handled naturally by FFTs)
!   - Semi-implicit split: backward-Euler for diffusion, explicit for reactions
!   - Supports Schnakenberg (model_flag=1) and Gierer–Meinhardt (model_flag=2)
!   - Uniform growth only: growth_rate = r (constant)
!   - Diffusion scaling matches your code: e^{-2 r t} factor in ∇^2 term
!   - Implements 2/3 de-aliasing for nonlinear terms (improves accuracy)
!   - Single step routine `rd_step_spectral_uniform_growth`
!   - Uses FFTW3 Fortran 2003 interface via include 'fftw3.f03' (no .mod needed)
!============================================================

module rd_spectral_fftw
  use, intrinsic :: iso_c_binding
  implicit none
  private
  public :: rd_spec_plan, init_rd_spec_plan, destroy_rd_spec_plan
  public :: rd_step_spectral_uniform_growth

  integer, parameter :: dp = selected_real_kind(15, 300)

  type rd_spec_plan
     integer :: Nx = 0, Ny = 0
     real(dp) :: Lx = 0.0_dp, Ly = 0.0_dp
     real(dp) :: dx = 0.0_dp, dy = 0.0_dp
     ! wave numbers (only kx >= 0 stored in r2c layout)
     real(dp), allocatable :: k2(:, :)         ! size: Ny x (Nx/2+1)
     real(dp), allocatable :: dealias_mask(:, :) ! same shape as k2, 1 or 0
     logical :: use_dealias = .true.

     ! FFTW pointers
     type(c_ptr) :: plan_r2c_u = c_null_ptr
     type(c_ptr) :: plan_c2r_u = c_null_ptr
     type(c_ptr) :: plan_r2c_v = c_null_ptr
     type(c_ptr) :: plan_c2r_v = c_null_ptr

     ! real-space work arrays
     real(dp), allocatable :: u(:, :), v(:, :), rhs_u(:, :), rhs_v(:, :)

     ! spectral work arrays (FFTW r2c shape)
     complex(c_double_complex), allocatable :: uhat(:, :), vhat(:, :), rhshat_u(:, :), rhshat_v(:, :)

     ! normalization factor after inverse FFT (FFTW is unnormalized by default)
     real(dp) :: invN = 1.0_dp
  end type rd_spec_plan

contains

  subroutine init_rd_spec_plan(P, Nx_in, Ny_in, Lx_in, Ly_in, use_dealias_in)
    ! Initialize spectral plan and precompute k^2 and de-alias mask
    include 'fftw3.f03'
    type(rd_spec_plan), intent(inout) :: P
    integer, intent(in) :: Nx_in, Ny_in
    real(dp), intent(in) :: Lx_in, Ly_in
    logical, intent(in), optional :: use_dealias_in
    integer :: nxc
    integer :: i, j
    real(dp) :: dkx, dky
    integer :: kx_idx, ky_idx

    if (present(use_dealias_in)) then
       P%use_dealias = use_dealias_in
    else
       P%use_dealias = .true.
    end if

    P%Nx = Nx_in; P%Ny = Ny_in
    P%Lx = Lx_in; P%Ly = Ly_in
    P%dx = Lx_in / real(Nx_in, dp)
    P%dy = Ly_in / real(Ny_in, dp)
    P%invN = 1.0_dp / real(Nx_in*Ny_in, dp)
    nxc = Nx_in/2 + 1

    allocate(P%u(P%Ny, P%Nx), P%v(P%Ny, P%Nx), P%rhs_u(P%Ny, P%Nx), P%rhs_v(P%Ny, P%Nx))
    allocate(P%uhat(P%Ny, nxc), P%vhat(P%Ny, nxc), P%rhshat_u(P%Ny, nxc), P%rhshat_v(P%Ny, nxc))
    allocate(P%k2(P%Ny, nxc))
    allocate(P%dealias_mask(P%Ny, nxc))

    ! Precompute squared wave numbers kx^2 + ky^2 on r2c grid
    dkx = 2.0_dp * acos(-1.0_dp) / P%Lx
    dky = 2.0_dp * acos(-1.0_dp) / P%Ly

    do j = 1, P%Ny
       ky_idx = j - 1
       if (ky_idx > P%Ny/2) ky_idx = ky_idx - P%Ny
       do i = 1, nxc
          kx_idx = i - 1
          P%k2(j, i) = (dkx*real(kx_idx, dp))**2 + (dky*real(ky_idx, dp))**2
       end do
    end do

    ! Build 2/3 de-alias mask (1 = keep, 0 = zero) on the r2c grid
    if (P%use_dealias) then
       do j = 1, P%Ny
          ky_idx = j - 1
          if (ky_idx > P%Ny/2) ky_idx = ky_idx - P%Ny
          do i = 1, nxc
             kx_idx = i - 1
             if ( abs(kx_idx) > P%Nx/3 .or. abs(ky_idx) > P%Ny/3 ) then
                P%dealias_mask(j,i) = 0.0_dp
             else
                P%dealias_mask(j,i) = 1.0_dp
             end if
          end do
       end do
    else
       P%dealias_mask = 1.0_dp
    end if

    ! Create FFTW plans. Using the allocated arrays here; executing plans with
    ! different in/out pointers at runtime is permitted by fftw (execute takes
    ! explicit pointers). We create separate plans for u and v transforms.
    P%plan_r2c_u = fftw_plan_dft_r2c_2d(P%Ny, P%Nx, P%u, P%uhat, FFTW_MEASURE)
    P%plan_c2r_u = fftw_plan_dft_c2r_2d(P%Ny, P%Nx, P%uhat, P%u, FFTW_MEASURE)
    P%plan_r2c_v = fftw_plan_dft_r2c_2d(P%Ny, P%Nx, P%v, P%vhat, FFTW_MEASURE)
    P%plan_c2r_v = fftw_plan_dft_c2r_2d(P%Ny, P%Nx, P%vhat, P%v, FFTW_MEASURE)
  end subroutine init_rd_spec_plan

  subroutine destroy_rd_spec_plan(P)
    include 'fftw3.f03'
    type(rd_spec_plan), intent(inout) :: P
    if (c_associated(P%plan_r2c_u)) call fftw_destroy_plan(P%plan_r2c_u)
    if (c_associated(P%plan_c2r_u)) call fftw_destroy_plan(P%plan_c2r_u)
    if (c_associated(P%plan_r2c_v)) call fftw_destroy_plan(P%plan_r2c_v)
    if (c_associated(P%plan_c2r_v)) call fftw_destroy_plan(P%plan_c2r_v)

    if (allocated(P%u)) deallocate(P%u)
    if (allocated(P%v)) deallocate(P%v)
    if (allocated(P%rhs_u)) deallocate(P%rhs_u)
    if (allocated(P%rhs_v)) deallocate(P%rhs_v)
    if (allocated(P%uhat)) deallocate(P%uhat)
    if (allocated(P%vhat)) deallocate(P%vhat)
    if (allocated(P%rhshat_u)) deallocate(P%rhshat_u)
    if (allocated(P%rhshat_v)) deallocate(P%rhshat_v)
    if (allocated(P%k2)) deallocate(P%k2)
    if (allocated(P%dealias_mask)) deallocate(P%dealias_mask)
  end subroutine destroy_rd_spec_plan

  !------------------------------------------------------------
  ! One time step using semi-implicit pseudo-spectral scheme
  !   u_t = gamm*f(u,v) - r*u + e^{-2 r t} ∇^2 u
  !   v_t = gamm*g(u,v) - r*v + d * e^{-2 r t} ∇^2 v
  ! where f,g implement Schnakenberg or Gierer–Meinhardt.
  ! Diffusion handled implicitly in spectral space at time t_n.
  ! Nonlinear terms are de-aliased using the 2/3-rule mask.
  !------------------------------------------------------------
  subroutine rd_step_spectral_uniform_growth(P, a,b,c,dcoef, r, gamm, dt, t, model_flag)
    include 'fftw3.f03'
    type(rd_spec_plan), intent(inout) :: P
    real(dp), intent(in) :: a, b, c, dcoef, r, gamm, dt, t
    integer, intent(in) :: model_flag
    integer :: Nx, Ny, i, j, nxc
    real(dp) :: e2rt, denom_u, denom_v

    Nx = P%Nx; Ny = P%Ny; nxc = Nx/2 + 1


    ! 1) Build reaction RHS in real space: rhs = gamm*f(u,v) - r*u
    select case(model_flag)
    case(1) ! Schnakenberg
       do j = 1, Ny
          do i = 1, Nx
             P%rhs_u(j,i) = gamm * ( a - P%u(j,i) + (P%u(j,i)**2)*P%v(j,i) ) - r*P%u(j,i)
             P%rhs_v(j,i) = gamm * ( b - (P%u(j,i)**2)*P%v(j,i) ) - r*P%v(j,i)
          end do
       end do
    case(2) ! Gierer–Meinhardt
       do j = 1, Ny
          do i = 1, Nx
             if (P%v(j,i) > 1.0e-12_dp) then
                P%rhs_u(j,i) = gamm * ( a - b*P%u(j,i) + (P%u(j,i)**2)/P%v(j,i) ) - r*P%u(j,i)
                P%rhs_v(j,i) = gamm * ( (P%u(j,i)**2) - c*P%v(j,i) ) - r*P%v(j,i)
             else
                P%rhs_u(j,i) = - r*P%u(j,i)
                P%rhs_v(j,i) = - r*P%v(j,i)
             end if
          end do
       end do
    case default
       stop 'Unknown model_flag in rd_step_spectral_uniform_growth'
    end select

    ! 2) Forward FFT of u, v, and RHS terms
    call fftw_execute_dft_r2c(P%plan_r2c_u, P%u, P%uhat)
    call fftw_execute_dft_r2c(P%plan_r2c_v, P%v, P%vhat)

    call fftw_execute_dft_r2c(P%plan_r2c_u, P%rhs_u, P%rhshat_u)
    call fftw_execute_dft_r2c(P%plan_r2c_v, P%rhs_v, P%rhshat_v)

    ! 2b) De-alias nonlinear RHS in spectral space (2/3 rule)
    if (P%use_dealias) then
       do j = 1, Ny
          do i = 1, nxc
             if (P%dealias_mask(j,i) == 0.0_dp) then
                P%rhshat_u(j,i) = (0.0_dp, 0.0_dp)
                P%rhshat_v(j,i) = (0.0_dp, 0.0_dp)
             end if
          end do
       end do
    end if

    ! 3) Implicit diffusion step in spectral space at time t_n
    !     (I - dt * e^{-2 r t_n} ∇^2) u^{n+1} = u^n + dt * rhs^n
    !     => û^{n+1} = (û^n + dt*rhŝ^n) / (1 + dt*e^{-2 r t_n} * k^2)
    e2rt = exp(-2.0_dp*r*t)

    do j = 1, Ny
       do i = 1, nxc
          denom_u = 1.0_dp + dt * e2rt * P%k2(j,i)
          denom_v = 1.0_dp + dt * dcoef * e2rt * P%k2(j,i)
          P%uhat(j,i) = (P%uhat(j,i) + dt * P%rhshat_u(j,i)) / denom_u
          P%vhat(j,i) = (P%vhat(j,i) + dt * P%rhshat_v(j,i)) / denom_v
       end do
    end do

    ! 4) Inverse FFT to real space and normalize
    call fftw_execute_dft_c2r(P%plan_c2r_u, P%uhat, P%u)
    call fftw_execute_dft_c2r(P%plan_c2r_v, P%vhat, P%v)

    P%u = P%u * P%invN
    P%v = P%v * P%invN
  end subroutine rd_step_spectral_uniform_growth

end module rd_spectral_fftw

