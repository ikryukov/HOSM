subroutine init_displacement(eta, x, y, nx, ny)
    use, intrinsic :: iso_c_binding
    implicit none
    integer, intent(in) :: nx, ny
    real(C_DOUBLE), dimension(nx, ny), intent(inout) :: eta
    real(C_DOUBLE), dimension(nx, ny), intent(in) :: x
    real(C_DOUBLE), dimension(nx, ny), intent(in) :: y
    real(C_DOUBLE), parameter :: A = 1.0_C_DOUBLE

    eta = A * cos(x)

end subroutine init_displacement

subroutine init_potential(phi, x, y, nx, ny)
    use, intrinsic :: iso_c_binding
    implicit none
    integer, intent(in) :: nx, ny
    real(C_DOUBLE), dimension(nx, ny), intent(inout) :: phi
    real(C_DOUBLE), dimension(nx, ny), intent(in) :: x
    real(C_DOUBLE), dimension(nx, ny), intent(in) :: y
    real(C_DOUBLE), parameter :: A = 1.0_C_DOUBLE
    real(C_DOUBLE), parameter :: g = 9.8_C_DOUBLE
    real(C_DOUBLE), parameter :: H = 1000.0_C_DOUBLE
    real(C_DOUBLE) :: omega, sigma

    omega = sqrt(g * tanh(H))
    sigma = 1.0_C_DOUBLE ! H inf <-> deep water
    phi = A * sin(x) * (-omega / sigma)

end subroutine init_potential

subroutine rhs_Phi(b, y, t, res)
    use, intrinsic :: iso_c_binding
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(in) :: b, y
    real(C_DOUBLE), intent(in) :: t
    complex(C_DOUBLE_COMPLEX), intent(out) :: res
    real(C_DOUBLE), parameter :: g = 9.8_C_DOUBLE

    res = -g * b

end subroutine rhs_Phi

subroutine runge_kutta4_Phi(yn, t, dt, b, res)
    use, intrinsic :: iso_c_binding
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(in) :: yn, b
    real(C_DOUBLE), intent(in) :: t, dt
    complex(C_DOUBLE_COMPLEX), intent(out) :: res
    complex(C_DOUBLE_COMPLEX) :: k1, k2, k3, k4
    complex(C_DOUBLE_COMPLEX) :: tmp

    call rhs_Phi(b, yn, t, tmp)
    k1 = dt * tmp
    call rhs_Phi(b, yn + 0.5_C_DOUBLE * k1, t + 0.5_C_DOUBLE * dt, tmp)
    k2 = dt * tmp
    call rhs_Phi(b, yn + 0.5_C_DOUBLE * k2, t + 0.5_C_DOUBLE * dt, tmp)
    k3 = dt * tmp
    call rhs_Phi(b, yn + k3, t + dt, tmp)
    k4 = dt * tmp

    res = yn + 1.0_C_DOUBLE / 6.0_C_DOUBLE * (k1 + 2.0_C_DOUBLE * k2 + 2.0_C_DOUBLE * k3 + k4)

end subroutine runge_kutta4_Phi

subroutine rhs_Eta(k_mod, sigma, a, y, t, res)
    use, intrinsic :: iso_c_binding
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(in) :: a, y
    real(C_DOUBLE), intent(in) :: k_mod, sigma, t
    complex(C_DOUBLE_COMPLEX), intent(out) :: res

    res = k_mod * sigma * a

end subroutine rhs_Eta

subroutine runge_kutta4_Eta(yn, t, dt, a, sigma, k_mod, res)
    use, intrinsic :: iso_c_binding
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(in) :: yn, a
    complex(C_DOUBLE_COMPLEX), intent(out) :: res
    real(C_DOUBLE), intent(in) :: t, dt, sigma, k_mod
    complex(C_DOUBLE_COMPLEX) :: k1, k2, k3, k4
    complex(C_DOUBLE_COMPLEX) :: tmp

    call rhs_Eta(k_mod, sigma, a, yn, t, tmp)
    k1 = dt * tmp
    call rhs_Eta(k_mod, sigma, a, yn + 0.5_C_DOUBLE * k1, t + 0.5_C_DOUBLE * dt, tmp)
    k2 = dt * tmp
    call rhs_Eta(k_mod, sigma, a, yn + 0.5_C_DOUBLE * k2, t + 0.5_C_DOUBLE * dt, tmp)
    k3 = dt * tmp
    call rhs_Eta(k_mod, sigma, a, yn + k3, t + dt, tmp)
    k4 = dt * tmp

    res = yn + 1.0_C_DOUBLE / 6.0_C_DOUBLE * (k1 + 2.0_C_DOUBLE * k2 + 2.0_C_DOUBLE * k3 + k4)

end subroutine runge_kutta4_Eta

subroutine print_matrix_real(matrix, size_i, size_j)
    use, intrinsic :: iso_c_binding
    implicit none
    integer, intent(in) :: size_i, size_j
    real(C_DOUBLE), dimension(size_i, size_j), intent(in) :: matrix
    integer :: i

    do i = 1, size_i
        write (*, *) matrix(i, :)
    end do

end subroutine print_matrix_real

subroutine print_matrix_cmplx(matrix, size_i, size_j)
    use, intrinsic :: iso_c_binding
    implicit none
    integer, intent(in) :: size_i, size_j
    complex(C_DOUBLE_COMPLEX), dimension(size_i, size_j), intent(in) :: matrix
    integer :: i

    do i = 1, size_i
        write (*, *) matrix(i, :)
    end do

end subroutine print_matrix_cmplx


program fftw_test
    ! C binding
    use, intrinsic :: iso_c_binding
    implicit none

    double precision, parameter :: pi = 4.0_C_DOUBLE * ATAN(1.0_C_DOUBLE)
    complex, parameter :: ii = (0.0, 1.0)

    integer(C_INT), parameter :: Nx = 8
    integer(C_INT), parameter :: Ny = Nx
    double precision, parameter :: Lx = 2.0 * pi, Ly = 2.0 * pi
    double precision, parameter :: dx = Lx / Nx, dy = Ly / Ny
    double precision, parameter :: H = 1000.0

    real(C_DOUBLE), dimension(Nx, Ny) :: x, y, eta, phi, in, dudx, dudxE, errdU
    real(C_DOUBLE), dimension(Nx / 2 + 1, Ny) :: kx, ky, k_mod, sigma
    real(C_DOUBLE) :: t, dt
    ! Fourier space variables
    complex(C_DOUBLE_COMPLEX), dimension(Nx / 2 + 1, Ny) :: fourier_eta, fourier_phi, fourier_eta_new, fourier_phi_new, out
    ! indices
    integer :: i, j
    !---FFTW plans
    type(C_PTR) :: pf, pb

    ! FFTW include
include 'fftw3.f03'

    write(*,'(A)') 'Starting...'

    ! Grid
    forall(i=1:Nx,j=1:Ny)
        x(i,j) = (i-1)*dx
        y(i,j) = (j-1)*dy
    end forall

    ! Compute the wavenumbers
    forall(i=1:Nx/2,j=1:Ny) kx(i,j)=2*pi*(i-1)/Lx
    kx(Nx/2+1,:) = 0.0
    forall(i=1:Nx/2+1,j=1:Ny/2) ky(i,j)=2*pi*(j-1)/Ly
    forall(i=1:Nx/2+1,j=Ny/2+1:Ny) ky(i,j)=2*pi*(j-Ny-1)/Ly

    forall(i = 1: Nx / 2 + 1, j = 1: Ny)
        k_mod(i, j) = sqrt(kx(i, j) ** 2 + ky(i, j) ** 2)
        sigma(i, j) = tanh(k_mod(i, j) * H)
    end forall

    ! Initial Condition
    call init_displacement(eta, x, y, Nx, Ny)
    write (*, *) "Eta matrix:"
    call print_matrix_real(eta, Nx, Ny)
    call init_potential(phi, x, y, Nx, Ny)
    write (*, *) "Phi matrix:"
    call print_matrix_real(phi, Nx, Ny)

    ! Prepare FFTW plans
    pf = fftw_plan_dft_r2c_2d(Nx, Ny, in, out, FFTW_ESTIMATE)

    ! Go to Fourier Space
    ! eta -> b
    in = eta
    call fftw_execute_dft_r2c(pf, in, out)
    fourier_eta = out
    write(*,'(A)') 'Fourier Eta...'
    call print_matrix_cmplx(fourier_eta, Nx / 2 + 1, Ny)
    ! phi -> a
    in = phi
    call fftw_execute_dft_r2c(pf, in, out)
    fourier_phi = out
    write(*,'(A)') 'Fourier Phi...'
    call print_matrix_cmplx(fourier_phi, Nx / 2 + 1, Ny)

    dt = 0.01_C_DOUBLE
    t = 0.0_C_DOUBLE

    forall(i = 1: Nx / 2 + 1, j = 1: Ny)
        fourier_phi_new(i, j) = 0.0_C_DOUBLE
        fourier_eta_new(i, j) = 0.0_C_DOUBLE
    end forall

    do i = 1, Nx / 2 + 1
        do j = 1, Ny
            ! d/dt a(n, p) = -g * b(n, p)
            ! ^
            ! v
            ! d/dt fourier_phi(i, j) = -g * fourier_eta(i, j)
            call runge_kutta4_Phi(fourier_phi(i, j), t, dt, fourier_eta(i, j), fourier_phi_new(i, j))
            ! d/dt b(n, p) = |k(n, p)| * sigma * a(n, p)
            ! ^
            ! v
            ! d/dt fourier_eta(i, j) = k_mod(i, j) * sigma * fourier_phi(i, j)
            call runge_kutta4_Eta(fourier_eta(i, j), t, dt, fourier_phi(i, j), sigma(i, j), k_mod(i, j), fourier_eta_new(i, j))
        end do
    end do

    ! Back to physical space
    out = fourier_phi_new
    in = phi
    pb = fftw_plan_dft_c2r_2d(Nx, Ny, out, in, FFTW_ESTIMATE)
    call fftw_execute_dft_c2r(pb, out, in)
    ! rescale
    phi = in / (Nx * Ny)

    out = fourier_eta_new
    in = eta
    pb = fftw_plan_dft_c2r_2d(Nx, Ny, out, in, FFTW_ESTIMATE)
    call fftw_execute_dft_c2r(pb, out, in)
    ! rescale
    eta = in / (Nx * Ny)

    write (*, *) "Final Eta matrix:"
    call print_matrix_real(eta, Nx, Ny)

    write (*, *) "Final Phi matrix:"
    call print_matrix_real(phi, Nx, Ny)

    call fftw_destroy_plan(pf)
    call fftw_destroy_plan(pb)

    write(*,'(A)') 'Done!'
end program fftw_test
