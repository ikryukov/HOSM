subroutine ddx(nx, ny, a, df)
    integer, intent(in) :: nx, ny
    real, dimension(:, :), intent(in) :: a
    real, dimension(size(a, 1), size(a, 2)), intent(out) :: df
    df = a
end subroutine ddx

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

subroutine print_matrix(matrix, size_i, size_j)
    use, intrinsic :: iso_c_binding
    implicit none
    integer, intent(in) :: size_i, size_j
    real(C_DOUBLE), dimension(size_i, size_j), intent(in) :: matrix
    integer :: i

    do i = 1, size_i
        write (*, *) matrix(i, :)
    end do

end subroutine print_matrix

program fftw_test
    ! C binding
    use, intrinsic :: iso_c_binding
    implicit none

    double precision, parameter :: pi = 4*ATAN(1.0_C_DOUBLE)
    complex, parameter :: ii = (0.0, 1.0)

    integer(C_INT), parameter :: Nx = 8
    integer(C_INT), parameter :: Ny = Nx
    double precision, parameter :: Lx = 2*pi, Ly = 2*pi
    ! Derived paramenter
    double precision, parameter :: dx = Lx/Nx, dy = Ly/Ny

    real(C_DOUBLE), dimension(Nx, Ny) :: x, y, u0, eta, phi, in, dudx, dudxE, errdU
    real(C_DOUBLE), dimension(Nx/2+1, Ny) :: kx, ky

    ! Fourier space variables
    complex(C_DOUBLE_COMPLEX), dimension(Nx / 2 + 1, Ny) :: fourier_eta, fourier_phi, out
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

    ! Initial Condition
    call init_displacement(eta, x, y, Nx, Ny)
    write (*, *) "Eta matrix:"
    call print_matrix(eta, Nx, Ny)
    call init_potential(phi, x, y, Nx, Ny)
    write (*, *) "Phi matrix:"
    call print_matrix(phi, Nx, Ny)

    ! Prepare FFTW plans
    pf = fftw_plan_dft_r2c_2d(Nx, Ny, in, out, FFTW_ESTIMATE)

    ! Go to Fourier Space
    ! eta
    in = eta
    call fftw_execute_dft_r2c(pf, in, out)
    fourier_eta = out
    write(*,'(A)') 'Fourier Eta...'
    call print_matrix(fourier_eta, Nx / 2 + 1, Ny)
    ! phi
    in = phi
    call fftw_execute_dft_r2c(pf, in, out)
    fourier_phi = out
    write(*,'(A)') 'Fourier Phi...'
    call print_matrix(fourier_phi, Nx / 2 + 1, Ny)

    ! Derivative
    out = ii*kx*out
!    out = ii*ky*out

    ! Back to physical space
    pb = fftw_plan_dft_c2r_2d(Nx, Ny, out, in, FFTW_ESTIMATE)
    call fftw_execute_dft_c2r(pb,out,in)

    ! rescale
    dudx = in / (Nx * Ny)

    ! check the derivative
    errdU = dudx - dudxE

    ! Write file
    write(*,*) 'Writing to files...'
    do i=1, Nx
        write (*, *) errdU(i, :)
    end do

    call fftw_destroy_plan(pf)
    call fftw_destroy_plan(pb)

    write(*,'(A)') 'Done!'
end program fftw_test
