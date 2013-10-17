subroutine ddx(nx, ny, a, df)
    integer, intent(in) :: nx, ny
    real, dimension(:, :), intent(in) :: a
    real, dimension(size(a, 1), size(a, 2)), intent(out) :: df
    df = a
end subroutine ddx

subroutine init_displacement(a, x, y, nx, ny)
    use, intrinsic :: iso_c_binding
    implicit none
    integer, intent(in) :: nx, ny
    real(C_DOUBLE), dimension(nx, ny), intent(inout) :: a
    real(C_DOUBLE), dimension(nx, ny), intent(in) :: x
    real(C_DOUBLE), dimension(nx, ny), intent(in) :: y

    a = sin(x)

end subroutine init_displacement

program fftw_test
    ! C binding
    use, intrinsic :: iso_c_binding
    implicit none

    double precision, parameter :: pi = 4*ATAN(1.0_C_DOUBLE)
    complex, parameter :: ii = (0.0, 1.0)

    integer(C_INT), parameter :: Nx = 16
    integer(C_INT), parameter :: Ny = Nx
    double precision, parameter :: Lx = 2*pi, Ly = 2*pi
    ! Derived paramenter
    double precision, parameter :: dx = Lx/Nx, dy = Ly/Ny

    real(C_DOUBLE), dimension(Nx, Ny) :: x,y, u0,in,dudx,dudxE, errdU
    real(C_DOUBLE), dimension(Nx/2+1, Ny) :: kx, ky

    ! Fourier space variables
    complex(C_DOUBLE_COMPLEX), dimension(Nx / 2 + 1, Ny) :: u_hat_x, out
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
!    u0 = sin(x)
    call init_displacement(u0, x, y, Nx, Ny)

    dudxE = cos(x)
    do i=1, Nx
        write (*, *) u0(i, :)
    end do
    ! Go to Fourier Space
    in = u0
    pf = fftw_plan_dft_r2c_2d(Nx, Ny, in, out, FFTW_ESTIMATE)
    call fftw_execute_dft_r2c(pf,in,out)

    write(*,'(A)') 'Fourier...'
    do i=1, Nx / 2 + 1
        write (*, *) out(i, :)
    end do

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
