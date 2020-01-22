program test
    use iso_c_binding
    use predictor_corrector
    implicit none

    integer(kind=c_int), parameter :: n_lambda = 32
    integer(kind=c_int), parameter :: n_time = 22
    real(kind=c_double), parameter :: radius = 3.15
    real(kind=c_double) :: weights(n_time, 4)
    real(kind=c_double) :: future_coeff
    integer :: i

    call weights_array(n_lambda, n_time, radius, weights, future_coeff)

    do i = 1, n_time
        write(*,*) weights(i,:)
    end do

    write(*,*) future_coeff
end program