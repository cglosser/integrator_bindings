module cfuncs
  interface
    subroutine linear_weights_array(n_lambda, n_time, radius, weights, future_coeff) bind(C, name="linear_weights_array")
        use iso_c_binding, only: C_INT, C_DOUBLE
        integer(kind=C_INT), intent(in) :: n_lambda, n_time
        real(kind=C_DOUBLE), intent(in) :: radius
        real(kind=C_DOUBLE), intent(out) :: weights(*), future_coeff
    end subroutine linear_weights_array
  end interface
end module

module predictor_corrector
    contains 
    subroutine weights_array(n_lambda, n_time, radius, weights, future_coeff)
        use iso_c_binding, only: C_INT, C_DOUBLE
        use cfuncs
        integer(kind=C_INT), intent(in) :: n_lambda, n_time
        real(kind=C_DOUBLE), intent(in) :: radius
        real(kind=C_DOUBLE), intent(out) :: weights(n_time, 4), future_coeff
        real(kind=C_DOUBLE) :: linear_weights(4 * n_time)

        future_coeff = 0
        linear_weights = 0

        call linear_weights_array(n_lambda, n_time, radius, linear_weights, future_coeff)
        weights = reshape(linear_weights, shape(weights))
    end
end