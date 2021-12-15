module const
    use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t

    implicit none

    ! Integer types
    integer, parameter :: int_32 = selected_int_kind(6)
    integer, parameter :: int_64 = selected_int_kind(15)

    ! Real types
    integer, parameter :: sp = selected_real_kind(6,37)
    integer, parameter :: dp = selected_real_kind(15,307)

    integer, parameter :: p = dp

    ! Highest precision pi
    real(p), parameter :: pi = 4.0_p * atan(1.0_p)

    real(p), parameter :: depsilon = 1.e-12_p

end module const
