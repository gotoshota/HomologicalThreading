module compute
    use omp_lib
    use ieee_arithmetic
    implicit none

contains

    subroutine threading(pd_i, pd_i_cup_j, threading_flags)
        implicit none

        double precision, intent(in) :: pd_i(:, :, :) ! shape: (nchains, npoints, 2)
        double precision, intent(in) :: pd_i_cup_j(:, :, :, :) ! shape: (nchains, nchains, npoints2, 2)
        logical, intent(inout) :: threading_flags(:, :, :) ! shape: (nchains, nchains, npoints) 

        integer :: nchains, npoints, npoints2, i, j, k, l
        double precision :: diff(2)
        double precision, parameter :: threshold = 1d-8

        nchains = size(pd_i, 1)
        npoints = size(pd_i, 2)
        npoints2 = size(pd_i_cup_j, 3)

        threading_flags = .true.

        !$omp parallel do private(i, j, k, l, diff)
        do i = 1, nchains  ! PD(i), i represents the passive ring
            do k = 1, npoints
                ! NaN チェック
                if (ieee_is_nan(pd_i(i, k, 1)) .or. ieee_is_nan(pd_i(i, k, 2))) then
                    threading_flags(i, :, k) = .false.
                    cycle
                end if

                do j = 1, nchains  ! j represents the active ring
                    if (i == j) then
                        threading_flags(i, j, k) = .false.
                        cycle
                    end if

                    ! Compute difference of set between PD(i) and PD(i cup j)
                    ! PD(i) / PD(i cup j) = threading_array(i, j, :, :)
                    do l = 1, npoints2
                        if (ieee_is_nan(pd_i_cup_j(i, j, l, 1)) .or. ieee_is_nan(pd_i_cup_j(i, j, l, 2))) cycle

                        diff = pd_i(i, k, :) - pd_i_cup_j(i, j, l, :)
                        if (all(abs(diff) < threshold)) then
                            threading_flags(i, j, k) = .false.
                        end if
                    end do
                end do
            end do
        end do
        !$omp end parallel do

    end subroutine threading
end module compute

