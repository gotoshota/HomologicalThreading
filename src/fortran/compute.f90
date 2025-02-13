module compute
    use omp_lib
    use ieee_arithmetic
    implicit none

contains

    subroutine threading(pd_i, pd_i_cup_j, threading_flags, threading_pd)
        implicit none
         
        double precision, intent(in) :: pd_i(:, :, :) ! shape: (2, npoints, nchains)
        double precision, intent(in) :: pd_i_cup_j(:, :, :, :) ! shape: (2, npoints2, nchains_a, nchains_p)
        logical, intent(inout) :: threading_flags(:, :) ! shape: (nchains, nchains)
        double precision, intent(inout) :: threading_pd(:, :, :, :) ! shape: (2, npoints, nchains_a, nchains_p)

        integer :: nchains, npoints, npoints2, i, j, k, l
        integer :: n
        logical, allocatable :: flags(:) ! shape: (npoitns)
        double precision :: diff(2)
        double precision :: target_point(2)
        double precision, parameter :: threshold = 1d-8

        nchains = size(pd_i, 3)
        npoints = size(pd_i, 2)
        npoints2 = size(pd_i_cup_j, 2)

        threading_flags = .false.
        threading_pd = -1d0 ! PD は 0 以上の値なので，-1 で初期化

        !$omp parallel do private(i, j, k, l, n, diff, flags, target_point) shared(pd_i, pd_i_cup_j, threading_flags, threading_pd)
        loop_passive_chain: do i = 1, nchains

            loop_active_chain: do j = 1, nchains
                if (i == j) cycle ! 次の j へ
                allocate(flags(npoints))
                flags = .true. ! 各点がthreading されているかどうかのフラグ， True で初期化

                loop_passive_point: do k = 1, npoints 
                    ! NaN なら，それ以降の要素も NaN なので，計算しない
                    if (ieee_is_nan(pd_i(1, k, i)) .or. ieee_is_nan(pd_i(2, k, i))) then
                        flags(k:) = .false.
                        exit loop_passive_point
                    end if
                    target_point = pd_i(:, k, i)

                    loop_active_point: do l = 1, npoints2 
                        ! NaN なら，それ以降の要素も NaN なので，計算しない
                        if (ieee_is_nan(pd_i_cup_j(1, l, j, i)) .or. ieee_is_nan(pd_i_cup_j(2, l, j, i))) exit loop_active_point
                        ! 点の距離が threshold 以下なら，同じ点とみなす → threading されていないループ
                        diff(:) = target_point(:) - pd_i_cup_j(:, l, j, i)
                        if (all(abs(diff) < threshold)) then
                            flags(k) = .false.
                            exit loop_active_point
                        end if
                    end do loop_active_point
                    
                end do loop_passive_point

                ! flag が True の点を threading_pd に格納
                n = 0
                loop_threading_point: do k = 1, npoints
                    if (flags(k)) then
                        n = n + 1
                        threading_pd(:, n, j, i) = pd_i(:, k, i)
                    end if
                end do loop_threading_point

                ! 1つでも true があれば，threading されている
                if (any(flags)) then
                    threading_flags(j, i) = .true.
                end if
                deallocate(flags)
            end do loop_active_chain

        end do loop_passive_chain
        !$omp end parallel do

    end subroutine threading
end module compute
