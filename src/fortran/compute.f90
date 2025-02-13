module compute
   use omp_lib
   use ieee_arithmetic
   implicit none

contains

   subroutine threading(pd_i, pd_i_cup_j, threading_flags, threading_pd)
      implicit none

      double precision, intent(in) :: pd_i(:, :, :) ! shape: (nchains, npoints, 2)
      double precision, intent(in) :: pd_i_cup_j(:, :, :, :) ! shape: (nchains, nchains, npoints2, 2)
      logical, intent(inout) :: threading_flags(:, :) ! shape: (nchains, nchains)
      double precision, intent(inout) :: threading_pd(:, :, :, :) ! shape: (nchains, nchains, npoints, 2)

      integer :: nchains, npoints, npoints2, i, j, k, l
      integer :: n
      logical, allocatable :: flags(:) ! shape: (npoitns)
      double precision :: diff(2)
      double precision, parameter :: threshold = 1d-8

      nchains = size(pd_i, 1)
      npoints = size(pd_i, 2)
      npoints2 = size(pd_i_cup_j, 3)

      allocate (flags(npoints))
      threading_flags = .false.
      threading_pd = -1d0

      !$omp parallel do private(i, j, k, l, n, diff) shared(pd_i, pd_i_cup_j, threading_flags, threading_pd) default(none)
      do i = 1, nchains  ! PD(i), i represents the passive ring
         do k = 1, npoints ! PD(i) の各点について，threading されているかどうかを調べる
            ! NaN なら，それ以降の要素も NaN なので，計算しない
            if (ieee_is_nan(pd_i(i, k, 1)) .or. ieee_is_nan(pd_i(i, k, 2))) continue

            do j = 1, nchains  ! j represents the active ring
               if (i == j) cycle ! 同じリング同士なら，計算しない

               n = 0 ! threading の数
               flags = .true. ! threading されているかどうかのフラグ， True で初期化

               do l = 1, npoints2
                  ! NaN なら，それ以降の要素も NaN なので，計算しない
                  if (ieee_is_nan(pd_i_cup_j(i, j, l, 1)) .or. ieee_is_nan(pd_i_cup_j(i, j, l, 2))) then
                     flags(l:npoints) = .false.
                     continue
                  end if

                  ! 点の距離が threshold 以下なら，同じ点とみなす → threading されていないループ
                  diff = pd_i(i, k, :) - pd_i_cup_j(i, j, l, :)
                  if (all(abs(diff) < threshold)) then
                     flags(k) = .false.
                  end if
               end do
               ! 1つでも true があれば，threading されている
               if (any(flags)) then
                  threading_flags(i, j) = .true.
                  ! flagがTrueの点を threading_pd に格納
                  n = 0
                  do l = 1, npoints
                     if (flags(l)) then
                        n = n + 1
                        threading_pd(i, j, n, :) = pd_i(i, l, :)
                     end if
                  end do
               end if
            end do
         end do
      end do
      !$omp end parallel do

   end subroutine threading
end module compute
