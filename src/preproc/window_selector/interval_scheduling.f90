module interval_scheduling
  use fwat_constants
  implicit none

  type :: interval_type
    integer :: original_index
    real(kind=dp) :: left
    real(kind=dp) :: right
    real(kind=dp) :: weight
  end type interval_type

contains

  subroutine schedule_weighted_intervals(intervals, selected_indices, n_selected)
    type(interval_type), dimension(:), intent(inout) :: intervals
    integer, allocatable, dimension(:), intent(out) :: selected_indices
    integer, intent(out) :: n_selected

    integer :: n, j
    integer, allocatable :: p(:)
    real(kind=dp), allocatable :: opt(:)
    
    n = size(intervals)
    if (n == 0) then
        n_selected = 0
        return
    end if

    ! Sort intervals by finish time (right)
    call sort_intervals(intervals)

    ! Compute previous intervals p[j]
    allocate(p(n))
    call compute_previous_intervals(intervals, p)

    ! Compute OPTs iteratively
    allocate(opt(0:n))
    opt(0) = 0.0_dp
    
    do j = 1, n
       ! p is 1-based index in Fortran, but 0 means no compatible interval.
       ! If p(j) == 0, then opt(p(j)) -> opt(0) which is 0.
       opt(j) = max(intervals(j)%weight + opt(p(j)), opt(j-1))
    end do

    ! Backtrack to find solution
    ! We need a list to store selected intervals.
    ! Let's allocate max size.
    allocate(selected_indices(n))
    n_selected = 0
    
    call compute_solution(n)

    ! The selected_indices will contain the original_index of the intervals.
    ! Reverse the selected_indices to be in chronological order
    block
      integer :: temp_idx, k
      do k = 1, n_selected / 2
         temp_idx = selected_indices(k)
         selected_indices(k) = selected_indices(n_selected - k + 1)
         selected_indices(n_selected - k + 1) = temp_idx
      end do
    end block
    
  contains
  
    recursive subroutine compute_solution(j)
      integer, intent(in) :: j
      if (j == 0) return
      
      if (intervals(j)%weight + opt(p(j)) > opt(j-1)) then
         n_selected = n_selected + 1
         selected_indices(n_selected) = intervals(j)%original_index
         call compute_solution(p(j))
      else
         call compute_solution(j-1)
      endif
    end subroutine compute_solution

  end subroutine schedule_weighted_intervals

  subroutine compute_previous_intervals(intervals, p)
    type(interval_type), dimension(:), intent(in) :: intervals
    integer, dimension(:), intent(out) :: p
    
    integer :: n, j, i
    real(kind=dp) :: start_val
    
    n = size(intervals)
    
    do j = 1, n
       start_val = intervals(j)%left
       ! Find the rightmost interval i such that intervals(i)%right <= start_val
       ! Since intervals are sorted by right, we can use binary search.
       ! We want largest i such that intervals(i)%right <= intervals(j)%left
       
       call bisect_right(intervals, start_val, j-1, i)
       p(j) = i
    end do
  end subroutine compute_previous_intervals

  subroutine bisect_right(intervals, val, high_idx, idx)
    type(interval_type), dimension(:), intent(in) :: intervals
    real(kind=dp), intent(in) :: val
    integer, intent(in) :: high_idx
    integer, intent(out) :: idx
    
    integer :: low, high, mid
    
    low = 1
    high = high_idx
    idx = 0
    
    do while (low <= high)
       mid = (low + high) / 2
       if (intervals(mid)%right <= val) then
          idx = mid
          low = mid + 1
       else
          high = mid - 1
       end if
    end do
  end subroutine bisect_right

  subroutine sort_intervals(arr)
    type(interval_type), dimension(:), intent(inout) :: arr
    call quicksort(arr, 1, size(arr))
  end subroutine sort_intervals

  recursive subroutine quicksort(arr, first, last)
    type(interval_type), dimension(:), intent(inout) :: arr
    integer, intent(in) :: first, last
    integer :: i, j, pivot_idx
    type(interval_type) :: temp, pivot
    
    if (first < last) then
       pivot_idx = (first + last) / 2
       pivot = arr(pivot_idx)
       
       ! Swap pivot to end
       temp = arr(last)
       arr(last) = arr(pivot_idx)
       arr(pivot_idx) = temp
       
       i = first - 1
       do j = first, last - 1
          ! Sort by right time (primary) and left time (secondary)
          ! This ensures deterministic order for intervals with same end time
          if (arr(j)%right < pivot%right - 1.0e-6_dp .or. &
             (abs(arr(j)%right - pivot%right) <= 1.0e-6_dp .and. &
              arr(j)%left <= pivot%left)) then
             i = i + 1
             temp = arr(i)
             arr(i) = arr(j)
             arr(j) = temp
          end if
       end do
       
       ! Swap pivot back
       temp = arr(i + 1)
       arr(i + 1) = arr(last)
       arr(last) = temp
       
       call quicksort(arr, first, i)
       call quicksort(arr, i + 2, last)
    end if
  end subroutine quicksort

end module interval_scheduling
