module heap_module
  !---------------------------------------------------------------------
  ! Binary min-heap with decrease-key support.
  !
  ! - val(k)   : value stored at heap position k
  ! - idx(k)   : index stored at heap position k
  ! - pos(i)   : heap position of index i (0 if not present)
  !
  ! The pos(:) array is what allows decrease_key in O(log n) instead
  ! of a linear search through the heap.
  !---------------------------------------------------------------------
  use precision_mod, only : wp                            !provides kind for eigther single or double precision (wp means working precision)
  implicit none

  type, public :: heap_t
     real(wp), allocatable :: val(:)
     integer,  allocatable :: idx(:)
     integer,  allocatable :: pos(:)   ! size npoints
     integer :: n = 0                  ! current number of elements in the heap
  end type heap_t

  public :: heap_init
  public :: heap_insert
  public :: heap_extract_min
  public :: heap_decrease_key
  public :: heap_is_empty
  public :: heap_contains
  public :: heap_destroy

contains

  subroutine heap_init(h, npoints, capacity)
    type(heap_t), intent(out) :: h
    integer, intent(in) :: npoints   ! total number of grid points
    integer, intent(in), optional :: capacity

    integer :: cap

    cap = npoints
    if (present(capacity)) cap = capacity

    allocate(h%val(cap))
    allocate(h%idx(cap))
    allocate(h%pos(npoints))
    h%pos = 0
    h%n = 0
  end subroutine heap_init

  logical function heap_is_empty(h)
    type(heap_t), intent(in) :: h
    heap_is_empty = (h%n == 0)
  end function heap_is_empty

  logical function heap_contains(h, i)
    ! Is grid point i currently in the heap?
    type(heap_t), intent(in) :: h
    integer, intent(in) :: i
    heap_contains = (h%pos(i) /= 0)
  end function heap_contains

  subroutine heap_insert(h, i, v)
    ! Insert grid point i with value v
    type(heap_t), intent(inout) :: h
    integer,  intent(in) :: i
    real(wp), intent(in) :: v

    h%n = h%n + 1
    if (h%n > size(h%val)) call heap_grow(h)

    h%val(h%n) = v
    h%idx(h%n) = i
    h%pos(i)   = h%n

    call sift_up(h, h%n)
  end subroutine heap_insert

  subroutine heap_decrease_key(h, i, v)
    ! Decrease the value associated with grid point i (already in the heap)
    ! If v is not smaller than the current value, does nothing.
    type(heap_t), intent(inout) :: h
    integer,  intent(in) :: i
    real(wp), intent(in) :: v

    integer :: k

    k = h%pos(i)
    if (k == 0) then
       ! not in the heap: insert directly
       call heap_insert(h, i, v)
       return
    end if

    if (v < h%val(k)) then
       h%val(k) = v
       call sift_up(h, k)
    end if
    ! if v >= h%val(k), nothing to do (this is a decrease-only operation)
  end subroutine heap_decrease_key

  subroutine heap_extract_min(h, i_min, v_min)
    ! Remove and return the grid point with the smallest value
    type(heap_t), intent(inout) :: h
    integer,  intent(out) :: i_min
    real(wp), intent(out) :: v_min

    i_min = h%idx(1)
    v_min = h%val(1)

    h%pos(i_min) = 0

    ! Move the last element to the root, then sift down
    h%val(1) = h%val(h%n)
    h%idx(1) = h%idx(h%n)
    h%pos(h%idx(1)) = 1
    h%n = h%n - 1

    if (h%n > 0) call sift_down(h, 1)
  end subroutine heap_extract_min

  subroutine sift_up(h, k0)
    type(heap_t), intent(inout) :: h
    integer, intent(in) :: k0

    integer :: k, parent

    k = k0
    do while (k > 1)
       parent = k / 2
       if (h%val(k) < h%val(parent)) then
          call swap(h, k, parent)
          k = parent
       else
          exit
       end if
    end do
  end subroutine sift_up

  subroutine sift_down(h, k0)
    type(heap_t), intent(inout) :: h
    integer, intent(in) :: k0

    integer :: k, left, right, smallest

    k = k0
    do
       left  = 2*k
       right = 2*k + 1
       smallest = k

       if (left  <= h%n) then
          if (h%val(left) < h%val(smallest)) smallest = left
       end if
       if (right <= h%n) then
          if (h%val(right) < h%val(smallest)) smallest = right
       end if

       if (smallest == k) exit

       call swap(h, k, smallest)
       k = smallest
    end do
  end subroutine sift_down

  subroutine swap(h, a, b)
    type(heap_t), intent(inout) :: h
    integer, intent(in) :: a, b

    real(wp) :: vtmp
    integer  :: itmp

    vtmp = h%val(a); h%val(a) = h%val(b); h%val(b) = vtmp
    itmp = h%idx(a); h%idx(a) = h%idx(b); h%idx(b) = itmp

    h%pos(h%idx(a)) = a
    h%pos(h%idx(b)) = b
  end subroutine swap

  subroutine heap_grow(h)
    ! Double the capacity if needed (rare if capacity = npoints)
    type(heap_t), intent(inout) :: h

    real(wp), allocatable :: val_tmp(:)
    integer,  allocatable :: idx_tmp(:)
    integer :: newcap

    newcap = max(2*size(h%val), h%n)

    allocate(val_tmp(newcap)); val_tmp(1:h%n-1) = h%val(1:h%n-1)
    allocate(idx_tmp(newcap)); idx_tmp(1:h%n-1) = h%idx(1:h%n-1)

    call move_alloc(val_tmp, h%val)
    call move_alloc(idx_tmp, h%idx)
  end subroutine heap_grow

  subroutine heap_destroy(h)
    ! Deallocate all components of h
    type(heap_t), intent(inout) :: h

    if (allocated(h%val)) deallocate(h%val)
    if (allocated(h%idx)) deallocate(h%idx)
    if (allocated(h%pos)) deallocate(h%pos)
    h%n = 0
  end subroutine heap_destroy
end module heap_module