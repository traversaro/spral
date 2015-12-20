integer(C_INT) function spral_core_analyse_basic_analyse(n, ptr, row, perm, &
      nnodes, csptr, csparent, crptr, crlist, nemin, nfact, nflops, base) &
      bind(C)
   use iso_c_binding
   use spral_core_analyse, only: basic_analyse
   implicit none

   integer(C_INT), value :: n
   integer(C_INT), dimension(n+1), intent(in) :: ptr
   integer(C_INT), dimension(*), intent(in) :: row
   integer(C_INT), dimension(n), intent(inout) :: perm
   integer(C_INT), intent(out) :: nnodes
   type(C_PTR), intent(out) :: csptr
   type(C_PTR), intent(out) :: csparent
   type(C_PTR), intent(out) :: crptr
   type(C_PTR), intent(out) :: crlist
   integer(C_INT), value :: nemin
   integer(C_LONG), intent(out) :: nfact
   integer(C_LONG), intent(out) :: nflops
   integer(C_INT), value :: base

   integer :: stat
   integer, dimension(:), allocatable :: ptr2, row2
   integer(C_INT), dimension(:), allocatable :: fsptr, fsparent, frlist
   integer(C_LONG), dimension(:), allocatable :: frptr
   integer(C_INT), dimension(:), pointer :: psptr, psparent, prlist
   integer(C_LONG), dimension(:), pointer :: prptr

   integer(C_INT) :: dummy_int, i
   integer(C_LONG) :: dummy_long

   interface
      type(C_PTR) function malloc(sz) bind(C)
         use iso_c_binding
         integer(C_SIZE_T), value :: sz
      end function malloc
   end interface

   if(base.eq.0) then
      allocate(ptr2(n+1), row2(ptr(n+1)))
      ptr2(1:n+1) = ptr(1:n+1) + 1
      row2(1:ptr2(n+1)-1) = row(1:ptr2(n+1)-1) + 1
      perm(1:n) = perm(1:n) + 1
      call basic_analyse(n, ptr2, row2, perm, nnodes, fsptr, fsparent, frptr, &
         frlist, nemin, spral_core_analyse_basic_analyse, stat, nfact, nflops)
      perm(1:n) = perm(1:n) - 1
      ! NB: other outputs are corrected for base as we copy them out below
   else
      call basic_analyse(n, ptr, row, perm, nnodes, fsptr, fsparent, frptr, &
         frlist, nemin, spral_core_analyse_basic_analyse, stat, nfact, nflops)
   endif

   if(allocated(fsptr)) then
      csptr = malloc(size(fsptr)*C_SIZEOF(dummy_int))
      call C_F_POINTER(csptr, psptr, shape(fsptr))
      psptr(:) = fsptr(:) - (1-base)
   endif
   if(allocated(fsparent)) then
      csparent = malloc(size(fsparent)*C_SIZEOF(dummy_int))
      call C_F_POINTER(csparent, psparent, shape(fsparent))
      psparent(:) = fsparent(:) - (1-base)
   endif
   if(allocated(frptr)) then
      crptr = malloc(size(frptr)*C_SIZEOF(dummy_long))
      call C_F_POINTER(crptr, prptr, shape(frptr))
      prptr(:) = frptr(:) - (1-base)
   endif
   if(allocated(frlist)) then
      crlist = malloc(size(frlist)*C_SIZEOF(dummy_int))
      call C_F_POINTER(crlist, prlist, shape(frlist))
      prlist(:) = frlist(:) - (1-base)
   endif
end function spral_core_analyse_basic_analyse
