integer(C_INT) function spral_core_analyse_basic_analyse(n, ptr, row, perm, &
      nnodes, csptr, csparent, crptr, crlist, nemin, nfact, nflops, base) &
      bind(C)
   use iso_c_binding
   use spral_core_analyse, only: basic_analyse
   implicit none

   integer(C_INT), value :: n
   integer(C_INT), dimension(*), intent(in) :: ptr
   integer(C_INT), dimension(*), intent(in) :: row
   integer(C_INT), dimension(*), intent(inout) :: perm
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

   interface
      type(C_PTR) function malloc(sz) bind(C)
         use iso_c_binding
         integer(C_SIZE_T) :: sz
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
      fsptr(:) = fsptr(:) - 1
      fsparent(:) = fsparent(:) - 1
      frptr(:) = frptr(:) - 1
      frlist(:) = frlist(:) - 1
   else
      call basic_analyse(n, ptr, row, perm, nnodes, fsptr, fsparent, frptr, &
         frlist, nemin, spral_core_analyse_basic_analyse, stat, nfact, nflops)
   endif

   if(allocated(fsptr)) then
      csptr = malloc(C_SIZEOF(fsptr))
      call C_F_POINTER(csptr, psptr, shape(fsptr))
      psptr(:) = fsptr(:)
   endif
   if(allocated(fsparent)) then
      csparent = malloc(C_SIZEOF(fsparent))
      call C_F_POINTER(csparent, psparent, shape(fsparent))
      psparent(:) = fsparent(:)
   endif
   if(allocated(frptr)) then
      crptr = malloc(C_SIZEOF(frptr))
      call C_F_POINTER(crptr, prptr, shape(frptr))
      prptr(:) = frptr(:)
   endif
   if(allocated(frlist)) then
      crlist = malloc(C_SIZEOF(frlist))
      call C_F_POINTER(crlist, prlist, shape(frlist))
      prlist(:) = frlist(:)
   endif
      

end function spral_core_analyse_basic_analyse
