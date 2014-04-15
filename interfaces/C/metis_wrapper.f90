integer(C_INT) function spral_metis_order(n, ptr, row, perm, invp, base) &
      bind(C)
   use iso_c_binding
   use spral_metis_wrapper, only: metis_order
   implicit none

   integer(C_INT), value :: n
   integer(C_INT), dimension(*), intent(in) :: ptr
   integer(C_INT), dimension(*), intent(in) :: row
   integer(C_INT), dimension(n), intent(out) :: perm
   integer(C_INT), dimension(n), intent(out) :: invp
   integer(C_INT), value :: base

   integer :: stat
   integer, dimension(:), allocatable :: ptr2, row2

   if(base.eq.0) then
      allocate(ptr2(n+1), row2(ptr(n+1)))
      ptr2(1:n+1) = ptr(1:n+1) + 1
      row2(1:ptr2(n+1)-1) = row(1:ptr2(n+1)-1)
      call metis_order(n, ptr2, row2, perm, invp, spral_metis_order, stat)
      perm(:) = perm(:) - 1
      invp(:) = invp(:) - 1
   else
      call metis_order(n, ptr, row, perm, invp, spral_metis_order, stat)
   endif
end function spral_metis_order
