module spral_ssids_cpu_solve
   use spral_ssids_datatypes
   implicit none

   private
   public :: solve_fwd_mult, &
             solve_fwd_one, &
             solve_diag_mult, &
             solve_diag_one, &
             solve_bwd_mult, &
             solve_bwd_one

contains

!*************************************************************************
!
! Forward substitution single rhs
!
subroutine solve_fwd_one(pos_def, rlist, x, blkm, blkn, nelim, nd, lcol, &
      lperm, xlocal, map)
   logical, intent(in) :: pos_def
   integer, dimension(*), intent(in) :: rlist
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(in) :: blkm
   integer, intent(in) :: blkn
   integer, intent(in) :: nelim
   integer, intent(in) :: nd
   real(wp), dimension(*), intent(in) :: lcol
   integer, dimension(*), intent(in) :: lperm
   real(wp), dimension(*), intent(out) :: xlocal
   integer, dimension(*), intent(out) :: map

   integer(long) :: ip, ip2
   integer :: i, j, k
   integer :: rp1
   real(wp) :: ri, ri2
   
   do i = 1, blkn
      map(i) = lperm(i)
   end do
   k = 1+blkn-nd
   do i = blkn+1, blkm
      map(i) = rlist(k)
      k = k + 1
   end do
   
   ! Copy eliminated variables into xlocal
   do i = 1, nelim
      rp1 = map(i)
      xlocal(i) = x(rp1)
   end do

   ! Perform the solve
   if (blkm.gt.10 .and. nelim.gt.4) then
      ! Work with xlocal

      if(pos_def) then
         call dtrsv('L','N','N', nelim, lcol, blkm, xlocal, 1)
      else
         call dtrsv('L','N','U', nelim, lcol, blkm, xlocal, 1)
      end if

      if (blkm-nelim.gt.0) then
         call dgemv('N', blkm-nelim, nelim, -one, lcol(nelim+1), blkm, &
            xlocal, 1, zero, xlocal(nelim+1), 1)
         ! Add contribution into x
         ! Delays first
         do i = nelim+1, blkm
            rp1 = map(i)
            x(rp1) = x(rp1) + xlocal(i)
         end do
      end if
   else
      do i = 1, nelim-1, 2
         ip = (i-1)*blkm
         ip2 = i*blkm
         if(pos_def) xlocal(i) = xlocal(i) / lcol(ip+i)
         ri = xlocal(i)
         xlocal(i+1) = xlocal(i+1) - ri * lcol(ip+i+1)
         if(pos_def) xlocal(i+1) = xlocal(i+1) / lcol(ip2+i+1)
         ri2 = xlocal(i+1)
         do j = i+2, nelim
            xlocal(j) = xlocal(j) - ri * lcol(ip+j) - ri2 * lcol(ip2+j)
         end do
         do j = nelim+1, blkm
            rp1 = map(j)
            x(rp1) = x(rp1) - ri * lcol(ip+j) - ri2 * lcol(ip2+j)
         end do
      end do
      if(mod(nelim,2).eq.1) then
         ip = (i-1)*blkm
         ip2 = i*blkm
         if(pos_def) xlocal(i) = xlocal(i) / lcol(ip+i)
         ri = xlocal(i)
         ! Commented loop redundant as i=nelim+1
         !do j = i+1, nelim
         !   xlocal(j) = xlocal(j) - ri * lcol(ip+j)
         !end do
         do j = nelim+1, blkm
            rp1 = map(j)
            x(rp1) = x(rp1) - ri * lcol(ip+j)
         end do
      end if
   end if

   ! Copy solution back from xlocal
   do i = 1, nelim
      rp1 = map(i)
      x(rp1) = xlocal(i)
   end do
end subroutine solve_fwd_one

!*************************************************************************
!
! Forward substitution multiple rhs
!
subroutine solve_fwd_mult(pos_def, rlist, nrhs, x, ldx, blkm, blkn, &
      nelim, nd, lcol, lperm, xlocal, map)
   logical, intent(in) :: pos_def
   integer, dimension(*), intent(in) :: rlist
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldx
   real(wp), dimension(ldx,*), intent(inout) :: x
   integer, intent(in) :: blkm
   integer, intent(in) :: blkn
   integer, intent(in) :: nelim
   integer, intent(in) :: nd
   real(wp), dimension(*), intent(in) :: lcol
   integer, dimension(*), intent(in) :: lperm
   real(wp), dimension(blkm,*), intent(out) :: xlocal
   integer, dimension(*), intent(out) :: map

   integer(long) :: ip
   integer :: i, j, k, r
   integer :: rp1
   real(wp) :: ri
   
   do i = 1, blkn
      map(i) = lperm(i)
   end do
   k = 1+blkn-nd
   do i = blkn+1, blkm
      map(i) = rlist(k)
      k = k + 1
   end do
   
   ! Copy eliminated variables into xlocal
   do r = 1, nrhs
      do i = 1, nelim
         rp1 = map(i)
         xlocal(i,r) = x(rp1, r)
      end do
   end do

   ! Perform the solve
   if (blkm.gt.10 .and. nelim.gt.4) then
      ! Work with xlocal

      if(pos_def) then
         call dtrsm('Left', 'Lower', 'Non-Trans', 'Non-Unit', nelim, nrhs, &
            one, lcol, blkm, xlocal, blkm)
      else
         call dtrsm('Left', 'Lower', 'Non-Trans', 'Unit', nelim, nrhs, &
            one, lcol, blkm, xlocal, blkm)
      end if

      if (blkm-nelim.gt.0) then
         call dgemm('N', 'N', blkm-nelim, nrhs, nelim, -one, &
            lcol(nelim+1), blkm, xlocal, blkm, zero, &
            xlocal(nelim+1,1), blkm)
         ! Add contribution into x
         do r = 1, nrhs
            ! Delays first
            do i = nelim+1, blkn
               rp1 = map(i)
               x(rp1,r) = x(rp1,r) + xlocal(i,r)
            end do
            ! Expected rows
            do j = blkn+1, blkm
               rp1 = map(j)
               x(rp1,r) = x(rp1,r) + xlocal(j,r)
            end do
         end do
      end if
   else
      do r = 1, nrhs
         do i = 1, nelim
            ip = (i-1)*blkm
            if(pos_def) xlocal(i,r) = xlocal(i,r) / lcol(ip+i)
            ri = xlocal(i,r)
            do j = i+1, nelim
               xlocal(j,r) = xlocal(j,r) - ri * lcol(ip+j)
            end do
            do j = nelim+1, blkm
               rp1 = map(j)
               x(rp1,r) = x(rp1,r) - ri * lcol(ip+j)
            end do
         end do
      end do
   end if

   ! Copy solution back from xlocal
   do r = 1, nrhs
      do i = 1, nelim
         rp1 = map(i)
         x(rp1,r) = xlocal(i,r)
      end do
   end do
end subroutine solve_fwd_mult

!*************************************************************************
!
! Back substitution (with diagonal solve) single rhs
!
subroutine solve_bwd_one(pos_def, job, rlist, x, blkm, blkn, nelim, &
      nd, lcol, d, lperm, xlocal, map)
   logical, intent(in) :: pos_def
   integer, intent(in) :: job  ! used to indicate whether diag. sol. required
      ! job = 3 : backsubs only ((PL)^Tx = b)
      ! job = 0 or 4 : diag and backsubs (D(PL)^Tx = b)
   integer, dimension(*), intent(in) :: rlist
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(in) :: blkm
   integer, intent(in) :: blkn
   integer, intent(in) :: nelim
   integer, intent(in) :: nd
   real(wp), dimension(*), intent(in) :: lcol
   real(wp), dimension(2*nelim) :: d
   integer, dimension(*), intent(in) :: lperm
   real(wp), dimension(*), intent(out) :: xlocal
   integer, dimension(*), intent(out) :: map

   integer(long) :: ip
   integer :: i, j, k
   integer :: rp1, rp2

   do i = 1, blkn
      map(i) = lperm(i)
   end do
   k = 1+blkn-nd
   do i = blkn+1, blkm
      map(i) = rlist(k)
      k = k + 1
   end do

   if (job.eq.SSIDS_SOLVE_JOB_BWD .or. pos_def) then
      ! No diagonal solve. Copy eliminated variables into xlocal.
      do i = 1, nelim
         rp1 = map(i)
         xlocal(i) = x(rp1)
      end do
   else 
      ! Copy eliminated vars into xlocal while performing diagonal solve
      j = 1
      do while(j .le. nelim)
         if (d(2*j).ne.0) then
            ! 2x2 pivot
            rp1 = map(j)
            rp2 = map(j+1)
            xlocal(j)   = d(2*j-1) * x(rp1) + &
                          d(2*j)   * x(rp2)
            xlocal(j+1) = d(2*j)   * x(rp1) + &
                          d(2*j+1) * x(rp2)
            j = j + 2
         else
            ! 1x1 pivot
            if (d(2*j-1).eq.0.0_wp) then
               ! Zero pivot column
               xlocal(j) = 0.0_wp
            else
               ! Proper pivot
               rp1 = map(j)
               xlocal(j) = x(rp1) * d(2*j-1)
            end if
            j = j + 1
         end if
      end do
   end if

   ! Perform the solve
   if (blkm.gt.10 .and. nelim.gt.4) then
      ! Delays
      do i = nelim+1, blkn
         rp1 = map(i)
         xlocal(i) = x(rp1)
      end do
      ! Expected rows
      do j = blkn+1, blkm
         rp1 = map(j)
         xlocal(j) = x(rp1)
      end do
      if (blkm-nelim.gt.0) then
         call dgemv('T', blkm-nelim, nelim, -one, lcol(nelim+1), blkm, &
            xlocal(nelim+1), 1, one, xlocal, 1)
      end if

      if(pos_def) then
         call dtrsv('L','T','N', nelim, lcol, blkm, xlocal, 1)
      else
         call dtrsv('L','T','U', nelim, lcol, blkm, xlocal, 1)
      end if

      ! Copy solution back from xlocal
      do i = 1, nelim
         rp1 = map(i)
         x(rp1) = xlocal(i)
      end do
   else
      ! Do update with indirect addressing
      do i = 1, nelim
         ip = (i-1)*blkm
         do j = nelim+1, blkm
            rp1 = map(j)
            xlocal(i) = xlocal(i) - x(rp1) * lcol(ip+j)
         end do
      end do

      ! Solve with direct addressing
      if(pos_def) then
         do i = nelim, 1, -1
            ip = (i-1)*blkm
            rp1 = map(i)
            xlocal(i) = xlocal(i) - &
               sum(xlocal(i+1:nelim) * lcol(ip+i+1:ip+nelim))
            xlocal(i) = xlocal(i) / lcol(ip+i)
            x(rp1) = xlocal(i)
         end do
      else
         do i = nelim, 1, -1
            ip = (i-1)*blkm
            rp1 = map(i)
            xlocal(i) = xlocal(i) - &
               sum(xlocal(i+1:nelim) * lcol(ip+i+1:ip+nelim))
            x(rp1) = xlocal(i)
         end do
      end if
   end if
end subroutine solve_bwd_one

!*************************************************************************
!
! Back substitution (with diagonal solve) multiple rhs
!
subroutine solve_bwd_mult(pos_def, job, rlist, nrhs, x, ldx, blkm, &
      blkn, nelim, nd, lcol, d, lperm, xlocal, map)
   logical, intent(in) :: pos_def
   integer, intent(in) :: job  ! used to indicate whether diag. sol. required
      ! job = 3 : backsubs only ((PL)^Tx = b)
      ! job = 0 or 4 : diag and backsubs (D(PL)^Tx = b)
   integer, dimension(*), intent(in) :: rlist
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldx
   real(wp), dimension(ldx,*), intent(inout) :: x
   integer, intent(in) :: blkm
   integer, intent(in) :: blkn
   integer, intent(in) :: nelim
   integer, intent(in) :: nd
   real(wp), dimension(*), intent(in) :: lcol
   real(wp), dimension(2*nelim) :: d
   integer, dimension(*), intent(in) :: lperm
   real(wp), dimension(blkm,*), intent(out) :: xlocal
   integer, dimension(*), intent(out) :: map

   integer(long) :: ip
   integer :: i, j, k, r
   integer :: rp1, rp2

   do i = 1, blkn
      map(i) = lperm(i)
   end do
   k = 1+blkn-nd
   do i = blkn+1, blkm
      map(i) = rlist(k)
      k = k + 1
   end do

   if (job == SSIDS_SOLVE_JOB_BWD .or. pos_def) then 
      ! no diagonal solve. Copy eliminated variables into xlocal
      do r = 1, nrhs
         do i = 1, nelim
            rp1 = map(i)
            xlocal(i,r) = x(rp1, r)
         end do
      end do
   else 
      ! Copy eliminated vars into xlocal while performing diagonal solve
      do r = 1, nrhs
         j = 1
         do while(j .le. nelim)
            if (d(2*j).ne.0) then
               ! 2x2 pivot
               rp1 = map(j)
               rp2 = map(j+1)
               xlocal(j,r)   = d(2*j-1) * x(rp1,r) + &
                               d(2*j)   * x(rp2,r)
               xlocal(j+1,r) = d(2*j)   * x(rp1,r) + &
                               d(2*j+1) * x(rp2,r)
               j = j + 2
            else
               ! 1x1 pivot
               if (d(2*j-1).eq.0.0_wp) then
                  ! Zero pivot column
                  xlocal(j,r) = 0.0_wp
               else
                  ! Proper pivot
                  rp1 = map(j)
                  xlocal(j,r) = x(rp1,r) * d(2*j-1)
               end if
               j = j + 1
            end if
         end do
      end do
   end if

   ! Perform the solve
   if (blkm.gt.10 .and. nelim.gt.4) then
      do r = 1, nrhs
         ! Delays
         do i = nelim+1, blkn
            rp1 = map(i)
            xlocal(i,r) = x(rp1,r)
         end do
         ! Expected rows
         do j = blkn+1, blkm
            rp1 = map(j)
            xlocal(j,r) = x(rp1,r)
         end do
      end do
      if (blkm-nelim.gt.0) then
         call dgemm('Trans', 'Non-trans', nelim, nrhs, blkm-nelim, -one, &
            lcol(nelim+1), blkm, xlocal(nelim+1,1), blkm, one, xlocal, &
            blkm)
      end if

      if(pos_def) then
         call dtrsm('Left','Lower','Trans','Non-Unit', nelim, nrhs, one, lcol, &
            blkm, xlocal, blkm)
      else
         call dtrsm('Left','Lower','Trans','Unit', nelim, nrhs, one, lcol, &
            blkm, xlocal, blkm)
      end if
      do r = 1, nrhs
         do i = 1, nelim
            rp1 = map(i)
            x(rp1,r) = xlocal(i,r)
         end do
      end do
   else
      ! Do update with indirect addressing
      do r = 1, nrhs
         do i = 1, nelim
            ip = (i-1)*blkm
            do j = nelim+1, blkm
               rp1 = map(j)
               xlocal(i,r) = xlocal(i,r) - x(rp1,r) * lcol(ip+j)
            end do
         end do

         ! Solve with direct addressing
         if(pos_def) then
            do i = nelim, 1, -1
               ip = (i-1)*blkm
               rp1 = map(i)
               xlocal(i,r) = xlocal(i,r) - &
                  sum(xlocal(i+1:nelim,r) * lcol(ip+i+1:ip+nelim))
               xlocal(i,r) = xlocal(i,r) / lcol(ip+i)
               x(rp1,r) = xlocal(i,r)
            end do
         else
            do i = nelim, 1, -1
               ip = (i-1)*blkm
               rp1 = map(i)
               xlocal(i,r) = xlocal(i,r) - &
                  sum(xlocal(i+1:nelim,r) * lcol(ip+i+1:ip+nelim))
               x(rp1,r) = xlocal(i,r)
            end do
         end if
      end do
   end if
end subroutine solve_bwd_mult

!*************************************************************************
!
! Diagonal solve one rhs
!
subroutine solve_diag_one(x, nelim, d, lperm)
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(in) :: nelim
   real(wp), dimension(2*nelim) :: d
   integer, dimension(*), intent(in) :: lperm

   integer :: j
   integer :: rp1, rp2
   real(wp) :: temp

   j = 1
   do while(j .le. nelim)
      if (d(2*j).ne.0) then
         ! 2x2 pivot
         rp1 = lperm(j)
         rp2 = lperm(j+1)
         temp   = d(2*j-1) * x(rp1) + &
                  d(2*j)   * x(rp2)
         x(rp2) = d(2*j)   * x(rp1) + &
                  d(2*j+1) * x(rp2)
         x(rp1) = temp
         j = j + 2
      else
         ! 1x1 pivot
         rp1 = lperm(j)
         x(rp1) = x(rp1) * d(2*j-1)
         j = j + 1
      end if
   end do
end subroutine solve_diag_one

!*************************************************************************
!
! Diagonal solve multiple rhs
!
subroutine solve_diag_mult(nrhs, x, ldx, nelim, d, lperm)
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldx
   real(wp), dimension(ldx,*), intent(inout) :: x
   integer, intent(in) :: nelim
   real(wp), dimension(2*nelim) :: d
   integer, dimension(*), intent(in) :: lperm

   integer :: j, r
   integer :: rp1, rp2
   real(wp) :: temp

   do r = 1, nrhs
      j = 1
      do while(j .le. nelim)
         if (d(2*j).ne.0) then
            ! 2x2 pivot
            rp1 = lperm(j)
            rp2 = lperm(j+1)
            temp     = d(2*j-1) * x(rp1,r) + &
                       d(2*j)   * x(rp2,r)
            x(rp2,r) = d(2*j)   * x(rp1,r) + &
                       d(2*j+1) * x(rp2,r)
            x(rp1,r) = temp
            j = j + 2
         else
            ! 1x1 pivot
            rp1 = lperm(j)
            x(rp1,r) = x(rp1,r) * d(2*j-1)
            j = j + 1
         end if
      end do
   end do
end subroutine solve_diag_mult

end module spral_ssids_cpu_solve
