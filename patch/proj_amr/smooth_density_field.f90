!GILEE
subroutine smooth_density_field(ilevel)
  use amr_commons
  use poisson_commons
  implicit none

  integer, intent(in) :: ilevel
  !-----------------------------------------------------------------------------
  ! This routine reads from 'rho' and computes new source term for Poisson 
  ! equation 'rho_gravity' by interpolating at a physical scale 'dr_smooth'
  ! to prevent artificial gravitational collapse below that scale. 
  ! The method is described in the Appendix of Gnedin (2016).
  !-----------------------------------------------------------------------------

  integer :: i, ind, iskip, igrid
  integer :: icell, icell_father
  integer :: l, level_d, level_u
  real(dp) :: w, dx_l, dx_level_u

  ! Do nothing if dr_smooth is not set
  if(dr_smooth <= 0.0d0) return

  ! Initialize rho_gravity with the original rho for current level
  do ind=1, twotondim
    iskip = ncoarse + (ind-1)*ngridmax
    do i=1, active(ilevel)%ngrid
      igrid = active(ilevel)%igrid(i)
      icell = igrid + iskip
      icell_father = father(icell)
      
      rho_father(icell) = rho(icell_father)
      rho_gravity(icell) = rho(icell) ! deep copy
    end do
  end do

  ! Find bracketing levels Ld and Lu
  level_d = 0
  do l=1, nlevelmax
    dx_l = 0.5d0**dble(l) * (boxlen_ini * 1.0d5 / h0) * aexp ! pkpc
    if(dx_l <= dr_smooth)then
      level_u = l
      level_d = l-1
      exit
    end if
  end do

  ! If dr_smooth is smaller than the finest cell sizse, do nothing
  if(level_d == 0) return

  ! Only perform interpolating for the relevant levels
  if(ilevel < level_u) return

  ! Calculate the weighting factor w
  dx_level_u = 0.5d0**level_u * (boxlen_ini * 1.0d5 / h0) * aexp ! pkpc
  w = log(dr_smooth / dx_level_u) / log(2.0d0)
  
  ! Smooth density fields
  do ind=1, twotondim 
    iskip = ncoarse + (ind-1)*ngridmax
    do i=1, active(ilevel)%ngrid
      igrid = active(ilevel)%igrid(i)
      icell = igrid + iskip
      icell_father = father(icell)

      rho_father(icell) = rho(icell_father)

      if(ilevel == level_u) then
        rho_gravity(icell) = (1.0d0 - w) * rho(icell) + w * rho_father(icell)
      elseif(ilevel > level_u) then 
        rho_gravity(icell) = rho_gravity(icell_father)
      end if
    end do
  end do

end subroutine smooth_density_field
!GILEE
