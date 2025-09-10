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

  integer :: i, ind, iskip, ncache, ngrid
  integer :: igrid
  integer :: icell, icell_father
  integer :: l, level_d, level_u
  real(dp) :: w, dx_l, dx_level_u

  ! Do nothing if smooth_gravity=.false. or dr_smooth is not set
  if(.not. smooth_gravity .or. dr_smooth <= 0.0d0) return

  ! 1. Initialize rho_for_gravity with the original rho for current level
  do ind=1, twotondim
    iskip = ncoarse + (ind-1)*ngridmax
    do i=1, active(ilevel)%ngrid
      igrid = active(ilevel)%igrid(i)
      icell = igrid + iskip
      rho_for_gravity(icell) = rho(icell)
    end do
  end do

  ! 2. Find bracketing levels Ld and Lu
  level_d = 0
  do l=1, nlevelmax
    dx_l = 0.5d0**l * (boxlen_ini * 1000 / h0 * 100) * aexp ! pkpc
    if(dx_l <= dr_smooth) then
      level_u = l
      level_d = l-1
      exit
    end if
  end do

  ! If dr_smooth is smaller than the finest cell size, do nothing
  if(level_d == 0) return

  ! Only perform smoothing for the relevant levels
  if(ilevel < level_u) return

  ! 3. Calculate the weighting factor w
  dx_level_u = 0.5d0**level_u * (boxlen_ini * 1000 / h0 * 100) * aexp ! pkpc
  w = log(dr_smooth / dx_level_u) / log(2.0d0)

  ! 4. Smooth density fields
  ncache = active(ilevel)%ngrid
  do i=1, ncache
    igrid = active(ilevel)%igrid(i)
    icell_father = father(igrid)

    if(ilevel == level_u) then
      do ind=1, twotondim
        iskip = ncoarse + (ind-1)*ngridmax
        icell = igrid + iskip
        rho_father(icell) = rho(icell_father)
        rho_for_gravity(icell) = (1.0d0 - w) * rho(icell) + w * rho_father(icell)
      end do
    elseif(ilevel > level_u) then
      do ind=1, twotondim
        iskip = ncoarse + (ind-1)*ngridmax
        icell = igrid + iskip
        rho_father(icell) = rho(icell_father)
        rho_for_gravity(icell) = rho_for_gravity(icell_father)
      end do
    end if
  end do

end subroutine smooth_density_field
!GILEE
