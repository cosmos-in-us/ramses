!GILEE
subroutine dump_density_fields(current_step)
  use amr_commons
  use poisson_commons
  implicit none

  integer, intent(in) :: current_step
  integer :: unit, i, l, ind, igrid, iskip
  integer :: icell, icell_father, igrid_son
  integer :: level_d, level_u
  real(dp):: w_val, dx_l, dx_level_u

  character(len=80) :: filename

  write(filename, '(A,I0.5,A,I0.5,A)') 'debug_density_', current_step, '_', myid, '.txt'
  open(newunit=unit, file=trim(filename), status='replace', recl=1024)
  write(unit, *) '# xg, yg, zg, igrid, iskip, icell_father, igrid_son, level, level_u, dx_level_u, w, rho, rho_father, rho_gravity'

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

  ! Loop over all levels
  do l=nlevelmax, levelmin, -1
    if(numbtot(1,l) > 0) then
      w_val = -1.0d0
      if(level_d > 0 .and. l>=level_u) then
        dx_level_u = 0.5d0**level_u * (boxlen_ini * 1.0d5 / h0) * aexp ! pkpc
        w_val = log(dr_smooth / dx_level_u) / log(2.0d0)
      end if

      do ind=1, twotondim
        iskip = ncoarse + (ind-1)*ngridmax
        do i=1, active(l)%ngrid
          igrid = active(l)%igrid(i)
          icell = igrid + iskip
          icell_father = father(icell)
          igrid_son = son(icell)

          write(unit, *) xg(igrid, 1), xg(igrid, 2), xg(igrid, 3), &
                               igrid, iskip, icell_father, igrid_son, &
                               l, level_u, dx_level_u, w_val, &
                               rho(icell), rho_father(icell), rho_gravity(icell)
        end do
      end do
    end if
  end do

  close(unit)

end subroutine dump_density_fields
!GILEE
