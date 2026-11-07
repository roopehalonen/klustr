!> KLUSTR rates for the conditions of the ACDC output, called by klustr_rates.sh.
!> Arguments: N and the number of conditions. Input, one condition per line: T, CS_A, C_A^mono, C_B, dG_1/kB T.
!> Settings of the comparison: sulfuric acid and dimethylamine, eta = 1, gamma = -1.7, open system, a cluster counted
!> as formed at its first collision with a molecule past N pairs, as in ACDC.
program klustr_rates
  use klustr
  implicit none
  character(len=32) :: arg
  integer :: n, ncell, i
  real(wp), allocatable :: x(:,:), j(:)

  call get_command_argument(1, arg); read(arg, *) n
  call get_command_argument(2, arg); read(arg, *) ncell
  allocate(x(5, ncell), j(ncell))
  read(*, *) x
  call klustr_setup(n, -1.7_wp, 98.08_wp, 45.08_wp, 1830.0_wp, 680.0_wp, 1.0_wp, closed=.false., first_collision=.true.)
  call klustr_rate_block(ncell, x(3, :), x(4, :), x(2, :), x(1, :), x(5, :), j)
  do i = 1, ncell
    write(arg, '(es15.8)') j(i)
    write(*, '(a)') trim(adjustl(arg))
  end do
end program klustr_rates
