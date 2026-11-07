!> Test of the Fortran version against ../python/test_reference.txt (both routines) and of J = 0 for empty inputs.
!> make test
program test_klustr
  use klustr
  implicit none
  character(len=512) :: path, line
  integer :: u, ios, n, iclosed, ifirst, ncase, nfail
  real(wp) :: eta, x(5), jref, j1, j2(1), worst, tol
  real(8) :: xd(6), etad
  real(wp) :: jb(3)

  ! J is zero before klustr_setup, and zero, not NaN, for non-positive concentrations in a block
  nfail = 0
  if (klustr_rate(1e6_wp, 1e7_wp, 1e-3_wp, 278.0_wp, -29.0_wp) /= 0) nfail = nfail + 1
  call klustr_rate_block(3, [1e6_wp, 1e6_wp, 1e6_wp], [1e7_wp, 1e7_wp, 1e7_wp], [1e-3_wp, 1e-3_wp, 1e-3_wp], &
                         [278.0_wp, 278.0_wp, 278.0_wp], [-29.0_wp, -29.0_wp, -29.0_wp], jb)
  if (any(jb /= 0)) nfail = nfail + 1
  call klustr_setup(8, -1.7_wp, 98.08_wp, 45.08_wp, 1830.0_wp, 680.0_wp, 3.0_wp)
  call klustr_rate_block(3, [-1e12_wp, 1e6_wp, 0.0_wp], [1e7_wp, -1e12_wp, 1e7_wp], [1e-3_wp, 1e-3_wp, 1e-3_wp], &
                         [278.0_wp, 278.0_wp, 278.0_wp], [-29.0_wp, -29.0_wp, -29.0_wp], jb)
  if (any(jb /= 0)) nfail = nfail + 1
  write(*, '(a, i0, a)') 'zero before setup and for non-positive inputs: ', nfail, ' failed'
  if (nfail > 0) error stop 1

  path = '../python/test_reference.txt'
  if (command_argument_count() > 0) call get_command_argument(1, path)
  tol = merge(1e-12_wp, 1e-4_wp, wp == kind(1.0d0))
  open(newunit=u, file=trim(path), status='old', action='read')
  ncase = 0; nfail = 0; worst = 0
  do
    read(u, '(a)', iostat=ios) line
    if (ios /= 0) exit
    if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
    read(line, *) n, iclosed, ifirst, etad, xd
    eta = real(etad, wp); x = real(xd(1:5), wp); jref = real(xd(6), wp)
    call klustr_setup(n, -1.7_wp, 98.08_wp, 45.08_wp, 1830.0_wp, 680.0_wp, eta, &
                      closed=(iclosed == 1), first_collision=(ifirst == 1))
    j1 = klustr_rate(x(1), x(2), x(3), x(4), x(5))
    call klustr_rate_block(1, x(1:1), x(2:2), x(3:3), x(4:4), x(5:5), j2)
    ncase = ncase + 1
    worst = max(worst, abs(j1/jref - 1), abs(j2(1)/jref - 1))
    if (abs(j1/jref - 1) > tol .or. abs(j2(1)/jref - 1) > tol) nfail = nfail + 1
  end do
  close(u)
  write(*, '(i0, a, es9.1, a, i0, a)') ncase, ' reference cases, largest relative difference', worst, ', ', &
       nfail, ' failed'
  if (nfail > 0 .or. ncase == 0) error stop 1
end program test_klustr
