!> KLUSTR v1.0: formation rate of new particles from two compounds A and B whose pair is the critical cluster or has no
!> barrier (Algorithm 1 of Halonen, R., Geosci. Model Dev., submitted, 2026). Fortran 2008, no dependencies.
!>   call klustr_setup(8, -1.7_wp, 98.08_wp, 45.08_wp, 1830.0_wp, 680.0_wp, 3.0_wp)   ! once
!>   J = klustr_rate(a_mono, b_free, cs_a, t, dg_kt)                                   ! cm-3 s-1, elemental
!> Copyright (c) 2026 Roope Halonen. MIT licence.
module klustr
  use, intrinsic :: iso_fortran_env, only: real32, real64
  implicit none
  private
#ifdef SINGLE
  integer, parameter, public :: wp = real32
#else
  integer, parameter, public :: wp = real64
#endif
  public :: klustr_setup, klustr_rate, klustr_rate_block, klustr_reduced

  integer, parameter :: blk = 64                       ! cells per block in klustr_rate_block
  real(wp), parameter :: kb = 1.380649e-23_wp, p0 = 101325.0_wp

  integer, save :: nmax = 0                            ! N, the boundary in pairs
  logical, save :: closed_sys = .true.                 ! pair limited to its saturation concentration
  logical, save :: first_coll = .false.                ! counted at the first collision with a molecule past N pairs
  real(wp), allocatable, save :: rho(:,:), rhoa(:), rhob(:), pab(:), sig(:)
  real(wp), save :: c0f = 0, babf = 0                  ! sat = c0f/T exp(dG_kT); beta_AB = babf sqrt(T)

contains

  !> Constants of one chemistry; call once. n: number of pairs past which a cluster is formed, gamma: exponent of the
  !> loss rate coefficient in the particle diameter (negative), ma, mb: molecular masses of A and B (amu), densa, densb:
  !> bulk densities of A and B (kg m-3), eta: factor on all collision coefficients (1: hard spheres, as in
  !> the verification of the paper; about 3 for neutral acid-base systems), closed: closed system (default true),
  !> first_collision: a cluster counted as formed at its first collision with a molecule past N pairs, as in ACDC
  !> (default false: when it has grown by a whole pair, Algorithm 1).
  subroutine klustr_setup(n, gamma, ma, mb, densa, densb, eta, closed, first_collision)
    integer, intent(in) :: n
    real(wp), intent(in) :: gamma, ma, mb, densa, densb, eta
    logical, intent(in), optional :: closed, first_collision
    real(real64) :: va, vb, ms(-1:n), rs(-1:n), gab, e   ! double precision: products of masses in kg underflow real32
    real(real64), parameter :: kb8 = 1.380649e-23_real64, amu8 = 1.66053907e-27_real64, na8 = 6.02214076e23_real64
    real(real64), parameter :: pi8 = 3.14159265358979324_real64
    integer :: i, j
    if (n < 1) error stop 'klustr_setup: n must be at least 1'
    if (ma <= 0 .or. mb <= 0 .or. densa <= 0 .or. densb <= 0) &
      error stop 'klustr_setup: masses and densities must be positive'
    e = real(eta, real64)
    closed_sys = .true.
    if (present(closed)) closed_sys = closed
    first_coll = .false.
    if (present(first_collision)) first_coll = first_collision
    nmax = n
    if (allocated(rho)) deallocate(rho, rhoa, rhob, pab, sig)
    allocate(rho(n, n), rhoa(n), rhob(n), pab(n), sig(n))
    va = real(ma, real64)*1.0e-3_real64/(real(densa, real64)*na8)
    vb = real(mb, real64)*1.0e-3_real64/(real(densb, real64)*na8)
    ms(-1) = real(ma, real64)*amu8; rs(-1) = rad(va)   ! index -1: A, 0: B, n >= 1: cluster of n pairs
    ms(0) = real(mb, real64)*amu8;  rs(0) = rad(vb)
    do i = 1, n
      ms(i) = i*real(ma + mb, real64)*amu8; rs(i) = rad(i*(va + vb))
    end do
    gab = g(-1, 0)
    do i = 1, n
      rhoa(i) = real(g(-1, i)/gab, wp); rhob(i) = real(g(0, i)/gab, wp); pab(i) = real(g(-1, i)*g(0, i)/gab**2, wp)
      sig(i) = real((i*(va + vb)/va)**(real(gamma, real64)/3), wp)
      do j = 1, n
        rho(i, j) = real(g(i, j)/gab, wp)
      end do
    end do
    c0f = p0/kb*1.0e-6_wp
    babf = real(e*1.0e6_real64*gab, wp)
  contains
    pure real(real64) function rad(v)
      real(real64), intent(in) :: v
      rad = (3*v/(4*pi8))**(1.0_real64/3)
    end function rad
    pure real(real64) function g(i, j)                 ! beta/sqrt(T) of species i and j (m3 s-1 K-1/2)
      integer, intent(in) :: i, j
      g = sqrt(8*kb8*(ms(i) + ms(j))/(pi8*ms(i)*ms(j)))*pi8*(rs(i) + rs(j))**2
    end function g
  end subroutine klustr_setup

  !> Reduced formation rate F = J/(beta_AB sat^2) for a = a_mono/sat, b = b_free/sat, s = CS_A/(beta_AB sat).
  pure real(wp) function klustr_reduced(a, b, s) result(f)
    real(wp), intent(in) :: a, b, s
    real(wp) :: k(nmax), y(nmax), a1, ab, mm, ss, d
    integer :: n, i, j, p
    f = 0
    if (nmax < 1 .or. a <= 0 .or. b <= 0) return     ! before klustr_setup, or no A or B
    a1 = a; ab = a*b
    do p = 1, 2                                        ! pair: quadratic in y1, kappa_1 from the previous pass
      mm = 1 + pab(1)*a1*b/(rhoa(1)*a1 + rhob(1)*b + s*sig(1)) + s*sig(1) + b
      y(1) = 2*ab/(mm + sqrt(mm*mm + 4*rho(1, 1)*ab))
      if (closed_sys) y(1) = min(y(1), 1.0_wp)
      a1 = a - y(1)
    end do
    do n = 1, nmax
      k(n) = pab(n)*a1*b/(rhoa(n)*a1 + rhob(n)*b + s*sig(n))
    end do
    if (first_coll) k(nmax) = rhoa(nmax)*a1 + rhob(nmax)*b
    do n = 2, nmax                                     ! larger clusters in order of size
      ss = k(n-1)*y(n-1); d = k(n) + s*sig(n)
      do i = 1, n - 1
        ss = ss + 0.5_wp*rho(i, n-i)*y(i)*y(n-i)
        d = d + rho(i, n)*y(i)
      end do
      y(n) = 2*ss/(d + sqrt(d*d + 4*rho(n, n)*ss))
    end do
    f = k(nmax)*y(nmax)                                ! growth past N and collisions that pass N
    do i = 1, nmax
      do j = max(i, nmax + 1 - i), nmax
        if (i == j) then
          f = f + 0.5_wp*rho(i, j)*y(i)*y(j)
        else
          f = f + rho(i, j)*y(i)*y(j)
        end if
      end do
    end do
  end function klustr_reduced

  !> Formation rate J (cm-3 s-1) past n pairs; elemental. a_mono: monomer concentration of A (all clusters with one
  !> A), b_free: free concentration of B (cm-3); cs_a: loss rate coefficient of A (s-1); t (K); dg_kt: standard
  !> formation free energy of the pair at 1 atm over kB T.
  elemental real(wp) function klustr_rate(a_mono, b_free, cs_a, t, dg_kt) result(j)
    real(wp), intent(in) :: a_mono, b_free, cs_a, t, dg_kt
    real(wp) :: sat, bab
    j = 0
    if (nmax < 1) return                               ! klustr_setup not called
    sat = c0f/t*exp(dg_kt)
    bab = babf*sqrt(t)
    j = bab*sat*sat*klustr_reduced(a_mono/sat, b_free/sat, cs_a/(bab*sat))
  end function klustr_rate

  !> The same for ncell cells, in blocks with the cell index innermost; a mask replaces the early return.
  subroutine klustr_rate_block(ncell, a_mono, b_free, cs_a, t, dg_kt, j)
    integer, intent(in) :: ncell
    real(wp), intent(in) :: a_mono(ncell), b_free(ncell), cs_a(ncell), t(ncell), dg_kt(ncell)
    real(wp), intent(out) :: j(ncell)
    real(wp) :: sat(blk), bab(blk), msk(blk), ra(blk), rb(blk), rs(blk), a1(blk), ab(blk), f(blk)
    real(wp) :: k(blk, nmax), y(blk, nmax), vs(blk), vd(blk), mm, w
    integer :: q0, q, m, n, i, i2, p
    if (nmax < 1) then                                 ! klustr_setup not called
      j = 0
      return
    end if
    do q0 = 1, ncell, blk
      m = min(blk, ncell - q0 + 1)
      do q = 1, m
        sat(q) = c0f/t(q0+q-1)*exp(dg_kt(q0+q-1))
        bab(q) = babf*sqrt(t(q0+q-1))
        msk(q) = merge(1.0_wp, 0.0_wp, a_mono(q0+q-1) > 0 .and. b_free(q0+q-1) > 0)
        ra(q) = merge(a_mono(q0+q-1)/sat(q), 1.0_wp, msk(q) > 0)   ! masked cells compute on safe values
        rb(q) = merge(b_free(q0+q-1)/sat(q), 1.0_wp, msk(q) > 0)
        rs(q) = cs_a(q0+q-1)/(bab(q)*sat(q))
        a1(q) = ra(q); ab(q) = ra(q)*rb(q)
      end do
      do p = 1, 2
        do q = 1, m
          mm = 1 + pab(1)*a1(q)*rb(q)/(rhoa(1)*a1(q) + rhob(1)*rb(q) + rs(q)*sig(1)) + rs(q)*sig(1) + rb(q)
          y(q, 1) = 2*ab(q)/(mm + sqrt(mm*mm + 4*rho(1, 1)*ab(q)))
          if (closed_sys) y(q, 1) = min(y(q, 1), 1.0_wp)
          a1(q) = ra(q) - y(q, 1)
        end do
      end do
      do n = 1, nmax
        do q = 1, m
          k(q, n) = pab(n)*a1(q)*rb(q)/(rhoa(n)*a1(q) + rhob(n)*rb(q) + rs(q)*sig(n))
        end do
      end do
      if (first_coll) then
        do q = 1, m
          k(q, nmax) = rhoa(nmax)*a1(q) + rhob(nmax)*rb(q)
        end do
      end if
      do n = 2, nmax                                   ! the cell loop innermost everywhere
        do q = 1, m
          vs(q) = k(q, n-1)*y(q, n-1); vd(q) = k(q, n) + rs(q)*sig(n)
        end do
        do i = 1, n - 1
          w = 0.5_wp*rho(i, n-i)
          do q = 1, m
            vs(q) = vs(q) + w*y(q, i)*y(q, n-i)
            vd(q) = vd(q) + rho(i, n)*y(q, i)
          end do
        end do
        do q = 1, m
          y(q, n) = 2*vs(q)/(vd(q) + sqrt(vd(q)*vd(q) + 4*rho(n, n)*vs(q)))
        end do
      end do
      do q = 1, m
        f(q) = k(q, nmax)*y(q, nmax)
      end do
      do i = 1, nmax
        do i2 = max(i, nmax + 1 - i), nmax
          w = rho(i, i2); if (i == i2) w = 0.5_wp*w
          do q = 1, m
            f(q) = f(q) + w*y(q, i)*y(q, i2)
          end do
        end do
      end do
      do q = 1, m
        j(q0+q-1) = msk(q)*bab(q)*sat(q)*sat(q)*f(q)
      end do
    end do
  end subroutine klustr_rate_block

end module klustr
