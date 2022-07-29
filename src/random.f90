
module random
    use precision

    implicit none
    private

    public :: RandomNumberGenerator

    type RandomNumberGenerator
        integer(di) :: state
    contains
        procedure :: constructor => rng_constructor
        procedure :: shuffle

        procedure, private :: random_uniform_0D
        procedure, private :: random_uniform_1D
        procedure, private :: random_uniform_2D
        procedure, private :: random_uniform_3D
        generic, public :: random_uniform => random_uniform_0D, random_uniform_1D, random_uniform_2D, random_uniform_3D
        procedure, private :: random_normal_0D
        procedure, private :: random_normal_1D
        procedure, private :: random_normal_2D
        procedure, private :: random_normal_3D
        generic, public :: random_normal => random_normal_0D, random_normal_1D, random_normal_2D, random_normal_3D
    end type RandomNumberGenerator
contains

    subroutine rng_constructor(self, seed)
        class(RandomNumberGenerator) :: self
        integer(di), intent(in) :: seed

        self%state = seed
    end subroutine


    ! Knuth shuffle for a list of integers
    subroutine shuffle(self, n, a)
        class(RandomNumberGenerator), intent(inout) :: self
        integer, intent(in) :: n
        integer, intent(inout) :: a(n)
        integer :: i, randpos, temp
        real(dp) :: r

        do i = n, 2, -1
            call self%random_uniform(r)
            randpos = int(r * i) + 1
            temp = a(randpos)
            a(randpos) = a(i)
            a(i) = temp
        end do

    end subroutine shuffle

    subroutine random_uniform_0D(self, x)
        class(RandomNumberGenerator), intent(inout) :: self
        real(dp), intent(out) :: x
        !$omp critical
        x = xorshift64star(self%state)
        !$omp end critical
    end subroutine random_uniform_0D

    subroutine random_uniform_1D(self, A)
        class(RandomNumberGenerator), intent(inout) :: self
        real(dp), intent(out) :: A(:)
        integer :: i
        integer :: n

        n = size(A)
        do i=1,n
            call random_uniform_0D(self, A(i))
        end do
    end subroutine

    subroutine random_uniform_2D(self, A)
        class(RandomNumberGenerator), intent(inout) :: self
        real(dp), intent(out) :: A(:,:)
        integer :: i, j
        integer :: n(2)

        n = shape(A)
        do j=1,n(2)
            do i=1,n(1)
                call random_uniform_0D(self, A(i,j))
            end do
        end do
    end subroutine

    subroutine random_uniform_3D(self, A)
        class(RandomNumberGenerator), intent(inout) :: self
        real(dp), intent(out) :: A(:,:,:)
        integer :: i, j, k
        integer :: n(3)

        n = shape(A)
        do k=1,n(3)
            do j=1,n(2)
                do i=1,n(1)
                    call random_uniform_0D(self, A(i,j,k))
                end do
            end do
        end do
    end subroutine

    ! Returns one normal distributed random number
    ! Implements the Box-Müller algorithm in a somewhat threadsafe fashion.
    ! However, the same rng state will be used by all threads, which makes this slow.
    ! For reproducable results, don't use OMP when calling the RNG.
    ! If performance is critical, use a seperate rng instance for each thread.
    subroutine random_normal_0D(self, rnor)
        class(RandomNumberGenerator), intent(inout) :: self
        real(dp), parameter :: PI = 4._dp * datan(1._dp)
        real(dp), intent(out) :: rnor
        real(dp) :: u(2)
        real(dp), save :: cache
        logical, save :: has_cache = .false.
        ! each thread has its own cache
        !$omp threadprivate(cache, has_cache)

        if (has_cache) then
            has_cache = .false.
            rnor = cache
        else
            call self%random_uniform_1D(u)
            rnor  = sqrt(-2._dp * log(u(1))) * sin(2._dp * PI * u(2))
            cache = sqrt(-2._dp * log(u(1))) * cos(2._dp * PI * u(2))
            has_cache = .true.
        end if
    end subroutine random_normal_0D

    subroutine random_normal_1D(self, A)
        class(RandomNumberGenerator), intent(inout) :: self
        real(dp), intent(out) :: A(:)
        integer :: i
        integer :: n

        n = size(A)
        do i=1,n
            call random_normal_0D(self, A(i))
        end do
    end subroutine

    subroutine random_normal_2D(self, A)
        class(RandomNumberGenerator), intent(inout) :: self
        real(dp), intent(out) :: A(:,:)
        integer :: i, j
        integer :: n(2)

        n = shape(A)
        do j=1,n(2)
            do i=1,n(1)
                call random_normal_0D(self, A(i,j))
            end do
        end do
    end subroutine

    subroutine random_normal_3D(self, A)
        class(RandomNumberGenerator), intent(inout) :: self
        real(dp), intent(out) :: A(:,:,:)
        integer :: i, j, k
        integer :: n(3)

        n = shape(A)
        do k=1,n(3)
            do j=1,n(2)
                do i=1,n(1)
                    call random_normal_0D(self, A(i,j,k))
                end do
            end do
        end do
    end subroutine

    ! Implementation of xorshift64*
    ! Took implementation of the algorithm in the Spire package as reference (MIT license)
    ! https://github.com/typelevel/spire/blob/master/extras/src/main/scala/spire/random/rng/XorShift64Star.scala
    ! has a period of 2**64-1
    function xorshift64star(state) result(rnd)
        implicit none
        real(dp) :: rnd
        integer(di), intent(inout) :: state
        integer(di) :: b(4)
        integer(di), parameter :: f(4) = [56605, 20332, 62609, 9541] ! 2685821657736338717 in base 2**16
        integer(di) :: m(4)

        state = ieor(state, ishft(state, -12))
        state = ieor(state, ishft(state,  25))
        state = ieor(state, ishft(state, -27))

        ! fortran has no unsigned ints.
        ! below would be how we convert signed to unsigned (twos complement)
        ! if (x < 0) then
        !   x = 2_16**64 + state
        ! else
        !   x = state
        ! end if
        ! doing the modulo is equivalent however
        ! BUT: no 128 bit int support from intel compiler
        ! we therefore use some trickery by transforming the number into base 2**16
        ! the version below would work with gnu compiler
        ! x = modulo(state * 2685821657736338717_16, 2_16**64)
        ! rnd = real(x, 8) / 2_16**64
        ! write(*,*) rnd
        ! this is the version for intel
        call toBase16(state, b)
        call multBase16(b, f, m)
        ! call fromBase16(m, xx)
        ! to return a double we divide by 2^64
        rnd = m(4) / 2.d0**16 + m(3) / 2.d0**32 + m(2) / 2.d0**48 + m(1) / 2.d0**64
        ! write(*,*) rnd

    end function

    ! becaus we cannot use 128 bit ints in ifort we represent our number using 4 ints in base 2**16
    ! x = sum_i b(i) * 2**(16*(i-1))
    ! treats x as unsigned int
    subroutine toBase16(x, b)
        implicit none
        integer(di), intent(in) :: x
        integer(di), intent(out) :: b(4)
        integer :: i
        integer(di) :: t

        b(:) = 0
        if (x<0) then
            t = (9223372036854775807_8 + x) + 1
            b(4) = 2_8**15
        else
          t = x
        end if

        do i=1,4
            b(i) = b(i) + modulo(t, 2_8**16)
            t = t / 2_8**16
        end do

    end subroutine

    ! does multiplication mod 2**64 of two numbers in the base of 2**16
    subroutine multBase16(a, b, c)
        implicit none
        integer(di), intent(in) :: a(4), b(4)
        integer(di), intent(out) :: c(4)

        integer :: i,j

        c(:) = 0
        do i=1,4
            do j=1,5-i
                c(i+j-1) = c(i+j-1) + a(i) * b(j)
            end do
        end do
        do i=1,3
            c(i+1) = c(i+1) + c(i) / 2_8**16
            c(i) = modulo(c(i), 2_8**16)
        end do
        c(4) = modulo(c(4), 2_8**16)

    end subroutine

end module random
