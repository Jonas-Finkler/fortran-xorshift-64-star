
program main
    use precision
    use random

    integer, parameter :: n = 2
    integer, parameter :: m = 10000

    type(RandomNumberGenerator) :: rng
    real(dp) :: x0
    real(dp) :: x1(n)
    real(dp) :: x2(n, n)
    real(dp) :: x3(n, n, n)
    real(dp) :: t(m)
    real(dp) :: x
    integer :: count, i, j

    ! initialize rng with seed
    rng = RandomNumberGenerator(123456_di)

    ! random uniform takes arrays with up to 3 dimensions
    call rng%random_uniform(x0)
    print*, x0
    call rng%random_uniform(x1)
    print*, x1
    call rng%random_uniform(x2)
    print*, x2
    call rng%random_uniform(x3)
    print*, x3

    ! The same is true for random_normal
    call rng%random_normal(x0)
    print*, x0
    call rng%random_normal(x1)
    print*, x1
    call rng%random_normal(x2)
    print*, x2
    call rng%random_normal(x3)
    print*, x3

    call rng%random_normal(t)
    print*, 'mean', sum(t) / m
    print*, 'std-dev', sqrt(sum((t-sum(t)/m)**2) / m)


    ! Now a simple test of the random_normal routine.
    ! Plot the numerical and analytical CDFs using:
    ! $./a.out | grep 'X' > cdf.txt
    ! $gnuplot
    ! $plot 'cdf.txt' u 2:3, 'cdf.txt' u 2:4 w lines
    do i=-300, 300
        x = real(i,dp) / 100._dp
        c = 0
        do j=1,m
            if (t(j) < x) then
                c = c + 1
            end if
        end do
        print*, 'X', x, real(c, dp) / m, 0.5_dp * (1 + erf(x / sqrt(2._dp)))
    end do


end program main