program project
    use projectsubs
    implicit none
    real(8),allocatable,dimension(:,:,:) :: omega
    real(8),allocatable,dimension(:,:)   :: r
    integer                              :: N, i, k, stride

    call init_random_seed()

    stride = 10
    N = 500
    k = N / stride

    allocate(r(N,3))
    allocate(omega(N,stride,3))




    call initializepolymer(r,omega,k,stride)


    do i = 1, N
        print*, i, r(i,:)
    end do

end program
