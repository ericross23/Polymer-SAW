program project
    use projectsubs
    implicit none
    real(8),allocatable,dimension(:,:,:) :: omega
    real(8),allocatable,dimension(:,:)   :: r, U, U_i, U_i_temp, d, d_temp
    real(8)                              :: kT_init, kT_final, kT_delta, sigma, epsil, kT, r_e
    integer                              :: N, i, j, k_stride, k, stride, nkT, nMC

    call init_random_seed()

    stride = 25
    N = 100
    kT_init  = 1.0
    kT_final = 2.0
    kT_delta = 1.0
    nMC = 1
    k_stride = N / stride
    nkT = int((kT_final - kT_init) / kT_delta) + 1
    kT = kT_init
    sigma = 1
    epsil = 1
    r_e = 2.891

    allocate(r(N,3))
    allocate(omega(N,stride,3))
    allocate(U(0:nMC,N))
    allocate(U_i(N,N))
    allocate(d(N,N))
    allocate(U_i_temp(N,N))
    allocate(d_temp(N,N))

    do i = 1, nkT
        call initializepolymer(r,omega,k_stride,stride,r_e)
        j = 0
        print*,'polymer initialized'
        call potential(U,U_i,U_i_temp,d,d_temp,j,k,r,sigma,epsil,N)
        print*, 'initial potential'
        do j = 1, nMC
            do k = 1, N
!                call moveparticle
                call potential(U,U_i,U_i_temp,d,d_temp,j,k,r,sigma,epsil,N)
!                call keepmove
!                if(j >= nMC/2 .and. mod(j,50) == 0) call gyration()
!                if(j >= nMC/2 .and. mod(j,50) == 0) call meandistance()
            end do
        end do
    end do

print*, U(1,:)

end program
