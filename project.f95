program project
    use projectsubs
    implicit none
    real(8),allocatable,dimension(:,:)   :: r, U, U_lj, U_m, U_b, U_lj_i, U_lj_i_temp, d, d_temp
    real(8),allocatable,dimension(:)     :: theta
    real(8)                              :: kT_init, kT_final, kT_delta, sigma, epsil, kT, r_e, D_morse, alpha, theta_0, k_theta
    integer                              :: N, i, j, k, nkT, nMC

    call init_random_seed()

    N = 100
    kT_init  = 1.0
    kT_final = 2.0
    kT_delta = 1.0
    nMC = 1
    nkT = int((kT_final - kT_init) / kT_delta) + 1
    kT = kT_init
    sigma = 7.559
    epsil = 0.00019
    r_e = 2.891
    D_morse = 0.1324
    alpha = 1.018
    theta_0 = 113.3
    k_theta = 0.04949


    allocate(r(N,3))
    allocate(U(0:nMC,N))
    allocate(U_lj_i(N,N))
    allocate(d(N,N))
    allocate(U_lj_i_temp(N,N))
    allocate(d_temp(N,N))
    allocate(U_lj(0:nMC,N))
    allocate(U_m(0:nMC,N))
    allocate(U_b(0:nMC,N))
    allocate(theta(N))

    do i = 1, nkT
        call initializepolymer(r,r_e,N)
        j = 0
        k = 1
        call ljpotential(U_lj,U_lj_i,U_lj_i_temp,d,d_temp,j,k,r,sigma,epsil,N)
        call morsepotential(U_m,d,D_morse,alpha,N,j,k,r_e)
        call bendingpotential(U_b,d,theta,theta_0,k_theta,N,j,k)
        do j = 1, nMC
            do k = 1, N
!                call moveparticle
                call ljpotential(U_lj,U_lj_i,U_lj_i_temp,d,d_temp,j,k,r,sigma,epsil,N)
                call morsepotential(U_m,d,D_morse,alpha,N,j,k,r_e)
                call bendingpotential(U_b,d,theta,theta_0,k_theta,N,j,k)
!                call keepmove
!                if(j >= nMC/2 .and. mod(j,50) == 0) call gyration()
!                if(j >= nMC/2 .and. mod(j,50) == 0) call meandistance()
            end do
        end do
    end do

    do i=1,N
        print*, r(i,:)
    end do

        print*, U_lj(0,1), U_m(0,1), U_b(0,1)
    do i=1,N
        print*, U_b(1,i), U_m(1,i)
    end do
end program
