module projectsubs

    contains

    subroutine init_random_seed()
        implicit none
        integer                          :: i, n, clock
        integer,allocatable,dimension(:) :: seed

        call random_seed(size = n)
        allocate(seed(n))

        call system_clock(count=clock)

        seed = clock + 37 * (/ (i-1, i = 1, n) /)
        call random_seed(put = seed)

        deallocate(seed)
    end subroutine

    subroutine initializepolymer(r,r_e,N)
        implicit none
        real(8),allocatable,dimension(:,:)   :: r
        real(8)                              :: randnum, r_e
        integer                              :: i, j, N

     10 r(1,:) = (/ 0, 0, 0 /)
        call random_number(randnum)
        if(randnum > 0 .and. randnum < (1./6)) then
            r(2,:) = (/ 1, 0, 0 /)
        else if(randnum > (1./6) .and. randnum < (2./6)) then
            r(2,:) = (/ -1, 0, 0 /)
        else if(randnum > (2./6) .and. randnum < (3./6)) then
            r(2,:) = (/ 0, 1, 0 /)
        else if(randnum > (3./6) .and. randnum < (4./6)) then
            r(2,:) = (/ 0, -1, 0 /)
        else if(randnum > (4./6) .and. randnum < (5./6)) then
            r(2,:) = (/ 0, 0, 1 /)
        else
            r(2,:) = (/ 0, 0, -1 /)
        end if

        do i = 3, N
            20 call random_number(randnum)
            if(randnum > 0 .and. randnum < (1./6)) then
                r(i,:) = r(i-1,:) + (/ 1, 0, 0 /)
            else if(randnum > (1./6) .and. randnum < (2./6)) then
                r(i,:) = r(i-1,:) + (/ -1, 0, 0 /)
            else if(randnum > (2./6) .and. randnum < (3./6)) then
                r(i,:) = r(i-1,:) + (/ 0, 1, 0 /)
            else if(randnum > (3./6) .and. randnum < (4./6)) then
                r(i,:) = r(i-1,:) + (/ 0, -1, 0 /)
            else if(randnum > (4./6) .and. randnum < (5./6)) then
                r(i,:) = r(i-1,:) + (/ 0, 0, 1 /)
            else
                r(i,:) = r(i-1,:) + (/ 0, 0, -1 /)
            end if

            if(r(i,1) - r(i-1,1) == - r(i-1,1) + r(i-2,1) .and. r(i,1) - r(i-1,1) /= 0) then
                go to 20
            else if(r(i,2) - r(i-1,2) == - r(i-1,2) + r(i-2,2) .and. r(i,2) - r(i-1,2) /= 0) then
                go to 20
            else if(r(i,3) - r(i-1,3) == - r(i-1,3) + r(i-2,3) .and. r(i,3) - r(i-1,3) /= 0) then
                go to 20
            end if

            do j=1, i
                if(r(i,1) == r(j,1) .and. r(i,2) == r(j,2) .and. r(i,3) == r(j,3) .and. i /= j) then
                    go to 10
                end if
            end do
        end do

        do j=1,N
            r(j,:) = r(j,:) * r_e
        end do


    end subroutine

!    subroutine dimrecursive(r,r_e,N)
!        implicit none
!        real(8),allocatable,dimension(:,:)   :: r
!        real(8)                              :: r_e
!        integer                              :: N, N1, N2
!
!        if(N <= N0) then
!            call initializepolymernr(r,r_e,N)
!        else
!            N1 = int(N/2)
!            N2 = N - N1
!            r(1,:) = dimrecursive(r,r_e,N)
!            r(2,:) = dimrecursive(r,r_e,N)
!
!    end subroutine


    subroutine ljpotential(U_lj,U_lj_i,U_lj_i_temp,d,d_temp,j,k,r,sigma,epsil,N)
        implicit none

        real(8),allocatable,dimension(:,:)  :: U_lj, U_lj_i, U_lj_i_temp, d, d_temp, r
        real(8)                             :: sigma, epsil
        integer                             :: j, k, l, m, N

        if(j == 0) then
            U_lj(0,1) = 0
            do l = 1, N-1
                do m = l+1, N
                    call distance(l,m,d,r,U_lj_i,sigma,epsil)
                    U_lj(0,1) = U_lj(0,1) + U_lj_i(l,m)
                end do
            end do
        else
            if(j==1 .and. k==1) then
                U_lj(j,k) = U_lj(0,1)
            else if(k==1) then
                U_lj(j,k) = U_lj(j-1,N)
            else
                U_lj(j,k) = U_lj(j,k-1)
            end if

            U_lj_i_temp = U_lj_i
            d_temp = d

            do m = k+1, N
                call distance(k,m,d,r,U_lj_i,sigma,epsil)
                U_lj(j,k) = U_lj(j,k) + U_lj_i(k,m) - U_lj_i_temp(k,m)
            end do

        end if
    end subroutine

    subroutine morsepotential(U_m,d,D_morse,alpha,N,j,k,r_e)
        implicit none
        real(8),allocatable,dimension(:,:) :: U_m, d
        real(8)                            :: D_morse, alpha, r_e
        integer                            :: i, j, k, N

        do i = 1, N-1
            U_m(j,k) = U_m(j,k) + D_morse * (1 - exp(-alpha*(d(i,i+1)-r_e))) ** 2
        end do

    end subroutine

    subroutine bendingpotential(U_b,d,theta,theta_0,k_theta,N,j,k)
        implicit none
        real(8),allocatable,dimension(:,:) :: U_b, d
        real(8),allocatable,dimension(:)   :: theta
        real(8),parameter                  :: pi = 4.0 * atan(1.0)
        real(8)                            :: theta_0, k_theta
        integer                            :: i, j, k, N

        do i = 1, N-2
            theta(i) = acos((d(i,i+1)**2 + d(i+1,i+2)**2 - d(i,i+2)**2) / (2*d(i,i+1)*d(i+1,i+2)))
            U_b(j,k) = U_b(j,k) + 0.5 * k_theta * (cos(theta(i)) - cos(theta_0)) ** 2.0
        end do

    end subroutine

    subroutine distance(l,m,d,r,U_lj_i,sigma,epsil)
        implicit none
        real(8),allocatable,dimension(:,:) :: U_lj_i, r, d
        real(8)                            :: sigma, epsil, dx, dy, dz
        integer                            :: l, m

        dx = r(l,1) - r(m,1)
        dy = r(l,2) - r(m,2)
        dz = r(l,3) - r(m,3)
        if(dx == 0 .and. dy == 0 .and. dz == 0) then
            print*, 'distance equals zero'
            print*, r(l,:), r(m,:), l, m
            stop
        end if

        d(l,m) = sqrt(dx**2 + dy**2 + dz**2)

        U_lj_i(l,m) = 2 * epsil  * (((sigma/d(l,m))**9)) - 1.5*(((sigma/d(l,m))**6))


    end subroutine

end module
