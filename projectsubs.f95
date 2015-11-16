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

    subroutine initializepolymer(r,omega,k_stride,stride,r_e)
        implicit none
        real(8),allocatable,dimension(:,:)   :: r
        real(8),allocatable,dimension(:,:,:) :: omega
        real(8)                              :: randnum, r_e
        integer                              :: i, j, k, k_stride, l, stride, randint



        10 do l = 1, k_stride
            r(1,:) = (/ 0, 0, 0 /)
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

            omega(l,1,:) = (/ 0, 0, 0 /)
            omega(l,2,:) = r(2,:)

            do i = 3, stride
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
                    if(r(i,1) == r(j,1) .and. i /= j) then
                        if(r(i,2) == r(j,2) .and. i /= j) then
                            if(r(i,3) == r(j,3) .and. i /= j) then
                                call init_random_seed()
                                go to 10
                            end if
                        end if
                    end if
                end do
                omega(l,i,:) = r(i,:)
            end do
        end do

        call random_number(randnum)
        randint = ceiling(randnum*k_stride)
        r(1:stride,:) = omega(randint,1:stride,:)
        call random_number(randnum)
        randint = ceiling(randnum*k_stride)
        r(stride+1:2*stride,:) = omega(randint,:,:)
        do i = 1, stride
            r(stride+l,:) = r(stride+l,:) + r(stride,:)
        end do
        call random_number(randnum)
        if(randnum > 0 .and. randnum < (1./6)) then
            do i = 1, stride
                r(stride+i,:) = r(stride+i,:) + (/ 1, 0, 0 /)
            end do
        else if(randnum > (1./6) .and. randnum < (2./6)) then
            do i = 1, stride
                r(stride+i,:) = r(stride+i,:) + (/ -1, 0, 0 /)
            end do
        else if(randnum > (2./6) .and. randnum < (3./6)) then
            do i = 1, stride
                r(stride+i,:) = r(stride+i,:) + (/ 0, 1, 0 /)
            end do
        else if(randnum > (3./6) .and. randnum < (4./6)) then
            do i = 1, stride
                r(stride+i,:) = r(stride+i,:) + (/ 0, -1, 0 /)
            end do
        else if(randnum > (4./6) .and. randnum < (5./6)) then
            do i = 1, stride
                r(stride+i,:) = r(stride+i,:) + (/ 0, 0, 1 /)
            end do
        else
            do i = 1, stride
                r(stride+i,:) = r(stride+i,:) + (/ 0, 0, -1 /)
            end do
        end if


        do i = 3, k_stride
            call random_number(randnum)
            randint = ceiling(randnum*k_stride)
            r(((i-1)*stride)+1:i*stride,:) = omega(randint,:,:)
            do j = 1, stride
                r(((i-1)*stride)+j,:) = r(((i-1)*stride)+j,:) + r((i-1)*stride,:)
            end do
            call random_number(randnum)
            if(randnum > 0 .and. randnum < (1./6)) then
                do l = 1, stride
                    r(((i-1)*stride)+l,:) = r(((i-1)*stride)+l,:) + (/ 1, 0, 0 /)
                end do
            else if(randnum > (1./6) .and. randnum < (2./6)) then
                do l = 1, stride
                    r(((i-1)*stride)+l,:) = r(((i-1)*stride)+l,:) + (/ -1, 0, 0 /)
                end do
            else if(randnum > (2./6) .and. randnum < (3./6)) then
                do l = 1, stride
                    r(((i-1)*stride)+l,:) = r(((i-1)*stride)+l,:) + (/ 0, 1, 0 /)
                end do
            else if(randnum > (3./6) .and. randnum < (4./6)) then
                do l = 1, stride
                    r(((i-1)*stride)+l,:) = r(((i-1)*stride)+l,:) + (/ 0, -1, 0 /)
                end do
            else if(randnum > (4./6) .and. randnum < (5./6)) then
                do l = 1, stride
                    r(((i-1)*stride)+l,:) = r(((i-1)*stride)+l,:) + (/ 0, 0, 1 /)
                end do
            else
                do l = 1, stride
                    r(((i-1)*stride)+l,:) = r(((i-1)*stride)+l,:) + (/ 0, 0, -1 /)
                end do
            end if

            do j=1, i*stride
                do k=j+1, i*stride
                    if(r(k,1) == r(j,1)) then
                        if(r(k,2) == r(j,2)) then
                            if(r(k,3) == r(j,3)) then
                                call init_random_seed()
                                go to 10
                            end if
                        end if
                    end if
                end do
            end do
        end do

        do i=1, k_stride * stride
            r(i,:) = r(i,:) * r_e
        end do

    end subroutine

    subroutine potential(U,U_i,U_i_temp,d,d_temp,j,k,r,sigma,epsil,N)
        implicit none
        real(8),allocatable,dimension(:,:) :: U, U_i, U_i_temp, d, d_temp, r
        real(8)                            :: sigma, epsil
        integer                            :: N, j, k

        call ljpotential(U,U_i,U_i_temp,d,d_temp,j,k,r,sigma,epsil,N)
!        call morsepotential()
!        call bendingpotential()

    end subroutine

    subroutine ljpotential(U,U_i,U_i_temp,d,d_temp,j,k,r,sigma,epsil,N)
        implicit none

        real(8),allocatable,dimension(:,:)  :: U, U_i, U_i_temp, d, d_temp, r
        real(8)                             :: sigma, epsil
        integer                             :: j, k, l, m, N

        if(j == 0) then
            U(0,1) = 0
            do l = 1, N-1
                do m = l+1, N
                    call distance(l,m,d,r,U_i,sigma,epsil)
                    U(0,1) = U(0,1) + U_i(l,m)
!                    print*, U(0,1),U_i(l,m), d(l,m)
                end do
            end do
        else
            if(j==1 .and. k==1) then
                U(j,k) = U(0,1)
            else if(k==1) then
                U(j,k) = U(j-1,N)
            else
                U(j,k) = U(j,k-1)
            end if

            U_i_temp = U_i
            d_temp = d

            do m = k+1, N
                call distance(k,m,d,r,U_i,sigma,epsil)
                U(j,k) = U(j,k) + U_i(k,m) - U_i_temp(k,m)
            end do

        end if
    end subroutine

!    subroutine morsepotential(N)
!        implicit none
!        integer :: i, N
!
!        do i = 1, N-1
!            U_morse = U_morse + D * (1 - exp(-alpha*(r()-r_e))) ** 2
!        end do
!
!    end subroutine

!    subroutine bendingpotential()
!        implicit none
!        integer :: i, N
!
!        do i = 1, N-2
!            U_theta = U_theta + 0.5 * k_theta * (cos(theta()) - cos(theta_0)) ** 2
!        end do
!    end subroutine

    subroutine distance(l,m,d,r,U_i,sigma,epsil)
        implicit none
        real(8),allocatable,dimension(:,:) :: U_i, r, d
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

        U_i(l,m) = 2 * epsil  * (((sigma/d(l,m))**9)) - 1.5*(((sigma/d(l,m))**6))


    end subroutine

end module
