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

    subroutine initializepolymer(r,omega,k,stride)
        implicit none
        real(8),allocatable,dimension(:,:)   :: r
        real(8),allocatable,dimension(:,:,:) :: omega
        real(8)                              :: randnum
        integer                              :: i, j, k, l, stride, randint



        do l = 1, k
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

            10 do i = 3, stride
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
                                go to 10
                            end if
                        end if
                    end if
                end do
                omega(l,i,:) = r(i,:)
            end do
        end do

        call random_number(randnum)
        randint = ceiling(randnum*k)
        r(1:stride,:) = omega(randint,1:stride,:)
        call random_number(randnum)
        randint = ceiling(randnum*k)
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


        30 do i = 3, k
            call random_number(randnum)
            randint = ceiling(randnum*k)
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

            do j=1, i
                if(r(i,1) == r(j,1) .and. i /= j) then
                    if(r(i,2) == r(j,2) .and. i /= j) then
                        if(r(i,3) == r(j,3) .and. i /= j) then
                            go to 30
                        end if
                    end if
                end if
            end do

    end subroutine
end module

