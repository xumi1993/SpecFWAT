module taper3d
    use utils, only: interp3
    implicit none

    type, public :: taper_cls
        integer, dimension(3) :: ndom
        real, dimension(:), allocatable :: xcrd, ycrd, zcrd
        real, dimension(:,:,:), allocatable :: taper
        contains
            procedure :: Create
            procedure :: Interp
            procedure :: Destroy
    end type taper_cls

contains
    subroutine Create(this, xmin, xmax, ymin, ymax, zmin, zmax, hx, hy, hz, &
                      taph1, taph2, tapz1, tapz2)
        implicit none
        class(taper_cls), intent(inout) :: this
        real :: xmin, xmax, ymin, ymax, zmin, zmax, hx, hy, hz, &
                taph1, taph2, tapz1, tapz2
        integer :: nx, ny, nz, i

        nx = floor((xmax-xmin)/hx+1)
        ny = floor((ymax-ymin)/hy+1)
        nz = floor((zmax-zmin)/hz+1)
        this%ndom(1) = nx
        this%ndom(2) = ny
        this%ndom(3) = nz

        allocate(this%xcrd(nx))
        allocate(this%ycrd(ny))
        allocate(this%zcrd(nz))
        allocate(this%taper(nx, ny, nz))

        do i = 1, nx
            this%xcrd(i) = xmin+(i-1)*hx
        enddo
        do i = 1, ny
            this%ycrd(i) = ymin+(i-1)*hy
        enddo
        do i = 1, nz
            this%zcrd(i) = zmin+(i-1)*hz
        enddo
        call create_taper(taph1, taph2, tapz1, tapz2, hx, hy, hz, nx, ny, nz, this%taper)
        
    end subroutine Create

    function Interp(this, x, y, z) result(val)
        class(taper_cls) :: this
        real, intent(in) :: x, y, z
        real :: val
        integer :: flag
        
        flag = 1
        if (x<this%xcrd(1) .or. x>this%xcrd(this%ndom(1))) flag = 0
        if (y<this%ycrd(1) .or. y>this%ycrd(this%ndom(2))) flag = 0
        if (z<this%zcrd(1) .or. z>this%zcrd(this%ndom(3))) flag = 0
        if (flag /= 0) then
            ! call trilin_interp(x,y,z,this%xcrd,this%ycrd,this%zcrd,&
                !  this%ndom(1),this%ndom(2),this%ndom(3),this%taper,val)
            val = interp3(this%xcrd,this%ycrd,this%zcrd,this%taper,x,y,z)
        else
            val = 0.
        endif
    end function Interp

    subroutine Destroy(this)
        class(taper_cls) :: this
        if (allocated(this%xcrd)) deallocate(this%xcrd)
        if (allocated(this%ycrd)) deallocate(this%ycrd)
        if (allocated(this%zcrd)) deallocate(this%zcrd)
        if (allocated(this%taper)) deallocate(this%taper)
    
    end subroutine Destroy
    
    subroutine create_taper(taph1, taph2, tapz1, tapz2, hx, hy, hz, nx, ny, nz, taper)
        integer :: nh1, nh2, x1, x2, x3, x4, y1, y2, y3, y4, z1, nx, ny, nz, i, j, k
        real :: hx, hy, hz, taph1, taph2, tapz1, tapz2
        real, dimension(nx) :: xtaper
        real, dimension(ny) :: ytaper
        real, dimension(nz) :: ztaper
        real, dimension(nx,ny,nz), intent(inout) :: taper

        xtaper = 0.
        ytaper = 0.
        ztaper = 0.

        nh1 = floor(taph1/hx)
        nh2 = floor(taph2/hx)
        x1 = nh1+1
        x2 = nh1 + nh2+1
        x3 = nx - nh1 - nh2
        x4 = nx - nh1
        call coswin(nh2, 1, xtaper(x1:x2-1))
        call coswin(nh2, 0, xtaper(x3+1:x4))
        xtaper(x2:x3) = 1.

        nh1 = floor(taph1/hy)
        nh2 = floor(taph2/hy)
        y1 = nh1+1
        y2 = nh1 + nh2+1
        y3 = ny - nh1 - nh2
        y4 = ny - nh1
        call coswin(nh2, 1, ytaper(y1:y2-1))
        call coswin(nh2, 0, ytaper(y3+1:y4))
        ytaper(y2:y3) = 1.

        nh1 = floor(tapz1/hz)
        nh2 = floor(tapz2/hz)
        z1 = nh1+nh2+1
        call coswin(nh2, 1, ztaper(nh1+1:z1-1))
        ztaper(z1:nz) = 1.

        do k=1,nz
            do j=1,ny
                do i=1,nx
                    taper(i,j,k) = xtaper(i) * ytaper(j) * ztaper(k)
                enddo
            enddo
        enddo
    
    end subroutine create_taper

    subroutine coswin(n, opt, winx)
        integer :: n, opt, i, nn
        real, parameter :: pi=3.1415926535
        real, dimension(:), allocatable :: win
        real, dimension(n) :: winx
    
        nn = 2 * n +1
        allocate(win(nn))
        do i = 1, nn
            win(i) = 0.5 - 0.5 *cos(2. * pi * (i-1) / (nn-1))
        enddo
        if (opt == 1) then
            winx = win(1:n)
        else
            winx = win(n+2:nn)
        endif
    end subroutine coswin

end module taper3d

