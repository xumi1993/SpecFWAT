!!! This program is to do a 2D rotation with respect to a reference point (x0,y0)
 program rotation_2d
   integer, parameter :: PI=3.14159265
    character(len=256) :: fname,arg
    integer :: i,j,ier,npts
    integer :: narg
    real(kind=4) x0,y0,angle,theta
    real(kind=4) , dimension(:,:), allocatable :: data
    real(kind=4) , dimension(:,:), allocatable :: data1

    narg=command_argument_count()
    if ( narg .ne. 5 ) then
      print *,'Usage: fsismo_bin2asc filename npts x0 y0 angle'
      print *,'fname- Input file'
      print *,'npts - number of points'
      print *,'x0   - x of reference point'
      print *,'y0   - y of reference point' 
      print *,'angle- Rotation angle' 
      stop
    endif
    call get_command_argument(1,fname) 
    call get_command_argument(2,arg) 
    read(arg,*,iostat=ier) npts
    call get_command_argument(3,arg) 
    read(arg,*,iostat=ier) x0
    call get_command_argument(4,arg) 
    read(arg,*,iostat=ier) y0
    call get_command_argument(5,arg) 
    read(arg,*,iostat=ier) angle
    !write(*,*)'fname,npts,x0,y0,ang=',trim(fname),npts,x0,y0,angle
    !!! read file
    allocate(data(npts,4))
    allocate(data1(npts,4))
    open(33,file=trim(fname),status='old')
    do i=1,npts
       read(33,*,iostat=ier) data(i,1),data(i,2),data(i,3),data(i,4)
       if (ier < 0 ) exit
       !write(*,*) data(i,1),data(i,2)
    enddo
    close(33)
    !!! do rotation
    ! shift data
    data(:,1)=data(:,1)-x0 
    data(:,2)=data(:,2)-y0 
    ! roatate
    theta=angle*PI/180.
    data1(:,1)=data(:,1)*cos(theta)-data(:,2)*sin(theta) 
    data1(:,2)=data(:,1)*sin(theta)+data(:,2)*cos(theta) 
    ! shift back
    data1(:,1)=data1(:,1)+x0
    data1(:,2)=data1(:,2)+y0
    do i=1,npts
       write(*,*) data1(i,1),data1(i,2),data(i,3),data(i,4)
    enddo

 end program rotation_2d
