program combine_h5
  use hdf5_interface
  use constants
  use utils

  implicit none

  type(hdf5_file) :: hp
  integer, parameter :: nargs = 3, ndata = 3
  character(len=MAX_STRING_LEN), dimension(ndata) :: data_names = ['vp ', 'vs ', 'rho']
  character(len=MAX_STRING_LEN), dimension(nargs) :: args
  character(len=MAX_STRING_LEN) :: indir, outdir, model_name
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: vp, vs, rho
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: x, y, z
  integer :: i

  do i = 1, nargs
    call get_command_argument(i, args(i))
  end do

  model_name = args(1)
  indir = args(2)
  outdir = args(3)

  do i = 1, ndata
    call h5read(trim(indir)//'/'//trim(data_names(i))//'_'//trim(model_name)//'.h5', '/'//trim(data_names(i)), vp)
    vp = transpose_3(vp)
    call h5read(trim(indir)//'/'//trim(data_names(i))//'_'//trim(model_name)//'.h5', '/'//trim(data_names(i)), vs)
    vs = transpose_3(vs)
    call h5read(trim(indir)//'/'//trim(data_names(i))//'_'//trim(model_name)//'.h5', '/'//trim(data_names(i)), rho)
    rho = transpose_3(rho)
  enddo
  call h5read(trim(indir)//'/'//trim(data_names(1))//'_'//trim(model_name)//'.h5', '/x', x)
  call h5read(trim(indir)//'/'//trim(data_names(1))//'_'//trim(model_name)//'.h5', '/y', y)
  call h5read(trim(indir)//'/'//trim(data_names(1))//'_'//trim(model_name)//'.h5', '/z', z)

  call hp%open(trim(outdir), status='new', action='write')
  call hp%add('/x', x)
  call hp%add('/y', y)
  call hp%add('/z', z)
  call hp%add('/vp', vp)
  call hp%add('/vs', vs)
  call hp%add('/rho', rho)
  call hp%close()


end program combine_h5