!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================


subroutine write_gradient_iso()

! file output for new model

  use tomography_kernels_iso
  implicit none
  character(len=MAX_STRING_LEN) :: m_file, fname

  ! user output
  if (myrank == 0) print *,'writing out gradient...'

  ! kernel updates
  fname = 'dbulk'
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) model_dbulk
  close(IOUT)

  fname = 'dbeta'
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) model_dbeta
  close(IOUT)

  fname = 'drho'
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) model_drho
  close(IOUT)

  if (myrank == 0) print *

end subroutine write_gradient_iso

!
!-------------------------------------------------------------------------------------------------
!

subroutine write_gradient_tiso()

! file output for new model

  use tomography_kernels_tiso
  implicit none

  character(len=MAX_STRING_LEN) :: m_file, fname

  ! user output
  if (myrank == 0) print *,'writing out gradient...'

  ! kernel updates
  if (USE_ALPHA_BETA_RHO_TISO) then
  fname = 'dalphav'
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) model_dalphav
  close(IOUT)

  fname = 'dalphah'
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) model_dalphah
  close(IOUT)
  else

  fname = 'dbulk_c'
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) model_dbulk
  close(IOUT)
  endif

  fname = 'dbetav'
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) model_dbetav
  close(IOUT)

  fname = 'dbetah'
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) model_dbetah
  close(IOUT)

  fname = 'deta'
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) model_deta
  close(IOUT)

  if (myrank == 0) print *

end subroutine write_gradient_tiso


