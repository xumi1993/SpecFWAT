!=====================================================================
!
!               Full Waveform Adjoint Tomography -v1.1
!               ---------------------------------------
!
!     Main historical authors: Kai Wang
!                              Macquarie Uni, Australia
!                            & University of Toronto, Canada
!                           (c) Martch 2020
!
!=====================================================================
!
subroutine run_precond()

  use fullwave_adjoint_tomo_par
  use fwat_input
  use my_mpi
  use specfem_par
  use specfem_interface
   
  implicit none

  real(kind=CUSTOM_REAL)                                                   :: zl
  integer                                                                  :: i,j,k,ispec, iglob

  if (PRECOND_TYPE=='z_precond')  then
    if(myrank==0)  write(OUT_FWAT_LOG,*) 'Use z_precond as hess ...'
    do ispec = 1, nspec_AB

      do k=1,NGLLZ; do j=1,NGLLY; do i=1,NGLLX

        iglob = ibool(i,j,k,ispec)
        if (zstore(iglob)>=0.) then
           zl =1.e-8 ! add a small vale to avoid zeros
        else
           zl = zstore(iglob) + 1.e-8 ! 
        endif 
        hess_kl(i,j,k,ispec)= 1./sqrt(zl*zl)/2. ! save_adjoint_kernel hess*2

      enddo; enddo; enddo ! NGLLZ,NGLLY,NGLLX
    enddo
  
  elseif (PRECOND_TYPE=='z_sq_precond') then
    if(myrank==0)  write(OUT_FWAT_LOG,*) 'Use z_sq_precond as hess ...'
    do ispec = 1, nspec_AB
      do k=1,NGLLZ; do j=1,NGLLY; do i=1,NGLLX

        iglob = ibool(i,j,k,ispec)
        if (zstore(iglob)>=0.) then
           zl =1.e-8 ! add a small vale to avoid zeros
        else
           zl = zstore(iglob) + 1.e-8 ! 
        endif 
        hess_kl(i,j,k,ispec)= 1./sqrt(sqrt(zl*zl))/2. ! save_adjoint_kernel hess*2

      enddo; enddo; enddo ! NGLLZ,NGLLY,NGLLX
    enddo

  else

    if(myrank==0)  write(OUT_FWAT_LOG,*) 'Use default hess ...'

  endif
 
end subroutine run_precond
