program test_io_netcdf
   use ftnunit
   use test_coordinate_reference_system, only: tests_coordinate_reference_system
   use test_netcdf_define_mode, only: tests_netcdf_define_mode

   implicit none

   ! Setup
   call prepareTests()
   call runtests_init()

   ! Run tests for modules
   call tests_coordinate_reference_system()
   call tests_netcdf_define_mode()

   ! Teardown
   call runtests_final()
   call showResult()

contains

!> Routine to start the testing
!! Note: This routine merely takes care that the unit tests are indeed run
   subroutine prepareTests

      integer :: lun !< LU-number

      open (newunit=lun, file='ftnunit.run')
      write (lun, '(a)') 'ALL'
      close (lun)

   end subroutine prepareTests

!> Start the browser to show the result
!!
   subroutine showResult()

      call system('ftnunit.html')

   end subroutine showResult

end program test_io_netcdf
