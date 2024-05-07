!!  Copyright (C)  Stichting Deltares, 2012-2024.
!!
!!  This program is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License version 3,
!!  as published by the Free Software Foundation.
!!
!!  This program is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program. If not, see <http://www.gnu.org/licenses/>.
!!
!!  contact: delft3d.support@deltares.nl
!!  Stichting Deltares
!!  P.O. Box 177
!!  2600 MH Delft, The Netherlands
!!
!!  All indications and logos of, and references to registered trademarks
!!  of Stichting Deltares remain the property of Stichting Deltares. All
!!  rights reserved.
module test_lateral
   use ftnunit
   use stdlib_kinds, only: dp
   use dfm_error, only: DFM_NOERR, DFM_GENERICERROR
   use m_lateral

   implicit none

   real(dp), parameter :: tolerance = 1.0e-10_dp

   contains
!
!
!==============================================================================
!   
subroutine tests_lateral
   ! initialization of global state variables for all tests in this module
   call setup_testcase()

   call test( test_get_lateral_discharge,     'Test computation of total discharge over laterals.')
   call test( test_add_lateral_load_and_sink, 'Test computation of constituents sinks and sources due to laterals.')

   ! deallocation of global state variables
   call finish_testcase()

end subroutine tests_lateral
!
!==============================================================================
!> Test computation of sinks and sources (discharge and transport load per cell) due to laterals
subroutine test_get_lateral_discharge
   use m_flow, only: vol1
   use m_transportdata, only: numconst
   use m_flowgeom, only: ndxi
   
   real(kind=dp), allocatable, dimension(:,:)   :: lateral_discharge_in                !< Lateral discharge going into the model (source)
   real(kind=dp), allocatable, dimension(:,:)   :: lateral_discharge_out               !< Lateral discharge extracted out of the model (sink)
   real(kind=dp), allocatable, dimension(:,:)   :: reference_lateral_discharge_in      !< Reference lateral discharge going into the model (source)
   real(kind=dp), allocatable, dimension(:,:)   :: reference_lateral_discharge_out     !< Reference lateral discharge extracted out of the model (sink)
   real(kind=dp), allocatable, dimension(:,:)   :: transport_load                      !< Load being transported into domain
   real(kind=dp), allocatable, dimension(:,:)   :: transport_sink                      !< Load being transported out 
   
   integer :: ierr                            ! error flag
   integer :: i_cell, i_const, i_lateral      ! loop counters

   ierr = 0
   allocate(lateral_discharge_in(numlatsg,ndxi),stat=ierr)
   allocate(lateral_discharge_out(numlatsg,ndxi),stat=ierr)
   allocate(reference_lateral_discharge_in(numlatsg,ndxi),stat=ierr)
   allocate(reference_lateral_discharge_out(numlatsg,ndxi),stat=ierr)
      
   reference_lateral_discharge_in = 0._dp
   reference_lateral_discharge_in(1,nnlat(1)) = 3._dp
   reference_lateral_discharge_in(1,nnlat(2)) = 3._dp
   reference_lateral_discharge_in(1,nnlat(3)) = 3._dp
   reference_lateral_discharge_out = 0._dp
   reference_lateral_discharge_out(2,nnlat(4)) = 5._dp
   reference_lateral_discharge_out(2,nnlat(5)) = 5._dp
   
   call get_lateral_discharge(lateral_discharge_in,lateral_discharge_out)
   
   do i_lateral = 1,numlatsg
      do i_cell=1,ndxi
         call assert_comparable(lateral_discharge_in(i_lateral,i_cell), reference_lateral_discharge_in(i_lateral,i_cell), tolerance, "get_lateral_discharge(): lateral_discharge_in is not correct" )
         call assert_comparable(lateral_discharge_out(i_lateral,i_cell), reference_lateral_discharge_out(i_lateral,i_cell), tolerance, "get_lateral_discharge(): lateral_discharge_out is not correct" )         
      enddo
   enddo
      
   deallocate(lateral_discharge_in)
   deallocate(lateral_discharge_out)
   deallocate(reference_lateral_discharge_in)
   deallocate(reference_lateral_discharge_out)
 
end subroutine test_get_lateral_discharge
!
!==============================================================================
!> Test computation of sinks and sources (discharge and transport load per cell) due to laterals
subroutine test_add_lateral_load_and_sink
   use m_flow, only: vol1
   use m_transportdata, only: numconst
   use m_flowgeom, only: ndxi
   
   real(kind=dp), allocatable, dimension(:,:)   :: discharge_in                !< Lateral discharge going into the model (source)
   real(kind=dp), allocatable, dimension(:,:)   :: discharge_out               !< Lateral discharge extracted out of the model (sink)
   real(kind=dp), allocatable, dimension(:,:)   :: transport_load              !< Load being transported into domain
   real(kind=dp), allocatable, dimension(:,:)   :: transport_sink              !< sink term due to transport into domain

   real(kind=dp), allocatable, dimension(:,:)   :: ref_load
   real(kind=dp) :: refval
   real(kind=dp) :: dvoli 
   integer :: ierr                            ! error flag
   integer :: i_cell, i_const, i_lateral      ! loop counters

   ierr = 0
   allocate(discharge_in(numlatsg,ndxi),stat=ierr)
   allocate(discharge_out(numlatsg,ndxi),stat=ierr)
   allocate(transport_load(numconst,ndxi),stat=ierr)
   allocate(transport_sink(numconst,ndxi),stat=ierr)
   allocate(ref_load(numconst,ndxi),stat=ierr)

   ! initialize transport to zero
   transport_load(:,:) = 0._dp
   transport_sink(:,:) = 0._dp

   ! first check that no discharge means no added transport
   discharge_in = 0._dp
   discharge_out = 0._dp
   call add_lateral_load_and_sink(transport_load,transport_sink,discharge_in,discharge_out,vol1,tolerance)

   call assert_comparable(sum(transport_load), 0._dp, tolerance, "todo")
   call assert_comparable(sum(transport_sink), 0._dp, tolerance, "todo")

   ! check transport into the domain
   i_lateral = 1 ! only the first lateral is a source
   discharge_in(i_lateral,nnlat(1)) = 5._dp
   discharge_in(i_lateral,nnlat(2)) = 5._dp
   discharge_in(i_lateral,nnlat(3)) = 5._dp
   call add_lateral_load_and_sink(transport_load,transport_sink,discharge_in,discharge_out,vol1,tolerance)

   do i_const = 1,numconst
      do i_cell=1,ndxi
         dvoli = 1/(vol1(i_cell))
         refval = dvoli*incoming_lat_concentration(1,i_const,i_lateral)*discharge_in(i_lateral,i_cell)
         call assert_comparable(transport_load(i_const,i_cell),refval,tolerance,"lateral_load value is not correct" )
      enddo
   enddo
   call assert_comparable(sum(transport_sink), 0._dp, tolerance, "todo")

   ! check transport out of the domain
   i_lateral = 2 ! only the second lateral is a sink
   discharge_in(:,:) = 0._dp
   discharge_out(i_lateral,nnlat(4)) = -5._dp
   discharge_out(i_lateral,nnlat(5)) = -5._dp
   ! copy values of transport_load
   ref_load(:,:) = transport_load(:,:)
   call add_lateral_load_and_sink(transport_load,transport_sink,discharge_in,discharge_out,vol1,tolerance)
   ! check that transport_load was not changed
   call assert_comparable(sum(transport_load), sum(ref_load), tolerance, "todo")
   do i_const = 1,numconst
      do i_cell=1,ndxi
         dvoli = 1/(vol1(i_cell))
         refval = dvoli*discharge_out(i_lateral,i_cell)
         call assert_comparable(transport_sink(i_const,i_cell), refval, tolerance, "lateral_sink value is not correct" )
      enddo
   enddo
      
   deallocate(discharge_in)
   deallocate(discharge_out)
   deallocate(transport_load)
   deallocate(transport_sink)
 
end subroutine test_add_lateral_load_and_sink
!
!> initialize a domain with an incoming and an outgoing lateral for three constituents
subroutine setup_testcase()
   use m_flow, only: vol1, hs
   use m_flowtimes, only: dts
   use m_partitioninfo, only: jampi
   use m_transportdata, only: numconst
   use m_cell_geometry, only: ba
   use m_flowgeom, only: ndxi

   integer :: ierr                  ! error flag
   integer :: i_cell, i_lateral, k1 ! loop counter
   integer :: k                     ! node number            

   jampi = 0
   dts = 1.0e-3_dp
   ! domain of 10 points, 2 laterals (1 incoming with 3 nodes, 1 outgoing with 2 nodes)
   ! and 3 constituents. Area (ba) and volume (vol1) of each cell are set to 1d0. 
   ! consider 3 constituents to represent salt, temperature and tracer transport. 
   ndxi = 10
   nlatnd = 5
   numlatsg = 2
   numconst = 3

   call initialize_lateraldata(numconst, ierr)
   allocate(n1latsg(numlatsg),stat=ierr)
   allocate(n2latsg(numlatsg),stat=ierr)
   allocate(nnlat(nlatnd),stat=ierr)
   allocate(qplat(numlatsg),stat=ierr)
   allocate(ba(ndxi),stat=ierr)
   allocate(balat(numlatsg),stat=ierr)
   allocate(vol1(ndxi),stat=ierr)
   allocate(hs(ndxi),stat=ierr)

   n1latsg(1) = 1
   n2latsg(1) = 3
   n1latsg(2) = 4
   n2latsg(2) = 5
   nnlat = (/1,2,3,5,8/)
   ba(:) = 0.1_dp
   vol1(:) = 0.1_dp
   hs(:) = 2_dp
   ! compute balat from ba
   do i_lateral = 1,numlatsg
      balat(i_lateral) = 0_dp
      do k1=n1latsg(i_lateral),n2latsg(i_lateral)
         k = nnlat(k1)
         balat(i_lateral) = balat(i_lateral) + ba(k)
      enddo
   enddo
   ! positive qplat is considered inflow (source), negative value outflow (sink) 
   qplat = (/9_dp,-10_dp/)
   ! top layer, per constituent, lateral 1, 
   incoming_lat_concentration(1,:,1) = (/31.0_dp,20.0_dp,0.23_dp/)  
   ! top layer, all constituents, lateral 2. 
   outgoing_lat_concentration(1,:,2) = 25_dp

end subroutine setup_testcase
!
!> reset to default values and deallocate arrays from other modules
subroutine finish_testcase()
   use m_flow, only: vol1
   use m_partitioninfo, only: jampi
   use m_transportdata, only: numconst
   use m_cell_geometry, only: ba
   use m_flowgeom, only: ndxi

   jampi = 1
   ndxi = 0
   numconst = 0

   call reset_lateral()
   call dealloc_lateraldata()
   deallocate(n1latsg)
   deallocate(n2latsg)
   deallocate(nnlat)
   deallocate(qplat)
   deallocate(ba)
   deallocate(balat)
   deallocate(vol1)

end subroutine finish_testcase

end module test_lateral
