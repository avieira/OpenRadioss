!Copyright>        OpenRadioss
!Copyright>        Copyright (C) 1986-2024 Altair Engineering Inc.
!Copyright>
!Copyright>        This program is free software: you can redistribute it and/or modify
!Copyright>        it under the terms of the GNU Affero General Public License as published by
!Copyright>        the Free Software Foundation, either version 3 of the License, or
!Copyright>        (at your option) any later version.
!Copyright>
!Copyright>        This program is distributed in the hope that it will be useful,
!Copyright>        but WITHOUT ANY WARRANTY; without even the implied warranty of
!Copyright>        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!Copyright>        GNU Affero General Public License for more details.
!Copyright>
!Copyright>        You should have received a copy of the GNU Affero General Public License
!Copyright>        along with this program.  If not, see <https://www.gnu.org/licenses/>.
!Copyright>
!Copyright>
!Copyright>        Commercial Alternative: Altair Radioss Software
!Copyright>
!Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
!Copyright>        software under a commercial license.  Contact Altair to discuss further if the
!Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
      !||====================================================================
      !||    polygon_mod            ../common_source/tools/clipping/polygon_mod.F90
      !||--- called by ------------------------------------------------------
      !||    polygon_clipping_mod   ../common_source/tools/clipping/polygon_clipping_mod.F90
      !||====================================================================
      module polygon_mod
        implicit none
#include  "my_real.inc"

        type polygon_point_
          my_real :: y
          my_real :: z
        end type polygon_point_

        type polygon_
          type(polygon_point_), allocatable, dimension(:) :: point
          integer :: numpoint ! defined points
          integer :: size ! allocated size numpoint <= size
          my_real :: area
          my_real :: diag ! maximum dimension along y and z
        end type polygon_

        type polygon_list_
          integer num_polygons
          type(polygon_),allocatable,dimension(:) :: polygon
        end type polygon_list_

        type polygon_linked_list_
          type(polygon_), pointer :: p
          type(polygon_linked_list_), pointer :: next
        end type polygon_linked_list_

      contains

! ======================================================================================================================
!                                                   FUNCTION
! ======================================================================================================================
!! \brief add 'point' in poly data structure.
!! \details  pre-condition, allocation must be correctly sized, otherwise an error message is displayed
      !||====================================================================
      !||    polygon_addpoint           ../common_source/tools/clipping/polygon_mod.F90
      !||--- called by ------------------------------------------------------
      !||    clipping_weiler_atherton   ../common_source/tools/clipping/polygon_clipping_mod.F90
      !||--- uses       -----------------------------------------------------
      !||    constant_mod               ../common_source/modules/constant_mod.F
      !||====================================================================
        function polygon_addpoint(poly, point) result(ierr)
          use constant_mod , only : zero
          implicit none
#include  "my_real.inc"
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
          type(polygon_), intent(inout) :: poly
          type(polygon_point_), intent(in) :: point
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local Variables
! ----------------------------------------------------------------------------------------------------------------------
          integer :: numpt, isize
          integer ierr
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
          numpt = poly%numpoint
          isize = poly%size
          ierr=1
          if(numpt >= isize)then
            write(*,*) "** ERROR : unexpected situation with polygon_addpoint"
            return
          end if
          numpt = numpt+1
          poly%numpoint = numpt
          poly%point(numpt)%y = point%y;
          poly%point(numpt)%z = point%z;
          ierr = 0
          !area not recomputed for performance reason. It has to be calculated once the polygon is fully built
        end function polygon_addpoint



! ======================================================================================================================
!                                                   PROCEDURES
! ======================================================================================================================
!! \brief allocate poly with size numnode and zeroing
!! \details
      !||====================================================================
      !||    polygon_create             ../common_source/tools/clipping/polygon_mod.F90
      !||--- called by ------------------------------------------------------
      !||    clipping_weiler_atherton   ../common_source/tools/clipping/polygon_clipping_mod.F90
      !||    init_inivol_2d_polygons    ../starter/source/initial_conditions/inivol/init_inivol_2D_polygons.F90
      !||--- uses       -----------------------------------------------------
      !||    constant_mod               ../common_source/modules/constant_mod.F
      !||====================================================================
        subroutine polygon_create(poly, numnodes)
          use constant_mod , only : zero
          implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
          type(polygon_), intent(out) :: poly
          integer, intent(in) :: numnodes
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
          poly%size = numnodes
          poly%numpoint = 0
          allocate(poly%point(numnodes));
          poly%point(1:numnodes)%y = zero;
          poly%point(1:numnodes)%z = zero;
          poly%area = zero
          poly%diag = zero
        end subroutine polygon_create



! ======================================================================================================================
!                                                   PROCEDURES
! ======================================================================================================================
      !||====================================================================
      !||    polygon_zeroing   ../common_source/tools/clipping/polygon_mod.F90
      !||--- uses       -----------------------------------------------------
      !||    constant_mod      ../common_source/modules/constant_mod.F
      !||====================================================================
        subroutine polygon_zeroing(poly)
          use constant_mod , only : zero
          implicit none
#include  "my_real.inc"
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
          type(polygon_) :: poly
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local Variables
! ----------------------------------------------------------------------------------------------------------------------
          integer ii
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
          poly%numpoint = 0
          do ii=1,poly%size
            poly%point(ii)%y = zero;
            poly%point(ii)%z = zero;
            poly%diag = zero
            poly%area = zero
          end do
        end subroutine polygon_zeroing



! ======================================================================================================================
!                                                   PROCEDURES
! ======================================================================================================================
      !||====================================================================
      !||    polygon_destroy           ../common_source/tools/clipping/polygon_mod.F90
      !||--- called by ------------------------------------------------------
      !||    init_inivol_2d_polygons   ../starter/source/initial_conditions/inivol/init_inivol_2D_polygons.F90
      !||    polygon_list_destroy      ../common_source/tools/clipping/polygon_mod.F90
      !||====================================================================
        subroutine polygon_destroy(poly)
          implicit none
#include  "my_real.inc"
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
          type(polygon_), intent(out) :: poly
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
          if (allocated(poly%point))deallocate(poly%point)
        end subroutine polygon_destroy



! ======================================================================================================================
!                                                   PROCEDURES
! ======================================================================================================================
      !||====================================================================
      !||    polygon_list_destroy   ../common_source/tools/clipping/polygon_mod.F90
      !||--- calls      -----------------------------------------------------
      !||    polygon_destroy        ../common_source/tools/clipping/polygon_mod.F90
      !||--- uses       -----------------------------------------------------
      !||    constant_mod           ../common_source/modules/constant_mod.F
      !||--- called by ------------------------------------------------------
      !||    polygon_linked_list_destroy   ../common_source/tools/clipping/polygon_mod.F90
      !||    polygon_linked_list_delete_element   ../common_source/tools/clipping/polygon_mod.F90
      !||====================================================================
        subroutine polygon_list_destroy(list)
          use constant_mod , only : zero
          implicit none
#include  "my_real.inc"
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
          type(polygon_list_), intent(out) :: list
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local Variables
! ----------------------------------------------------------------------------------------------------------------------
          integer :: ii
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
          do ii=1,list%num_polygons
            call polygon_destroy(list%polygon(ii))
          end do
        end subroutine polygon_list_destroy



! ======================================================================================================================
!                                                   PROCEDURES
! ======================================================================================================================
!! \brief copy Base_polygon into Target_polygon (allocated inside this subroutine)
!! \details
      !||====================================================================
      !||    polygon_copy   ../common_source/tools/clipping/polygon_mod.F90
      !||--- called by ------------------------------------------------------
      !||    polygon_linked_list_insert_before   ../common_source/tools/clipping/polygon_mod.F90
      !||    polygon_linked_list_insert_after   ../common_source/tools/clipping/polygon_mod.F90
      !||====================================================================
        subroutine polygon_copy(Base_polygon, Target_polygon)
          implicit none
#include  "my_real.inc"
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
          type(polygon_), intent(in) :: Base_polygon
          type(polygon_), intent(out) :: Target_polygon
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local Variables
! ----------------------------------------------------------------------------------------------------------------------
          integer Base_size
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
          Base_size = Base_polygon%size
          if(Base_size > 0)then
            allocate(Target_polygon%point(Base_polygon%size))
            Target_polygon%point(1:Base_size)%y = Base_polygon%point(1:Base_size)%y
            Target_polygon%point(1:Base_size)%z = Base_polygon%point(1:Base_size)%z
            Target_polygon%numpoint = Base_polygon%numpoint
            Target_polygon%size = Base_polygon%size
            Target_polygon%area = Base_polygon%area
          else
            ! not expected
            stop 220582
          end if
        end subroutine polygon_copy


! ======================================================================================================================
!                                                   PROCEDURES
! ======================================================================================================================
      !||====================================================================
      !||    polygon_linked_list_insert_before   ../common_source/tools/clipping/polygon_mod.F90
      !||--- calls      -----------------------------------------------------
      !||    polygon_copy   ../common_source/tools/clipping/polygon_mod.F90
      !||====================================================================
        function polygon_linked_list_insert_before(list, p) result(new_elem)
          implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
          type(polygon_linked_list_), intent(in), pointer :: list
          type(polygon_), intent(in) :: p
          type(polygon_linked_list_), pointer :: new_elem

! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local Variables
! ----------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
          allocate(new_elem)
          polygon_copy(p, new_elem%p)
          new_elem%next => list
        end function polygon_linked_list_insert_before

! ======================================================================================================================
!                                                   PROCEDURES
! ======================================================================================================================
      !||====================================================================
      !||    polygon_linked_list_insert_after      ../common_source/tools/clipping/polygon_mod.F90
      !||--- called by ------------------------------------------------------
      !||    polygon_linked_list_copy              ../common_source/tools/clipping/polygon_mod.F90
      !||    Clipping_Sutherland_Hodgman           ../common_source/tools/clipping/polygon_clipping_mod.F90
      !||--- calls      -----------------------------------------------------
      !||    polygon_copy                          ../common_source/tools/clipping/polygon_mod.F90
      !||====================================================================
        subroutine polygon_linked_list_insert_after(list, p)
          implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
          type(polygon_linked_list_), intent(inout) :: list
          type(polygon_), intent(in) :: p

! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local Variables
! ----------------------------------------------------------------------------------------------------------------------
          type(polygon_), pointer :: old_next
          type(polygon_linked_list_), target :: new_elem
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
          if (.not. associated(list)) then 
            allocate(list)
            polygon_copy(p, list%p)
            nullify(list%next)
          else 
            allocate(new_elem)
            polygon_copy(p, new_elem%p)
            new_elem%next => list%next
            list%next => new_elem
          end if
        end subroutine polygon_linked_list_insert_after

! ======================================================================================================================
!                                                   PROCEDURES
! ======================================================================================================================
      !||====================================================================
      !||    polygon_linked_list_copy   ../common_source/tools/clipping/polygon_mod.F90
      !||--- calls      -----------------------------------------------------
      !||    polygon_linked_list_insert_after   ../common_source/tools/clipping/polygon_mod.F90
      !||====================================================================
        subroutine polygon_linked_list_copy(original, targ)
          implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
          type(polygon_linked_list_), intent(in) :: original
          type(polygon_linked_list_), intent(out) :: targ

! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local Variables
! ----------------------------------------------------------------------------------------------------------------------
          type(polygon_linked_list_) :: current_targ, current_original
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
          allocate(targ)
          if (associated(original)) then
            nullify(targ%next)
            current_original => original
            current_targ => targ
            do while ( associated(current_original) )
              polygon_linked_list_insert_after(current_targ, current%p)
              current_original => current_original%next
              current_targ => current_targ%next
            end do
          else
            nullify(targ)
          end if
        end subroutine polygon_linked_list_copy
    
! ======================================================================================================================
!                                                   PROCEDURES
! ======================================================================================================================
      !||====================================================================
      !||    polygon_linked_list_destroy   ../common_source/tools/clipping/polygon_mod.F90
      !||--- calls      -----------------------------------------------------
      !||    polygon_destroy   ../common_source/tools/clipping/polygon_mod.F90
      !||====================================================================
        subroutine polygon_linked_list_destroy( list )
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
          type(polygon_linked_list_), pointer  :: list
        
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local Variables
! ----------------------------------------------------------------------------------------------------------------------
          type(polygon_linked_list_), pointer  :: current
          type(polygon_linked_list_), pointer  :: next
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
          current => list
          do while ( associated(current%next) )
              next => current%next
              polygon_destroy(current%p)
              deallocate( current )
              current => next
          enddo
        end subroutine polygon_linked_list_destroy

! ======================================================================================================================
!                                                   PROCEDURES
! ======================================================================================================================
      !||====================================================================
      !||    polygon_linked_list_delete_element   ../common_source/tools/clipping/polygon_mod.F90
      !||--- calls      -----------------------------------------------------
      !||    polygon_destroy   ../common_source/tools/clipping/polygon_mod.F90
      !||====================================================================
        subroutine polygon_linked_list_delete_element( list, elem )
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
          type(polygon_linked_list_), pointer  :: list
          type(polygon_linked_list_), pointer  :: elem

! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local Variables
! ----------------------------------------------------------------------------------------------------------------------
          type(polygon_linked_list_), pointer  :: current
          type(polygon_linked_list_), pointer  :: prev

! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
          if ( associated(list,elem) ) then
              list => elem%next
              polygon_destroy(elem%p)
              deallocate( elem )
          else
              current => list
              prev    => list
              do while ( associated(current) )
                  if ( associated(current,elem) ) then
                      prev%next => current%next
                      polygon_destroy(elem%p)
                      deallocate( current ) ! Is also "elem"
                      exit
                  endif
                  prev    => current
                  current => current%next
              enddo
          endif
        end subroutine polygon_linked_list_delete_element
    end module polygon_mod
