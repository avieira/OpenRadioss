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
!||    compute_lambdas_mod      ../common_source/tools/clipping/SutherlandHodgmanUtil_mod.F90
!||--- called by ------------------------------------------------------
!||--- uses       -----------------------------------------------------
!||    polygon_mod            ../common_source/tools/clipping/polygon_mod.F90
!||====================================================================
module compute_lambdas_mod
  use polygon_mod
  implicit none
#include "my_real.inc"
  
  contains 
  function scalProd(v1, v2)
    implicit none 
#include "my_real.inc"
    type(polygon_point_), intent(in) :: v1, v2
    my_real:: scalProd

    scalProd = v1%y*v2%y + v1%z*v2%z
  end function scalProd


    !\brief Compute the length of the points stored in pts_on_edge.
    !\details The signed legnth of the segments in pts_on_edge are summed.
    !         All the points in pts_on_edge must be in the segment [p1, p2]
    !         The sign of the length is determined by the orientation given by [p1, p2]
    !         Given the points i1, i2 stored in the list,
    !         if the vectors p1p2 and i1i2 are in the same direction, the length of [p1, p2] is positive
    !         otherwise, it is counted negatively.
    !||====================================================================
    !||    compute_effective_surface             ../engine/source/cutcell.F90
    !||--- calls      -----------------------------------------------------
    !||    scalProd                              ../engine/source/cutcell.F90
    !||--- uses       -----------------------------------------------------
    !||    constant_mod                          ../common_source/modules/constant_mod.F
    !||====================================================================
  function compute_effective_surface(pts_on_edge, p1, p2) result(lambda)
    use constant_mod, only : zero
    implicit none
#include  "my_real.inc"
    type(pair_points_linked_list_), pointer, intent(in) :: pts_on_edge !<list of points on the edge
    type(polygon_point_), intent(in) :: p1, p2 !<points defining the edge, needed for orientation
    my_real :: lambda !< effective surface

    type(pair_points_linked_list_), pointer :: curr_pts_on_edge
    type(polygon_point_) :: i1, i2, v1, v2
    my_real :: dy, dz, a1, a2

    lambda = zero
    curr_pts_on_edge => pts_on_edge
    do while(associated(curr_pts_on_edge))
        i1 = curr_pts_on_edge%pts(1)
        i2 = curr_pts_on_edge%pts(2)
        dy = i2%y - i1%y ; dz = i2%z - i1%z ;

        !We have i1 and i2 between p1 and p2. We want to know if i1->i2 is in the same direction as p1->p2.
        !We necessarily have p1->ij = aj * p1->p2 for some 0<aj<1 and j=1,2.
        !Thus, doing the scalar product with p1->p2 on each side, we have:
        !aj = (p1->ij).(p1->p2) / (|p1->p2|^2)
        !if p1->p2 and i1->i2 are in the same direction, then i1 is closer to p1 than i2, thus a1<a2.
        v1%y = i1%y - p1%y; v1%z = i1%z - p1%z;
        v2%y = p2%y - p1%y; v2%z = p2%z - p1%z;
        a1 = scalProd(v1, v2)/scalProd(v2, v2)

        v1%y = i2%y - p1%y; v1%z = i2%z - p1%z;
        a2 = scalProd(v1, v2)/scalProd(v2, v2)

        if (a1 <= a2) then !edge and points in the same direction
            lambda = lambda + sqrt(dy*dy + dz*dz)
        else
            lambda = lambda - sqrt(dy*dy + dz*dz)
        end if
        curr_pts_on_edge => curr_pts_on_edge%next
    end do
  end function compute_effective_surface

  function shoelace_area(poly)
    use constant_mod, only : zero, half
    implicit none
#include  "my_real.inc"
    type(polygon_), intent(in) :: poly
    my_real :: shoelace_area
    my_real :: left_sum, right_sum
    integer :: i, j

    left_sum = zero
    right_sum = zero
    do i=1,poly%numpoint-1
        left_sum  =  left_sum + poly%point(i)%y   * poly%point(i+1)%z
        right_sum = right_sum + poly%point(i+1)%y * poly%point(i)%z
    end do
    left_sum = left_sum + poly%point(poly%numpoint)%y * poly%point(1)%z
    right_sum = right_sum + poly%point(1)%y * poly%point(poly%numpoint)%z

    shoelace_area = half*(left_sum - right_sum)
  end function shoelace_area

    !\brief Compute the area of the polygon.
    !||====================================================================
    !||    compute_effective_volume             ../engine/source/cutcell.F90
    !||--- calls      -----------------------------------------------------
    !||    shoelace_area                              ../engine/source/cutcell.F90
    !||--- uses       -----------------------------------------------------
    !||    constant_mod                          ../common_source/modules/constant_mod.F
    !||====================================================================
  function compute_effective_volume(Clipped_in) result(Lambda)
    use constant_mod, only : zero
    implicit none
#include  "my_real.inc"
    type(polygon_linked_list_), pointer, intent(in) :: Clipped_in 
    my_real :: Lambda !< effective volume

    type(polygon_linked_list_), pointer :: curr_poly

    Lambda = zero
    curr_poly => Clipped_in
    do while(associated(curr_poly))
        Lambda = Lambda + shoelace_area(curr_poly%p)
        curr_poly => curr_poly%next
    end do
  end function compute_effective_volume
end module compute_lambdas_mod