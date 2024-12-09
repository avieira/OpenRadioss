      !||====================================================================
      !||    SutherlandHodgmanUtil_mod      ../common_source/tools/clipping/SutherlandHodgmanUtil_mod.F90
      !||--- called by ------------------------------------------------------
      !||    Clipping_Sutherland_Hodgman   ../common_source/tools/clipping/polygon_clipping_mod.F90
      !||--- uses       -----------------------------------------------------
      !||    polygon_mod               ../common_source/tools/clipping/polygon_mod.F90
      !||    constant_mod              ../common_source/modules/constant_mod.F
      !||====================================================================
      module SutherlandHodgmanUtil_mod
        use polygon_mod
        implicit none
#include  "my_real.inc"
        ! functions and type needed for Sutherland-Hodgman algorithm
      
        contains 
        
        ! -------------------------------------------------------- !
          ! 
! ======================================================================================================================
!                                                   PROCEDURES
! ======================================================================================================================
! \brief Make the clipping of the polygons in `poly_list` by the line (y1,y2).
! \details This produces two lists of polygon (`clipped_in` and `clipped_out`), a list of points `pts_on_edge` that are
!          extreme points of segments that are included in the segment  [y1, y2], the list of points clipped which are
!          inside `clipped_in` but not on the boundary.
!          `last_edge_clipping` serves as a flag for a special treatment if (y1, y2) is the last line with which we are clipping.
        subroutine edgeClipping( y1, y2, poly_list, clipped_in, clipped_out, pts_on_edge, pts_inside, last_edge_clipping )
          use constant_mod , only : ep20
          implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
          type(polygon_linked_list_), pointer, intent(in) :: poly_list 
          type(polygon_linked_list_), pointer, intent(out) :: clipped_in, clipped_out
          type(polygon_point_), intent(in) :: y1, y2
          type(pair_points_linked_list_), pointer, intent(out) :: pts_on_edge, pts_inside
          logical, value, intent(in):: last_edge_clipping
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local Variables
! ----------------------------------------------------------------------------------------------------------------------
          type(polygon_linked_list_), pointer :: curr_poly, curr_clip_in, curr_clip_out
          type(polygon_), pointer :: poly
          type(polygon_) :: workin, workout
          type(polygon_point_) :: x1, x2, intersecPt, firstpt, p1, newedge1, newedge2
          integer ::  i, ierr
          logical :: positive_side, face_started, d1, d2, isin
          
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
          nullify(clipped_in)
          nullify(clipped_out)
          curr_poly => poly_list
          curr_clip_in => clipped_in
          curr_clip_out => clipped_out
          workin%numpoint = 0
          workout%numpoint = 0

          do while(associated(curr_poly))
            poly => curr_poly%p
            firstpt%y = ep20; firstpt%z = ep20;
            face_started = .false.; positive_side = .false.;
            call polygon_create(workin,  2*poly%numpoint)
            call polygon_create(workout, 2*poly%numpoint)
            
            if (poly%numpoint>1) then
              x2 = poly%point(1)
              d2 = SHUinside(x2, y1, y2)
              if (d2) then 
                ierr = polygon_addpoint(workin, x2)
              else
                ierr = polygon_addpoint(workout, x2)
              end if
            end if

            do i=1,poly%numpoint-1 ! for each edge i of poly
              if (firstpt%y == ep20 .and. firstpt%z == ep20) then 
                firstpt = poly%point(i)
                x2 = firstpt
                d2 = SHUinside(firstpt, y1, y2)
              end if
              x1 = x2   ! vertex 1 of edge i
              d1 = d2
              x2 = poly%point(i+1) ! vertex 2 of edge i
              d2 = SHUinside(x2, y1, y2)
              
              if (d1 .and. d2) then !both points are inside
                ierr = polygon_addpoint(workin, x2)
                positive_side = .true.
                if (last_edge_clipping) then 
                  call pair_points_linked_list_insert_after(pts_inside, [x1, x2])
                end if
              elseif (.not. d1 .and. .not. d2) then !both points are outside
                ierr = polygon_addpoint(workout, x2)
                positive_side = .false.
              else !One point inside and one point outside
                intersecPt = SHUintersection(x1, x2, y1, y2)
                if (.not. face_started) then !opening a new face
                  if (d1) then !starting inside, ending outside
                    ierr = polygon_addpoint(workin, intersecPt)
                    ierr = polygon_addpoint(workout, intersecPt)
                    ierr = polygon_addpoint(workout, x2)
                  else !starting outside, ending inside
                    ierr = polygon_addpoint(workout, intersecPt)
                    ierr = polygon_addpoint(workin, intersecPt)
                    ierr = polygon_addpoint(workin, x2)
                  end if
                  p1 = firstpt
                else !closing a face
                  if (d1) then !starting inside, ending outside
                    ierr = polygon_addpoint(workin, intersecPt)
                    ierr = polygon_addpoint(workin, p1) ! we close a polygon here
                    call polygon_linked_list_insert_after(curr_clip_in, workin)
                    if (.not. associated(clipped_in)) then 
                      clipped_in => curr_clip_in
                    end if
                    curr_clip_in => curr_clip_in%next
                    workin%numpoint = 0
                    isin = SHUinside_edge(intersecPt, p1, y1, y2, newedge1, newedge2)
                    if (isin) then 
                      call pair_points_linked_list_insert_after(pts_on_edge, [newedge1, newedge2])
                    end if
                    if (last_edge_clipping) then 
                      call pair_points_linked_list_insert_after(pts_inside, [x1, intersecPt])
                    end if

                    ierr = polygon_addpoint(workout, intersecPt)
                    ierr = polygon_addpoint(workout, x2)
                  else !starting outside, ending inside
                    ierr = polygon_addpoint(workout, intersecPt)
                    ierr = polygon_addpoint(workout, p1) ! we close a polygon here
                    call polygon_linked_list_insert_after(curr_clip_out, workout)
                    if (.not. associated(clipped_out)) then 
                      clipped_out => curr_clip_out
                    end if
                    curr_clip_out => curr_clip_out%next
                    workout%numpoint = 0

                    ierr = polygon_addpoint(workin, intersecPt)
                    ierr = polygon_addpoint(workin, x2)
                    isin = SHUinside_edge(p1, intersecPt, y1, y2, newedge1, newedge2)
                    if (isin) then 
                      call pair_points_linked_list_insert_after(pts_on_edge, [newedge1, newedge2])
                    end if
                    if (last_edge_clipping) then 
                      call pair_points_linked_list_insert_after(pts_inside, [intersecPt, x2])
                    end if
                  end if
                end if
                positive_side = .not. positive_side
                face_started = .not. face_started
              end if

              if (firstpt%y == x2%y .and. firstpt%z == x2%z) then 
                if (workin%numpoint>0) then 
                  if ((workin%point(1)%y == workin%point(workin%numpoint)%y) .and. & 
                     (workin%point(1)%z == workin%point(workin%numpoint)%z)) then
                    call polygon_linked_list_insert_after(curr_clip_in, workin)
                    if (.not. associated(clipped_in)) then 
                      clipped_in => curr_clip_in
                    end if
                    curr_clip_in => curr_clip_in%next
                  end if
                  workin%numpoint = 0
                end if
                if (workout%numpoint>0) then 
                  if ((workout%point(1)%y == workout%point(workout%numpoint)%y) .and. &
                  (workout%point(1)%z == workout%point(workout%numpoint)%z)) then
                    call polygon_linked_list_insert_after(curr_clip_out, workout)
                    if (.not. associated(clipped_out)) then 
                      clipped_out => curr_clip_out
                    end if
                    curr_clip_out => curr_clip_out%next
                  end if
                  workout%numpoint = 0
                end if
              end if
            end do

            if (workin%numpoint>0) then 
              if ((workin%point(1)%y == workin%point(workout%numpoint)%y) .and. &
                  (workin%point(1)%z == workin%point(workout%numpoint)%z)) then
                call polygon_linked_list_insert_after(curr_clip_in, workin)
                if (.not. associated(clipped_in)) then 
                  clipped_in => curr_clip_in
                end if
                curr_clip_in => curr_clip_in%next
              end if
              workin%numpoint = 0
            end if
            if (workout%numpoint>0) then 
              if ((workout%point(1)%y == workout%point(workout%numpoint)%y) .and. & 
                 (workout%point(1)%z == workout%point(workout%numpoint)%z)) then
                call polygon_linked_list_insert_after(curr_clip_out, workout)
                if (.not. associated(clipped_out)) then 
                  clipped_out => curr_clip_out
                end if
                curr_clip_out => curr_clip_out%next
              end if
              workout%numpoint = 0
            end if
          
            !next polygon
            curr_poly => curr_poly%next
            call polygon_destroy(workin)
            call polygon_destroy(workout)
          end do
        end subroutine edgeClipping
        
        ! -------------------------------------------------------- !
        function SHUintersection( x1, x2, y1, y2)
          use constant_mod , only : zero, one
          implicit none
#include "my_real.inc"
          ! computes the intersection between segment [x1x2] 
          ! and line the line (y1y2) 
      
          ! -- parameters of the function --
          type(polygon_point_) :: x1, x2, &  ! points of the segment
                                  y1, y2     ! points of the line
          
          type(polygon_point_) :: SHUintersection, vx, vy, x1y1 
          my_real :: a
        
          vx%y = x2%y - x1%y;  vx%z = x2%z - x1%z;
          vy%y = y2%y - y1%y;  vy%z = y2%z - y1%z; 
      
          ! if the vectors are colinear
          if ( SHUcrossProduct(vx,vy) .eq. zero) then
            x1y1%y = y1%y - x1%y; x1y1%z = y1%z - x1%z; 
            ! if the the segment [x1x2] is included in the line (y1y2)
            if ( SHUcrossProduct(x1y1,vx) .eq. zero) then
              ! the intersection is the last point of the segment
              SHUintersection = x2
            end if
          else ! the vectors are not colinear
            ! we want to find the inersection between [x1x2]
            ! and (y1,y2).
            ! mathematically, we want to find a in [0;1] such
            ! that :
            !     x1 + a vx = y1 + b vy        
            ! <=> a vx = x1y1 + b vy
            ! <=> a vx^vy = x1y1^vy      , ^ is cross product
            ! <=> a = x1y1^vy / vx^vy
           
            x1y1%y = y1%y - x1%y; x1y1%z = y1%z - x1%z; 
            ! we compute a
            a = SHUcrossProduct(x1y1,vy)/SHUcrossProduct(vx,vy)
            ! if a is not in [0;1]
            if ( (a .lt. one) .and. (a .gt. zero)) then
              SHUintersection%y = x1%y + a*vx%y; SHUintersection%z = x1%z + a*vx%z; 
            end if
          end if
      
        end function SHUintersection
        
        
        ! -------------------------------------------------------- !
        function SHUinside( p, y1, y2)
          use constant_mod , only : zero 
          implicit none
          ! function that tells if the point p is at left of the line (y1y2)
          
          type(polygon_point_) :: p, y1, y2, v1, v2
          logical :: SHUinside

          v1%y = y2%y -  y1%y; v1%z = y2%z -  y1%z;
          v2%y =  p%y -  y1%y; v2%z =  p%z -  y1%z;  
          if ( SHUcrossProduct(v1,v2) .ge. zero ) then
            SHUinside = .true.
          else 
            SHUinside = .false.
          end if
        end function SHUinside
      
        function SHUinside_edge( p1, p2, y1, y2, newp1, newp2)
          use constant_mod , only : zero 
          implicit none
#include "my_real.inc"
          ! function that tells if the intersection of [p1, p2] and [y1, y2] is empty (normally, both segments are aligned)
          ! if [p1, p2] is inluded in [y1, y2], newp1 = p1, newp2 = p2 and returns true
          ! if only p1 (resp. p2) is outside of [y1, y2], then p1 (resp. p2) is projected on [y1, y2] and newp1 = projection of p1 and newp2 = p2 (resp. if p2 is outside), and the function returns true
          ! if both points are outside of [y1, y2], then newpj = pj, j=1,2, and returns false.
          
          type(polygon_point_), intent(in) :: p1, p2, y1, y2 
          type(polygon_point_), intent(out) :: newp1, newp2
          type(polygon_point_) :: v1, v2, p1y1, p1y2, p2y1, p2y2
          my_real :: disty1, disty2
          logical :: SHUinside_edge

          p1y1%y = p1%y-y1%y ; p1y1%z = p1%z-y1%z ; 
          p1y2%y = p1%y-y2%y ; p1y2%z = p1%z-y2%z ; 
          p2y1%y = p2%y-y1%y ; p2y1%z = p2%z-y1%z ; 
          p2y2%y = p2%y-y2%y ; p2y2%z = p2%z-y2%z ; 
        
          if ( p1y1%y * p1y2%y .le. zero .and. p1y1%z * p1y2%z .le. zero ) then !p1 is in [y1, y2]
            if ( p2y1%y * p2y2%y .le. zero .and. p2y1%z * p2y2%z .le. zero ) then !p2 is in [y1, y2]
              newp1 = p1
              newp2 = p2
              SHUinside_edge = .true.
            else !p2 is not in [y1, y2] : we have to project p2 on [y1, y2]
              disty1 = p2y1%y**2 + p2y1%z**2
              disty2 = p2y2%y**2 + p2y2%z**2
              if (disty1 .le. disty2) then 
                newp2 = y1
              else
                newp2 = y2
              end if
              newp1 = p1
              SHUinside_edge = .true.
            end if
          else !p1 is not in [y1, y2]
            if ( p2y1%y * p2y2%y .le. zero .and. p2y1%z * p2y2%z .le. zero ) then !p2 is in [y1, y2] : we have to project p1 on [y1, y2]
              disty1 = p1y1%y**2 + p1y1%z**2
              disty2 = p1y2%y**2 + p1y2%z**2
              if (disty1 .le. disty2) then 
                newp1 = y1
              else
                newp1 = y2
              end if
              newp2 = p2
              SHUinside_edge = .true.
            else !p2 is not in [y1, y2] : no intersection
              newp1 = p1
              newp2 = p2
              SHUinside_edge = .false.
            end if
          end if
        end function SHUinside_edge

        ! -------------------------------------------------------- !
        function SHUcrossProduct( v1, v2)
          implicit none
#include "my_real.inc"
          ! compute the crossproduct of vectors v1 and v2
          type(polygon_point_) :: v1
          type(polygon_point_) :: v2
          my_real :: SHUcrossProduct

          SHUcrossProduct = v1%y*v2%z - v1%z*v2%y
        end function SHUcrossProduct
      
      end module SutherlandHodgmanUtil_mod
