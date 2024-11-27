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
! \brief Make the clipping  of the polygon by the line (y1y2)
      !||====================================================================
      !||    edgeClipping                          ../common_source/tools/clipping/SutherlandHodgmanUtil_mod.F90
      !||--- calls      -----------------------------------------------------
      !||    SHUintersection                       ../common_source/tools/clipping/SutherlandHodgmanUtil_mod.F90
      !||    SHUinside                             ../common_source/tools/clipping/SutherlandHodgmanUtil_mod.F90
      !||--- uses       -----------------------------------------------------
      !||    constant_mod                          ../common_source/modules/constant_mod.F
      !||====================================================================
        subroutine edgeClipping( y1, y2, poly_list, clipped_in, clipped_out )
          use constant_mod , only : ep20
          implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
          type(polygon_linked_list_), pointer :: poly_list, clipped_in, clipped_out
          type(polygon_point_) :: y1, y2
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local Variables
! ----------------------------------------------------------------------------------------------------------------------
          type(polygon_linked_list_), pointer :: curr_poly, curr_clip_in, curr_clip_out
          type(polygon_), pointer :: poly
          type(polygon_) :: workin, workout
          type(polygon_point_) :: x1, x2, intersecPt, firstpt, p1
          integer ::  i
          logical :: positive_side, face_started, d1, d2
          
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
            polygon_create(workin, curr_poly%p%size)
            polygon_create(workout, curr_poly%p%size)
            
            if (poly%numpoint>1) then
              x2 = poly%point(1)
              d2 = SHUinside(x2, y1, y2)
              if (d2) then 
                polygon_addpoint(workin, x2)
              else
                polygon_addpoint(workout, x2)
              end if
            end if

            do i=1,poly%numpoint-1 ! for each edge i of poly
              if (firstpt%y == ep20 .and. firstpt%z == ep20) then 
                firstpt = poly%point(i)
                d2 = SHUinside(firstpt, y1, y2)
              end if
              x1 = x2   ! vertex 1 of edge i
              d1 = d2
              x2 = poly%point(i+1) ! vertex 2 of edge i
              d2 = SHUinside(x2, y1, y2)
              
              if (d1 .and. d2) then !both points are inside
                polygon_addpoint(workin, x2)
                positive_side = .true.
              elseif (.not. d1 .and. .not. d2) !both points are outside
                polygon_addpoint(workout, x2)
                positive_side = .false.
              else !One point inside and one point outside
                intersecPt = SHUintersection(x1, x2, y1, y2)
                if (.not. face_started) then !opening a new face
                  if (d1) then !starting inside, ending outside
                    polygon_addpoint(workin, intersecPt)
                    polygon_addpoint(workout, intersecPt)
                    polygon_addpoint(workout, x2)
                  else !starting outside, ending inside
                    polygon_addpoint(workout, intersecPt)
                    polygon_addpoint(workin, intersecPt)
                    polygon_addpoint(workin, x2)
                  end if
                  p1 = firstpt
                else !closing a face
                  if (d1) then !starting inside, ending outside
                    polygon_addpoint(workin, intersecPt)
                    polygon_addpoint(workin, p1) ! we close a polygon here
                    polygon_linked_list_insert_after(curr_clip_in, workin)
                    if (.not. associated(clipped_in)) then 
                      clipped_in => curr_clip_in
                    end if
                    curr_clip_in => curr_clip_in%next
                    workin%numpoint = 0

                    polygon_addpoint(workout, intersecPt)
                    polygon_addpoint(workout, x2)
                  else !starting outside, ending inside
                    polygon_addpoint(workout, intersecPt)
                    polygon_addpoint(workout, p1) ! we close a polygon here
                    polygon_linked_list_insert_after(curr_clip_out, workout)
                    if (.not. associated(clipped_out)) then 
                      clipped_out => curr_clip_out
                    end if
                    curr_clip_out => curr_clip_out%next
                    workout%numpoint = 0

                    polygon_addpoint(workin, intersecPt)
                    polygon_addpoint(workin, x2)
                  end if
                end if
                positive_side = .not. positive_side
                face_started = .not. face_started
              end if

              if (firstpt%y == x2%y .and. firstpt%z == x2%z) then 
                if (workin%numpoint>0) then 
                  if ((workin%point(1)%y == workin%point(numpoint)%y) .and. (workin%point(1)%z == workin%point(numpoint)%z)) then
                    polygon_linked_list_insert_after(curr_clip_in, workin)
                    if (.not. associated(clipped_in)) then 
                      clipped_in => curr_clip_in
                    end if
                    curr_clip_in => curr_clip_in%next
                  end if
                  workin%numpoint = 0
                end if
                if (workout%numpoint>0) then 
                  if ((workout%point(1)%y == workout%point(numpoint)%y) .and. (workout%point(1)%z == workout%point(numpoint)%z)) then
                    polygon_linked_list_insert_after(curr_clip_out, workout)
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
              if ((workin%point(1)%y == workin%point(numpoint)%y) .and. (workin%point(1)%z == workin%point(numpoint)%z)) then
                polygon_linked_list_insert_after(curr_clip_in, workin)
                if (.not. associated(clipped_in)) then 
                  clipped_in => curr_clip_in
                end if
                curr_clip_in => curr_clip_in%next
              end if
              workin%numpoint = 0
            end if
            if (workout%numpoint>0) then 
              if ((workout%point(1)%y == workout%point(numpoint)%y) .and. (workout%point(1)%z == workout%point(numpoint)%z)) then
                polygon_linked_list_insert_after(curr_clip_out, workout)
                if (.not. associated(clipped_out)) then 
                  clipped_out => curr_clip_out
                end if
                curr_clip_out => curr_clip_out%next
              end if
              workout%numpoint = 0
            end if
          
            !next polygon
            curr_poly => curr_poly%next
            polygon_destroy(workin)
            polygon_destroy(workout)
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
      
      end module SutherlandHodgmanUtil
