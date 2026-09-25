module geometric_rebuilder_mod
  use precision_mod, only : wp                            !provides kind for eigther single or double precision (wp means working precision)
  use level_set_definition_mod

  implicit none
 
  integer, parameter :: dp = kind(1.0d0)
 
  ! Tolerances 
  real(dp), parameter :: TOL_ZERO  = 1.0d-10   
  real(dp), parameter :: TOL_DENOM = 1.0d-10   
  real(dp), parameter :: ETA_MIN   = 1.0d-3    
 
contains
 
  subroutine find_0pt_Q1(level_set_tn, level_set_tnp1, xi, eta, zeta, ierr)
    !--------------------------------------------------------------------------
    ! Local Q1 coordinated = (xi, eta, zeta) on [0,1]^3.
    !
    ! We are looking for a point on the 0 level set of the Q1 function defined
    ! by level_set_tn at zeta=0 and level_set_tnp1 at zeta=1. We are only looking
    ! for it at zeta = 0.5.
    !
    ! level_set_tn, level_set_tnp1 : 4 values arrays (values at the 4 vertices)
    ! xi, eta, zeta                : coordinates of the point found (output)
    ! ierr                         : 0 if a point has been found
    !                                1 if no point has been found
    !--------------------------------------------------------------------------
    real(dp), dimension(4), intent(in)  :: level_set_tn
    real(dp), dimension(4), intent(in)  :: level_set_tnp1
    real(dp),                intent(out) :: xi
    real(dp),                intent(out) :: eta
    real(dp),                intent(out) :: zeta
    integer,                 intent(out) :: ierr
 
    real(dp), dimension(4) :: level
    real(dp) :: denom, num
    real(dp) :: eta_loc, zeta_loc
 
    ierr = 0
 
    level = 0.5d0 * (level_set_tn + level_set_tnp1)
 
    if (sqrt(sum(level*level)) < TOL_ZERO) then
      xi   = 0.5d0
      eta  = 0.5d0
      zeta = 0.5d0
      return
    end if
 
    !--------------------------------------------------------------------
    ! First trial : looking for a correct xi given eta, diminishing
    ! eta staarting from 1/2.
    !--------------------------------------------------------------------
    eta_loc = 0.5d0
    do while (eta_loc > ETA_MIN)
      denom = eta_loc*(level(1) - level(2) + level(3) - level(4)) - (level(1) - level(2))
      num   = eta_loc*(level(1) - level(4)) - level(1)
      if (abs(denom) > TOL_DENOM) then
        zeta_loc = num/denom
        if (zeta_loc > 0.0d0 .and. zeta_loc < 1.0d0) then
          xi   = zeta_loc
          eta  = eta_loc
          zeta = 0.5d0
          return
        end if
      end if
      eta_loc = eta_loc / 2.0d0
    end do
 
    !--------------------------------------------------------------------
    ! If we get here, no admissible eta has been found.
    ! Let us increase eta from 0.5 to 1 instead.
    !--------------------------------------------------------------------
    eta_loc = 0.5d0
    do while (eta_loc > ETA_MIN)
      denom = (1.0d0-eta_loc)*(level(1) - level(2) + level(3) - level(4)) - (level(1) - level(2))
      num   = (1.0d0-eta_loc)*(level(1) - level(4)) - level(1)
      if (abs(denom) > TOL_DENOM) then
        zeta_loc = num/denom
        if (zeta_loc > 0.0d0 .and. zeta_loc < 1.0d0) then
          xi   = zeta_loc
          eta  = 1.0d0 - eta_loc
          zeta = 0.5d0
          return
        end if
      end if
      eta_loc = eta_loc / 2.0d0
    end do
 
    !--------------------------------------------------------------------
    ! If nothing had been found yet, it could be a straight line.
    ! We thus look for a point using zeta = 0.5.
    !--------------------------------------------------------------------
    zeta_loc = 0.5d0
    denom = zeta_loc*(level(1) - level(2) + level(3) - level(4)) - (level(1) - level(4))
    num   = zeta_loc*(level(1) - level(2)) - level(1)
    if (abs(denom) > TOL_DENOM) then
      eta_loc = num/denom
      if (eta_loc > 0.0d0 .and. eta_loc < 1.0d0) then
        xi   = zeta_loc
        eta  = eta_loc
        zeta = 0.5d0
        return
      end if
    end if
 
    !--------------------------------------------------------------------
    ! No point found: ierr=1.
    !--------------------------------------------------------------------
    ierr = 1
    xi   = 0.0d0
    eta  = 0.0d0
    zeta = 0.0d0
 
  end subroutine find_0pt_Q1
 
  !! \brief Computes the normal(s) to the 0 level set in cell #i_cell.
  !!        The level-set is defined at the nodes of the grid, stored in phi_full.
  subroutine build_polygon_from_level_set(IXQ, X, i_cell, phi, normals, nb_normals)
    use polygon_cutcell_mod, only : Point2D

    IMPLICIT NONE

    integer, dimension(:,:), intent(in) :: IXQ
    real(kind=wp), dimension(:,:), intent(in) :: X
    integer, intent(in) :: i_cell
    real(kind=wp), dimension(:), intent(in) :: phi
    type(Point2D), dimension(2), intent(out) :: normals
    integer, intent(out) :: nb_normals

    integer, dimension(4) :: ips
    integer :: nb_negative, i_negphi, i_negphi2
    integer :: i
    real(kind=wp) :: lambda
    type(Point2D):: pt1, pt2

    ips(1) = IXQ(2, i_cell)
    ips(2) = IXQ(3, i_cell)
    ips(3) = IXQ(4, i_cell)
    ips(4) = IXQ(5, i_cell)

    nb_negative = 0
    do i=1,4
      if (phi(i) <= 0) then
        nb_negative = nb_negative + 1
      endif
    enddo


    if ((nb_negative == 0) .or. (nb_negative == 4)) then
      nb_normals = 0
      return
    elseif ((nb_negative == 1) .or. (nb_negative == 3)) then
      nb_normals = 1

      i_negphi = 0
      if (nb_negative == 1) then
        do i=1,4
          if (phi(i) <= 0.0) then
            i_negphi = i 
            stop
          endif
        enddo
      else
        do i=1,4
          if (phi(i) > 0.0) then
            i_negphi = i 
            stop
          endif
        enddo
      endif

      if (i_negphi == 4) then
        lambda = phi(i_negphi)/(phi(i_negphi)-phi(1))
        pt1%y = lambda*X(2, ips(i_negphi)) + (1-lambda)*(X(2, ips(1)))
        pt1%z = lambda*X(3, ips(i_negphi)) + (1-lambda)*(X(3, ips(1)))
      else
        lambda = phi(i_negphi)/(phi(i_negphi)-phi(i_negphi+1))
        pt1%y = lambda*X(2, ips(i_negphi)) + (1-lambda)*(X(2, ips(i_negphi+1)))
        pt1%z = lambda*X(3, ips(i_negphi)) + (1-lambda)*(X(3, ips(i_negphi+1)))
      endif
      if (i_negphi == 1) then
        lambda = phi(i_negphi)/(phi(i_negphi)-phi(4))
        pt2%y = lambda*X(2, ips(i_negphi)) + (1-lambda)*(X(2, ips(4)))
        pt2%z = lambda*X(3, ips(i_negphi)) + (1-lambda)*(X(3, ips(4)))
      else
        lambda = phi(i_negphi)/(phi(i_negphi)-phi(i_negphi-1))
        pt2%y = lambda*X(2, ips(i_negphi)) + (1-lambda)*(X(2, ips(i_negphi-1)))
        pt2%z = lambda*X(3, ips(i_negphi)) + (1-lambda)*(X(3, ips(i_negphi-1)))
      endif

      normals(1)%y = pt1%z - pt2%z
      normals(1)%z = pt2%y - pt1%y
      lambda = sqrt(normals(1)%y*normals(1)%y + normals(1)%z*normals(1)%z)
      if (nb_negative == 1) then
        normals(1)%y = normals(1)%y/lambda
        normals(1)%z = normals(1)%z/lambda
      else
        normals(1)%y = -normals(1)%y/lambda
        normals(1)%z = -normals(1)%z/lambda
      endif
    elseif (nb_negative == 2) then
      i_negphi = 0
      i_negphi2 = 0
      do i=1,4
        if (phi(i) <= 0.0) then
          if (i_negphi == 0) then
            i_negphi = i
          else
            i_negphi2 = i
            stop
          endif
        endif
      enddo

      !if ((i_negphi2 == 4) .and. (i_negphi == 1)) then
      !  i = i_negphi2
      !  i_negphi2 = i_negphi
      !  i_negphi = i 
      !endif

      !Given the loop built above, i_negphi < i_negphi2
      !Two possible cases : negative points are next to each other
      ! or they are diagonally distributed.
      if (i_negphi == i_negphi2 - 1) then
        nb_normals = 1

        if (i_negphi == 1) then
          lambda = phi(i_negphi)/(phi(i_negphi)-phi(4))
          pt1%y = lambda*X(2, ips(i_negphi)) + (1-lambda)*(X(2, ips(4)))
          pt1%z = lambda*X(3, ips(i_negphi)) + (1-lambda)*(X(3, ips(4)))
        else
          lambda = phi(i_negphi)/(phi(i_negphi)-phi(i_negphi-1))
          pt1%y = lambda*X(2, ips(i_negphi)) + (1-lambda)*(X(2, ips(i_negphi-1)))
          pt1%z = lambda*X(3, ips(i_negphi)) + (1-lambda)*(X(3, ips(i_negphi-1)))
        endif
        if (i_negphi2 == 4) then
          lambda = phi(i_negphi2)/(phi(i_negphi2)-phi(1))
          pt2%y = lambda*X(2, ips(i_negphi2)) + (1-lambda)*(X(2, ips(1)))
          pt2%z = lambda*X(3, ips(i_negphi2)) + (1-lambda)*(X(3, ips(1)))
        else
          lambda = phi(i_negphi2)/(phi(i_negphi2)-phi(i_negphi2-1))
          pt2%y = lambda*X(2, ips(i_negphi2)) + (1-lambda)*(X(2, ips(i_negphi2-1)))
          pt2%z = lambda*X(3, ips(i_negphi2)) + (1-lambda)*(X(3, ips(i_negphi2-1)))
        endif

        normals(1)%y = pt1%z - pt2%z
        normals(1)%z = pt2%y - pt1%y
        lambda = sqrt(normals(1)%y*normals(1)%y + normals(1)%z*normals(1)%z)
        normals(1)%y = -normals(1)%y/lambda
        normals(1)%z = -normals(1)%z/lambda
      elseif ((i_negphi == 1) .and. (i_negphi2 == 4)) then
        nb_normals = 1

        lambda = phi(1)/(phi(1)-phi(2))
        pt2%y = lambda*X(2, ips(1)) + (1-lambda)*(X(2, ips(2)))
        pt2%z = lambda*X(3, ips(1)) + (1-lambda)*(X(3, ips(2)))
        lambda = phi(4)/(phi(4)-phi(3))
        pt1%y = lambda*X(2, ips(4)) + (1-lambda)*(X(2, ips(3)))
        pt1%z = lambda*X(3, ips(4)) + (1-lambda)*(X(3, ips(3)))

        normals(1)%y = pt1%z - pt2%z
        normals(1)%z = pt2%y - pt1%y
        lambda = sqrt(normals(1)%y*normals(1)%y + normals(1)%z*normals(1)%z)
        normals(1)%y = -normals(1)%y/lambda
        normals(1)%z = -normals(1)%z/lambda
      elseif ((i_negphi == 1) .and. (i_negphi2 == 3)) then
        nb_normals = 2

        lambda = phi(1)/(phi(1)-phi(2))
        pt1%y = lambda*X(2, ips(1)) + (1-lambda)*(X(2, ips(2)))
        pt1%z = lambda*X(3, ips(1)) + (1-lambda)*(X(3, ips(2)))
        lambda = phi(1)/(phi(1)-phi(4))
        pt2%y = lambda*X(2, ips(1)) + (1-lambda)*(X(2, ips(4)))
        pt2%z = lambda*X(3, ips(1)) + (1-lambda)*(X(3, ips(4)))

        normals(1)%y = pt1%z - pt2%z
        normals(1)%z = pt2%y - pt1%y
        lambda = sqrt(normals(1)%y*normals(1)%y + normals(1)%z*normals(1)%z)
        normals(1)%y = normals(1)%y/lambda
        normals(1)%z = normals(1)%z/lambda

        lambda = phi(3)/(phi(3)-phi(4))
        pt1%y = lambda*X(2, ips(3)) + (1-lambda)*(X(2, ips(4)))
        pt1%z = lambda*X(3, ips(3)) + (1-lambda)*(X(3, ips(4)))
        lambda = phi(3)/(phi(3)-phi(2))
        pt2%y = lambda*X(2, ips(3)) + (1-lambda)*(X(2, ips(2)))
        pt2%z = lambda*X(3, ips(3)) + (1-lambda)*(X(3, ips(2)))

        normals(2)%y = pt1%z - pt2%z
        normals(2)%z = pt2%y - pt1%y
        lambda = sqrt(normals(2)%y*normals(2)%y + normals(2)%z*normals(2)%z)
        normals(2)%y = normals(2)%y/lambda
        normals(2)%z = normals(2)%z/lambda
      else !((i_negphi == 2) .and. (i_negphi2 == 4))
        nb_normals = 2

        lambda = phi(2)/(phi(2)-phi(3))
        pt1%y = lambda*X(2, ips(2)) + (1-lambda)*(X(2, ips(3)))
        pt1%z = lambda*X(3, ips(2)) + (1-lambda)*(X(3, ips(3)))
        lambda = phi(2)/(phi(2)-phi(1))
        pt2%y = lambda*X(2, ips(2)) + (1-lambda)*(X(2, ips(1)))
        pt2%z = lambda*X(3, ips(2)) + (1-lambda)*(X(3, ips(1)))

        normals(1)%y = pt1%z - pt2%z
        normals(1)%z = pt2%y - pt1%y
        lambda = sqrt(normals(1)%y*normals(1)%y + normals(1)%z*normals(1)%z)
        normals(1)%y = normals(1)%y/lambda
        normals(1)%z = normals(1)%z/lambda

        lambda = phi(4)/(phi(4)-phi(1))
        pt1%y = lambda*X(2, ips(4)) + (1-lambda)*(X(2, ips(1)))
        pt1%z = lambda*X(3, ips(4)) + (1-lambda)*(X(3, ips(1)))
        lambda = phi(4)/(phi(4)-phi(3))
        pt2%y = lambda*X(2, ips(4)) + (1-lambda)*(X(2, ips(3)))
        pt2%z = lambda*X(3, ips(4)) + (1-lambda)*(X(3, ips(3)))

        normals(2)%y = pt1%z - pt2%z
        normals(2)%z = pt2%y - pt1%y
        lambda = sqrt(normals(2)%y*normals(2)%y + normals(2)%z*normals(2)%z)
        normals(2)%y = normals(2)%y/lambda
        normals(2)%z = normals(2)%z/lambda
      endif
    endif

  end subroutine build_polygon_from_level_set

  !!Builds a polyhedron clipped in space-cell i_cell of grid from the information given in level-set at time tn and tn + dt.
  !!The boolean says if the polyhedron encapsulates the second fluid (if true) or the first one (if false)
  subroutine rebuild_polyhedron(level_set_tn, level_set_tnp1, dt, is_reversed)
    use polygon_cutcell_mod, only : Point2D

    IMPLICIT NONE

    real(kind=wp), dimension(:), intent(in) :: level_set_tn
    real(kind=wp), dimension(:), intent(in) :: level_set_tnp1
    real(kind=wp), intent(in) :: dt
    logical, intent(out) :: is_reversed

    integer :: nb_edges, i
    integer(kind=8) :: nb_tn, nb_tnp1, is_reversed_c

    nb_tn = 0
    nb_tnp1 = 0
    do i=1,4
      if (level_set_tn(i) <= 0) then
        nb_tn = nb_tn + 1
      endif
      if (level_set_tnp1(i) <= 0) then
        nb_tnp1 = nb_tnp1 + 1
      endif
    enddo

    call rebuild_polyhedron_fortran(level_set_tn, level_set_tnp1, dt, nb_tn, nb_tnp1, is_reversed_c)
    is_reversed = (is_reversed_c == 0)
  end subroutine rebuild_polyhedron

  !! \brief Computes the values of phi at the corner nodes of cell i_cell.
  subroutine compute_nodes_level_set(NUMELQ, NUMELTG, IXQ, IXTG, X, grid, dg_info_ref, phi_full, i_cell, phi_corner)
    use grid2D_struct_multicutcell_mod
    use polygon_cutcell_mod, only : Point2D

    implicit none 

    integer :: NUMELQ, NUMELTG
    integer, dimension(:,:) :: IXQ, IXTG
    real(kind=wp), dimension(:,:) :: X
    type(grid2D_struct_multicutcell), dimension(:, :) :: grid
    type(DGCell_GaussLegendre_ref), intent(in) :: dg_info_ref
    type(DG_coeffs), dimension(:), intent(in) :: phi_full
    integer, intent(in) :: i_cell
    real(kind=wp), dimension(4), intent(out) :: phi_corner

    integer :: e, ips
    type(Point2D) :: pt

    do e = 1,4
      ips = IXQ(e+1, i_cell)
      pt = Point2D(X(2, ips), X(3, ips))
      phi_corner(e) = evaluate(NUMELQ, NUMELTG, IXQ, IXTG, X, grid, pt, dg_info_ref, phi_full, i_cell)
    end do
  end subroutine compute_nodes_level_set

!! \brief Compute cell occupancies for each phase in grid.
  subroutine multicutcell_compute_lambdas_LS(NUMELQ, NUMELTG, IXQ, IXTG, X, grid, dt, &
                                          dg_info_ref, level_set_tn_full, level_set_tnp1_full, &
                                          rho, vely, velz, p, gamma, wave_type, flag_conservative)
    use polygon_cutcell_mod
    use grid2D_struct_multicutcell_mod
    use riemann_solver_mod
  
    IMPLICIT NONE
  
    ! INPUT argument
    integer, intent(in) :: NUMELQ, NUMELTG
    integer, dimension(:,:), intent(in) :: IXQ, IXTG
    real(kind=wp), dimension(:,:), intent(in) :: X
    real(kind=wp), intent(in) :: dt
    type(DGCell_GaussLegendre_ref), intent(in) :: dg_info_ref
    type(DG_coeffs), dimension(:), intent(in) :: level_set_tn_full
    type(DG_coeffs), dimension(:), intent(in) :: level_set_tnp1_full
    real(kind=wp), dimension(:,:), intent(in) :: rho
    real(kind=wp), dimension(:,:), intent(in) :: vely
    real(kind=wp), dimension(:,:), intent(in) :: velz
    real(kind=wp), dimension(:,:), intent(in) :: p
    real(kind=wp), dimension(:), intent(in) :: gamma
    integer, intent(in) :: wave_type
    logical, intent(in) :: flag_conservative
    ! IN/OUTPUT argument
    type(grid2D_struct_multicutcell), dimension(:, :), intent(inout) :: grid

    !Local variables
    integer(kind=8), parameter:: max_length_array=1000 !I hope 1000 is enough... But it is reasonable.
    integer(kind=8) :: nb_cell, nb_regions, nb_edges
    integer(kind=8) :: j, k
    integer :: i
    real(kind=wp), dimension(:), allocatable :: ptr_lambdas_arr
    real(kind=wp), dimension(:), allocatable :: ptr_big_lambda_n, ptr_big_lambda_np1
    real(kind=wp), dimension(max_length_array) :: normals_y 
    real(kind=wp), dimension(max_length_array) :: normals_z
    real(kind=wp), dimension(max_length_array) :: normals_t
    integer(kind=8) :: nb_normals
    integer(kind=8) :: is_narrowband
    real(kind=wp), dimension(4) :: ys, zs
    integer(kind=8), dimension(4) :: pt_indices
    real(kind=wp), dimension(4) :: level_set_tn, level_set_tnp1
    integer :: print_nb_cell
    real(kind=wp) :: us, vsL, vsR, ps, norm
    real(kind=wp) :: ny, nz, nt
    real(kind=wp) :: psy, psz, pst
    integer(kind=8), dimension(max_length_array) :: local_index_edge
    logical :: is_reversed

    nb_cell = size(grid, 1)
    nb_regions = size(grid, 2)
    if (NUMELQ>0) then
      nb_edges = 4
    else
      nb_edges = 3
    end if

    allocate(ptr_lambdas_arr(4*nb_regions))
    allocate(ptr_big_lambda_n(nb_regions))
    allocate(ptr_big_lambda_np1(nb_regions))

    print_nb_cell = nb_cell/10
    write(*,*) "Doing cell number ", 1, "/", nb_cell
    do i = 1,nb_cell
      if (i>print_nb_cell-1) then 
        write(*,*) "Doing cell number ", i, "/", nb_cell
        print_nb_cell = print_nb_cell + nb_cell/10
        call system('sync')
      end if

      if (grid(i,1)%close_cells) then
        pt_indices = IXQ(2:2+nb_edges-1, i)
        ys(1:nb_edges) = X(2, pt_indices)
        zs(1:nb_edges) = X(3, pt_indices)
        call compute_nodes_level_set(NUMELQ, NUMELTG, IXQ, IXTG, X, grid, dg_info_ref, level_set_tn_full, i, level_set_tn)
        call compute_nodes_level_set(NUMELQ, NUMELTG, IXQ, IXTG, X, grid, dg_info_ref, level_set_tnp1_full, i, level_set_tnp1)
        call build_grid_from_points_fortran(ys, zs, pt_indices, nb_edges) 
        call rebuild_polyhedron(level_set_tn, level_set_tnp1, dt, is_reversed)
        call compute_lambdas2d_fortran(dt, ptr_lambdas_arr, ptr_big_lambda_n, ptr_big_lambda_np1, &
                                      normals_y, normals_z, normals_t, local_index_edge, max_length_array, nb_normals,&
                                      is_narrowband)

        nb_normals = min(nb_normals, max_length_array) !if there is more than max_length_array normals, we only keep the first max_length_array of them...
        ny = 0.0
        nz = 0.0
        nt = 0.0
        psy = 0.0
        psz = 0.0
        pst = 0.0
        do j=1,nb_normals
          ny = ny + normals_y(j)
          nz = nz + normals_z(j)
          nt = nt + normals_t(j)

          if (flag_conservative) then
            norm = sqrt(ny*ny + nz*nz)
            if (norm > 0) then
              call solve_riemann_problem(gamma(1), gamma(2), &
                                          rho(i,1), rho(i,2), &
                                          vely(i,1), vely(i,2), &
                                          velz(i,1), velz(i,2), &
                                          p(i,1), p(i,2), wave_type, &
                                          ny / norm, nz / norm, &
                                          us, vsL, vsR, ps)
              psy = psy + ps * normals_y(j)
              psz = psz + ps * normals_z(j)
              pst = pst + ps * normals_t(j)
            end if
          end if
        end do

        if (is_reversed) then
          grid(i,1)%lambdan_per_cell = ptr_big_lambda_n(1)
          grid(i,1)%lambdanp1_per_cell = ptr_big_lambda_np1(1)
          grid(i,1)%is_narrowband = (is_narrowband>0)
          do k=1,nb_edges
            grid(i,1)%lambda_per_edge(k) = ptr_lambdas_arr((k-1)*nb_regions + 1)
          end do
          grid(i,2)%lambdan_per_cell = ptr_big_lambda_n(2)
          grid(i,2)%lambdanp1_per_cell = ptr_big_lambda_np1(2)
          grid(i,2)%is_narrowband = (is_narrowband>0)
          do k=1,nb_edges
            grid(i,2)%lambda_per_edge(k) = ptr_lambdas_arr((k-1)*nb_regions + 2)
          end do
          grid(i,1)%normal_intern_face_space%y   = ny
          grid(i,1)%normal_intern_face_space%z   = nz
          grid(i,1)%normal_intern_face_time      = nt
          grid(i,2)%normal_intern_face_space%y   = ny
          grid(i,2)%normal_intern_face_space%z   = nz
          grid(i,2)%normal_intern_face_time      = nt

          if (flag_conservative) then
            grid(i,1)%p_normal_intern_face_space%y = psy
            grid(i,1)%p_normal_intern_face_space%z = psz
            grid(i,1)%p_normal_intern_face_time    = pst
            grid(i,2)%p_normal_intern_face_space%y = psy
            grid(i,2)%p_normal_intern_face_space%z = psz
            grid(i,2)%p_normal_intern_face_time    = pst
          end if
        else 
          grid(i,2)%lambdan_per_cell = ptr_big_lambda_n(1)
          grid(i,2)%lambdanp1_per_cell = ptr_big_lambda_np1(1)
          grid(i,2)%is_narrowband = (is_narrowband>0)
          do k=1,nb_edges
            grid(i,2)%lambda_per_edge(k) = ptr_lambdas_arr((k-1)*nb_regions + 1)
          end do
          grid(i,1)%lambdan_per_cell = ptr_big_lambda_n(2)
          grid(i,1)%lambdanp1_per_cell = ptr_big_lambda_np1(2)
          grid(i,1)%is_narrowband = (is_narrowband>0)
          do k=1,nb_edges
            grid(i,1)%lambda_per_edge(k) = ptr_lambdas_arr((k-1)*nb_regions + 2)
          end do
          grid(i,1)%normal_intern_face_space%y   = -ny
          grid(i,1)%normal_intern_face_space%z   = -nz
          grid(i,1)%normal_intern_face_time      = -nt
          grid(i,2)%normal_intern_face_space%y   = -ny
          grid(i,2)%normal_intern_face_space%z   = -nz
          grid(i,2)%normal_intern_face_time      = -nt

          if (flag_conservative) then
            grid(i,1)%p_normal_intern_face_space%y = -psy
            grid(i,1)%p_normal_intern_face_space%z = -psz
            grid(i,1)%p_normal_intern_face_time    = -pst
            grid(i,2)%p_normal_intern_face_space%y = -psy
            grid(i,2)%p_normal_intern_face_space%z = -psz
            grid(i,2)%p_normal_intern_face_time    = -pst
          end if
        end if
      end if
    end do

    deallocate(ptr_lambdas_arr)
    deallocate(ptr_big_lambda_n)
    deallocate(ptr_big_lambda_np1)

    !call compute_close_cells(NUMELQ, NUMELTG, NUMNOD, IXQ, IXTG, grid)

  end subroutine multicutcell_compute_lambdas_LS

 
end module geometric_rebuilder_mod
 
