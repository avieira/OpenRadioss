module level_set_definition_mod
  use precision_mod, only : wp                            !provides kind for eigther single or double precision (wp means working precision)
  use polygon_cutcell_mod, only : Point2D
  use grid2D_struct_multicutcell_mod, only : grid2D_struct_multicutcell
  use quadrature_DG_mod, only : quadrature, basis1d, basis1d_derivative, compute_matrices_proj
  use constant_mod, only : INFINITY, ONE

  implicit none 
  real(wp), parameter :: INF   = INFINITY 
  integer, parameter :: FAR    = 0
  integer, parameter :: NARROW = 1
  integer, parameter :: FROZEN = 2

  type DG_neighbour
    integer :: up, down, left, right
  end type DG_neighbour

  type DGCell_GaussLegendre_ref
    integer :: degree, ndof
    real(wp), dimension(:), allocatable :: xigp
    type(Point2D), dimension(:), allocatable :: vertices
    real(wp), dimension(:, :), allocatable :: Minv, Mprimeinv, W, Wx, Wy
  end type DGCell_GaussLegendre_ref

  type DG_coeffs
    real(wp), dimension(:), allocatable :: coeff
  end type DG_coeffs

  type DG_local_global_index
    integer :: i_cell, local_ind
  end type DG_local_global_index

contains
  subroutine DGCell_GaussLegendre_create(this, degree)
    type(DGCell_GaussLegendre_ref), intent(inout) :: this
    integer, intent(in) :: degree

    integer :: ndof, ngp
    integer :: iq, jq, ix, iy, i, j
    real(wp), dimension(:), allocatable :: xigp, wgp
    real(wp), dimension(:), allocatable :: psi_tmp, psi_x_tmp
    real(wp), dimension(:, :), allocatable :: psix, psiy, psix_prime, psiy_prime

    ndof = degree + 1
    ngp  = ndof + 1

    ! --- quadrature points/weights on the reference element ---
    allocate(xigp(ngp), wgp(ngp))
    call quadrature(ngp, xigp, wgp)

    ! --- 1D basis functions (and derivatives) at each quadrature point ---
    allocate(psix(ndof, ngp), psiy(ndof, ngp))
    allocate(psix_prime(ndof, ngp), psiy_prime(ndof, ngp))
    allocate(psi_tmp(ndof), psi_x_tmp(ndof))
    allocate(this%vertices(ngp*ngp))

    do jq = 1, ngp
        call basis1d(degree, ndof, xigp(jq), psi_tmp)
        call basis1d_derivative(degree, ndof, xigp(jq), psi_x_tmp)
        psix(:, jq)       = psi_tmp
        psix_prime(:, jq) = psi_x_tmp
        psiy(:, jq)       = psi_tmp
        psiy_prime(:, jq) = psi_x_tmp

        do iq = 1, ngp
            i = iq + (jq - 1) * ngp 
            this%vertices(i)%y = xigp(iq) 
            this%vertices(i)%z = xigp(jq) 
        end do
    end do

    ! --- weight matrices W, Wx, Wy ---
    allocate(this%W(ndof*ndof, ngp*ngp), this%Wx(ndof*ndof, ngp*ngp), this%Wy(ndof*ndof, ngp*ngp))
    this%W  = 0.0_wp
    this%Wx = 0.0_wp
    this%Wy = 0.0_wp

    do jq = 1, ngp
        do iq = 1, ngp
            i = iq + (jq - 1) * ngp 
            do ix = 1, ndof
                do iy = 1, ndof
                    j = ix + (iy - 1) * ndof
                    this%W(j, i)  = this%W(j, i)  + wgp(iq) * wgp(jq) * psix(ix, iq)       * psiy(iy, jq)
                    this%Wx(j, i) = this%Wx(j, i) + wgp(iq) * wgp(jq) * psix_prime(ix, iq) * psiy(iy, jq)
                    this%Wy(j, i) = this%Wy(j, i) + wgp(iq) * wgp(jq) * psix(ix, iq)       * psiy_prime(iy, jq)
                end do
            end do
        end do
    end do

    ! --- mass matrix and its "prime" counterpart, already inverted ---
    allocate(this%Minv(ndof*ndof, ndof*ndof), this%Mprimeinv(ndof*ndof, ndof*ndof))
    call compute_matrices_proj(degree, this%Minv, this%Mprimeinv)

    ! --- remaining scalar/array fields ---
    this%degree = degree
    this%ndof   = ndof
    this%xigp   = xigp

    deallocate(xigp, wgp, psix, psiy, psix_prime, psiy_prime, psi_tmp, psi_x_tmp)
  end subroutine DGCell_GaussLegendre_create

  subroutine DGCell_GaussLegendre_destroy(this)
    type(DGCell_GaussLegendre_ref), intent(inout) :: this

    deallocate(this%xigp, this%vertices, this%Minv, this%Mprimeinv)
    deallocate(this%W, this%Wx, this%Wy)
  end subroutine DGCell_GaussLegendre_destroy
 
  subroutine create_dg_coeff(dg_info_ref, nb_cells, coeffs)
    implicit none 
    type(DGCell_GaussLegendre_ref), intent(in) :: dg_info_ref
    integer, intent(in) :: nb_cells
    type(DG_coeffs), dimension(nb_cells), intent(out) :: coeffs 

    integer :: i_cell, ndof

    ndof = (dg_info_ref%ndof)**2
    do i_cell = 1,nb_cells
      allocate(coeffs(i_cell)%coeff(ndof))
    end do
  end subroutine create_dg_coeff

  subroutine destroy_dg_coeff(nb_cells, coeffs)
    implicit none 
    integer, intent(in) :: nb_cells
    type(DG_coeffs), dimension(nb_cells), intent(inout) :: coeffs 

    integer :: i_cell

    do i_cell = 1,nb_cells
      deallocate(coeffs(i_cell)%coeff)
    end do
  end subroutine destroy_dg_coeff

  subroutine compute_normals_LS(NUMELQ, NUMELTG, IXQ, IXTG, X, i_cell, normals, nb_normals)
    use polygon_cutcell_mod

    IMPLICIT NONE

    integer, intent(in) :: i_cell
    integer, intent(in) :: NUMELQ, NUMELTG
    integer, dimension(:,:), intent(in) :: IXQ, IXTG
    real(kind=wp), dimension(:,:), intent(in) :: X
    type(Point2D), dimension(4) :: normals
    integer, intent(out) :: nb_normals

    integer :: j 
    type(Point2D) :: P1, P2, P3, P4
    real(kind=wp) :: norm

    if (NUMELQ > 0) then
      nb_normals = 4

      P1%y = X(2, IXQ(2, i_cell))
      P1%z = X(3, IXQ(2, i_cell))

      P2%y = X(2, IXQ(3, i_cell))
      P2%z = X(3, IXQ(3, i_cell))

      P3%y = X(2, IXQ(4, i_cell))
      P3%z = X(3, IXQ(4, i_cell))

      P4%y = X(2, IXQ(5, i_cell))
      P4%z = X(3, IXQ(5, i_cell))

      normals(1)%y = (P2%z-P1%z)
      normals(1)%z =-(P2%y-P1%y)
      normals(2)%y = (P3%z-P2%z)
      normals(2)%z =-(P3%y-P2%y)
      normals(3)%y = (P4%z-P3%z)
      normals(3)%z =-(P4%y-P3%y)
      normals(4)%y = (P1%z-P4%z)
      normals(4)%z =-(P1%y-P4%y)
    elseif (NUMELTG > 0) then
      nb_normals = 3

      P1%y = X(2, IXTG(2, i_cell))
      P1%z = X(3, IXTG(2, i_cell))

      P2%y = X(2, IXTG(3, i_cell))
      P2%z = X(3, IXTG(3, i_cell))

      P3%y = X(2, IXTG(4, i_cell))
      P3%z = X(3, IXTG(4, i_cell))

      normals(1)%y = (P2%z-P1%z)
      normals(1)%z =-(P2%y-P1%y)
      normals(2)%y = (P3%z-P2%z)
      normals(2)%z =-(P3%y-P2%y)
      normals(3)%y = (P1%z-P3%z)
      normals(3)%z =-(P1%y-P3%y)
    end if

    do j=1,nb_normals
      !Normalize
      norm = sqrt(normals(j)%y*normals(j)%y + normals(j)%z*normals(j)%z)
      normals(j)%y = normals(j)%y / norm
      normals(j)%z = normals(j)%z / norm
    end do
  end subroutine compute_normals_LS

  !! \brief Find the Q1 coordinates (xi, eta) of point pt inside a convex quad defined by vertices.
  subroutine inverse_Q1(pt, vertices, xi, eta, tol, max_iter, ierr) 
    implicit none 
    type(Point2D), intent(in) :: pt ! Target physical point 
    type(Point2D), intent(in) :: vertices(4) ! Quadrilateral vertices 
    real(wp), intent(in) :: tol ! Residual tolerance 
    integer, intent(in) :: max_iter 
    real(wp), intent(out) :: xi, eta 
    integer, intent(out) :: ierr 

    !Local variables 
    integer :: iter 
    real(wp) :: bx, by 
    real(wp) :: cx, cy 
    real(wp) :: dx, dy 
    real(wp) :: pz, py 
    real(wp) :: fx, fy 
    real(wp) :: res2 
    real(wp) :: j11, j12 
    real(wp) :: j21, j22 
    real(wp) :: det 
    real(wp) :: dxi, deta 

    !Bilinear mapping coefficients 
    !p(xi,eta) = x1 + b*xi + c*eta + d*xi*eta 

    bx = vertices(2)%y - vertices(1)%y
    by = vertices(2)%z - vertices(1)%z
    cx = vertices(4)%y - vertices(1)%y
    cy = vertices(4)%z - vertices(1)%z 
    dx = vertices(1)%y - vertices(2)%y - vertices(4)%y + vertices(3)%y 
    dy = vertices(1)%z - vertices(2)%z - vertices(4)%z + vertices(3)%z 

    !Initial guess for Newton method
    xi = 0.5 
    eta = 0.5 
    ierr = 1 

    !Newton iterations 
    do iter = 1, max_iter 
      !Evaluate bilinear mapping 
      py = vertices(1)%y + bx*xi + cx*eta + dx*xi*eta 
      pz = vertices(1)%z + by*xi + cy*eta + dy*xi*eta 

      ! Residual 
      fx = pt%y - py 
      fy = pt%z - pz 

      ! Squared residual norm 
      res2 = fx*fx + fy*fy 
      if (res2 < tol*tol) then 
        ierr = 0 
        return 
      end if 

      ! Jacobian 
      ! J = [ dp_x/dxi dp_x/deta ] 
      !     [ dp_y/dxi dp_y/deta ] 
      j11 = bx + dx*eta 
      j12 = cx + dx*xi 
      j21 = by + dy*eta 
      j22 = cy + dy*xi 

      ! Determinant  
      det = j11*j22 - j12*j21 
      if (abs(det) < 1.d-14) then 
        ierr = 2 
        return 
      end if 

      ! Solve J * [dxi, deta] = [fx, fy] 
      dxi = ( j22*fx - j12*fy ) / det 
      deta = (-j21*fx + j11*fy ) / det 

      ! Newton update 
      xi = xi + dxi 
      eta = eta + deta 
    end do 
  end subroutine inverse_Q1

  !! \brief Evaluate the DG function defined by dg_info_ref, dgc_coeff at point pt.
  function evaluate(NUMELQ, NUMELTG, IXQ, IXTG, X, grid, pt, dg_info_ref, dgc_coeff, i_cell)
    implicit none 
    integer :: NUMELQ, NUMELTG
    integer, dimension(:,:) :: IXQ, IXTG
    real(kind=wp), dimension(:,:) :: X
    type(grid2D_struct_multicutcell), dimension(:, :) :: grid
    type(Point2D) :: pt
    type(DGCell_GaussLegendre_ref) :: dg_info_ref
    type(DG_coeffs), dimension(:) :: dgc_coeff 
    integer:: i_cell
    real(wp) :: evaluate

    integer :: j
    type(Point2D), dimension(4) :: normals
    type(Point2D) :: pt_grid
    real(kind=wp) :: d
    integer :: nb_normals
    logical:: is_inside
    type(Point2D), dimension(4) :: vertices
    real(wp) :: xi, eta
    real(wp), parameter :: tol = 1e-10
    integer, parameter :: max_iter = 20
    integer :: k, nb_cell
    integer :: ierr

    nb_cell = size(grid, 1)
    if (i_cell < 0) then
        do j = 1,nb_cell
          if (any(grid(j, :)%close_cells)) then
            call compute_normals_LS(NUMELQ, NUMELTG, IXQ, IXTG, X, j, normals, nb_normals)

            is_inside = .true.
            do k=1,nb_normals
              pt_grid%y = X(2, IXQ(k+1, j))
              pt_grid%z = X(3, IXQ(k+1, j))
              d = normals(k)%y*(pt_grid%y-pt%y) + normals(k)%z*(pt_grid%z-pt%z)
              if (d<0) then
                is_inside = .false.
                exit
              endif
            enddo

            if (is_inside) then
              i_cell = j
              exit
            endif
          endif
        enddo
    end if

    do j=1,4
      vertices(j)%y = X(2, IXQ(j+1, i_cell))
      vertices(j)%z = X(3, IXQ(j+1, i_cell))
    end do

    call inverse_Q1(pt, vertices, xi, eta, tol, max_iter, ierr)

    if (ierr == 0) then
      evaluate = evaluate_in_ref(xi, eta, dg_info_ref, dgc_coeff(i_cell))
    else
      evaluate = 0.0
    end if
    contains  
    function evaluate_in_ref(xi, eta, dg_info_ref, dgc_coeff) result(phi_value)
      implicit none 
      real(wp) :: xi, eta
      type(DGCell_GaussLegendre_ref) :: dg_info_ref
      type(DG_coeffs) :: dgc_coeff 
      real(wp) :: phi_value

      real(wp), dimension(:), allocatable :: psix, psiy
      integer :: ndof
      integer :: kx, ky
      
      ndof = dg_info_ref%ndof
      allocate(psix(ndof), psiy(ndof))
      call basis1d(ndof-1, ndof, xi, psix)
      call basis1d(ndof-1, ndof, eta, psiy)
      do ky=1,ndof
        do kx=1,ndof
            phi_value = phi_value + dgc_coeff%coeff(kx+(ky-1)*ndof) * psix(kx) * psiy(ky)
        end do
      end do

      deallocate(psix, psiy)
    end function evaluate_in_ref
  end function evaluate

  !!\brief Project function f on DG basis
  subroutine DG_project(NUMELQ, NUMELTG, IXQ, X, dg_info_ref, dgc_coeff, f)
    implicit none 
    interface
      pure function f_real(y,z) result(res)
      use precision_mod, only : wp
      implicit none
      real(kind=wp), intent(in) :: y,z
      real(kind=wp)             :: res
      end function f_real
    end interface

    integer, intent(in) :: NUMELQ, NUMELTG
    integer, dimension(:,:), intent(in) :: IXQ
    real(kind=wp), dimension(:,:), intent(in) :: X
    type(DGCell_GaussLegendre_ref), intent(in) :: dg_info_ref
    type(DG_coeffs), dimension(:), intent(out) :: dgc_coeff 
    procedure(f_real) :: f

    !Local variables
    integer :: nb_cell, ndof, ngp
    real(wp), dimension(:), allocatable :: RHS
    integer :: jq, iq, kq, i_cell
    real(wp) :: pty, ptz, eta, xi
    type(Point2D), dimension(4) :: vertices

    nb_cell = NUMELQ+NUMELTG 
    ndof = dg_info_ref%ndof 
    ngp  = ndof+1

    allocate(RHS(ngp*ngp))

    do i_cell = 1,nb_cell
      do jq=1,4
        vertices(jq)%y = X(2, IXQ(jq+1, i_cell))
        vertices(jq)%z = X(3, IXQ(jq+1, i_cell))
      end do

      do jq=1,ngp
        eta = dg_info_ref%xigp(jq)
        do iq=1,ngp
          xi = dg_info_ref%xigp(iq)
          kq = iq + (jq-1)*ngp

          pty = (1-eta)*(1-xi)*vertices(1)%y + xi*(1-eta)*vertices(2)%y + &
                eta*xi*vertices(3)%y + eta*(1-xi)*vertices(4)%y
          ptz = (1-eta)*(1-xi)*vertices(1)%z + xi*(1-eta)*vertices(2)%z + &
                  eta*xi*vertices(3)%z + eta*(1-xi)*vertices(4)%z
          RHS(kq) = f(pty, ptz)
        end do
      end do

      dgc_coeff(i_cell)%coeff = matmul(dg_info_ref%Minv, matmul(dg_info_ref%W, RHS))
    end do
    
    deallocate(RHS)
  end subroutine DG_project

  !! \brief Lists all cells neighbours of i_cell.
  subroutine all_adjacency_face(ALE_CONNECT, i_cell, other_cells)
    use ALE_CONNECTIVITY_MOD

    implicit none
    TYPE(t_ale_connectivity), INTENT(IN) :: ALE_CONNECT
    integer, intent(in) :: i_cell
    integer, dimension(4), intent(out) :: other_cells

    integer :: IAD2

    IAD2 = ALE_CONNECT%ee_connect%iad_connect(i_cell)
    other_cells = ALE_CONNECT%ee_connect%connected(IAD2 : IAD2 + 3)
  end subroutine all_adjacency_face
 
  !! \brief Defines two corresponding arrays which maps the local gaussian points in a cell to a global numbering of these points in the domain.
  !! \details For point j in cell i, its global numerotation will be local_to_global(i,j).
  !!          For point with global numerotation k, global_to_local(k)%i_cell defines the cell number in which it is located,
  !!                                                global_to_local(k)%local_ind defines its local index in i_cell.
  subroutine build_list_local_to_global(NUMELQ, NUMELTG, IXQ, X, dg_info_ref, &
                                        global_to_local, local_to_global, global_point_coord)
    implicit none 
    integer, intent(in) :: NUMELQ, NUMELTG
    integer, dimension(:,:), intent(in) :: IXQ
    real(kind=wp), dimension(:,:), intent(in) :: X
    type(DGCell_GaussLegendre_ref), intent(in) :: dg_info_ref
    integer, dimension(:,:), allocatable, intent(out) :: local_to_global
    type(Point2D), dimension(:), allocatable, intent(out) :: global_point_coord
    type(DG_local_global_index), dimension(:), allocatable, intent(out) :: global_to_local

    !Local variables
    integer :: ndof, nb_cells, nb_pts, i_cell
    integer :: ind_local, ind_global, i, j
    real(wp) :: eta, xi
    type(Point2D), dimension(4) :: vertices

    ndof = dg_info_ref%ndof + 1
    nb_cells = NUMELQ + NUMELTG
    nb_pts = nb_cells * ndof * ndof

    if (.not. allocated(local_to_global)) then
      allocate(local_to_global(nb_cells, ndof*ndof))
    end if
    if (.not. allocated(global_to_local)) then
      allocate(global_to_local(nb_pts))
    end if
    if (.not. allocated(global_point_coord)) then
      allocate(global_point_coord(nb_pts))
    end if

    do i_cell = 1,nb_cells
      do j=1,4
        vertices(j)%y = X(2, IXQ(j+1, i_cell))
        vertices(j)%z = X(3, IXQ(j+1, i_cell))
      end do

      do j=1,ndof
        eta = dg_info_ref%xigp(j)
        do i=1,ndof
          xi = dg_info_ref%xigp(i)
          ind_local = i + (j-1)*ndof
          ind_global = (i_cell-1) * ndof * ndof + ind_local

          global_to_local(ind_global)%i_cell = i_cell 
          global_to_local(ind_global)%local_ind = ind_local
          local_to_global(i_cell, ind_local) = ind_global
          global_point_coord(ind_global)%y = (1-eta)*(1-xi)*vertices(1)%y + xi*(1-eta)*vertices(2)%y + &
                                                eta*xi*vertices(3)%y + eta*(1-xi)*vertices(4)%y
          global_point_coord(ind_global)%z = (1-eta)*(1-xi)*vertices(1)%z + xi*(1-eta)*vertices(2)%z + &
                                                eta*xi*vertices(3)%z + eta*(1-xi)*vertices(4)%z
        end do
      end do
    end do
  end subroutine build_list_local_to_global

  !! \brief For each point used for the intgration of the level-set, defines its neighbours in the whole domain.
  subroutine build_DG_neighbours(NUMELQ, NUMELTG, ALE_CONNECT, dg_info_ref, local_to_global, global_neighs)
    use ALE_CONNECTIVITY_MOD
    implicit none 
    integer, intent(in) :: NUMELQ, NUMELTG
    TYPE(t_ale_connectivity), INTENT(IN) :: ALE_CONNECT
    type(DGCell_GaussLegendre_ref), intent(in) :: dg_info_ref
    integer, dimension(:,:), intent(in) :: local_to_global
    type(DG_neighbour), dimension(:), allocatable, intent(out) :: global_neighs

    !Local variables
    integer :: ndof, nb_cells, nb_nodes, global_ind, i, i_cell, ind, j
    integer, dimension(4) :: other_cells
    integer :: left_cell, right_cell, down_cell, up_cell    
    integer :: dgc_neighs_right, dgc_neighs_up, dgc_neighs_down ,dgc_neighs_left

    ndof = dg_info_ref%ndof + 1
    nb_cells = NUMELQ + NUMELTG
    nb_nodes = nb_cells * ndof * ndof

    if (.not. allocated(global_neighs)) then
      allocate(global_neighs(nb_nodes))
    end if

    do i_cell=1,nb_cells !for each cell
      call all_adjacency_face(ALE_CONNECT, i_cell, other_cells)

      down_cell  = other_cells(1)
      right_cell = other_cells(2)
      up_cell    = other_cells(3)
      left_cell  = other_cells(4)

      !~~~~~~~~~~~~~~~~~~~~~BOTTOMMOST EDGE~~~~~~~~~~~~~~~~~~~~~~~~~
      j=1
      !---------------------LEFT CORNER-----------------------
      i=1
      ind = 1
      dgc_neighs_right = local_to_global(i_cell, ind + 1)
      dgc_neighs_up    = local_to_global(i_cell, ind + ndof)
      dgc_neighs_down  = -1
      if (down_cell > 0) then
        dgc_neighs_down = local_to_global(down_cell, ndof*(ndof-1) + 1)
      end if
      dgc_neighs_left = -1
      if (left_cell > 0) then
        dgc_neighs_left = local_to_global(left_cell, ndof)
      end if
      global_ind = local_to_global(i_cell, ind)
      global_neighs(global_ind) = DG_neighbour(dgc_neighs_up,&
      dgc_neighs_down,&
      dgc_neighs_left,&
      dgc_neighs_right)
      !---------------------MIDDLE OF EDGE-----------------------
      do i=2,ndof-1
        ind = i
        dgc_neighs_right = local_to_global(i_cell, ind + 1)
        dgc_neighs_up = local_to_global(i_cell, ind + ndof)
        dgc_neighs_left = local_to_global(i_cell, ind - 1)
        dgc_neighs_down = -1
        if (down_cell > 0) then
          dgc_neighs_down = local_to_global(down_cell, i + ndof*(ndof-1))
        end if
        global_ind = local_to_global(i_cell, ind)
        global_neighs(global_ind) = DG_neighbour(dgc_neighs_up,&
        dgc_neighs_down,&
        dgc_neighs_left,&
        dgc_neighs_right)
      end do
      !---------------------RIGHT CORNER-----------------------
      i = ndof
      ind = i
      dgc_neighs_right = -1
      if (right_cell > 0) then
        dgc_neighs_right = local_to_global(right_cell, 1)
      end if
      dgc_neighs_up = local_to_global(i_cell, ind + ndof)
      dgc_neighs_left = local_to_global(i_cell, ind - 1)
      dgc_neighs_down = -1
      if (down_cell > 0) then
        dgc_neighs_down = local_to_global(down_cell, ndof*ndof)
      end if
      global_ind = local_to_global(i_cell, ind)
      global_neighs(global_ind) = DG_neighbour(dgc_neighs_up,&
      dgc_neighs_down,&
      dgc_neighs_left,&
      dgc_neighs_right)

      !~~~~~~~~~~~~~~~~~~~~~NOT EXTREME Y POINTS~~~~~~~~~~~~~~~~~~~
      do j=2,ndof-1
        !---------------------LEFT PART-----------------------
        i=1
        ind = i + (j-1)*ndof
        dgc_neighs_right = local_to_global(i_cell, ind + 1)
        dgc_neighs_up = local_to_global(i_cell, ind + ndof)
        dgc_neighs_left = -1
        if (left_cell > 0) then
          dgc_neighs_left = local_to_global(left_cell, ndof + (j-1)*ndof)
        end if
        dgc_neighs_down = local_to_global(i_cell, ind - ndof)
        global_ind = local_to_global(i_cell, ind)
        global_neighs(global_ind) = DG_neighbour(dgc_neighs_up,&
        dgc_neighs_down,&
        dgc_neighs_left,&
        dgc_neighs_right)
        !---------------------MIDDLE OF CELL-----------------------
        do i=2,ndof-1
          ind = i + (j-1)*ndof
          dgc_neighs_right = local_to_global(i_cell, ind + 1)
          dgc_neighs_up = local_to_global(i_cell, ind + ndof)
          dgc_neighs_left = local_to_global(i_cell, ind - 1)
          dgc_neighs_down = local_to_global(i_cell, ind - ndof)
          global_ind = local_to_global(i_cell, ind)
          global_neighs(global_ind) = DG_neighbour(dgc_neighs_up,&
          dgc_neighs_down,&
          dgc_neighs_left,&
          dgc_neighs_right)
        end do
        !---------------------RIGHT PART-----------------------
        i=ndof
        ind = i + (j-1)*ndof
        dgc_neighs_right = -1
        if (right_cell > 0) then 
          dgc_neighs_right = local_to_global(right_cell, 1 + (j-1)*ndof)
        end if
        dgc_neighs_up = local_to_global(i_cell, ind + ndof)
        dgc_neighs_left = local_to_global(i_cell, ind - 1)
        dgc_neighs_down = local_to_global(i_cell, ind - ndof)
        global_ind = local_to_global(i_cell, ind)
        global_neighs(global_ind) = DG_neighbour(dgc_neighs_up,&
        dgc_neighs_down,&
        dgc_neighs_left,&
        dgc_neighs_right)
      end do

      !~~~~~~~~~~~~~~~~~~~~~TOPMOST EDGE~~~~~~~~~~~~~~~~~~~~~~~~~
      j=ndof
      !---------------------LEFT CORNER-----------------------
      i=1
      ind = i + (j-1)*ndof
      dgc_neighs_right = local_to_global(i_cell, ind + 1)
      dgc_neighs_up = -1
      if (up_cell > 0) then
        dgc_neighs_up = local_to_global(up_cell, i)
      end if
      dgc_neighs_left = -1
      if (left_cell > 0) then
        dgc_neighs_left = local_to_global(left_cell, ndof*ndof)
      end if
      dgc_neighs_down = local_to_global(i_cell, ind - ndof)
      global_ind = local_to_global(i_cell, ind)
      global_neighs(global_ind) = DG_neighbour(dgc_neighs_up,&
      dgc_neighs_down,&
      dgc_neighs_left,&
      dgc_neighs_right)
      !---------------------MIDDLE OF EDGE-----------------------
      do i=2,ndof-1
        ind = i + (j-1)*ndof
        dgc_neighs_right = local_to_global(i_cell, ind + 1)
        dgc_neighs_up = -1
        if (up_cell > 0) then
          dgc_neighs_up = local_to_global(up_cell, i)
        end if
        dgc_neighs_left = local_to_global(i_cell, ind - 1)
        dgc_neighs_down = local_to_global(i_cell, ind - ndof)
        global_ind = local_to_global(i_cell, ind)
        global_neighs(global_ind) = DG_neighbour(dgc_neighs_up,&
        dgc_neighs_down,&
        dgc_neighs_left,&
        dgc_neighs_right)
      end do
      !---------------------RIGHT CORNER-----------------------
      i = ndof
      ind = i + (j-1)*ndof
      dgc_neighs_up = -1
      if (up_cell > 0) then
        dgc_neighs_up = local_to_global(up_cell, i)
      end if
      dgc_neighs_right = -1
      if (right_cell > 0) then
        dgc_neighs_right = local_to_global(right_cell, 1 + (j-1)*ndof)
      end if
      dgc_neighs_left = local_to_global(i_cell, ind - 1)
      dgc_neighs_down = local_to_global(i_cell, ind - ndof)
      global_ind = local_to_global(i_cell, ind)
      global_neighs(global_ind) = DG_neighbour(dgc_neighs_up,&
      dgc_neighs_down,&
      dgc_neighs_left,&
      dgc_neighs_right)
    end do
  end subroutine build_DG_neighbours

  !! \brief Reinitialize dgc_coeff to define a distance function (ie |grad(phi)|=1), which is the distance to
  !!        the 0-level curve.
  !! \details Uses the fast marching algorithm of Sethian et al.
  subroutine reinitialize_signed_distance_2d(NUMELQ, NUMELTG, IXQ, IXTG, X, grid, dg_info_ref, dgc_coeff,&
    global_neighs, local_to_global, global_to_local, global_point_coord, bandwidth)
    use heap_module
    implicit none
    integer, intent(in) :: NUMELQ, NUMELTG
    integer, dimension(:,:), intent(in) :: IXQ, IXTG
    real(kind=wp), dimension(:,:), intent(in) :: X
    type(grid2D_struct_multicutcell), dimension(:, :) :: grid
    type(DGCell_GaussLegendre_ref), intent(in) :: dg_info_ref
    type(DG_coeffs), dimension(:), intent(inout) :: dgc_coeff 
    type(DG_neighbour), dimension(:), intent(in) :: global_neighs
    type(Point2D), dimension(:), intent(in) :: global_point_coord
    integer, dimension(:,:), intent(in) :: local_to_global
    type(DG_local_global_index), dimension(:), intent(in) :: global_to_local
    real(wp), intent(in) :: bandwidth
    
    !Local variables
    integer :: ndof, nb_cells, n, ngp
    integer :: i_cell, i_loc, i, j, ci
    real(wp), dimension(:), allocatable :: sign_phi, phi, phi_vals, T
    integer, dimension(:), allocatable :: state
    integer, dimension(4) :: neighs
    integer :: size_neighs, nb
    type(heap_t) :: pq
    type(Point2D) :: vi, vj
    real(wp) :: phi_i, phi_j, abs_i, abs_j, denom, d_i, d_j, f_min, t_new
    real(wp), dimension(:), allocatable :: RHS
    integer :: iq, jq, locq, globq
    real(wp) :: h
    
    ndof = dg_info_ref%ndof + 1
    nb_cells = NUMELQ + NUMELTG
    n = nb_cells * ndof * ndof
    ngp = dg_info_ref%ndof + 1
    
    allocate(sign_phi(n), phi(n), phi_vals(n))
    allocate(T(n), state(n))
    
    T        = INF
    state    = FAR
    call heap_init(pq, n)
    
    !Save old values
    do i_cell = 1,nb_cells
      do i_loc = 1,ngp*ngp
        i = local_to_global(i_cell, i_loc)
        if (i < 1) then
          continue
        end if
        vi = global_point_coord(i)
        phi_vals(i) = evaluate(NUMELQ, NUMELTG, IXQ, IXTG, X, grid, vi, dg_info_ref, dgc_coeff, i_cell)
        sign_phi(i) = sign(ONE, phi_vals(i))
      end do
    end do
    
    !First: initialize points closest to 0 level line.
    do i_cell = 1,nb_cells
      do i_loc = 1,ngp*ngp
        i = local_to_global(i_cell, i_loc)
        if (i < 1) then
          continue
        end if
        phi_i = phi_vals(i)
        
        ! Up point
        j = global_neighs(i)%up
        if (j > 0) then
          phi_j = phi_vals(j)
          if (phi_i * phi_j <= 0.0) then
            vi = global_point_coord(i)
            vj = global_point_coord(j)
            h = sqrt((vi%z - vj%z)**2 + (vi%y - vj%y)**2)
            
            abs_i = abs(phi_i)
            abs_j = abs(phi_j)
            denom = abs_i + abs_j
            
            if (denom == 0.0) then
              d_i = 0.0
              d_j = 0.0
            else
              d_i = h * abs_i / denom
              d_j = h * abs_j / denom
            end if
            
            if (d_i < T(i)) then
              T(i)     = d_i
              state(i) = NARROW
              call heap_insert(pq, i, d_i)
            end if
            if (d_j < T(j)) then
              T(j)     = d_j
              state(j) = NARROW
              call heap_insert(pq, j, d_j)
            end if
          end if
        end if
        
        ! Bottom point
        j = global_neighs(i)%down
        if (j > 0) then
          phi_j = phi_vals(j)
          if (phi_i * phi_j <= 0.0) then
            vi = global_point_coord(i)
            vj = global_point_coord(j)
            h = sqrt((vi%z - vj%z)**2 + (vi%y - vj%y)**2)
            
            abs_i = abs(phi_i)
            abs_j = abs(phi_j)
            denom = abs_i + abs_j
            
            if (denom == 0.0) then
              d_i = 0.0
              d_j = 0.0
            else
              d_i = h * abs_i / denom
              d_j = h * abs_j / denom
            end if
            
            if (d_i < T(i)) then
              T(i)     = d_i
              state(i) = NARROW
              call heap_insert(pq, i, d_i)
            end if
            if (d_j < T(j)) then
              T(j)     = d_j
              state(j) = NARROW
              call heap_insert(pq, j, d_j)
            end if
          end if
        end if
        
        ! Left point
        j = global_neighs(i)%left
        if (j > 0) then
          phi_j = phi_vals(j)
          if (phi_i * phi_j <= 0.0) then
            vi = global_point_coord(i)
            vj = global_point_coord(j)
            h = sqrt((vi%z - vj%z)**2 + (vi%y - vj%y)**2)
            
            abs_i = abs(phi_i)
            abs_j = abs(phi_j)
            denom = abs_i + abs_j
            
            if (denom == 0.0) then
              d_i = 0.0
              d_j = 0.0
            else
              d_i = h * abs_i / denom
              d_j = h * abs_j / denom
            end if
            
            if (d_i < T(i)) then
              T(i)     = d_i
              state(i) = NARROW
              call heap_insert(pq, i, d_i)
            end if
            if (d_j < T(j)) then
              T(j)     = d_j
              state(j) = NARROW
              call heap_insert(pq, j, d_j)
            end if
          end if
        end if
        
        ! Right point
        j = global_neighs(i)%right
        if (j > 0) then
          phi_j = phi_vals(j)
          if (phi_i * phi_j <= 0.0) then
            vi = global_point_coord(i)
            vj = global_point_coord(j)
            h = sqrt((vi%z - vj%z)**2 + (vi%y - vj%y)**2)
            
            abs_i = abs(phi_i)
            abs_j = abs(phi_j)
            denom = abs_i + abs_j
            
            if (denom == 0.0) then
              d_i = 0.0
              d_j = 0.0
            else
              d_i = h * abs_i / denom
              d_j = h * abs_j / denom
            end if
            
            if (d_i < T(i)) then
              T(i)     = d_i
              state(i) = NARROW
              call heap_insert(pq, i, d_i)
            end if
            if (d_j < T(j)) then
              T(j)     = d_j
              state(j) = NARROW
              call heap_insert(pq, j, d_j)
            end if
          end if
        end if
      end do
    end do
    
    ! FMM propagation restricted to the narrow band
    do while (.not. (heap_is_empty(pq)))
      call heap_extract_min(pq, ci, f_min)
      if (state(ci) == FROZEN) then
        continue
      end if
      state(ci) = FROZEN
      
      size_neighs = 0
      i = global_neighs(ci)%up
      if (i>0) then
        size_neighs = size_neighs + 1
        neighs(size_neighs) = i
      end if
      i = global_neighs(ci)%down
      if (i>0) then
        size_neighs = size_neighs + 1
        neighs(size_neighs) = i
      end if
      i = global_neighs(ci)%left
      if (i>0) then
        size_neighs = size_neighs + 1
        neighs(size_neighs) = i
      end if
      i = global_neighs(ci)%right
      if (i>0) then
        size_neighs = size_neighs + 1
        neighs(size_neighs) = i
      end if
      
      do i=1,size_neighs
        nb = neighs(i)
        if (state(nb) == FROZEN) then
          continue
        end if
        t_new = solve_eikonal_2d(T, global_neighs, global_point_coord, global_to_local, nb)
        
        if ((t_new < T(nb)) .and. (t_new < bandwidth)) then
          T(nb)     = t_new
          state(nb) = NARROW
          call heap_insert(pq, nb, t_new)
        end if
      end do
    end do
    
    ! Outside the band: keep the original |φ| (sign restored below)
    do i_cell = 1,nb_cells
      do i_loc = 1,ngp*ngp
        i = local_to_global(i_cell, i_loc)
        if (i < 1) then
          continue
        end if
        if (state(i) == FAR) then
          vi = global_point_coord(i)
          T(i) = abs(evaluate(NUMELQ, NUMELTG, IXQ, IXTG, X, grid, vi, dg_info_ref, dgc_coeff, i_cell))
        end if
      end do
    end do
    
    phi = T*sign_phi
    
    
    !Reproject it on DG
    allocate(RHS(ngp*ngp))
    do i_cell = 1,nb_cells
      RHS = 0.0
      do iq=1,ngp
        do jq=1,ngp
          locq = iq + (jq-1)*ngp
          globq = local_to_global(i_cell, locq)
          RHS(locq) = phi(globq)
        end do
      end do
      
      dgc_coeff(i_cell)%coeff = matmul(dg_info_ref%Minv, matmul(dg_info_ref%W, RHS))
    end do

    deallocate(sign_phi, phi, phi_vals)
    deallocate(T, state)
    deallocate(RHS)
    call heap_destroy(pq)
    
    contains
    function solve_eikonal_2d(T, global_neighs, global_point_coord, global_to_local, nb_mark) result(t_min)
      implicit none
      real(wp), dimension(:) :: T
      type(DG_neighbour), dimension(:) :: global_neighs
      type(Point2D), dimension(:) :: global_point_coord
      type(DG_local_global_index), dimension(:) :: global_to_local
      integer :: nb_mark
      real(wp) :: t_min
      
      !Local variables
      integer :: i_cell, i_loc, i, size_nbs, size_frozen, nb, j
      integer, dimension(4) :: nbs
      real(wp), dimension(4) :: frozen1, frozen2
      type(Point2D) :: v_mark, vn
      real(wp) :: T1, T2, h1, h2, a, b, c, disc, h, t_cand
      
      t_min  = INF
      i_cell = global_to_local(nb_mark)%i_cell
      i_loc = global_to_local(nb_mark)%local_ind
      
      size_nbs = 0
      i = global_neighs(nb_mark)%up
      if (i>0) then
        size_nbs = size_nbs + 1
        nbs(size_nbs) = i 
      end if
      i = global_neighs(nb_mark)%down
      if (i>0) then
        size_nbs = size_nbs + 1
        nbs(size_nbs) = i 
      end if
      i = global_neighs(nb_mark)%left
      if (i>0) then
        size_nbs = size_nbs + 1
        nbs(size_nbs) = i 
      end if
      i = global_neighs(nb_mark)%right
      if (i>0) then
        size_nbs = size_nbs + 1
        nbs(size_nbs) = i 
      end if
      v_mark = global_point_coord(nb_mark)
      
      size_frozen = 0
      do i = 1,size_nbs
        nb = nbs(i)
        if (T(nb) < INF) then
          vn = global_point_coord(nb)
          h  = sqrt((v_mark%y - vn%y)**2 + (v_mark%z - vn%z)**2)
          size_frozen = size_frozen + 1
          frozen1(size_frozen) = T(nb)
          frozen1(size_frozen) = h
          t_min = min(t_min, T(nb) + h)
        end if
      end do
      
      ! 2-D quadratic update
      do i = 1, size_frozen
        do j = i+1,size_frozen
          T1  = frozen1(i)
          h1  = frozen2(i)
          T2  = frozen1(j)
          h2  = frozen2(j)
          a    =  1.0/h1**2 + 1.0/h2**2
          b    = -2.0*(T1/h1**2 + T2/h2**2)
          c    =  T1**2/h1**2 + T2**2/h2**2 - 1.0
          disc = b**2 - 4.0*a*c
          if (disc >= 0.0) then
            t_cand = (-b + sqrt(disc)) / (2.0 * a)
            if (t_cand >= max(T1, T2)) then
              t_min = min(t_min, t_cand)
            end if
          end if
        end do
      end do
    end function solve_eikonal_2d
  end subroutine reinitialize_signed_distance_2d



end module level_set_definition_mod
