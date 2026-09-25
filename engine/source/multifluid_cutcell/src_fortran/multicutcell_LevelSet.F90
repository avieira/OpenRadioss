module level_set_mod
  use precision_mod, only : wp                            !provides kind for eigther single or double precision (wp means working precision)
  use polygon_cutcell_mod, only : Point2D
  use quadrature_DG_mod, only : basis1d, basis1d_derivative
  use grid2D_struct_multicutcell_mod, only : grid2D_struct_multicutcell
  use constant_mod, only : INFINITY
  use level_set_definition_mod

  implicit none
 
  ! Tolerances 
  real(wp), parameter :: TOL_ZERO  = 1.0d-10   
  real(wp), parameter :: TOL_DENOM = 1.0d-10   
  real(wp), parameter :: ETA_MIN   = 1.0d-3    

contains


  !! \brief Solve the advection level-set equation using an advection then projection of DG Gaussian nodes.
  subroutine solve_levelset_projection(NUMELQ, NUMELTG, IXQ, IXTG, X, ALE_CONNECT, grid, &
                                    dg_info_ref, phi, phi_new,&
                                    rho, vely, velz, p, gamma, dt,&
                                    reinit, full_reinit, bandwidth, h_max&!,
                                    !!local_to_global::Array{Int}, global_to_local::Vector{Set{Tuple{Int, Int}}}, 
                                    !global_point_coord::Vector{Point2D}, 
                                    )
    use ALE_CONNECTIVITY_MOD
    
    implicit none 

    integer, intent(in) :: NUMELQ, NUMELTG
    integer, dimension(:,:), intent(in) :: IXQ, IXTG
    real(kind=wp), dimension(:,:), intent(in) :: X
    TYPE(t_ale_connectivity), INTENT(IN) :: ALE_CONNECT
    type(grid2D_struct_multicutcell), dimension(:, :), intent(in) :: grid
    type(DGCell_GaussLegendre_ref), intent(in) :: dg_info_ref
    type(DG_coeffs), dimension(:), intent(in) :: phi
    type(DG_coeffs), dimension(:), intent(out) :: phi_new
    real(kind=wp), dimension(:,:), intent(in) :: vely
    real(kind=wp), dimension(:,:), intent(in) :: velz
    real(kind=wp), dimension(:,:), intent(in) :: rho
    real(kind=wp), dimension(:,:), intent(in) :: p
    real(kind=wp), dimension(:), intent(in) :: gamma
    real(kind=wp), intent(in) :: dt
    logical, intent(in) :: reinit, full_reinit
    real(kind=wp), intent(inout) :: bandwidth, h_max

    !Local variables
    type(Point2D), dimension(:), allocatable :: middle_pts
    real(kind=wp), dimension(:), allocatable :: phi_vals
    logical, dimension(:), allocatable :: in_band
    integer :: i_cell, jq, iq, kq
    type(Point2D) :: pt, ptp, ptm, ptn, middle_pt
    type(Point2D), dimension(4) :: vertices
    integer :: ndof, ngp, i_cellm, i_cellp, nb_cell
    real(kind=wp), dimension(:), allocatable :: RHS, RHSx, RHSy
    real(kind=wp) :: eta, xi, phi_i
    type(Point2D), dimension(:), allocatable :: vel_LS

    if (bandwidth < 0.0) then
      bandwidth = 10.0 * h_max
    end if

    nb_cell = size(vely, 1)
    allocate(middle_pts(nb_cell))
    allocate(phi_vals(nb_cell))
    allocate(in_band(nb_cell))
    allocate(vel_LS(nb_cell))

    do i_cell = 1,nb_cell
      pt = Point2D(0.0,0.0)
      do kq = 1,4
          pt%y = pt%y + 0.25*X(2, IXQ(kq+1, i_cell))
          pt%z = pt%z + 0.25*X(3, IXQ(kq+1, i_cell))
      end do
      middle_pts(i_cell) = pt
      phi_vals(i_cell) = evaluate(NUMELQ, NUMELTG, IXQ, IXTG, X, grid, pt, dg_info_ref, phi, i_cell)
      in_band(i_cell) = (abs(phi_vals(i_cell)) < bandwidth) 
    end do

    call extend_velocity_2d(NUMELQ, NUMELTG, IXQ, IXTG, X, ALE_CONNECT, grid, &
                             phi, dg_info_ref, middle_pts, phi_vals, &
                             rho, vely, velz, p, gamma,&
                             in_band, vel_LS)

    ndof = dg_info_ref%ndof
    ngp = ndof + 1
    allocate(RHS(ngp*ngp))
    allocate(RHSx(ngp*ngp))
    allocate(RHSy(ngp*ngp))
    
    !On each cell, pull back the quadrature points on the velocity field
    do i_cell=1,nb_cell
      if (in_band(i_cell)) then
        RHS = 0.0
        middle_pt = middle_pts(i_cell)
        phi_i = phi_vals(i_cell)
        
        do jq=1,4
          vertices(jq)%y = X(2, IXQ(jq+1, i_cell))
          vertices(jq)%z = X(3, IXQ(jq+1, i_cell))
        end do

        do jq=1,ngp
          eta = dg_info_ref%xigp(jq)
          do iq=1,ngp
            xi = dg_info_ref%xigp(iq)
            kq = iq + (jq-1)*ngp
            ptn%y = (1-eta)*(1-xi)*vertices(1)%y + xi*(1-eta)*vertices(2)%y + &
                    eta*xi*vertices(3)%y + eta*(1-xi)*vertices(4)%y
            ptn%z = (1-eta)*(1-xi)*vertices(1)%z + xi*(1-eta)*vertices(2)%z + &
                    eta*xi*vertices(3)%z + eta*(1-xi)*vertices(4)%z
            !RK2
            ptm%y = ptn%y - 0.5*dt*vel_LS(i_cell)%y
            ptm%z = ptn%z - 0.5*dt*vel_LS(i_cell)%z
            i_cellm = find_ind_encompassing_cell(NUMELQ, NUMELTG, IXQ, IXTG, X, grid, ptm)
            if (i_cellm > 0) then
              ptp%y = ptn%y - dt * vel_LS(i_cellm)%y
              ptp%z = ptn%y - dt * vel_LS(i_cellm)%z
              i_cellp = find_ind_encompassing_cell(NUMELQ, NUMELTG, IXQ, IXTG, X, grid, ptp)

              if (i_cellp > 0) then
                RHS(kq) = evaluate(NUMELQ, NUMELTG, IXQ, IXTG, X, grid, ptp, dg_info_ref, phi, i_cellp)
                ptm = evaluate_grad(NUMELQ, NUMELTG, IXQ, IXTG, X, grid, ptp, dg_info_ref, phi, i_cellp)
                RHSx(kq) = ptm%y
                RHSy(kq) = ptm%z
              else
                RHS(kq) = evaluate(NUMELQ, NUMELTG, IXQ, IXTG, X, grid, pt, dg_info_ref, phi, i_cell)
                ptm = evaluate_grad(NUMELQ, NUMELTG, IXQ, IXTG, X, grid, pt, dg_info_ref, phi, i_cell)
                RHSx(kq) = ptm%y
                RHSy(kq) = ptm%z
              end if
            else 
              RHS(kq) = evaluate(NUMELQ, NUMELTG, IXQ, IXTG, X, grid, pt, dg_info_ref, phi, i_cell)
              ptm = evaluate_grad(NUMELQ, NUMELTG, IXQ, IXTG, X, grid, pt, dg_info_ref, phi, i_cell)
              RHSx(kq) = ptm%y
              RHSy(kq) = ptm%z
            end if
          end do
        end do
        phi_new(i_cell)%coeff = matmul(dg_info_ref%Mprimeinv,  (&
                                        matmul(dg_info_ref%W, RHS) +&
                                        matmul(dg_info_ref%Wx, RHSx) +&
                                        matmul(dg_info_ref%Wy, RHSy)&
                                      ))
      end if
    end do

    if ((reinit) .and. (.not. (full_reinit))) then
    !    println("Reinitializing...")
    !    @time reinitialize_signed_distance_2d(grid, dg_info, phi_new, local_to_global, global_to_local, global_point_coord, 2*bandwidth)
    end if
    if (full_reinit) then
    !    println("Reinitializing all level set...")
    !    @time reinitialize_signed_distance_2d(grid, dg_info, phi_new, local_to_global, global_to_local, global_point_coord, Inf)
    end if

    deallocate(middle_pts)
    deallocate(phi_vals)
    deallocate(in_band)
    deallocate(RHS)
    deallocate(RHSx)
    deallocate(RHSy)
    deallocate(vel_LS)

    contains      
    !! \brief Given pt, finds in what cell this point is.
    !! \details If pt is out of the domain, is_pt_cell == -1.
    function find_ind_encompassing_cell(NUMELQ, NUMELTG, IXQ, IXTG, X, grid, pt) result(id_pt_cell)
      use grid2D_struct_multicutcell_mod
      use polygon_cutcell_mod
  
      IMPLICIT NONE
  
      ! INPUT argument
      integer :: NUMELQ, NUMELTG
      integer, dimension(:,:) :: IXQ, IXTG
      real(kind=wp), dimension(:,:) :: X
      type(grid2D_struct_multicutcell), dimension(:, :) :: grid
      type(Point2D) :: pt
      ! OUTPUT argument
      integer :: id_pt_cell

      !Dummy arguments
      type(Point2D), dimension(4) :: normals
      type(Point2D) :: pt_grid
      real(kind=wp) :: d
      integer :: nb_normals
      integer :: nb_cell
      integer :: i, j, k
      logical:: is_inside

      nb_cell = NUMELQ+NUMELTG !size(vely, 1)
      id_pt_cell = -1
      call get_clipped_ith_vertex_fortran(i, pt) 
      is_inside = .true.

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
            id_pt_cell = j
            exit
          endif
        endif
      enddo
    end function find_ind_encompassing_cell

    !! \brief Evaluate the gradient of the DG function (defined with dg_info_ref and dgc_coeff) in the reference cell, at the reference coordinates (xi, eta).
    function evaluate_grad_in_ref(xi, eta, dg_info_ref, dgc_coeff) result(phi_value)
      implicit none 
      real(wp) :: xi, eta
      type(DGCell_GaussLegendre_ref) :: dg_info_ref
      type(DG_coeffs) :: dgc_coeff 
      type(Point2D) :: phi_value

      real(wp), dimension(:), allocatable :: psix, psiy
      real(wp), dimension(:), allocatable :: psix_derivative, psiy_derivative
      integer :: ndof
      integer :: kx, ky

      ndof = dg_info_ref%ndof
      allocate(psix(ndof), psiy(ndof))
      allocate(psix_derivative(ndof), psiy_derivative(ndof))
      call basis1d(ndof-1, ndof, xi, psix)
      call basis1d(ndof-1, ndof, eta, psiy)
      call basis1d_derivative(ndof-1, ndof, xi, psix_derivative)
      call basis1d_derivative(ndof-1, ndof, eta, psiy_derivative)
      do ky=1,ndof
        do kx=1,ndof
            phi_value%y = phi_value%y + dgc_coeff%coeff(kx+(ky-1)*ndof) * psix_derivative(kx) * psiy(ky)
            phi_value%z = phi_value%z + dgc_coeff%coeff(kx+(ky-1)*ndof) * psix(kx) * psiy_derivative(ky)
        end do
      end do

      deallocate(psix, psiy)
      deallocate(psix_derivative, psiy_derivative)
    end function evaluate_grad_in_ref

    !! \brief Evaluate the gradient of the DG function (defined with dg_info_ref and dgc_coeff) at point pt in cell i_cell.
    function evaluate_grad(NUMELQ, NUMELTG, IXQ, IXTG, X, grid, pt, dg_info_ref, dgc_coeff, i_cell)
      implicit none 
      integer :: NUMELQ, NUMELTG
      integer, dimension(:,:) :: IXQ, IXTG
      real(kind=wp), dimension(:,:) :: X
      type(grid2D_struct_multicutcell), dimension(:, :) :: grid
      type(Point2D) :: pt
      type(DGCell_GaussLegendre_ref) :: dg_info_ref
      type(DG_coeffs), dimension(:) :: dgc_coeff 
      integer:: i_cell
      type(Point2D) :: evaluate_grad

      integer :: j
      type(Point2D), dimension(4) :: normals
      type(Point2D) :: pt_grid, grad_ref
      real(kind=wp) :: d, det
      integer :: nb_normals
      logical:: is_inside
      type(Point2D), dimension(4) :: vertices
      real(wp) :: xi, eta
      real(wp), parameter :: tol = 1e-10
      integer, parameter :: max_iter = 20
      integer :: ierr
      real(wp) :: bx, by 
      real(wp) :: cx, cy 
      real(wp) :: dx, dy 
      real(wp) :: j11, j12, j21, j22
      integer(kind=8) :: k

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
        grad_ref = evaluate_grad_in_ref(xi, eta, dg_info_ref, dgc_coeff(i_cell))

        bx = vertices(2)%y - vertices(1)%y
        by = vertices(2)%z - vertices(1)%z
        cx = vertices(4)%y - vertices(1)%y
        cy = vertices(4)%z - vertices(1)%z 
        dx = vertices(1)%y - vertices(2)%y - vertices(4)%y + vertices(3)%y 
        dy = vertices(1)%z - vertices(2)%z - vertices(4)%z + vertices(3)%z 

        j11 = bx + dx*eta 
        j12 = cx + dx*xi 
        j21 = by + dy*eta 
        j22 = cy + dy*xi 

        det = j11*j22 - j12*j21 

        evaluate_grad%y = (j22*grad_ref%y - j21*grad_ref%z)/det
        evaluate_grad%z = (-j12*grad_ref%y + j22*grad_ref%z)/det
      else
        evaluate_grad%y = 0.0
        evaluate_grad%z = 0.0
      end if
    end function evaluate_grad

    !! \brief Computes the interface velocity at the 0 level of the level-set function phi and extends it
    !!        in a bandwidth around the 0 level using a fast marching tiles algorithm.
    subroutine extend_velocity_2d(NUMELQ, NUMELTG, IXQ, IXTG, X, ALE_CONNECT, grid, &
                               phi, dg_info_ref, middle_pts, phi_vals, &
                               rho, vely, velz, p, gamma,&
                               in_band, F)
      use heap_module
      use riemann_solver_mod, only : solve_riemann_problem
      use ALE_CONNECTIVITY_MOD
      use geometric_rebuilder_mod

      implicit none 

      integer :: NUMELQ, NUMELTG
      integer, dimension(:,:) :: IXQ, IXTG
      real(kind=wp), dimension(:,:) :: X
      TYPE(t_ale_connectivity), INTENT(IN) :: ALE_CONNECT
      type(grid2D_struct_multicutcell), dimension(:, :) :: grid
      type(DG_coeffs), dimension(:), intent(in) :: phi
      type(DGCell_GaussLegendre_ref), intent(in) :: dg_info_ref
      type(Point2D), dimension(:), intent(in) :: middle_pts
      real(kind=wp), dimension(:) :: phi_vals
      real(kind=wp), dimension(:,:), intent(in) :: vely
      real(kind=wp), dimension(:,:), intent(in) :: velz
      real(kind=wp), dimension(:,:), intent(in) :: rho
      real(kind=wp), dimension(:,:), intent(in) :: p
      real(kind=wp), dimension(:), intent(in) :: gamma
      logical, dimension(:), intent(in) :: in_band
      type(Point2D), dimension(:), intent(out) :: F

      !Locaal variables
      integer :: nb_cell, i, e, k1, k2, nb, nb2
      real(wp), dimension(:), allocatable :: T
      integer, dimension(:), allocatable :: state
      logical :: is_interface
      real(wp) :: phi_middle, t_new
      type(heap_t) :: pq
      type(Point2D), dimension(2) :: normalVecEdge
      integer :: nb_normals, ci
      type(Point2D) :: tang_edge, F_star
      integer, parameter :: wave_type = 2
      real(wp) :: us, vsL, vsR, ps, f_min, a, b, c, disc
      real(wp) :: T1, h1, T2, h2, h, t_cand, normVec
      integer, dimension(4) :: neighs, nbs
      real(wp), dimension(4) :: frozen1, frozen2, phi_corner
      integer :: size_frozen, ips
      type(Point2D) :: vi, vn

      nb_cell = size(vely, 1)

      allocate(T(nb_cell))
      allocate(state(nb_cell))

      T       = INF
      state   = FAR
      F       = Point2D(0.0, 0.0)
      call heap_init(pq, nb_cell)

      ! Seed from zero-crossings
      do i = 1,nb_cell
        if (any(grid(i,:)%is_narrowband)) then
          is_interface = .false.
          phi_middle = phi_vals(i)

          do e = 1,4
            ips = IXQ(e+1, i)
            vi = Point2D(X(2, ips), X(3, ips))
            phi_corner(e) = evaluate(NUMELQ, NUMELTG, IXQ, IXTG, X, grid, vi, dg_info_ref, phi, i)
          end do

          call build_polygon_from_level_set(IXQ, X, i, phi_corner, normalVecEdge, nb_normals)

          F_star%y = 0.0
          F_star%z = 0.0
          do e = 1,nb_normals
            is_interface = .true.
            normVec = sqrt(normalVecEdge(e)%y*normalVecEdge(e)%y + normalVecEdge(e)%z*normalVecEdge(e)%z)
            if (normVec > 0.) then
              normalVecEdge(e)%y = normalVecEdge(e)%y / normVec
              normalVecEdge(e)%z = normalVecEdge(e)%z / normVec
              tang_edge%y = -normalVecEdge(e)%z 
              tang_edge%z =  normalVecEdge(e)%y
              call solve_riemann_problem(gamma(1), gamma(2), rho(i, 1), rho(i, 2), &
                                    vely(i, 1), vely(i, 2), velz(i, 1), velz(i, 2), p(i, 1), p(i, 2), wave_type, &
                                    normalVecEdge(e)%y, normalVecEdge(e)%z, us, vsL, vsR, ps)
            else
              us = 0.0 
              vsL = 0.0 
              vsR = 0.0
              tang_edge%y = 0.0
              tang_edge%z = 0.0
            end if
            F_star%y = F_star%y + us * normalVecEdge(e)%y + (rho(i,1)*vsL + rho(i,2)*vsR)/(rho(i,1) + rho(i,2)) * tang_edge%y
            F_star%z = F_star%z + us * normalVecEdge(e)%z + (rho(i,1)*vsL + rho(i,2)*vsR)/(rho(i,1) + rho(i,2)) * tang_edge%z
          end do

          if (is_interface) then
            T(i)     = abs(phi_middle)
            F(i)%y   = F_star%y
            F(i)%z   = F_star%z
            state(i) = FROZEN
            call heap_insert(pq, i, abs(phi_middle))
          end if
        end if
      end do

      ! FMM propagation restricted to the narrow band 
      do while (.not. (heap_is_empty(pq)))
        call heap_extract_min(pq, ci, f_min)
        state(ci) = FROZEN

        call all_adjacency_face(ALE_CONNECT, ci, neighs)
        do i=1,4
          nb = neighs(i)
          if ((nb < 1) .or. (state(nb) == FROZEN) .or. (.not. (in_band(nb)))) then
            continue
          end if

          t_new  = INF
          call all_adjacency_face(ALE_CONNECT, nb, nbs)
          vi = middle_pts(nb)

          size_frozen = 0
          do k1=1,4
            nb2 = nbs(k1)
            if ((nb2 > 0) .and. (T(nb2) < INF)) then
              vn = middle_pts(nb2)
              h  = sqrt((vi%y - vn%y)**2 + (vi%z - vn%z)**2)
              size_frozen = size_frozen + 1
              frozen1(size_frozen) = T(nb2)
              frozen2(size_frozen) = h
              t_new = min(t_new, T(nb2) + h)
            end if
          end do

          !2-D quadratic update
          do k1=1,size_frozen
            do k2 =k1 + 1,size_frozen
              T1 = frozen1(k1)
              h1 = frozen2(k1)
              T2 = frozen1(k2)
              h2 = frozen2(k2)
              a    =  1.0/h1**2 + 1.0/h2**2
              b    = -2.0*(T1/h1**2 + T2/h2**2)
              c    =  T1**2/h1**2 + T2**2/h2**2 - 1.0
              disc = b**2 - 4.0*a*c
              if (disc >= 0.0) then
                t_cand = (-b + sqrt(disc)) / (2.0 * a)
                if (t_cand >= max(T1, T2)) then
                  t_new = min(t_new, t_cand)
                end if
              end if
            end do
          end do

          if (t_new < T(nb)) then
              T(nb)     = t_new
              F(nb)     = F(ci)
              state(nb) = NARROW
              call heap_insert(pq, nb, t_new)
          end if
        end do
      end do

      deallocate(T)
      deallocate(state)
      call heap_destroy(pq)

    end subroutine extend_velocity_2d
  end subroutine solve_levelset_projection


end module level_set_mod