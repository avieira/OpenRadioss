#include "multicutcell_geometric_rebuilder.h"
#include "array_Points.h"
#include <math.h>

void vertex_neighbors(int i, int *neigh1, int *neigh2, int *e_neigh1, int *e_neigh2){
    if (i==0){
        *neigh1 = 3;
        *neigh2 = 1;
        *e_neigh1 = 3;
        *e_neigh2 = 0;
    } else if (i==1){
        *neigh1 = 0;
        *neigh2 = 2;
        *e_neigh1 = 0;
        *e_neigh2 = 1;
    } else if (i==2){
        *neigh1 = 1;
        *neigh2 = 3;
        *e_neigh1 = 1;
        *e_neigh2 = 2;
    } else if (i==3){
        *neigh1 = 2;
        *neigh2 = 0;
        *e_neigh1 = 2;
        *e_neigh2 = 3;
    }
}

/// @brief Builds the Polyhedron3D corresponding to a cell.
///         Based on the number of vertices (in space-time) with negative level-set value,
///         the polyhedron will reconstruct the first phase (positive part of the level-set, is_reversed will be true)
//          or the second phase (negative part, is_revsered will be false).
void rebuild_polyhedron_level_set(const Polygon2D* grid, Polyhedron3D* clipped3D, my_real_c dt,\
                                const my_real_c* level_set_tn, const my_real_c* level_set_tnp1, \
                                const long long int nb_tn, const long long int nb_tnp1, long long int *is_reversed){
    if ((nb_tn +  nb_tnp1 == 0) || (nb_tn + nb_tnp1 == 8)){
        *is_reversed = (nb_tn +  nb_tnp1 == 0);
    } else if ((nb_tn + nb_tnp1 == 1) || (nb_tn + nb_tnp1 == 7)){
        *is_reversed = (nb_tn +  nb_tnp1 == 1);
        polygon_from_level_set_1_pt(grid, clipped3D, dt, level_set_tn, level_set_tnp1, nb_tn, nb_tnp1);
    } else if ((nb_tn + nb_tnp1 == 2) || (nb_tn + nb_tnp1 == 6)) {
        *is_reversed = (nb_tn +  nb_tnp1 == 2);
        polygon_from_level_set_2_pts(grid, clipped3D, dt, level_set_tn, level_set_tnp1, nb_tn, nb_tnp1);
    } else if ((nb_tn + nb_tnp1 == 3) || (nb_tn + nb_tnp1 == 5)) {
        *is_reversed = (nb_tn +  nb_tnp1 == 3);
        polygon_from_level_set_3_pts(grid, clipped3D, dt, level_set_tn, level_set_tnp1, nb_tn, nb_tnp1);
    } else if (nb_tn + nb_tnp1 == 4) {
        *is_reversed = true;
        polygon_from_level_set_4_pts(grid, clipped3D, dt, level_set_tn, level_set_tnp1, nb_tn, nb_tnp1);
    }
}

/// @brief Builds the Polyhedron3D corresponding to a cell where exactly one vertex
///        (in space-time) is on the "inside" side of the level-set, at t^n or t^{n+1}.
void polygon_from_level_set_1_pt(const Polygon2D* grid, Polyhedron3D* clipped3D, my_real_c dt,
                                  const my_real_c* level_set_tn, const my_real_c* level_set_tnp1,
                                  const long long int nb_tn, const long long int nb_tnp1){
    GrB_Info infogrb;
    GrB_Matrix *edges   = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    GrB_Matrix *faces   = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    GrB_Matrix *volumes = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    Vector_points3D *vertices;
    Vector_int *status_face;
    Polyhedron3D *built;

    int i_pt, i_b, i_a, e_b, e_a;
    my_real_c lambda_b, lambda_a, lambda_t;
    Point2D *pt, *ptb, *pta;
    Point3D v1, v2, v3, v4;
    long int s1, s2, s3, s4;
    int8_t face_sign; // +1 if the faces below are correctly oriented, -1 if they must be flipped

    vertices = alloc_with_capacity_vec_pts3D(4);
    infogrb = GrB_Matrix_new(edges,   GrB_INT8, 4, 6);
    infogrb = GrB_Matrix_new(faces,   GrB_INT8, 6, 4);
    infogrb = GrB_Matrix_new(volumes, GrB_INT8, 4, 1);
    status_face = alloc_with_capacity_vec_int(4);

    if (nb_tn == 1 || nb_tn == 3){
        // Isolated vertex at t^n
        i_pt = 0;
        if (nb_tn == 1){
            while (level_set_tn[i_pt] > 0) i_pt++;
        } else {
            while (level_set_tn[i_pt] <= 0) i_pt++;
        }

        vertex_neighbors(i_pt, &i_b, &i_a, &e_b, &e_a);

        lambda_b = level_set_tn[i_pt] / (level_set_tn[i_pt] - level_set_tn[i_b]);
        lambda_a = level_set_tn[i_pt] / (level_set_tn[i_pt] - level_set_tn[i_a]);
        lambda_t = level_set_tn[i_pt] / (level_set_tn[i_pt] - level_set_tnp1[i_pt]);

        pt  = get_ith_elem_vec_pts2D(grid->vertices, i_pt);
        ptb = get_ith_elem_vec_pts2D(grid->vertices, i_b);
        pta = get_ith_elem_vec_pts2D(grid->vertices, i_a);

        v1 = (Point3D){pt->x + lambda_b*(ptb->x - pt->x), pt->y + lambda_b*(ptb->y - pt->y), 0};
        v2 = (Point3D){pt->x, pt->y, 0};
        v3 = (Point3D){pt->x + lambda_a*(pta->x - pt->x), pt->y + lambda_a*(pta->y - pt->y), 0};
        v4 = (Point3D){pt->x, pt->y, lambda_t * dt};

        s1 = 1;          // status_face[1] = 1
        face_sign = (int8_t) -1;
    } else { // nb_tnp1 == 1 || nb_tnp1 == 3
        i_pt = 0;
        if (nb_tnp1 == 1){
            while (level_set_tnp1[i_pt] > 0) i_pt++;
        } else {
            while (level_set_tnp1[i_pt] <= 0) i_pt++;
        }

        vertex_neighbors(i_pt, &i_b, &i_a, &e_b, &e_a);

        lambda_b = level_set_tnp1[i_pt] / (level_set_tnp1[i_pt] - level_set_tnp1[i_b]);
        lambda_a = level_set_tnp1[i_pt] / (level_set_tnp1[i_pt] - level_set_tnp1[i_a]);
        lambda_t = level_set_tnp1[i_pt] / (level_set_tnp1[i_pt] - level_set_tn[i_pt]);

        pt  = get_ith_elem_vec_pts2D(grid->vertices, i_pt);
        ptb = get_ith_elem_vec_pts2D(grid->vertices, i_b);
        pta = get_ith_elem_vec_pts2D(grid->vertices, i_a);

        v1 = (Point3D){pt->x + lambda_b*(ptb->x - pt->x), pt->y + lambda_b*(ptb->y - pt->y), dt};
        v2 = (Point3D){pt->x, pt->y, dt};
        v3 = (Point3D){pt->x + lambda_a*(pta->x - pt->x), pt->y + lambda_a*(pta->y - pt->y), dt};
        v4 = (Point3D){pt->x, pt->y, (1 - lambda_t) * dt};

        s1 = 2;         // status_face[1] = 2
        face_sign = (int8_t) 1;
    }

    push_back_vec_pts3D(&vertices, &v1);
    push_back_vec_pts3D(&vertices, &v2);
    push_back_vec_pts3D(&vertices, &v3);
    push_back_vec_pts3D(&vertices, &v4);

    // edges[vertex, edge] = +-1  (4x6, 0-based indices)
    infogrb = GrB_Matrix_setElement(*edges, -1, 0, 0); infogrb = GrB_Matrix_setElement(*edges,  1, 1, 0);
    infogrb = GrB_Matrix_setElement(*edges, -1, 1, 1); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 1);
    infogrb = GrB_Matrix_setElement(*edges, -1, 0, 2); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 2);
    infogrb = GrB_Matrix_setElement(*edges, -1, 0, 3); infogrb = GrB_Matrix_setElement(*edges,  1, 3, 3);
    infogrb = GrB_Matrix_setElement(*edges, -1, 1, 4); infogrb = GrB_Matrix_setElement(*edges,  1, 3, 4);
    infogrb = GrB_Matrix_setElement(*edges, -1, 2, 5); infogrb = GrB_Matrix_setElement(*edges,  1, 3, 5);

    // faces[edge, face] = +-1 * face_sign  (6x4, 0-based indices)
    infogrb = GrB_Matrix_setElement(*faces, -face_sign, 0, 0);
    infogrb = GrB_Matrix_setElement(*faces, -face_sign, 1, 0);
    infogrb = GrB_Matrix_setElement(*faces,  face_sign, 2, 0);

    infogrb = GrB_Matrix_setElement(*faces,  face_sign, 0, 1);
    infogrb = GrB_Matrix_setElement(*faces,  face_sign, 4, 1);
    infogrb = GrB_Matrix_setElement(*faces, -face_sign, 3, 1);

    infogrb = GrB_Matrix_setElement(*faces,  face_sign, 1, 2);
    infogrb = GrB_Matrix_setElement(*faces,  face_sign, 5, 2);
    infogrb = GrB_Matrix_setElement(*faces, -face_sign, 4, 2);

    infogrb = GrB_Matrix_setElement(*faces, -face_sign, 2, 3);
    infogrb = GrB_Matrix_setElement(*faces, -face_sign, 5, 3);
    infogrb = GrB_Matrix_setElement(*faces,  face_sign, 3, 3);

    // status_face[2]/[3]: status of edges e_b / e_a of the reference cell
    s2 = *get_ith_elem_vec_int(grid->status_edge, e_b) + 2;
    s3 = *get_ith_elem_vec_int(grid->status_edge, e_a) + 2;
    s4 = -1;

    push_back_vec_int(&status_face, &s1);
    push_back_vec_int(&status_face, &s2);
    push_back_vec_int(&status_face, &s3);
    push_back_vec_int(&status_face, &s4);

    // volumes .= 1
    infogrb = GrB_Matrix_setElement(*volumes, 1, 0, 0);
    infogrb = GrB_Matrix_setElement(*volumes, 1, 1, 0);
    infogrb = GrB_Matrix_setElement(*volumes, 1, 2, 0);
    infogrb = GrB_Matrix_setElement(*volumes, 1, 3, 0);

    built = new_Polyhedron3D_vefvs(vertices, edges, faces, volumes, status_face);
    copy_Polyhedron3D(built, clipped3D);

    dealloc_Polyhedron3D(built); free(built);
    dealloc_vec_pts3D(vertices); free(vertices);
    dealloc_vec_int(status_face); free(status_face);
    GrB_free(edges);   free(edges);
    GrB_free(faces);   free(faces);
    GrB_free(volumes); free(volumes);
}

/// @brief Linear interpolation between two space-time points: (1-lambda)*p0 + lambda*p1.
static Point3D lerp_pt3D(my_real_c lambda, Point3D p0, Point3D p1){
    return (Point3D){ (1 - lambda) * p0.x + lambda * p1.x,
                       (1 - lambda) * p0.y + lambda * p1.y,
                       (1 - lambda) * p0.t + lambda * p1.t };
}


/// @brief Builds the Polyhedron3D corresponding to a cell where exactly two vertices
///        (in space-time) are on the "inside" side of the level-set.
void polygon_from_level_set_2_pts(const Polygon2D* grid, Polyhedron3D* clipped3D, my_real_c dt,
                                   const my_real_c* level_set_tn, const my_real_c* level_set_tnp1,
                                   const long long int nb_tn, const long long int nb_tnp1){
    GrB_Info infogrb;
    int i, i_pt1, i_pt2, case_1_face;
    Vector_points3D *vertices;
    Vector_int *status_face;
    GrB_Matrix *edges   = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    GrB_Matrix *faces   = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    GrB_Matrix *volumes = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    Polyhedron3D *built;

    // --- Find the two relevant vertices (0-based) ---
    i_pt1 = -1; i_pt2 = -1;
    if (nb_tn == 2 && nb_tnp1 == 0){
        for (i = 0; i < 4; i++){
            if (level_set_tn[i] <= 0){ if (i_pt1 == -1) i_pt1 = i; else i_pt2 = i; }
        }
    } else if (nb_tnp1 == 4 && nb_tn == 2){
        for (i = 0; i < 4; i++){
            if (level_set_tn[i] > 0){ if (i_pt1 == -1) i_pt1 = i; else i_pt2 = i; }
        }
    } else if (nb_tnp1 == 2 && nb_tn == 0){
        for (i = 0; i < 4; i++){
            if (level_set_tnp1[i] <= 0){ if (i_pt1 == -1) i_pt1 = i; else i_pt2 = i; }
        }
    } else if (nb_tnp1 == 2 && nb_tn == 4){
        for (i = 0; i < 4; i++){
            if (level_set_tnp1[i] > 0){ if (i_pt1 == -1) i_pt1 = i; else i_pt2 = i; }
        }
    } else if (nb_tn == 1 && nb_tnp1 == 1){
        for (i = 0; i < 4; i++){ if (level_set_tn[i]   <= 0){ i_pt1 = i; break; } }
        for (i = 0; i < 4; i++){ if (level_set_tnp1[i]  <= 0){ i_pt2 = i; break; } }
    } else if (nb_tn == 3 && nb_tnp1 == 3){
        for (i = 0; i < 4; i++){ if (level_set_tn[i]   > 0){ i_pt1 = i; break; } }
        for (i = 0; i < 4; i++){ if (level_set_tnp1[i]  > 0){ i_pt2 = i; break; } }
    } else {
        printf("It should not happen: nb_tn + nb_tnp1 should equal 2 or 6, but here it equals %lld\n", nb_tn + nb_tnp1);
        return;
    }

    // --- Decide whether this is the "one internal face" case (1 volume) or the "two tetrahedra" case (2 volumes) ---
    // Vertex labels below are 0-based versions of the original 1-based adjacency test.
    case_1_face = (i_pt1 == i_pt2) ||
        ((nb_tn == 2 || nb_tn == 6 || nb_tnp1 == 2 || nb_tnp1 == 6) &&
         ((i_pt1 == i_pt2 + 1) || (i_pt1 == 3 && i_pt2 == 0) ||
          (i_pt1 == i_pt2 - 1) || (i_pt1 == 0 && i_pt2 == 3)));

    if (case_1_face){
        // In this case, there is only one volume with an internal face inside the time-space volume
        // (that we actually cut in two to approximate the inner curvature).
        // Two cases may happen: both points have negative level-set at the same time, or not.
        Point3D pt1, pt2, pta1, pta2, ptb1, ptb2;
        my_real_c lambda_a1, lambda_b1, lambda_a2, lambda_b2;
        long int s1, s2, s3, s4, s5, s6;
        int8_t flip_faces;

        vertices    = alloc_with_capacity_vec_pts3D(6);
        infogrb = GrB_Matrix_new(edges,   GrB_INT8, 6, 10);
        infogrb = GrB_Matrix_new(faces,   GrB_INT8, 10, 6);
        infogrb = GrB_Matrix_new(volumes, GrB_INT8, 6, 1);
        status_face = alloc_with_capacity_vec_int(6);

        s3 = -1; s4 = -1;

        if (i_pt1 == i_pt2){
            int i_b1, i_a1, i_b2, i_a2;
            int neigh1, neigh2, e_neigh1, e_neigh2;
            Point2D *pt2D, *pta2D, *ptb2D;

            i_b1 = (i_pt1 == 0 ? 3 : i_pt1 - 1);
            i_a1 = (i_pt1 == 3 ? 0 : i_pt1 + 1);
            i_b2 = i_b1;
            i_a2 = i_a1;

            pt2D  = get_ith_elem_vec_pts2D(grid->vertices, i_pt1);
            pta2D = get_ith_elem_vec_pts2D(grid->vertices, i_a1);
            ptb2D = get_ith_elem_vec_pts2D(grid->vertices, i_b1);
            pt1  = (Point3D){pt2D->x, pt2D->y, 0.0};
            pt2  = (Point3D){pt2D->x, pt2D->y, dt};
            pta1 = (Point3D){pta2D->x, pta2D->y, 0.0};
            pta2 = (Point3D){pta2D->x, pta2D->y, dt};
            ptb1 = (Point3D){ptb2D->x, ptb2D->y, 0.0};
            ptb2 = (Point3D){ptb2D->x, ptb2D->y, dt};

            vertex_neighbors(i_pt1, &neigh1, &neigh2, &e_neigh1, &e_neigh2);
            s1 = 1;
            s2 = 2;

            s5 = 2 + *get_ith_elem_vec_int(grid->status_edge, e_neigh1);
            s6 = 2 + *get_ith_elem_vec_int(grid->status_edge, e_neigh2);

            lambda_a1 = level_set_tn[i_pt1]   / (level_set_tn[i_pt1]   - level_set_tn[i_a1]);
            lambda_b1 = level_set_tn[i_pt1]   / (level_set_tn[i_pt1]   - level_set_tn[i_b1]);
            lambda_a2 = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tnp1[i_a2]);
            lambda_b2 = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tnp1[i_b2]);
        } else {
            int i_b1, i_b2, i_a1, i_a2;
            int neigh1, neigh2, e_neigh1, e_neigh2;
            Point2D *pt2D_1, *pt2D_2, *pta2D, *ptb2D;

            i_a1 = i_pt1;
            i_a2 = i_pt2;
            if (i_pt1 == 0){
                i_b1 = (i_pt2 == 1 ? 3 : 1);
                i_b2 = 2;
            } else if (i_pt1 == 1){
                i_b1 = (i_pt2 == 2 ? 0 : 2);
                i_b2 = 3;
            } else if (i_pt1 == 2){
                i_b1 = (i_pt2 == 1 ? 3 : 1);
                i_b2 = 0;
            } else { // i_pt1 == 3
                i_b1 = (i_pt2 == 2 ? 0 : 2);
                i_b2 = 1;
            }

            pt2D_1 = get_ith_elem_vec_pts2D(grid->vertices, i_pt1);
            pt2D_2 = get_ith_elem_vec_pts2D(grid->vertices, i_pt2);
            pta2D  = get_ith_elem_vec_pts2D(grid->vertices, i_b1);
            ptb2D  = get_ith_elem_vec_pts2D(grid->vertices, i_b2);

            if (nb_tn == 2){
                pt1  = (Point3D){pt2D_1->x, pt2D_1->y, 0.0};
                pt2  = (Point3D){pt2D_2->x, pt2D_2->y, 0.0};
                pta1 = (Point3D){pt2D_1->x, pt2D_1->y, dt};
                pta2 = (Point3D){pt2D_2->x, pt2D_2->y, dt};
                ptb1 = (Point3D){pta2D->x, pta2D->y, 0.0};
                ptb2 = (Point3D){ptb2D->x, ptb2D->y, 0.0};

                s5 = 1;
                lambda_a1 = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tnp1[i_a1]);
                lambda_b1 = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tn[i_b1]);
                lambda_a2 = level_set_tn[i_pt2] / (level_set_tn[i_pt2] - level_set_tnp1[i_a2]);
                lambda_b2 = level_set_tn[i_pt2] / (level_set_tn[i_pt2] - level_set_tn[i_b2]);
            } else {
                pt1  = (Point3D){pt2D_1->x, pt2D_1->y, dt};
                pt2  = (Point3D){pt2D_2->x, pt2D_2->y, dt};
                pta1 = (Point3D){pt2D_1->x, pt2D_1->y, 0.0};
                pta2 = (Point3D){pt2D_2->x, pt2D_2->y, 0.0};
                ptb1 = (Point3D){pta2D->x, pta2D->y, dt};
                ptb2 = (Point3D){ptb2D->x, ptb2D->y, dt};

                s5 = 2;
                lambda_a1 = level_set_tnp1[i_pt1] / (level_set_tnp1[i_pt1] - level_set_tn[i_a1]);
                lambda_b1 = level_set_tnp1[i_pt1] / (level_set_tnp1[i_pt1] - level_set_tnp1[i_b1]);
                lambda_a2 = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tn[i_a2]);
                lambda_b2 = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tnp1[i_b2]);
            }

            vertex_neighbors(i_pt1, &neigh1, &neigh2, &e_neigh1, &e_neigh2);
            if (neigh1 == i_pt2){
                s1 = 2 + *get_ith_elem_vec_int(grid->status_edge, e_neigh2);
                s6 = 2 + *get_ith_elem_vec_int(grid->status_edge, e_neigh1);
            } else {
                s1 = 2 + *get_ith_elem_vec_int(grid->status_edge, e_neigh1);
                s6 = 2 + *get_ith_elem_vec_int(grid->status_edge, e_neigh2);
            }
            vertex_neighbors(i_pt2, &neigh1, &neigh2, &e_neigh1, &e_neigh2);
            if (neigh1 == i_pt1){
                s2 = 2 + *get_ith_elem_vec_int(grid->status_edge, e_neigh2);
            } else {
                s2 = 2 + *get_ith_elem_vec_int(grid->status_edge, e_neigh1);
            }
        }

        // net sign to apply to the "faces" values listed below, folding in the
        // conditional flip and the final unconditional negation 
        flip_faces = (int8_t)((i_pt1 == 0 && i_pt2 == 3) || (i_pt2 == 0 && i_pt1 == 3));
        flip_faces = ((nb_tn == 2 && flip_faces == 0) || (nb_tnp1 == 2 && flip_faces == 1)) ? (int8_t) 1 : (int8_t) -1;

        vertices->size = 0; // (alloc_with_capacity leaves size at 0; push_back below fills it)
        push_back_vec_pts3D(&vertices, &pt1);
        push_back_vec_pts3D(&vertices, &pt2);
        { Point3D v = lerp_pt3D(lambda_a1, pt1, pta1); push_back_vec_pts3D(&vertices, &v); }
        { Point3D v = lerp_pt3D(lambda_a2, pt2, pta2); push_back_vec_pts3D(&vertices, &v); }
        { Point3D v = lerp_pt3D(lambda_b1, pt1, ptb1); push_back_vec_pts3D(&vertices, &v); }
        { Point3D v = lerp_pt3D(lambda_b2, pt2, ptb2); push_back_vec_pts3D(&vertices, &v); }

        // edges[vertex, edge] = +-1  (6x10, 0-based indices)
        infogrb = GrB_Matrix_setElement(*edges, -1, 0, 0); infogrb = GrB_Matrix_setElement(*edges,  1, 1, 0);
        infogrb = GrB_Matrix_setElement(*edges, -1, 0, 1); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 1);
        infogrb = GrB_Matrix_setElement(*edges, -1, 0, 2); infogrb = GrB_Matrix_setElement(*edges,  1, 4, 2);
        infogrb = GrB_Matrix_setElement(*edges, -1, 1, 3); infogrb = GrB_Matrix_setElement(*edges,  1, 3, 3);
        infogrb = GrB_Matrix_setElement(*edges, -1, 1, 4); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 4);
        infogrb = GrB_Matrix_setElement(*edges, -1, 2, 5); infogrb = GrB_Matrix_setElement(*edges,  1, 4, 5);
        infogrb = GrB_Matrix_setElement(*edges, -1, 3, 6); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 6);
        infogrb = GrB_Matrix_setElement(*edges, -1, 2, 7); infogrb = GrB_Matrix_setElement(*edges,  1, 3, 7);
        infogrb = GrB_Matrix_setElement(*edges, -1, 4, 8); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 8);
        infogrb = GrB_Matrix_setElement(*edges, -1, 3, 9); infogrb = GrB_Matrix_setElement(*edges,  1, 4, 9);

        // faces[edge, face] = +-1 * flip_faces  (10x6, 0-based indices)
        infogrb = GrB_Matrix_setElement(*faces, -flip_faces, 1, 0); infogrb = GrB_Matrix_setElement(*faces, -flip_faces, 5, 0); infogrb = GrB_Matrix_setElement(*faces,  flip_faces, 2, 0);
        infogrb = GrB_Matrix_setElement(*faces,  flip_faces, 3, 1); infogrb = GrB_Matrix_setElement(*faces,  flip_faces, 6, 1); infogrb = GrB_Matrix_setElement(*faces, -flip_faces, 4, 1);
        infogrb = GrB_Matrix_setElement(*faces, -flip_faces, 7, 2); infogrb = GrB_Matrix_setElement(*faces,  flip_faces, 5, 2); infogrb = GrB_Matrix_setElement(*faces, -flip_faces, 9, 2);
        infogrb = GrB_Matrix_setElement(*faces,  flip_faces, 8, 3); infogrb = GrB_Matrix_setElement(*faces, -flip_faces, 6, 3); infogrb = GrB_Matrix_setElement(*faces,  flip_faces, 9, 3);
        infogrb = GrB_Matrix_setElement(*faces, -flip_faces, 2, 4); infogrb = GrB_Matrix_setElement(*faces, -flip_faces, 8, 4); infogrb = GrB_Matrix_setElement(*faces,  flip_faces, 4, 4); infogrb = GrB_Matrix_setElement(*faces,  flip_faces, 0, 4);
        infogrb = GrB_Matrix_setElement(*faces, -flip_faces, 0, 5); infogrb = GrB_Matrix_setElement(*faces,  flip_faces, 1, 5); infogrb = GrB_Matrix_setElement(*faces, -flip_faces, 3, 5); infogrb = GrB_Matrix_setElement(*faces,  flip_faces, 7, 5);

        push_back_vec_int(&status_face, &s1);
        push_back_vec_int(&status_face, &s2);
        push_back_vec_int(&status_face, &s3);
        push_back_vec_int(&status_face, &s4);
        push_back_vec_int(&status_face, &s5);
        push_back_vec_int(&status_face, &s6);

        // volumes .= 1  (dense 6x1 fill)
        for (i = 0; i < 6; i++) infogrb = GrB_Matrix_setElement(*volumes, 1, i, 0);

    } else {
        // In this case, there are two volumes: two tetrahedra pointing towards two corners of the time-space cube.
        Point3D pt1, pt2, pta1, pta2, ptb, ptc;
        my_real_c lambda_a1, lambda_b1, lambda_c1, lambda_a2, lambda_b2, lambda_c2;
        long int s1, s2, s3, s4, s5, s6, s7, s8;
        int neigh1_1, neigh1_2, e1_1, e1_2;
        int neigh2_1, neigh2_2, e2_1, e2_2;
        int e_ptb1, e_ptc1, e_ptb2, e_ptc2;
        int8_t edge_sign;
        int face_sign_1to4, face_sign_5to8;
        Point2D *pt12D, *pt22D;

        vertices    = alloc_with_capacity_vec_pts3D(8);
        infogrb = GrB_Matrix_new(edges,   GrB_INT8, 8, 12);
        infogrb = GrB_Matrix_new(faces,   GrB_INT8, 12, 8);
        infogrb = GrB_Matrix_new(volumes, GrB_INT8, 8, 2);
        status_face = alloc_with_capacity_vec_int(8);

        vertex_neighbors(i_pt1, &neigh1_1, &neigh1_2, &e1_1, &e1_2);
        vertex_neighbors(i_pt2, &neigh2_1, &neigh2_2, &e2_1, &e2_2);

        pt12D = get_ith_elem_vec_pts2D(grid->vertices, i_pt1);
        pt22D = get_ith_elem_vec_pts2D(grid->vertices, i_pt2);

        if ((i_pt1 == i_pt2 + 1) || (i_pt1 == 3 && i_pt2 == 0) ||
            (i_pt1 == i_pt2 - 1) || (i_pt1 == 0 && i_pt2 == 3)){
            // One point at each time, diagonal on the same space-time face.
            int i1_ptb, i1_ptc, i2_ptb, i2_ptc;
            Point2D *pb2D, *pc2D;

            infogrb = GrB_Matrix_extractElement(&edge_sign, *(grid->edges), neigh1_1, e1_1);
            if (edge_sign < 0){
                i1_ptb = neigh1_1; i1_ptc = neigh1_2; e_ptb1 = e1_1; e_ptc1 = e1_2;
            } else {
                i1_ptc = neigh1_1; i1_ptb = neigh1_2; e_ptc1 = e1_1; e_ptb1 = e1_2;
            }

            lambda_a1 = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tnp1[i_pt1]);
            lambda_b1 = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tn[i1_ptb]);
            lambda_c1 = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tn[i1_ptc]);

            pb2D = get_ith_elem_vec_pts2D(grid->vertices, i1_ptb);
            ptb  = (Point3D){pb2D->x, pb2D->y, 0.0};
            pc2D = get_ith_elem_vec_pts2D(grid->vertices, i1_ptc);
            ptc  = (Point3D){pc2D->x, pc2D->y, 0.0};
            pt1  = (Point3D){pt12D->x, pt12D->y, 0.0};
            pta1 = (Point3D){pt12D->x, pt12D->y, dt};

            vertices    = alloc_with_capacity_vec_pts3D(8); // (re-declared above; kept here for clarity of order)
            push_back_vec_pts3D(&vertices, &pt1);
            { Point3D v = lerp_pt3D(lambda_a1, pt1, pta1); push_back_vec_pts3D(&vertices, &v); }
            { Point3D v = lerp_pt3D(lambda_b1, pt1, ptb);  push_back_vec_pts3D(&vertices, &v); }
            { Point3D v = lerp_pt3D(lambda_c1, pt1, ptc);  push_back_vec_pts3D(&vertices, &v); }

            infogrb = GrB_Matrix_extractElement(&edge_sign, *(grid->edges), neigh2_1, e2_1);
            if (edge_sign < 0){
                i2_ptc = neigh2_1; i2_ptb = neigh2_2; e_ptc2 = e2_1; e_ptb2 = e2_2;
            } else {
                i2_ptb = neigh2_1; i2_ptc = neigh2_2; e_ptb2 = e2_1; e_ptc2 = e2_2;
            }

            pb2D = get_ith_elem_vec_pts2D(grid->vertices, i2_ptb);
            ptb  = (Point3D){pb2D->x, pb2D->y, dt};
            pc2D = get_ith_elem_vec_pts2D(grid->vertices, i2_ptc);
            ptc  = (Point3D){pc2D->x, pc2D->y, dt};
            pt2  = (Point3D){pt22D->x, pt22D->y, dt};
            pta2 = (Point3D){pt22D->x, pt22D->y, 0.0};

            lambda_a2 = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tn[i_pt2]);
            lambda_b2 = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tnp1[i2_ptb]);
            lambda_c2 = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tnp1[i2_ptc]);

            push_back_vec_pts3D(&vertices, &pt2);
            { Point3D v = lerp_pt3D(lambda_a2, pt2, pta2); push_back_vec_pts3D(&vertices, &v); }
            { Point3D v = lerp_pt3D(lambda_b2, pt2, ptb);  push_back_vec_pts3D(&vertices, &v); }
            { Point3D v = lerp_pt3D(lambda_c2, pt2, ptc);  push_back_vec_pts3D(&vertices, &v); }

        } else {
            // Opposite (diagonal) points of the reference cell.
            int i_ptb, i_ptc;
            Point2D *pb2D, *pc2D;

            if (neigh2_1 != neigh1_1){
                int tmp = e2_1; e2_1 = e2_2; e2_2 = tmp;
            }
            infogrb = GrB_Matrix_extractElement(&edge_sign, *(grid->edges), neigh1_1, e1_1);
            if (edge_sign < 0){
                i_ptb = neigh1_1; i_ptc = neigh1_2;
                e_ptb1 = e1_1; e_ptc1 = e1_2;
                e_ptb2 = e2_1; e_ptc2 = e2_2;
            } else {
                i_ptc = neigh1_1; i_ptb = neigh1_2;
                e_ptb1 = e1_2; e_ptc1 = e1_1;
                e_ptb2 = e2_2; e_ptc2 = e2_1;
            }

            vertices = alloc_with_capacity_vec_pts3D(8);

            if (nb_tn == 2){
                pb2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptb);
                pc2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptc);
                ptb = (Point3D){pb2D->x, pb2D->y, 0.0};
                ptc = (Point3D){pc2D->x, pc2D->y, 0.0};

                pt1  = (Point3D){pt12D->x, pt12D->y, 0.0};
                pta1 = (Point3D){pt12D->x, pt12D->y, dt};
                pt2  = (Point3D){pt22D->x, pt22D->y, 0.0};
                pta2 = (Point3D){pt22D->x, pt22D->y, dt};

                lambda_a1 = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tnp1[i_pt1]);
                lambda_a2 = level_set_tn[i_pt2] / (level_set_tn[i_pt2] - level_set_tnp1[i_pt2]);
                lambda_b1 = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tn[i_ptb]);
                lambda_b2 = level_set_tn[i_pt2] / (level_set_tn[i_pt2] - level_set_tn[i_ptb]);
                lambda_c1 = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tn[i_ptc]);
                lambda_c2 = level_set_tn[i_pt2] / (level_set_tn[i_pt2] - level_set_tn[i_ptc]);

                push_back_vec_pts3D(&vertices, &pt1);
                { Point3D v = lerp_pt3D(lambda_a1, pt1, pta1); push_back_vec_pts3D(&vertices, &v); }
                { Point3D v = lerp_pt3D(lambda_b1, pt1, ptb);  push_back_vec_pts3D(&vertices, &v); }
                { Point3D v = lerp_pt3D(lambda_c1, pt1, ptc);  push_back_vec_pts3D(&vertices, &v); }
                push_back_vec_pts3D(&vertices, &pt2);
                { Point3D v = lerp_pt3D(lambda_a2, pt2, pta2); push_back_vec_pts3D(&vertices, &v); }
                { Point3D v = lerp_pt3D(lambda_b2, pt2, ptb);  push_back_vec_pts3D(&vertices, &v); }
                { Point3D v = lerp_pt3D(lambda_c2, pt2, ptc);  push_back_vec_pts3D(&vertices, &v); }

            } else if (nb_tnp1 == 2){
                pb2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptb);
                pc2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptc);
                ptb = (Point3D){pb2D->x, pb2D->y, dt};
                ptc = (Point3D){pc2D->x, pc2D->y, dt};

                pt1  = (Point3D){pt12D->x, pt12D->y, dt};
                pta1 = (Point3D){pt12D->x, pt12D->y, 0.0};
                pt2  = (Point3D){pt22D->x, pt22D->y, dt};
                pta2 = (Point3D){pt22D->x, pt22D->y, 0.0};

                lambda_a1 = level_set_tnp1[i_pt1] / (level_set_tnp1[i_pt1] - level_set_tn[i_pt1]);
                lambda_a2 = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tn[i_pt2]);
                lambda_b1 = level_set_tnp1[i_pt1] / (level_set_tnp1[i_pt1] - level_set_tnp1[i_ptb]);
                lambda_b2 = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tnp1[i_ptb]);
                lambda_c1 = level_set_tnp1[i_pt1] / (level_set_tnp1[i_pt1] - level_set_tnp1[i_ptc]);
                lambda_c2 = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tnp1[i_ptc]);

                push_back_vec_pts3D(&vertices, &pt1);
                { Point3D v = lerp_pt3D(lambda_a1, pt1, pta1); push_back_vec_pts3D(&vertices, &v); }
                { Point3D v = lerp_pt3D(lambda_b1, pt1, ptb);  push_back_vec_pts3D(&vertices, &v); }
                { Point3D v = lerp_pt3D(lambda_c1, pt1, ptc);  push_back_vec_pts3D(&vertices, &v); }
                push_back_vec_pts3D(&vertices, &pt2);
                { Point3D v = lerp_pt3D(lambda_a2, pt2, pta2); push_back_vec_pts3D(&vertices, &v); }
                { Point3D v = lerp_pt3D(lambda_b2, pt2, ptb);  push_back_vec_pts3D(&vertices, &v); }
                { Point3D v = lerp_pt3D(lambda_c2, pt2, ptc);  push_back_vec_pts3D(&vertices, &v); }

            } else {
                // nb_tn == 1 && nb_tnp1 == 1  (or nb_tn == 3 && nb_tnp1 == 3)
                pb2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptb);
                pc2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptc);
                ptb = (Point3D){pb2D->x, pb2D->y, 0.0};
                ptc = (Point3D){pc2D->x, pc2D->y, 0.0};
                pt1  = (Point3D){pt12D->x, pt12D->y, 0.0};
                pta1 = (Point3D){pt12D->x, pt12D->y, dt};

                lambda_a1 = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tnp1[i_pt1]);
                lambda_b1 = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tn[i_ptb]);
                lambda_c1 = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tn[i_ptc]);

                push_back_vec_pts3D(&vertices, &pt1);
                { Point3D v = lerp_pt3D(lambda_a1, pt1, pta1); push_back_vec_pts3D(&vertices, &v); }
                { Point3D v = lerp_pt3D(lambda_b1, pt1, ptb);  push_back_vec_pts3D(&vertices, &v); }
                { Point3D v = lerp_pt3D(lambda_c1, pt1, ptc);  push_back_vec_pts3D(&vertices, &v); }

                ptb.t = dt;
                ptc.t = dt;
                pt2  = (Point3D){pt22D->x, pt22D->y, dt};
                pta2 = (Point3D){pt22D->x, pt22D->y, 0.0};

                lambda_a2 = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tn[i_pt2]);
                lambda_b2 = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tnp1[i_ptb]);
                lambda_c2 = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tnp1[i_ptc]);

                push_back_vec_pts3D(&vertices, &pt2);
                { Point3D v = lerp_pt3D(lambda_a2, pt2, pta2); push_back_vec_pts3D(&vertices, &v); }
                { Point3D v = lerp_pt3D(lambda_b2, pt2, ptb);  push_back_vec_pts3D(&vertices, &v); }
                { Point3D v = lerp_pt3D(lambda_c2, pt2, ptc);  push_back_vec_pts3D(&vertices, &v); }
            }
        }

        // edges[vertex, edge] = +-1  (8x12, 0-based indices)
        infogrb = GrB_Matrix_setElement(*edges, -1, 0, 0);  infogrb = GrB_Matrix_setElement(*edges,  1, 1, 0);
        infogrb = GrB_Matrix_setElement(*edges, -1, 0, 1);  infogrb = GrB_Matrix_setElement(*edges,  1, 2, 1);
        infogrb = GrB_Matrix_setElement(*edges, -1, 0, 2);  infogrb = GrB_Matrix_setElement(*edges,  1, 3, 2);
        infogrb = GrB_Matrix_setElement(*edges, -1, 1, 3);  infogrb = GrB_Matrix_setElement(*edges,  1, 2, 3);
        infogrb = GrB_Matrix_setElement(*edges, -1, 1, 4);  infogrb = GrB_Matrix_setElement(*edges,  1, 3, 4);
        infogrb = GrB_Matrix_setElement(*edges, -1, 2, 5);  infogrb = GrB_Matrix_setElement(*edges,  1, 3, 5);
        infogrb = GrB_Matrix_setElement(*edges, -1, 4, 6);  infogrb = GrB_Matrix_setElement(*edges,  1, 5, 6);
        infogrb = GrB_Matrix_setElement(*edges, -1, 4, 7);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 7);
        infogrb = GrB_Matrix_setElement(*edges, -1, 4, 8);  infogrb = GrB_Matrix_setElement(*edges,  1, 7, 8);
        infogrb = GrB_Matrix_setElement(*edges, -1, 5, 9);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 9);
        infogrb = GrB_Matrix_setElement(*edges, -1, 5, 10); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 10);
        infogrb = GrB_Matrix_setElement(*edges, -1, 6, 11); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 11);

        // net sign for faces columns 1-4 (pt1's tetrahedron) and 5-8 (pt2's tetrahedron),
        // folding in both conditional flips and the final unconditional negation
        face_sign_1to4 = (nb_tnp1 == 2) ? 1 : -1;
        face_sign_5to8 = (nb_tnp1 == 0 || nb_tnp1 == 4) ? -1 : 1;

        // faces[edge, face] = +-1 * face_sign  (12x8, 0-based indices)
        infogrb = GrB_Matrix_setElement(*faces,  face_sign_1to4, 0, 0); infogrb = GrB_Matrix_setElement(*faces, -face_sign_1to4, 1, 0); infogrb = GrB_Matrix_setElement(*faces,  face_sign_1to4, 3, 0);
        infogrb = GrB_Matrix_setElement(*faces,  face_sign_1to4, 0, 1); infogrb = GrB_Matrix_setElement(*faces, -face_sign_1to4, 4, 1); infogrb = GrB_Matrix_setElement(*faces, -face_sign_1to4, 2, 1);
        infogrb = GrB_Matrix_setElement(*faces, -face_sign_1to4, 3, 2); infogrb = GrB_Matrix_setElement(*faces,  face_sign_1to4, 4, 2); infogrb = GrB_Matrix_setElement(*faces, -face_sign_1to4, 5, 2);
        infogrb = GrB_Matrix_setElement(*faces, -face_sign_1to4, 1, 3); infogrb = GrB_Matrix_setElement(*faces,  face_sign_1to4, 2, 3); infogrb = GrB_Matrix_setElement(*faces,  face_sign_1to4, 5, 3);
        infogrb = GrB_Matrix_setElement(*faces, -face_sign_5to8, 6, 4); infogrb = GrB_Matrix_setElement(*faces,  face_sign_5to8, 7, 4); infogrb = GrB_Matrix_setElement(*faces, -face_sign_5to8, 9, 4);
        infogrb = GrB_Matrix_setElement(*faces, -face_sign_5to8, 6, 5); infogrb = GrB_Matrix_setElement(*faces,  face_sign_5to8, 10, 5); infogrb = GrB_Matrix_setElement(*faces, -face_sign_5to8, 8, 5);
        infogrb = GrB_Matrix_setElement(*faces,  face_sign_5to8, 9, 6); infogrb = GrB_Matrix_setElement(*faces,  face_sign_5to8, 10, 6); infogrb = GrB_Matrix_setElement(*faces,  face_sign_5to8, 11, 6);
        infogrb = GrB_Matrix_setElement(*faces, -face_sign_5to8, 7, 7); infogrb = GrB_Matrix_setElement(*faces,  face_sign_5to8, 8, 7); infogrb = GrB_Matrix_setElement(*faces, -face_sign_5to8, 11, 7);

        s1 = 2 + *get_ith_elem_vec_int(grid->status_edge, e_ptb1);
        s2 = 2 + *get_ith_elem_vec_int(grid->status_edge, e_ptc1);
        s3 = -1;
        s4 = ((nb_tnp1 == 2) || (nb_tnp1 == 6)) ? 2 : 1;
        s5 = 2 + *get_ith_elem_vec_int(grid->status_edge, e_ptb2);
        s6 = 2 + *get_ith_elem_vec_int(grid->status_edge, e_ptc2);
        s7 = -1;
        s8 = ((nb_tn == 2) || (nb_tn == 6)) ? 1 : 2;

        push_back_vec_int(&status_face, &s1);
        push_back_vec_int(&status_face, &s2);
        push_back_vec_int(&status_face, &s3);
        push_back_vec_int(&status_face, &s4);
        push_back_vec_int(&status_face, &s5);
        push_back_vec_int(&status_face, &s6);
        push_back_vec_int(&status_face, &s7);
        push_back_vec_int(&status_face, &s8);

        // volumes[1:4, 1] .= 1 ; volumes[5:8, 2] .= 1  (dense fill, 0-based)
        for (i = 0; i < 4; i++) infogrb = GrB_Matrix_setElement(*volumes, 1, i,     0);
        for (i = 4; i < 8; i++) infogrb = GrB_Matrix_setElement(*volumes, 1, i,     1);
    }

    built = new_Polyhedron3D_vefvs(vertices, edges, faces, volumes, status_face);
    copy_Polyhedron3D(built, clipped3D);

    dealloc_Polyhedron3D(built); free(built);
    dealloc_vec_pts3D(vertices); free(vertices);
    dealloc_vec_int(status_face); free(status_face);
    GrB_free(edges);   free(edges);
    GrB_free(faces);   free(faces);
    GrB_free(volumes); free(volumes);
}

// --- shared helpers (add once to the file) ---

static int is_adjacent_ref_square(int a, int b){
    return (a == b + 1) || (a == 3 && b == 0) || (a == b - 1) || (a == 0 && b == 3);
}

static void find_indices_le(const my_real_c* ls, int n_expected, int* idx){
    int cnt = 0;
    for (int i = 0; i < 4 && cnt < n_expected; i++) if (ls[i] <= 0) idx[cnt++] = i;
}
static void find_indices_gt(const my_real_c* ls, int n_expected, int* idx){
    int cnt = 0;
    for (int i = 0; i < 4 && cnt < n_expected; i++) if (ls[i] > 0) idx[cnt++] = i;
}
static void find_indices_ge(const my_real_c* ls, int n_expected, int* idx){
    int cnt = 0;
    for (int i = 0; i < 4 && cnt < n_expected; i++) if (ls[i] >= 0) idx[cnt++] = i;
}

static int8_t grid_edge_sign(const Polygon2D* grid, int vertex, int edge){
    int8_t v;
    GrB_Matrix_extractElement(&v, *(grid->edges), vertex, edge);
    return v;
}

/// @brief Builds the Polyhedron3D corresponding to a cell where exactly three vertices
///        (in space-time) are on the "inside" side of the level-set.
void polygon_from_level_set_3_pts(const Polygon2D* grid, Polyhedron3D* clipped3D, my_real_c dt,
                                   const my_real_c* level_set_tn, const my_real_c* level_set_tnp1,
                                   const long long int nb_tn, const long long int nb_tnp1){
    GrB_Info infogrb;
    int i_pt1, i_pt2, i_pt3;
    int neighs_12, neighs_13, neighs_23;
    int one_face, three_faces;
    int idx[3];
    Vector_points3D *vertices;
    Vector_int *status_face;
    GrB_Matrix *edges   = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    GrB_Matrix *faces   = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    GrB_Matrix *volumes = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    Polyhedron3D *built;

    // --- find i_pt1, i_pt2, i_pt3 (0-based) ---
    if (nb_tn == 3 && nb_tnp1 == 0){
        find_indices_le(level_set_tn, 3, idx); i_pt1=idx[0]; i_pt2=idx[1]; i_pt3=idx[2];
    } else if (nb_tnp1 == 4 && nb_tn == 1){
        find_indices_gt(level_set_tn, 3, idx); i_pt1=idx[0]; i_pt2=idx[1]; i_pt3=idx[2];
    } else if (nb_tnp1 == 3 && nb_tn == 0){
        find_indices_le(level_set_tnp1, 3, idx); i_pt1=idx[0]; i_pt2=idx[1]; i_pt3=idx[2];
    } else if (nb_tnp1 == 1 && nb_tn == 4){
        find_indices_ge(level_set_tnp1, 3, idx); i_pt1=idx[0]; i_pt2=idx[1]; i_pt3=idx[2];
    } else if (nb_tn == 1 && nb_tnp1 == 2){
        find_indices_le(level_set_tn, 1, idx); i_pt1 = idx[0];
        find_indices_le(level_set_tnp1, 2, idx); i_pt2=idx[0]; i_pt3=idx[1];
    } else if (nb_tn == 3 && nb_tnp1 == 2){
        find_indices_gt(level_set_tn, 1, idx); i_pt1 = idx[0];
        find_indices_gt(level_set_tnp1, 2, idx); i_pt2=idx[0]; i_pt3=idx[1];
    } else if (nb_tn == 2 && nb_tnp1 == 1){
        find_indices_le(level_set_tn, 2, idx); i_pt1=idx[0]; i_pt2=idx[1];
        find_indices_le(level_set_tnp1, 1, idx); i_pt3=idx[0];
    } else if (nb_tn == 2 && nb_tnp1 == 3){
        find_indices_gt(level_set_tn, 2, idx); i_pt1=idx[0]; i_pt2=idx[1];
        find_indices_gt(level_set_tnp1, 1, idx); i_pt3=idx[0];
    } else {
        printf("It should not happen: nb_tn + nb_tnp1 should equal 3 or 5, but here it equals %lld\n", nb_tn + nb_tnp1);
        return;
    }

    neighs_12 = is_adjacent_ref_square(i_pt1, i_pt2);
    neighs_13 = is_adjacent_ref_square(i_pt1, i_pt3);
    neighs_23 = is_adjacent_ref_square(i_pt2, i_pt3);

    one_face = (nb_tn == 3 && nb_tnp1 == 0) || (nb_tnp1 == 3 && nb_tn == 0) ||
               (nb_tnp1 == 4 && nb_tn == 1) || (nb_tnp1 == 1 && nb_tn == 4) ||
               (((nb_tn == 1 && nb_tnp1 == 2) || (nb_tn == 3 && nb_tnp1 == 2)) && neighs_23 && ((i_pt1 == i_pt2) || (i_pt1 == i_pt3))) ||
               (((nb_tn == 2 && nb_tnp1 == 1) || (nb_tn == 2 && nb_tnp1 == 3)) && neighs_12 && ((i_pt1 == i_pt3) || (i_pt2 == i_pt3)));
    three_faces = !one_face && ((((nb_tn == 1 && nb_tnp1 == 2) || (nb_tn == 3 && nb_tnp1 == 2)) && !neighs_23 && (i_pt1 != i_pt2) && (i_pt1 != i_pt3)) ||
                  (((nb_tn == 2 && nb_tnp1 == 1) || (nb_tn == 2 && nb_tnp1 == 3)) && !neighs_12 && (i_pt1 != i_pt3) && (i_pt2 != i_pt3)));

    if (one_face){
        Point3D verts[8];
        long int sf[8];
        int k;

        vertices    = alloc_with_capacity_vec_pts3D(8);
        infogrb = GrB_Matrix_new(edges,   GrB_INT8, 8, 14);
        infogrb = GrB_Matrix_new(faces,   GrB_INT8, 14, 8);
        infogrb = GrB_Matrix_new(volumes, GrB_INT8, 8, 1);
        status_face = alloc_with_capacity_vec_int(8);

        for (k = 0; k < 8; k++) infogrb = GrB_Matrix_setElement(*volumes, 1, k, 0);

        if ((nb_tn == 3 && nb_tnp1 == 0) || (nb_tnp1 == 4 && nb_tn == 1)){
            // I want point 2 to be in the middle, neighbour to 1 and 3, but 1 and 3 not neighbours.
            int i_pt4, e_14, e_12, e_34, e_23;
            int n1, n2, e1, e2;
            Point2D *p2D;
            Point3D pt13D_n, pt13D_np1, pt23D_n, pt23D_np1, pt33D_n, pt33D_np1, pt4;
            my_real_c lam1, lam2, lam3, lam14, lam34;

            if (neighs_13 && !neighs_12){
                int tmp = i_pt3; i_pt3 = i_pt2; i_pt2 = tmp;
                neighs_12 = is_adjacent_ref_square(i_pt1, i_pt2);
                neighs_13 = is_adjacent_ref_square(i_pt1, i_pt3);
                neighs_23 = is_adjacent_ref_square(i_pt2, i_pt3);
            }
            if (neighs_13 && !neighs_23){
                int tmp = i_pt1; i_pt1 = i_pt2; i_pt2 = tmp;
                neighs_12 = is_adjacent_ref_square(i_pt1, i_pt2);
                neighs_13 = is_adjacent_ref_square(i_pt1, i_pt3);
                neighs_23 = is_adjacent_ref_square(i_pt2, i_pt3);
            }

            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt1);
            pt13D_n   = (Point3D){p2D->x, p2D->y, 0.0};
            pt13D_np1 = (Point3D){p2D->x, p2D->y, dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt2);
            pt23D_n   = (Point3D){p2D->x, p2D->y, 0.0};
            pt23D_np1 = (Point3D){p2D->x, p2D->y, dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt3);
            pt33D_n   = (Point3D){p2D->x, p2D->y, 0.0};
            pt33D_np1 = (Point3D){p2D->x, p2D->y, dt};

            vertex_neighbors(i_pt1, &n1, &n2, &e1, &e2);
            if (n1 != i_pt2){ i_pt4 = n1; e_14 = e1; e_12 = e2; }
            else            { i_pt4 = n2; e_14 = e2; e_12 = e1; }
            vertex_neighbors(i_pt3, &n1, &n2, &e1, &e2);
            if (n1 != i_pt2){ e_34 = e1; e_23 = e2; }
            else            { e_34 = e2; e_23 = e1; }

            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt4);
            pt4 = (Point3D){p2D->x, p2D->y, 0.0};

            lam1  = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tnp1[i_pt1]);
            lam2  = level_set_tn[i_pt2] / (level_set_tn[i_pt2] - level_set_tnp1[i_pt2]);
            lam3  = level_set_tn[i_pt3] / (level_set_tn[i_pt3] - level_set_tnp1[i_pt3]);
            lam14 = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tn[i_pt4]);
            lam34 = level_set_tn[i_pt3] / (level_set_tn[i_pt3] - level_set_tn[i_pt4]);

            verts[0] = pt13D_n;
            verts[1] = pt23D_n;
            verts[2] = pt33D_n;
            verts[3] = lerp_pt3D(lam34, pt33D_n, pt4);
            verts[4] = lerp_pt3D(lam14, pt13D_n, pt4);
            verts[5] = lerp_pt3D(lam1,  pt13D_n, pt13D_np1);
            verts[6] = lerp_pt3D(lam2,  pt23D_n, pt23D_np1);
            verts[7] = lerp_pt3D(lam3,  pt33D_n, pt33D_np1);
            for (k = 0; k < 8; k++) push_back_vec_pts3D(&vertices, verts + k);

            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 0);  infogrb = GrB_Matrix_setElement(*edges,  1, 1, 0);
            infogrb = GrB_Matrix_setElement(*edges, -1, 1, 1);  infogrb = GrB_Matrix_setElement(*edges,  1, 2, 1);
            infogrb = GrB_Matrix_setElement(*edges, -1, 2, 2);  infogrb = GrB_Matrix_setElement(*edges,  1, 3, 2);
            infogrb = GrB_Matrix_setElement(*edges, -1, 3, 3);  infogrb = GrB_Matrix_setElement(*edges,  1, 4, 3);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 4);  infogrb = GrB_Matrix_setElement(*edges,  1, 0, 4);
            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 5);  infogrb = GrB_Matrix_setElement(*edges,  1, 5, 5);
            infogrb = GrB_Matrix_setElement(*edges, -1, 1, 6);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 6);
            infogrb = GrB_Matrix_setElement(*edges, -1, 2, 7);  infogrb = GrB_Matrix_setElement(*edges,  1, 7, 7);
            infogrb = GrB_Matrix_setElement(*edges, -1, 5, 8);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 8);
            infogrb = GrB_Matrix_setElement(*edges, -1, 6, 9);  infogrb = GrB_Matrix_setElement(*edges,  1, 7, 9);
            infogrb = GrB_Matrix_setElement(*edges, -1, 7, 10); infogrb = GrB_Matrix_setElement(*edges,  1, 3, 10);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 11); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 11);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 12); infogrb = GrB_Matrix_setElement(*edges,  1, 6, 12);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 13); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 13);

            // faces: net sign = -1 (literal value negated by the final unconditional flip)
            infogrb = GrB_Matrix_setElement(*faces,  1, 0, 0); infogrb = GrB_Matrix_setElement(*faces,  1, 1, 0); infogrb = GrB_Matrix_setElement(*faces,  1, 2, 0); infogrb = GrB_Matrix_setElement(*faces,  1, 3, 0); infogrb = GrB_Matrix_setElement(*faces,  1, 4, 0);
            infogrb = GrB_Matrix_setElement(*faces, -1, 8, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 11, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 12, 1);
            infogrb = GrB_Matrix_setElement(*faces, -1, 9, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 12, 2); infogrb = GrB_Matrix_setElement(*faces,  1, 13, 2);
            infogrb = GrB_Matrix_setElement(*faces, -1, 3, 3); infogrb = GrB_Matrix_setElement(*faces, -1, 10, 3); infogrb = GrB_Matrix_setElement(*faces, -1, 13, 3);
            infogrb = GrB_Matrix_setElement(*faces, -1, 4, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 5, 4); infogrb = GrB_Matrix_setElement(*faces,  1, 11, 4);
            infogrb = GrB_Matrix_setElement(*faces, -1, 0, 5); infogrb = GrB_Matrix_setElement(*faces,  1, 5, 5); infogrb = GrB_Matrix_setElement(*faces, -1, 6, 5); infogrb = GrB_Matrix_setElement(*faces,  1, 8, 5);
            infogrb = GrB_Matrix_setElement(*faces, -1, 1, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 6, 6); infogrb = GrB_Matrix_setElement(*faces, -1, 7, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 9, 6);
            infogrb = GrB_Matrix_setElement(*faces, -1, 2, 7); infogrb = GrB_Matrix_setElement(*faces,  1, 7, 7); infogrb = GrB_Matrix_setElement(*faces,  1, 10, 7);

            sf[0] = 1; sf[1] = -1; sf[2] = -1; sf[3] = -1;
            sf[4] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_14);
            sf[5] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_12);
            sf[6] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_23);
            sf[7] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_34);
            for (k = 0; k < 8; k++) push_back_vec_int(&status_face, &sf[k]);

        } else if ((nb_tnp1 == 3 && nb_tn == 0) || (nb_tn == 4 && nb_tnp1 == 1)){
            int i_pt4, e_14, e_12, e_34, e_23;
            int n1, n2, e1, e2;
            Point2D *p2D;
            Point3D pt13D_n, pt13D_np1, pt23D_n, pt23D_np1, pt33D_n, pt33D_np1, pt4;
            my_real_c lam1, lam2, lam3, lam14, lam34;

            if (neighs_13 && !neighs_12){
                int tmp = i_pt3; i_pt3 = i_pt2; i_pt2 = tmp;
                neighs_12 = is_adjacent_ref_square(i_pt1, i_pt2);
                neighs_13 = is_adjacent_ref_square(i_pt1, i_pt3);
                neighs_23 = is_adjacent_ref_square(i_pt2, i_pt3);
            }
            if (neighs_13 && !neighs_23){
                int tmp = i_pt1; i_pt1 = i_pt2; i_pt2 = tmp;
                neighs_12 = is_adjacent_ref_square(i_pt1, i_pt2);
                neighs_13 = is_adjacent_ref_square(i_pt1, i_pt3);
                neighs_23 = is_adjacent_ref_square(i_pt2, i_pt3);
            }

            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt1);
            pt13D_n   = (Point3D){p2D->x, p2D->y, 0.0};
            pt13D_np1 = (Point3D){p2D->x, p2D->y, dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt2);
            pt23D_n   = (Point3D){p2D->x, p2D->y, 0.0};
            pt23D_np1 = (Point3D){p2D->x, p2D->y, dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt3);
            pt33D_n   = (Point3D){p2D->x, p2D->y, 0.0};
            pt33D_np1 = (Point3D){p2D->x, p2D->y, dt};

            vertex_neighbors(i_pt1, &n1, &n2, &e1, &e2);
            if (n1 != i_pt2){ i_pt4 = n1; e_14 = e1; e_12 = e2; }
            else            { i_pt4 = n2; e_14 = e2; e_12 = e1; }
            vertex_neighbors(i_pt3, &n1, &n2, &e1, &e2);
            if (n1 != i_pt2){ e_34 = e1; e_23 = e2; }
            else            { e_34 = e2; e_23 = e1; }

            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt4);
            pt4 = (Point3D){p2D->x, p2D->y, dt};

            lam1  = level_set_tnp1[i_pt1] / (level_set_tnp1[i_pt1] - level_set_tn[i_pt1]);
            lam2  = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tn[i_pt2]);
            lam3  = level_set_tnp1[i_pt3] / (level_set_tnp1[i_pt3] - level_set_tn[i_pt3]);
            lam14 = level_set_tnp1[i_pt1] / (level_set_tnp1[i_pt1] - level_set_tnp1[i_pt4]);
            lam34 = level_set_tnp1[i_pt3] / (level_set_tnp1[i_pt3] - level_set_tnp1[i_pt4]);

            verts[0] = pt13D_np1;
            verts[1] = pt23D_np1;
            verts[2] = pt33D_np1;
            verts[3] = lerp_pt3D(lam34, pt33D_np1, pt4);
            verts[4] = lerp_pt3D(lam14, pt13D_np1, pt4);
            verts[5] = lerp_pt3D(lam1,  pt13D_np1, pt13D_n);
            verts[6] = lerp_pt3D(lam2,  pt23D_np1, pt23D_n);
            verts[7] = lerp_pt3D(lam3,  pt33D_np1, pt33D_n);
            for (k = 0; k < 8; k++) push_back_vec_pts3D(&vertices, &verts[k]);

            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 0);  infogrb = GrB_Matrix_setElement(*edges,  1, 1, 0);
            infogrb = GrB_Matrix_setElement(*edges, -1, 1, 1);  infogrb = GrB_Matrix_setElement(*edges,  1, 2, 1);
            infogrb = GrB_Matrix_setElement(*edges, -1, 2, 2);  infogrb = GrB_Matrix_setElement(*edges,  1, 3, 2);
            infogrb = GrB_Matrix_setElement(*edges, -1, 3, 3);  infogrb = GrB_Matrix_setElement(*edges,  1, 4, 3);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 4);  infogrb = GrB_Matrix_setElement(*edges,  1, 0, 4);
            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 5);  infogrb = GrB_Matrix_setElement(*edges,  1, 5, 5);
            infogrb = GrB_Matrix_setElement(*edges, -1, 1, 6);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 6);
            infogrb = GrB_Matrix_setElement(*edges, -1, 2, 7);  infogrb = GrB_Matrix_setElement(*edges,  1, 7, 7);
            infogrb = GrB_Matrix_setElement(*edges, -1, 5, 8);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 8);
            infogrb = GrB_Matrix_setElement(*edges, -1, 6, 9);  infogrb = GrB_Matrix_setElement(*edges,  1, 7, 9);
            infogrb = GrB_Matrix_setElement(*edges, -1, 7, 10); infogrb = GrB_Matrix_setElement(*edges,  1, 3, 10);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 11); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 11);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 12); infogrb = GrB_Matrix_setElement(*edges,  1, 6, 12);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 13); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 13);

            infogrb = GrB_Matrix_setElement(*faces, -1, 0, 0); infogrb = GrB_Matrix_setElement(*faces, -1, 1, 0); infogrb = GrB_Matrix_setElement(*faces, -1, 2, 0); infogrb = GrB_Matrix_setElement(*faces, -1, 3, 0); infogrb = GrB_Matrix_setElement(*faces, -1, 4, 0);
            infogrb = GrB_Matrix_setElement(*faces,  1, 8, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 11, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 12, 1);
            infogrb = GrB_Matrix_setElement(*faces,  1, 9, 2); infogrb = GrB_Matrix_setElement(*faces,  1, 12, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 13, 2);
            infogrb = GrB_Matrix_setElement(*faces,  1, 3, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 10, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 13, 3);
            infogrb = GrB_Matrix_setElement(*faces,  1, 4, 4); infogrb = GrB_Matrix_setElement(*faces,  1, 5, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 11, 4);
            infogrb = GrB_Matrix_setElement(*faces,  1, 0, 5); infogrb = GrB_Matrix_setElement(*faces, -1, 5, 5); infogrb = GrB_Matrix_setElement(*faces,  1, 6, 5); infogrb = GrB_Matrix_setElement(*faces, -1, 8, 5);
            infogrb = GrB_Matrix_setElement(*faces,  1, 1, 6); infogrb = GrB_Matrix_setElement(*faces, -1, 6, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 7, 6); infogrb = GrB_Matrix_setElement(*faces, -1, 9, 6);
            infogrb = GrB_Matrix_setElement(*faces,  1, 2, 7); infogrb = GrB_Matrix_setElement(*faces, -1, 7, 7); infogrb = GrB_Matrix_setElement(*faces, -1, 10, 7);

            sf[0] = 2; sf[1] = -1; sf[2] = -1; sf[3] = -1;
            sf[4] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_14);
            sf[5] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_12);
            sf[6] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_23);
            sf[7] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_34);
            for (k = 0; k < 8; k++) push_back_vec_int(&status_face, &sf[k]);

        } else {
            // "I want point 2 to be the point apart, and 1 and 3 to be stacked."
            int i_pt4, e_13, e_12, e_24;
            int n1, n2, e1, e2;
            Point2D *p2D;
            Point3D pt13D_n, pt13D_np1, pt23D_n, pt23D_np1, pt33D_n, pt33D_np1, pt4;
            my_real_c lam12, lam24, lam13_n, lam13_np1, lam2;
            int face_sign;

            if (i_pt1 != i_pt3){
                if (nb_tn == 1 || nb_tn == 3){ int tmp = i_pt3; i_pt3 = i_pt2; i_pt2 = tmp; }
                else                          { int tmp = i_pt1; i_pt1 = i_pt2; i_pt2 = tmp; }
            }

            vertex_neighbors(i_pt1, &n1, &n2, &e1, &e2);
            if (n1 == i_pt2){ i_pt3 = n2; e_13 = e2; e_12 = e1; }
            else            { i_pt3 = n1; e_13 = e1; e_12 = e2; }
            vertex_neighbors(i_pt2, &n1, &n2, &e1, &e2);
            if (n1 == i_pt1){ i_pt4 = n2; e_24 = e2; }
            else            { i_pt4 = n1; e_24 = e1; }

            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt1);
            pt13D_n   = (Point3D){p2D->x, p2D->y, 0.0};
            pt13D_np1 = (Point3D){p2D->x, p2D->y, dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt2);
            pt23D_n   = (Point3D){p2D->x, p2D->y, 0.0};
            pt23D_np1 = (Point3D){p2D->x, p2D->y, dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt3);
            pt33D_n   = (Point3D){p2D->x, p2D->y, 0.0};
            pt33D_np1 = (Point3D){p2D->x, p2D->y, dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt4);

            if (nb_tn == 1 || nb_tn == 3){
                pt4 = (Point3D){p2D->x, p2D->y, dt};
                lam12 = level_set_tn[i_pt1]   / (level_set_tn[i_pt1]   - level_set_tn[i_pt2]);
                lam24 = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tnp1[i_pt4]);
                verts[0] = pt23D_np1;
                verts[1] = lerp_pt3D(lam12, pt13D_n, pt23D_n);
                verts[2] = lerp_pt3D(lam24, pt23D_np1, pt4);

                lam13_n   = level_set_tn[i_pt1]   / (level_set_tn[i_pt1]   - level_set_tn[i_pt3]);
                lam13_np1 = level_set_tnp1[i_pt1] / (level_set_tnp1[i_pt1] - level_set_tnp1[i_pt3]);
                lam2      = level_set_tn[i_pt2]   / (level_set_tn[i_pt2]   - level_set_tnp1[i_pt2]);

                verts[3] = lerp_pt3D(lam2, pt23D_n, pt23D_np1);
                verts[4] = pt13D_n;
                verts[5] = lerp_pt3D(lam13_n, pt13D_n, pt33D_n);
                verts[6] = pt13D_np1;
                verts[7] = lerp_pt3D(lam13_np1, pt13D_np1, pt33D_np1);
            } else {
                pt4 = (Point3D){p2D->x, p2D->y, 0.0};
                lam12 = level_set_tnp1[i_pt1] / (level_set_tnp1[i_pt1] - level_set_tnp1[i_pt2]);
                lam24 = level_set_tn[i_pt2]   / (level_set_tn[i_pt2]   - level_set_tn[i_pt4]);
                verts[0] = pt23D_n;
                verts[1] = lerp_pt3D(lam12, pt13D_np1, pt23D_np1);
                verts[2] = lerp_pt3D(lam24, pt23D_n, pt4);

                lam13_n   = level_set_tn[i_pt1]   / (level_set_tn[i_pt1]   - level_set_tn[i_pt3]);
                lam13_np1 = level_set_tnp1[i_pt1] / (level_set_tnp1[i_pt1] - level_set_tnp1[i_pt3]);
                lam2      = level_set_tn[i_pt2]   / (level_set_tn[i_pt2]   - level_set_tnp1[i_pt2]);

                verts[3] = lerp_pt3D(lam2, pt23D_n, pt23D_np1);
                verts[4] = pt13D_np1;
                verts[5] = lerp_pt3D(lam13_np1, pt13D_np1, pt33D_np1);
                verts[6] = pt13D_n;
                verts[7] = lerp_pt3D(lam13_n, pt13D_n, pt33D_n);
            }
            for (k = 0; k < 8; k++) push_back_vec_pts3D(&vertices, &verts[k]);

            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 0);  infogrb = GrB_Matrix_setElement(*edges,  1, 3, 0);
            infogrb = GrB_Matrix_setElement(*edges, -1, 3, 1);  infogrb = GrB_Matrix_setElement(*edges,  1, 1, 1);
            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 2);  infogrb = GrB_Matrix_setElement(*edges,  1, 2, 2);
            infogrb = GrB_Matrix_setElement(*edges, -1, 2, 3);  infogrb = GrB_Matrix_setElement(*edges,  1, 3, 3);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 4);  infogrb = GrB_Matrix_setElement(*edges,  1, 1, 4);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 5);  infogrb = GrB_Matrix_setElement(*edges,  1, 5, 5);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 6);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 6);
            infogrb = GrB_Matrix_setElement(*edges, -1, 6, 7);  infogrb = GrB_Matrix_setElement(*edges,  1, 7, 7);
            infogrb = GrB_Matrix_setElement(*edges, -1, 7, 8);  infogrb = GrB_Matrix_setElement(*edges,  1, 2, 8);
            infogrb = GrB_Matrix_setElement(*edges, -1, 6, 9);  infogrb = GrB_Matrix_setElement(*edges,  1, 0, 9);
            infogrb = GrB_Matrix_setElement(*edges, -1, 5, 10); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 10);
            infogrb = GrB_Matrix_setElement(*edges, -1, 5, 11); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 11);
            infogrb = GrB_Matrix_setElement(*edges, -1, 5, 12); infogrb = GrB_Matrix_setElement(*edges,  1, 3, 12);
            infogrb = GrB_Matrix_setElement(*edges, -1, 5, 13); infogrb = GrB_Matrix_setElement(*edges,  1, 1, 13);

            if (((nb_tnp1 == 1) || (nb_tnp1 == 3)) && (grid_edge_sign(grid, i_pt2, e_12) > 0)) face_sign = 1;
            else if ((nb_tnp1 == 2) && (grid_edge_sign(grid, i_pt2, e_12) < 0))                face_sign = 1;
            else                                                                                face_sign = -1;

            infogrb = GrB_Matrix_setElement(*faces, -face_sign, 4, 0); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 5, 0); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 13, 0);
            infogrb = GrB_Matrix_setElement(*faces,  face_sign, 2, 1); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 7, 1); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 8, 1); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 9, 1);
            infogrb = GrB_Matrix_setElement(*faces, -face_sign, 5, 2); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 6, 2); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 7, 2); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 10, 2);
            infogrb = GrB_Matrix_setElement(*faces, -face_sign, 0, 3); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 1, 3); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 4, 3); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 6, 3); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 9, 3);
            infogrb = GrB_Matrix_setElement(*faces,  face_sign, 0, 4); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 2, 4); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 3, 4);
            infogrb = GrB_Matrix_setElement(*faces,  face_sign, 8, 5); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 10, 5); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 11, 5);
            infogrb = GrB_Matrix_setElement(*faces,  face_sign, 3, 6); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 11, 6); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 12, 6);
            infogrb = GrB_Matrix_setElement(*faces,  face_sign, 1, 7); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 12, 7); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 13, 7);

            sf[0] = (nb_tn == 1 || nb_tn == 3) ? 1 : 2;
            sf[1] = (nb_tn == 1 || nb_tn == 3) ? 2 : 1;
            sf[2] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_13);
            sf[3] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_12);
            sf[4] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_24);
            sf[5] = -1; sf[6] = -1; sf[7] = -1;
            for (k = 0; k < 8; k++) push_back_vec_int(&status_face, &sf[k]);
        }

    } else if (three_faces){
        int n1, n2, e1, e2;
        int i1_ptb, i1_ptc, e_ptb1, e_ptc1;
        int i2_ptb, i2_ptc, e_ptb2, e_ptc2;
        int i3_ptb, i3_ptc, e_ptb3, e_ptc3;
        Point2D *p2D;
        Point3D pt1, pta1, ptb, ptc, pt2, pta2, pt3, pta3;
        my_real_c lam_a, lam_b, lam_c;
        Point3D verts[12];
        long int sf[12];
        int k;

        vertices    = alloc_with_capacity_vec_pts3D(12);
        infogrb = GrB_Matrix_new(edges,   GrB_INT8, 12, 18);
        infogrb = GrB_Matrix_new(faces,   GrB_INT8, 18, 12);
        infogrb = GrB_Matrix_new(volumes, GrB_INT8, 12, 3);
        status_face = alloc_with_capacity_vec_int(12);

        vertex_neighbors(i_pt1, &n1, &n2, &e1, &e2);
        if (grid_edge_sign(grid, n1, e1) < 0){ i1_ptb = n1; i1_ptc = n2; e_ptb1 = e1; e_ptc1 = e2; }
        else                                  { i1_ptc = n1; i1_ptb = n2; e_ptc1 = e1; e_ptb1 = e2; }

        lam_a = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tnp1[i_pt1]);
        lam_b = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tn[i1_ptb]);
        lam_c = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tn[i1_ptc]);

        p2D = get_ith_elem_vec_pts2D(grid->vertices, i1_ptb); ptb = (Point3D){p2D->x, p2D->y, 0.0};
        p2D = get_ith_elem_vec_pts2D(grid->vertices, i1_ptc); ptc = (Point3D){p2D->x, p2D->y, 0.0};
        p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt1);
        pt1  = (Point3D){p2D->x, p2D->y, 0.0};
        pta1 = (Point3D){p2D->x, p2D->y, dt};

        verts[0] = pt1;
        verts[1] = lerp_pt3D(lam_a, pt1, pta1);
        verts[2] = lerp_pt3D(lam_b, pt1, ptb);
        verts[3] = lerp_pt3D(lam_c, pt1, ptc);

        vertex_neighbors(i_pt2, &n1, &n2, &e1, &e2);
        if (nb_tn == 1 || nb_tn == 3){ // second point is at t^{n+1}
            if (grid_edge_sign(grid, n1, e1) > 0){ i2_ptc = n1; i2_ptb = n2; e_ptc2 = e1; e_ptb2 = e2; }
            else                                  { i2_ptb = n1; i2_ptc = n2; e_ptb2 = e1; e_ptc2 = e2; }
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i2_ptb); ptb = (Point3D){p2D->x, p2D->y, dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i2_ptc); ptc = (Point3D){p2D->x, p2D->y, dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt2);
            pt2  = (Point3D){p2D->x, p2D->y, dt};
            pta2 = (Point3D){p2D->x, p2D->y, 0.0};
            lam_a = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tn[i_pt2]);
            lam_b = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tnp1[i2_ptb]);
            lam_c = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tnp1[i2_ptc]);
        } else { // second point is at t^n
            if (grid_edge_sign(grid, n1, e1) < 0){ i2_ptc = n1; i2_ptb = n2; e_ptc2 = e1; e_ptb2 = e2; }
            else                                  { i2_ptb = n1; i2_ptc = n2; e_ptb2 = e1; e_ptc2 = e2; }
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i2_ptb); ptb = (Point3D){p2D->x, p2D->y, 0.0};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i2_ptc); ptc = (Point3D){p2D->x, p2D->y, 0.0};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt2);
            pt2  = (Point3D){p2D->x, p2D->y, 0.0};
            pta2 = (Point3D){p2D->x, p2D->y, dt};
            lam_a = level_set_tn[i_pt2] / (level_set_tn[i_pt2] - level_set_tnp1[i_pt2]);
            lam_b = level_set_tn[i_pt2] / (level_set_tn[i_pt2] - level_set_tn[i2_ptb]);
            lam_c = level_set_tn[i_pt2] / (level_set_tn[i_pt2] - level_set_tn[i2_ptc]);
        }
        verts[4] = pt2;
        verts[5] = lerp_pt3D(lam_a, pt2, pta2);
        verts[6] = lerp_pt3D(lam_b, pt2, ptb);
        verts[7] = lerp_pt3D(lam_c, pt2, ptc);

        vertex_neighbors(i_pt3, &n1, &n2, &e1, &e2);
        if (grid_edge_sign(grid, n1, e1) > 0){ i3_ptc = n1; i3_ptb = n2; e_ptc3 = e1; e_ptb3 = e2; }
        else                                  { i3_ptb = n1; i3_ptc = n2; e_ptb3 = e1; e_ptc3 = e2; }
        p2D = get_ith_elem_vec_pts2D(grid->vertices, i3_ptb); ptb = (Point3D){p2D->x, p2D->y, dt};
        p2D = get_ith_elem_vec_pts2D(grid->vertices, i3_ptc); ptc = (Point3D){p2D->x, p2D->y, dt};
        p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt3);
        pt3  = (Point3D){p2D->x, p2D->y, dt};
        pta3 = (Point3D){p2D->x, p2D->y, 0.0};
        lam_a = level_set_tnp1[i_pt3] / (level_set_tnp1[i_pt3] - level_set_tn[i_pt3]);
        lam_b = level_set_tnp1[i_pt3] / (level_set_tnp1[i_pt3] - level_set_tnp1[i3_ptb]);
        lam_c = level_set_tnp1[i_pt3] / (level_set_tnp1[i_pt3] - level_set_tnp1[i3_ptc]);

        verts[8]  = pt3;
        verts[9]  = lerp_pt3D(lam_a, pt3, pta3);
        verts[10] = lerp_pt3D(lam_b, pt3, ptb);
        verts[11] = lerp_pt3D(lam_c, pt3, ptc);

        for (k = 0; k < 12; k++) push_back_vec_pts3D(&vertices, &verts[k]);

        infogrb = GrB_Matrix_setElement(*edges, -1, 0, 0);  infogrb = GrB_Matrix_setElement(*edges,  1, 1, 0);
        infogrb = GrB_Matrix_setElement(*edges, -1, 0, 1);  infogrb = GrB_Matrix_setElement(*edges,  1, 2, 1);
        infogrb = GrB_Matrix_setElement(*edges, -1, 0, 2);  infogrb = GrB_Matrix_setElement(*edges,  1, 3, 2);
        infogrb = GrB_Matrix_setElement(*edges, -1, 1, 3);  infogrb = GrB_Matrix_setElement(*edges,  1, 2, 3);
        infogrb = GrB_Matrix_setElement(*edges, -1, 1, 4);  infogrb = GrB_Matrix_setElement(*edges,  1, 3, 4);
        infogrb = GrB_Matrix_setElement(*edges, -1, 2, 5);  infogrb = GrB_Matrix_setElement(*edges,  1, 3, 5);
        infogrb = GrB_Matrix_setElement(*edges, -1, 4, 6);  infogrb = GrB_Matrix_setElement(*edges,  1, 5, 6);
        infogrb = GrB_Matrix_setElement(*edges, -1, 4, 7);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 7);
        infogrb = GrB_Matrix_setElement(*edges, -1, 4, 8);  infogrb = GrB_Matrix_setElement(*edges,  1, 7, 8);
        infogrb = GrB_Matrix_setElement(*edges, -1, 5, 9);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 9);
        infogrb = GrB_Matrix_setElement(*edges, -1, 5, 10); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 10);
        infogrb = GrB_Matrix_setElement(*edges, -1, 6, 11); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 11);
        infogrb = GrB_Matrix_setElement(*edges, -1, 8, 12); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 12);
        infogrb = GrB_Matrix_setElement(*edges, -1, 8, 13); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 13);
        infogrb = GrB_Matrix_setElement(*edges, -1, 8, 14); infogrb = GrB_Matrix_setElement(*edges,  1, 11, 14);
        infogrb = GrB_Matrix_setElement(*edges, -1, 9, 15); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 15);
        infogrb = GrB_Matrix_setElement(*edges, -1, 9, 16); infogrb = GrB_Matrix_setElement(*edges,  1, 11, 16);
        infogrb = GrB_Matrix_setElement(*edges, -1, 10, 17);infogrb = GrB_Matrix_setElement(*edges,  1, 11, 17);

        // faces: net sign = -1
        infogrb = GrB_Matrix_setElement(*faces, -1, 0, 0); infogrb = GrB_Matrix_setElement(*faces,  1, 1, 0); infogrb = GrB_Matrix_setElement(*faces, -1, 3, 0);
        infogrb = GrB_Matrix_setElement(*faces,  1, 0, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 4, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 2, 1);
        infogrb = GrB_Matrix_setElement(*faces,  1, 3, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 4, 2); infogrb = GrB_Matrix_setElement(*faces,  1, 5, 2);
        infogrb = GrB_Matrix_setElement(*faces, -1, 1, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 2, 3); infogrb = GrB_Matrix_setElement(*faces, -1, 5, 3);
        infogrb = GrB_Matrix_setElement(*faces,  1, 6, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 7, 4); infogrb = GrB_Matrix_setElement(*faces,  1, 9, 4);
        infogrb = GrB_Matrix_setElement(*faces, -1, 6, 5); infogrb = GrB_Matrix_setElement(*faces, -1, 10, 5); infogrb = GrB_Matrix_setElement(*faces,  1, 8, 5);
        infogrb = GrB_Matrix_setElement(*faces, -1, 9, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 10, 6); infogrb = GrB_Matrix_setElement(*faces, -1, 11, 6);
        infogrb = GrB_Matrix_setElement(*faces,  1, 7, 7); infogrb = GrB_Matrix_setElement(*faces, -1, 8, 7); infogrb = GrB_Matrix_setElement(*faces,  1, 11, 7);
        infogrb = GrB_Matrix_setElement(*faces,  1, 12, 8); infogrb = GrB_Matrix_setElement(*faces, -1, 13, 8); infogrb = GrB_Matrix_setElement(*faces,  1, 15, 8);
        infogrb = GrB_Matrix_setElement(*faces, -1, 12, 9); infogrb = GrB_Matrix_setElement(*faces, -1, 16, 9); infogrb = GrB_Matrix_setElement(*faces,  1, 14, 9);
        infogrb = GrB_Matrix_setElement(*faces, -1, 15, 10); infogrb = GrB_Matrix_setElement(*faces,  1, 16, 10); infogrb = GrB_Matrix_setElement(*faces, -1, 17, 10);
        infogrb = GrB_Matrix_setElement(*faces,  1, 13, 11); infogrb = GrB_Matrix_setElement(*faces, -1, 14, 11); infogrb = GrB_Matrix_setElement(*faces,  1, 17, 11);

        sf[0] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_ptb1);
        sf[1] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_ptc1);
        sf[2] = -1;
        sf[3] = 1;
        sf[4] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_ptb2);
        sf[5] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_ptc2);
        sf[6] = -1;
        sf[7] = ((nb_tn == 1) || (nb_tn == 3)) ? 2 : 1;
        sf[8] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_ptb3);
        sf[9] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_ptc3);
        sf[10] = -1;
        sf[11] = 2;
        for (k = 0; k < 12; k++) push_back_vec_int(&status_face, &sf[k]);

        // volumes[1:4,1]/[5:8,2]/[9:12,3] .= 1  (dense fill, 0-based).
        // The third fill was missing from the original Julia code; added here as confirmed.
        for (k = 0; k < 4; k++)  infogrb = GrB_Matrix_setElement(*volumes, 1, k, 0);
        for (k = 4; k < 8; k++)  infogrb = GrB_Matrix_setElement(*volumes, 1, k, 1);
        for (k = 8; k < 12; k++) infogrb = GrB_Matrix_setElement(*volumes, 1, k, 2);

    } else { // two_faces
        int i_pta, i_ptb, i_ptc, corner_tn;
        int n1, n2, e1, e2;
        int i_ptc1, i_ptc2, e_ptc1, e_ptc2;
        Point2D *p2D;
        Point3D ptc3D_n, ptc3D_np1, ptc13D, ptc23D;
        my_real_c lam_a, lam_b, lam_c;
        Point3D verts[10];
        long int sf[10];
        int k;

        vertices    = alloc_with_capacity_vec_pts3D(10);
        infogrb = GrB_Matrix_new(edges,   GrB_INT8, 10, 16);
        infogrb = GrB_Matrix_new(faces,   GrB_INT8, 16, 10);
        infogrb = GrB_Matrix_new(volumes, GrB_INT8, 10, 2);
        status_face = alloc_with_capacity_vec_int(10);

        // --- find the isolated corner point ---
        if (nb_tn == 1 || nb_tn == 3){
            if (neighs_23){ i_ptc = i_pt1; i_pta = i_pt2; i_ptb = i_pt3; corner_tn = 1; }
            else if (i_pt1 == i_pt2){ i_ptc = i_pt3; i_pta = i_pt1; i_ptb = i_pt2; corner_tn = 0; }
            else { i_ptc = i_pt2; i_pta = i_pt1; i_ptb = i_pt3; corner_tn = 0; }
        } else {
            if (neighs_12){ i_ptc = i_pt3; i_pta = i_pt1; i_ptb = i_pt2; corner_tn = 0; }
            else if (i_pt3 == i_pt1){ i_ptc = i_pt2; i_pta = i_pt1; i_ptb = i_pt3; corner_tn = 1; }
            else { i_ptc = i_pt1; i_pta = i_pt3; i_ptb = i_pt2; corner_tn = 1; }
        }

        // --- build the corner tetrahedron first ---
        p2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptc);
        ptc3D_n   = (Point3D){p2D->x, p2D->y, 0.0};
        ptc3D_np1 = (Point3D){p2D->x, p2D->y, dt};

        vertex_neighbors(i_ptc, &n1, &n2, &e1, &e2);
        if (grid_edge_sign(grid, n1, e1) < 0){ i_ptc1 = n1; i_ptc2 = n2; e_ptc1 = e1; e_ptc2 = e2; }
        else                                  { i_ptc2 = n1; i_ptc1 = n2; e_ptc2 = e1; e_ptc1 = e2; }

        if (corner_tn){
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptc1); ptc13D = (Point3D){p2D->x, p2D->y, 0.0};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptc2); ptc23D = (Point3D){p2D->x, p2D->y, 0.0};

            lam_a = level_set_tn[i_ptc] / (level_set_tn[i_ptc] - level_set_tn[i_ptc1]);
            lam_b = level_set_tn[i_ptc] / (level_set_tn[i_ptc] - level_set_tn[i_ptc2]);
            lam_c = level_set_tn[i_ptc] / (level_set_tn[i_ptc] - level_set_tnp1[i_ptc]);

            verts[0] = ptc3D_n;
            verts[1] = lerp_pt3D(lam_a, ptc3D_n, ptc13D);
            verts[2] = lerp_pt3D(lam_b, ptc3D_n, ptc23D);
            verts[3] = lerp_pt3D(lam_c, ptc3D_n, ptc3D_np1);
        } else {
            // Since we are at t^{n+1}, swap i_ptc1/i_ptc2 (and their edges) for a correct orientation.
            int tmp = i_ptc1; i_ptc1 = i_ptc2; i_ptc2 = tmp;
            tmp = e_ptc1; e_ptc1 = e_ptc2; e_ptc2 = tmp;

            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptc1); ptc13D = (Point3D){p2D->x, p2D->y, dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptc2); ptc23D = (Point3D){p2D->x, p2D->y, dt};

            lam_a = level_set_tnp1[i_ptc] / (level_set_tnp1[i_ptc] - level_set_tnp1[i_ptc1]);
            lam_b = level_set_tnp1[i_ptc] / (level_set_tnp1[i_ptc] - level_set_tnp1[i_ptc2]);
            lam_c = level_set_tn[i_ptc]   / (level_set_tn[i_ptc]   - level_set_tnp1[i_ptc]);

            verts[0] = ptc3D_np1;
            verts[1] = lerp_pt3D(lam_a, ptc3D_np1, ptc13D);
            verts[2] = lerp_pt3D(lam_b, ptc3D_np1, ptc23D);
            verts[3] = lerp_pt3D(lam_c, ptc3D_n, ptc3D_np1);
        }

        infogrb = GrB_Matrix_setElement(*edges, -1, 0, 0); infogrb = GrB_Matrix_setElement(*edges,  1, 1, 0);
        infogrb = GrB_Matrix_setElement(*edges, -1, 0, 1); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 1);
        infogrb = GrB_Matrix_setElement(*edges, -1, 0, 2); infogrb = GrB_Matrix_setElement(*edges,  1, 3, 2);
        infogrb = GrB_Matrix_setElement(*edges, -1, 3, 3); infogrb = GrB_Matrix_setElement(*edges,  1, 1, 3);
        infogrb = GrB_Matrix_setElement(*edges, -1, 3, 4); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 4);
        infogrb = GrB_Matrix_setElement(*edges, -1, 1, 5); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 5);

        infogrb = GrB_Matrix_setElement(*faces, -1, 0, 0); infogrb = GrB_Matrix_setElement(*faces,  1, 1, 0); infogrb = GrB_Matrix_setElement(*faces, -1, 5, 0);
        infogrb = GrB_Matrix_setElement(*faces,  1, 0, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 2, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 3, 1);
        infogrb = GrB_Matrix_setElement(*faces, -1, 1, 2); infogrb = GrB_Matrix_setElement(*faces,  1, 2, 2); infogrb = GrB_Matrix_setElement(*faces,  1, 4, 2);
        infogrb = GrB_Matrix_setElement(*faces,  1, 3, 3); infogrb = GrB_Matrix_setElement(*faces, -1, 4, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 5, 3);

        sf[0] = corner_tn ? 1 : 2;
        sf[1] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_ptc1);
        sf[2] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_ptc2);
        sf[3] = -1;

        for (k = 0; k < 4; k++) infogrb = GrB_Matrix_setElement(*volumes, 1, k, 0);

        // --- build the second volume ---
        if (i_pta == i_ptb){
            int i_ptb2, i_ptc2b, e_ptb, e_ptc;
            Point3D pta3D_n, pta3D_np1, ptb3D_n, ptb3D_np1, ptcc3D_n, ptcc3D_np1;

            vertex_neighbors(i_pta, &n1, &n2, &e1, &e2);
            if (grid_edge_sign(grid, n1, e1) < 0){ i_ptb2 = n1; i_ptc2b = n2; e_ptb = e1; e_ptc = e2; }
            else                                  { i_ptc2b = n1; i_ptb2 = n2; e_ptc = e1; e_ptb = e2; }

            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pta);
            pta3D_n = (Point3D){p2D->x, p2D->y, 0.0}; pta3D_np1 = (Point3D){p2D->x, p2D->y, dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptb2);
            ptb3D_n = (Point3D){p2D->x, p2D->y, 0.0}; ptb3D_np1 = (Point3D){p2D->x, p2D->y, dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptc2b);
            ptcc3D_n = (Point3D){p2D->x, p2D->y, 0.0}; ptcc3D_np1 = (Point3D){p2D->x, p2D->y, dt};

            lam_b = level_set_tn[i_pta] / (level_set_tn[i_pta] - level_set_tn[i_ptb2]);
            lam_c = level_set_tn[i_pta] / (level_set_tn[i_pta] - level_set_tn[i_ptc2b]);
            verts[4] = pta3D_n;
            verts[5] = lerp_pt3D(lam_b, pta3D_n, ptb3D_n);
            verts[6] = lerp_pt3D(lam_c, pta3D_n, ptcc3D_n);

            lam_b = level_set_tnp1[i_pta] / (level_set_tnp1[i_pta] - level_set_tnp1[i_ptb2]);
            lam_c = level_set_tnp1[i_pta] / (level_set_tnp1[i_pta] - level_set_tnp1[i_ptc2b]);
            verts[7] = pta3D_np1;
            verts[8] = lerp_pt3D(lam_b, pta3D_np1, ptb3D_np1);
            verts[9] = lerp_pt3D(lam_c, pta3D_np1, ptcc3D_np1);

            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 6);  infogrb = GrB_Matrix_setElement(*edges,  1, 7, 6);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 7);  infogrb = GrB_Matrix_setElement(*edges,  1, 5, 7);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 8);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 8);
            infogrb = GrB_Matrix_setElement(*edges, -1, 7, 9);  infogrb = GrB_Matrix_setElement(*edges,  1, 8, 9);
            infogrb = GrB_Matrix_setElement(*edges, -1, 7, 10); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 10);
            infogrb = GrB_Matrix_setElement(*edges, -1, 5, 11); infogrb = GrB_Matrix_setElement(*edges,  1, 8, 11);
            infogrb = GrB_Matrix_setElement(*edges, -1, 6, 12); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 12);
            infogrb = GrB_Matrix_setElement(*edges, -1, 5, 13); infogrb = GrB_Matrix_setElement(*edges,  1, 6, 13);
            infogrb = GrB_Matrix_setElement(*edges, -1, 8, 14); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 14);
            infogrb = GrB_Matrix_setElement(*edges, -1, 6, 15); infogrb = GrB_Matrix_setElement(*edges,  1, 8, 15);

            infogrb = GrB_Matrix_setElement(*faces, -1, 7, 4); infogrb = GrB_Matrix_setElement(*faces,  1, 8, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 13, 4);
            infogrb = GrB_Matrix_setElement(*faces,  1, 9, 5); infogrb = GrB_Matrix_setElement(*faces, -1, 10, 5); infogrb = GrB_Matrix_setElement(*faces,  1, 14, 5);
            infogrb = GrB_Matrix_setElement(*faces,  1, 6, 6); infogrb = GrB_Matrix_setElement(*faces, -1, 8, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 10, 6); infogrb = GrB_Matrix_setElement(*faces, -1, 12, 6);
            infogrb = GrB_Matrix_setElement(*faces, -1, 6, 7); infogrb = GrB_Matrix_setElement(*faces,  1, 7, 7); infogrb = GrB_Matrix_setElement(*faces, -1, 9, 7); infogrb = GrB_Matrix_setElement(*faces,  1, 11, 7);
            infogrb = GrB_Matrix_setElement(*faces, -1, 11, 8); infogrb = GrB_Matrix_setElement(*faces,  1, 13, 8); infogrb = GrB_Matrix_setElement(*faces,  1, 15, 8);
            infogrb = GrB_Matrix_setElement(*faces,  1, 12, 9); infogrb = GrB_Matrix_setElement(*faces, -1, 14, 9); infogrb = GrB_Matrix_setElement(*faces, -1, 15, 9);

            sf[4] = 1; sf[5] = 2;
            sf[6] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_ptc);
            sf[7] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_ptb);
            sf[8] = -1; sf[9] = -1;

        } else if (nb_tn == 2){
            int i_pt1l, i_pt2l, i_pt3l, i_pt4l, e_13, e_12, e_24;
            Point3D pt13D_n, pt13D_np1, pt23D_n, pt23D_np1, pt3, pt4;
            my_real_c lam_a2, lam_b2;

            i_pt1l = i_pta; i_pt2l = i_ptb;
            vertex_neighbors(i_pt1l, &n1, &n2, &e1, &e2);
            if (n1 == i_pt2l){
                if (grid_edge_sign(grid, i_pt2l, e1) > 0){
                    int tmp = i_pt2l; i_pt2l = i_pt1l; i_pt1l = tmp;
                    vertex_neighbors(i_pt1l, &n1, &n2, &e1, &e2);
                }
            } else {
                if (grid_edge_sign(grid, i_pt2l, e2) > 0){
                    int tmp = i_pt2l; i_pt2l = i_pt1l; i_pt1l = tmp;
                    vertex_neighbors(i_pt1l, &n1, &n2, &e1, &e2);
                }
            }
            if (n1 != i_pt2l){ i_pt3l = n1; e_13 = e1; e_12 = e2; }
            else              { i_pt3l = n2; e_13 = e2; e_12 = e1; }
            vertex_neighbors(i_pt2l, &n1, &n2, &e1, &e2);
            if (n1 != i_pt1l){ i_pt4l = n1; e_24 = e1; }
            else              { i_pt4l = n2; e_24 = e2; }

            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt1l);
            pt13D_n = (Point3D){p2D->x, p2D->y, 0.0}; pt13D_np1 = (Point3D){p2D->x, p2D->y, dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt2l);
            pt23D_n = (Point3D){p2D->x, p2D->y, 0.0}; pt23D_np1 = (Point3D){p2D->x, p2D->y, dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt3l);
            pt3 = (Point3D){p2D->x, p2D->y, 0.0};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt4l);
            pt4 = (Point3D){p2D->x, p2D->y, 0.0};

            lam_a2 = level_set_tn[i_pt1l] / (level_set_tn[i_pt1l] - level_set_tnp1[i_pt1l]);
            lam_b2 = level_set_tn[i_pt1l] / (level_set_tn[i_pt1l] - level_set_tn[i_pt3l]);
            verts[4] = pt13D_n;
            verts[5] = lerp_pt3D(lam_a2, pt13D_n, pt13D_np1);
            verts[6] = lerp_pt3D(lam_b2, pt13D_n, pt3);

            lam_a2 = level_set_tn[i_pt2l] / (level_set_tn[i_pt2l] - level_set_tnp1[i_pt2l]);
            lam_b2 = level_set_tn[i_pt2l] / (level_set_tn[i_pt2l] - level_set_tn[i_pt4l]);
            verts[7] = pt23D_n;
            verts[8] = lerp_pt3D(lam_a2, pt23D_n, pt23D_np1);
            verts[9] = lerp_pt3D(lam_b2, pt23D_n, pt4);

            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 6);  infogrb = GrB_Matrix_setElement(*edges,  1, 7, 6);
            infogrb = GrB_Matrix_setElement(*edges, -1, 5, 7);  infogrb = GrB_Matrix_setElement(*edges,  1, 8, 7);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 8);  infogrb = GrB_Matrix_setElement(*edges,  1, 5, 8);
            infogrb = GrB_Matrix_setElement(*edges, -1, 7, 9);  infogrb = GrB_Matrix_setElement(*edges,  1, 8, 9);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 10); infogrb = GrB_Matrix_setElement(*edges,  1, 6, 10);
            infogrb = GrB_Matrix_setElement(*edges, -1, 7, 11); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 11);
            infogrb = GrB_Matrix_setElement(*edges, -1, 5, 12); infogrb = GrB_Matrix_setElement(*edges,  1, 6, 12);
            infogrb = GrB_Matrix_setElement(*edges, -1, 8, 13); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 13);
            infogrb = GrB_Matrix_setElement(*edges, -1, 6, 14); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 14);
            infogrb = GrB_Matrix_setElement(*edges, -1, 5, 15); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 15);

            // faces (corrected)
            infogrb = GrB_Matrix_setElement(*faces, -1, 6, 4);  infogrb = GrB_Matrix_setElement(*faces,  1, 10, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 11, 4); infogrb = GrB_Matrix_setElement(*faces,  1, 14, 4);
            infogrb = GrB_Matrix_setElement(*faces,  1, 6, 5);  infogrb = GrB_Matrix_setElement(*faces, -1, 7, 5);  infogrb = GrB_Matrix_setElement(*faces, -1, 8, 5);  infogrb = GrB_Matrix_setElement(*faces,  1, 9, 5);
            infogrb = GrB_Matrix_setElement(*faces,  1, 8, 6);  infogrb = GrB_Matrix_setElement(*faces, -1, 10, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 12, 6);
            infogrb = GrB_Matrix_setElement(*faces, -1, 9, 7);  infogrb = GrB_Matrix_setElement(*faces,  1, 11, 7); infogrb = GrB_Matrix_setElement(*faces, -1, 13, 7);
            infogrb = GrB_Matrix_setElement(*faces, -1, 12, 8); infogrb = GrB_Matrix_setElement(*faces, -1, 14, 8); infogrb = GrB_Matrix_setElement(*faces,  1, 15, 8);
            infogrb = GrB_Matrix_setElement(*faces,  1, 7, 9);  infogrb = GrB_Matrix_setElement(*faces,  1, 13, 9); infogrb = GrB_Matrix_setElement(*faces, -1, 15, 9);

            sf[4] = 1;
            sf[5] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_12);
            sf[6] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_13);
            sf[7] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_24);
            sf[8] = -1; sf[9] = -1;

        } else { // nb_tnp1 == 2
            int i_pt1l, i_pt2l, i_pt3l, i_pt4l, e_13, e_12, e_24;
            Point3D pt13D_n, pt13D_np1, pt23D_n, pt23D_np1, pt3, pt4;
            my_real_c lam_a2, lam_b2;

            i_pt1l = i_pta; i_pt2l = i_ptb;
            vertex_neighbors(i_pt1l, &n1, &n2, &e1, &e2);
            if (n1 == i_pt2l){
                if (grid_edge_sign(grid, i_pt2l, e1) < 0){
                    int tmp = i_pt2l; i_pt2l = i_pt1l; i_pt1l = tmp;
                    vertex_neighbors(i_pt1l, &n1, &n2, &e1, &e2);
                }
            } else {
                if (grid_edge_sign(grid, i_pt2l, e2) < 0){
                    int tmp = i_pt2l; i_pt2l = i_pt1l; i_pt1l = tmp;
                    vertex_neighbors(i_pt1l, &n1, &n2, &e1, &e2);
                }
            }
            if (n1 != i_pt2l){ i_pt3l = n1; e_13 = e1; e_12 = e2; }
            else              { i_pt3l = n2; e_13 = e2; e_12 = e1; }
            vertex_neighbors(i_pt2l, &n1, &n2, &e1, &e2);
            if (n1 != i_pt1l){ i_pt4l = n1; e_24 = e1; }
            else              { i_pt4l = n2; e_24 = e2; }

            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt1l);
            pt13D_n = (Point3D){p2D->x, p2D->y, 0.0}; pt13D_np1 = (Point3D){p2D->x, p2D->y, dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt2l);
            pt23D_n = (Point3D){p2D->x, p2D->y, 0.0}; pt23D_np1 = (Point3D){p2D->x, p2D->y, dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt3l);
            pt3 = (Point3D){p2D->x, p2D->y, dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt4l);
            pt4 = (Point3D){p2D->x, p2D->y, dt};

            lam_a2 = level_set_tn[i_pt1l]   / (level_set_tn[i_pt1l]   - level_set_tnp1[i_pt1l]);
            lam_b2 = level_set_tnp1[i_pt1l] / (level_set_tnp1[i_pt1l] - level_set_tnp1[i_pt3l]);
            verts[4] = pt13D_np1;
            verts[5] = lerp_pt3D(lam_a2, pt13D_n, pt13D_np1);
            verts[6] = lerp_pt3D(lam_b2, pt13D_np1, pt3);

            lam_a2 = level_set_tn[i_pt2l]   / (level_set_tn[i_pt2l]   - level_set_tnp1[i_pt2l]);
            lam_b2 = level_set_tnp1[i_pt2l] / (level_set_tnp1[i_pt2l] - level_set_tnp1[i_pt4l]);
            verts[7] = pt23D_np1;
            verts[8] = lerp_pt3D(lam_a2, pt23D_n, pt23D_np1);
            verts[9] = lerp_pt3D(lam_b2, pt23D_np1, pt4);

            // same edge/face topology as the nb_tn == 2 sub-case above
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 6);  infogrb = GrB_Matrix_setElement(*edges,  1, 7, 6);
            infogrb = GrB_Matrix_setElement(*edges, -1, 5, 7);  infogrb = GrB_Matrix_setElement(*edges,  1, 8, 7);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 8);  infogrb = GrB_Matrix_setElement(*edges,  1, 5, 8);
            infogrb = GrB_Matrix_setElement(*edges, -1, 7, 9);  infogrb = GrB_Matrix_setElement(*edges,  1, 8, 9);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 10); infogrb = GrB_Matrix_setElement(*edges,  1, 6, 10);
            infogrb = GrB_Matrix_setElement(*edges, -1, 7, 11); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 11);
            infogrb = GrB_Matrix_setElement(*edges, -1, 5, 12); infogrb = GrB_Matrix_setElement(*edges,  1, 6, 12);
            infogrb = GrB_Matrix_setElement(*edges, -1, 8, 13); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 13);
            infogrb = GrB_Matrix_setElement(*edges, -1, 6, 14); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 14);
            infogrb = GrB_Matrix_setElement(*edges, -1, 5, 15); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 15);

            infogrb = GrB_Matrix_setElement(*faces, -1, 6, 4);  infogrb = GrB_Matrix_setElement(*faces,  1, 10, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 11, 4); infogrb = GrB_Matrix_setElement(*faces,  1, 14, 4);
            infogrb = GrB_Matrix_setElement(*faces,  1, 6, 5);  infogrb = GrB_Matrix_setElement(*faces, -1, 7, 5);  infogrb = GrB_Matrix_setElement(*faces, -1, 8, 5);  infogrb = GrB_Matrix_setElement(*faces,  1, 9, 5);
            infogrb = GrB_Matrix_setElement(*faces,  1, 8, 6);  infogrb = GrB_Matrix_setElement(*faces, -1, 10, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 12, 6);
            infogrb = GrB_Matrix_setElement(*faces, -1, 9, 7);  infogrb = GrB_Matrix_setElement(*faces,  1, 11, 7); infogrb = GrB_Matrix_setElement(*faces, -1, 13, 7);
            infogrb = GrB_Matrix_setElement(*faces, -1, 12, 8); infogrb = GrB_Matrix_setElement(*faces, -1, 14, 8); infogrb = GrB_Matrix_setElement(*faces,  1, 15, 8);
            infogrb = GrB_Matrix_setElement(*faces,  1, 7, 9);  infogrb = GrB_Matrix_setElement(*faces,  1, 13, 9); infogrb = GrB_Matrix_setElement(*faces, -1, 15, 9);

            sf[4] = 2;
            sf[5] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_12);
            sf[6] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_13);
            sf[7] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_24);
            sf[8] = -1; sf[9] = -1;
        }

        for (k = 0; k < 10; k++) push_back_vec_pts3D(&vertices, &verts[k]);
        for (k = 0; k < 10; k++) push_back_vec_int(&status_face, &sf[k]);
        for (k = 4; k < 10; k++) infogrb = GrB_Matrix_setElement(*volumes, 1, k, 1);
    }

    built = new_Polyhedron3D_vefvs(vertices, edges, faces, volumes, status_face);
    copy_Polyhedron3D(built, clipped3D);

    dealloc_Polyhedron3D(built); free(built);
    dealloc_vec_pts3D(vertices); free(vertices);
    dealloc_vec_int(status_face); free(status_face);
    GrB_free(edges);   free(edges);
    GrB_free(faces);   free(faces);
    GrB_free(volumes); free(volumes);
}

static Point2D bilinear_pt2D(const Polygon2D* grid, my_real_c xi, my_real_c eta){
    Point2D *v0 = get_ith_elem_vec_pts2D(grid->vertices, 0);
    Point2D *v1 = get_ith_elem_vec_pts2D(grid->vertices, 1);
    Point2D *v2 = get_ith_elem_vec_pts2D(grid->vertices, 2);
    Point2D *v3 = get_ith_elem_vec_pts2D(grid->vertices, 3);
    my_real_c w0 = (1-xi)*(1-eta), w1 = xi*(1-eta), w2 = xi*eta, w3 = (1-xi)*eta;
    Point2D r;
    r.x = w0*v0->x + w1*v1->x + w2*v2->x + w3*v3->x;
    r.y = w0*v0->y + w1*v1->y + w2*v2->y + w3*v3->y;
    return r;
}


#include <math.h>  // pour sqrt/fabs, si pas déjà inclus dans le fichier

/// @brief Finds a point (xi, eta, zeta) in [0,1]^3 (Q1 local coordinates) where the
///        bilinearly/trilinearly interpolated level-set vanishes, searching only at zeta = 0.5.
static void find_0pt_Q1(const my_real_c* level_set_tn, const my_real_c* level_set_tnp1,
                  my_real_c* xi, my_real_c* eta, my_real_c* zeta){
    my_real_c level[4];
    my_real_c sum_sq, denom, num, eta_loc, xi_loc;
    int i;

    for (i = 0; i < 4; i++) level[i] = 0.5 * (level_set_tn[i] + level_set_tnp1[i]);

    sum_sq = 0.0;
    for (i = 0; i < 4; i++) sum_sq += level[i]*level[i];
    if (sqrt(sum_sq) < 1e-10){
        *xi = 0.5; *eta = 0.5; *zeta = 0.5;
        return;
    }

    eta_loc = 0.5;
    while (eta_loc > 1e-3){
        denom = eta_loc*(level[0] - level[1] + level[2] - level[3]) - (level[0] - level[1]);
        num   = eta_loc*(level[0] - level[3]) - level[0];
        if (fabs(denom) > 1e-10){
            xi_loc = num / denom;
            if (xi_loc > 0 && xi_loc < 1){
                *xi = xi_loc; *eta = eta_loc; *zeta = 0.5;
                return;
            }
        }
        eta_loc /= 2;
    }

    // Haven't found a suitable eta decreasing from 0.5: try increasing it instead.
    eta_loc = 0.5;
    while (eta_loc > 1e-3){
        denom = (1-eta_loc)*(level[0] - level[1] + level[2] - level[3]) - (level[0] - level[1]);
        num   = (1-eta_loc)*(level[0] - level[3]) - level[0];
        if (fabs(denom) > 1e-10){
            xi_loc = num / denom;
            if (xi_loc > 0 && xi_loc < 1){
                *xi = xi_loc; *eta = (1 - eta_loc); *zeta = 0.5;
                return;
            }
        }
        eta_loc /= 2;
    }

    // Could still be a straight line: search for a point given xi = zeta = 0.5 instead.
    {
        my_real_c xi_fixed = 0.5;
        denom = xi_fixed*(level[0] - level[1] + level[2] - level[3]) - (level[0] - level[3]);
        num   = xi_fixed*(level[0] - level[1]) - level[0];
        if (fabs(denom) > 1e-10){
            my_real_c eta_val = num / denom;
            if (eta_val > 0 && eta_val < 1){
                *xi = xi_fixed; *eta = eta_val; *zeta = 0.5;
                return;
            }
        }
    }

    printf("Error: could not find a suitable point on the 0 level curve.\n");
    *xi = 0.5; *eta = 0.5; *zeta = 0.5; // fallback so callers don't read uninitialized values
}

/// @brief Builds the Polyhedron3D corresponding to a cell where exactly four vertices
///        (in space-time) are on the "inside" side of the level-set.
void polygon_from_level_set_4_pts(const Polygon2D* grid, Polyhedron3D* clipped3D, my_real_c dt,
                                   const my_real_c* level_set_tn, const my_real_c* level_set_tnp1,
                                   const long long int nb_tn, const long long int nb_tnp1){
    GrB_Info infogrb;
    int i_pt1, i_pt2, i_pt3, i_pt4;
    int neighs_12, neighs_34;
    int case1, case2, case3, case4, case5;
    int idx[4];
    Vector_points3D *vertices;
    Vector_int *status_face;
    GrB_Matrix *edges   = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    GrB_Matrix *faces   = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    GrB_Matrix *volumes = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    Polyhedron3D *built;

    // --- find i_pt1..i_pt4 (0-based) ---
    if (nb_tn == 4){
        find_indices_le(level_set_tn, 4, idx);
        i_pt1=idx[0]; i_pt2=idx[1]; i_pt3=idx[2]; i_pt4=idx[3];
    } else if (nb_tnp1 == 4){
        find_indices_le(level_set_tnp1, 4, idx);
        i_pt1=idx[0]; i_pt2=idx[1]; i_pt3=idx[2]; i_pt4=idx[3];
    } else if (nb_tn == 1 && nb_tnp1 == 3){
        find_indices_le(level_set_tn, 1, idx); i_pt1 = idx[0];
        find_indices_le(level_set_tnp1, 3, idx); i_pt2=idx[0]; i_pt3=idx[1]; i_pt4=idx[2];
    } else if (nb_tn == 3 && nb_tnp1 == 1){
        find_indices_le(level_set_tn, 3, idx); i_pt1=idx[0]; i_pt2=idx[1]; i_pt3=idx[2];
        find_indices_le(level_set_tnp1, 1, idx); i_pt4 = idx[0];
    } else if (nb_tn == 2 && nb_tnp1 == 2){
        find_indices_le(level_set_tn, 2, idx); i_pt1=idx[0]; i_pt2=idx[1];
        find_indices_le(level_set_tnp1, 2, idx); i_pt3=idx[0]; i_pt4=idx[1];
    } else {
        printf("It should not happen: nb_tn + nb_tnp1 should equal 4, but here it equals %lld\n", nb_tn + nb_tnp1);
        return;
    }

    neighs_12 = is_adjacent_ref_square(i_pt1, i_pt2);
    neighs_34 = is_adjacent_ref_square(i_pt3, i_pt4);

    case1 = (nb_tn == 4) || (nb_tnp1 == 4) ||
            ((nb_tn == 2) && neighs_12 && neighs_34 && i_pt1 == i_pt3 && i_pt2 == i_pt4);

    case2 = !case1 && (
                (nb_tn   == 3 && (i_pt4 != i_pt1 && i_pt4 != i_pt2 && i_pt4 != i_pt3)) ||
                (nb_tnp1 == 3 && (i_pt1 != i_pt4 && i_pt1 != i_pt2 && i_pt1 != i_pt3)) ||
                (nb_tn   == 2 && ((!neighs_12 && neighs_34) || (!neighs_34 && neighs_12)))
            );

    case3 = !case1 && !case2 && (
                (nb_tn == 2 && ((neighs_12 && neighs_34 && !(
                                    (i_pt1 == i_pt3 && i_pt2 != i_pt4) ||
                                    (i_pt1 == i_pt4 && i_pt2 != i_pt3) ||
                                    (i_pt2 == i_pt3 && i_pt1 != i_pt4) ||
                                    (i_pt2 == i_pt4 && i_pt1 != i_pt3)
                                    )) ||
                                (i_pt1 == i_pt3 && i_pt2 == i_pt4)))
            );

    case4 = !case1 && !case2 && !case3 && (
                (nb_tn   == 1 && (i_pt1 == i_pt4 || i_pt1 == i_pt2 || i_pt1 == i_pt3)) ||
                (nb_tnp1 == 1 && (i_pt4 == i_pt1 || i_pt4 == i_pt2 || i_pt4 == i_pt3)) ||
                (nb_tn == 2 && neighs_12 && neighs_34 && (
                        (i_pt1 == i_pt3 && i_pt2 != i_pt4) ||
                        (i_pt1 == i_pt4 && i_pt2 != i_pt3) ||
                        (i_pt2 == i_pt3 && i_pt1 != i_pt4) ||
                        (i_pt2 == i_pt4 && i_pt1 != i_pt3)
                        ))
            );

    case5 = !case1 && !case2 && !case3 && !case4 && (
                nb_tn == 2 && !neighs_12 && !neighs_34 &&
                i_pt1 != i_pt3 && i_pt1 != i_pt4 && i_pt2 != i_pt3 && i_pt2 != i_pt4
            );
    
    // --- case1 (8 vertices, 8x13 edges, 13x7 faces, 7 status_face, 7x1 volumes) ---
    if (case1){
        Point3D verts[8];
        long int sf[7];
        int k;

        vertices    = alloc_with_capacity_vec_pts3D(8);
        infogrb = GrB_Matrix_new(edges,   GrB_INT8, 8, 13);
        infogrb = GrB_Matrix_new(faces,   GrB_INT8, 13, 7);
        infogrb = GrB_Matrix_new(volumes, GrB_INT8, 7, 1);
        status_face = alloc_with_capacity_vec_int(7);

        if ((nb_tn == 4) || (nb_tnp1 == 4)){
            Point2D *p2D;
            Point3D pt13D_n, pt13D_np1, pt23D_n, pt23D_np1, pt33D_n, pt33D_np1, pt43D_n, pt43D_np1;
            my_real_c lam1, lam2, lam3, lam4;
            int face_sign;

            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt1);
            pt13D_n = (Point3D){p2D->x,p2D->y,0.0}; pt13D_np1 = (Point3D){p2D->x,p2D->y,dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt2);
            pt23D_n = (Point3D){p2D->x,p2D->y,0.0}; pt23D_np1 = (Point3D){p2D->x,p2D->y,dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt3);
            pt33D_n = (Point3D){p2D->x,p2D->y,0.0}; pt33D_np1 = (Point3D){p2D->x,p2D->y,dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt4);
            pt43D_n = (Point3D){p2D->x,p2D->y,0.0}; pt43D_np1 = (Point3D){p2D->x,p2D->y,dt};

            lam1 = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tnp1[i_pt1]);
            lam2 = level_set_tn[i_pt2] / (level_set_tn[i_pt2] - level_set_tnp1[i_pt2]);
            lam3 = level_set_tn[i_pt3] / (level_set_tn[i_pt3] - level_set_tnp1[i_pt3]);
            lam4 = level_set_tn[i_pt4] / (level_set_tn[i_pt4] - level_set_tnp1[i_pt4]);

            if (nb_tn == 4){
                verts[0]=pt13D_n; verts[1]=pt23D_n; verts[2]=pt33D_n; verts[3]=pt43D_n;
            } else {
                verts[0]=pt13D_np1; verts[1]=pt23D_np1; verts[2]=pt33D_np1; verts[3]=pt43D_np1;
            }
            verts[4] = lerp_pt3D(lam1, pt13D_n, pt13D_np1);
            verts[5] = lerp_pt3D(lam2, pt23D_n, pt23D_np1);
            verts[6] = lerp_pt3D(lam3, pt33D_n, pt33D_np1);
            verts[7] = lerp_pt3D(lam4, pt43D_n, pt43D_np1);
            for (k=0;k<8;k++) push_back_vec_pts3D(&vertices, &verts[k]);

            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 0);  infogrb = GrB_Matrix_setElement(*edges,  1, 4, 0);
            infogrb = GrB_Matrix_setElement(*edges, -1, 1, 1);  infogrb = GrB_Matrix_setElement(*edges,  1, 5, 1);
            infogrb = GrB_Matrix_setElement(*edges, -1, 2, 2);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 2);
            infogrb = GrB_Matrix_setElement(*edges, -1, 3, 3);  infogrb = GrB_Matrix_setElement(*edges,  1, 7, 3);
            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 4);  infogrb = GrB_Matrix_setElement(*edges,  1, 1, 4);
            infogrb = GrB_Matrix_setElement(*edges, -1, 1, 5);  infogrb = GrB_Matrix_setElement(*edges,  1, 2, 5);
            infogrb = GrB_Matrix_setElement(*edges, -1, 2, 6);  infogrb = GrB_Matrix_setElement(*edges,  1, 3, 6);
            infogrb = GrB_Matrix_setElement(*edges, -1, 3, 7);  infogrb = GrB_Matrix_setElement(*edges,  1, 0, 7);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 8);  infogrb = GrB_Matrix_setElement(*edges,  1, 5, 8);
            infogrb = GrB_Matrix_setElement(*edges, -1, 5, 9);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 9);
            infogrb = GrB_Matrix_setElement(*edges, -1, 6, 10); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 10);
            infogrb = GrB_Matrix_setElement(*edges, -1, 7, 11); infogrb = GrB_Matrix_setElement(*edges,  1, 4, 11);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 12); infogrb = GrB_Matrix_setElement(*edges,  1, 6, 12);

            face_sign = (nb_tnp1 == 4) ? 1 : -1;

            infogrb = GrB_Matrix_setElement(*faces, -face_sign, 4, 0); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 5, 0); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 6, 0); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 7, 0);
            infogrb = GrB_Matrix_setElement(*faces,  face_sign, 8, 1); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 9, 1); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 12, 1);
            infogrb = GrB_Matrix_setElement(*faces,  face_sign, 10, 2); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 11, 2); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 12, 2);
            infogrb = GrB_Matrix_setElement(*faces, -face_sign, 0, 3); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 1, 3); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 4, 3); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 8, 3);
            infogrb = GrB_Matrix_setElement(*faces, -face_sign, 1, 4); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 2, 4); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 5, 4); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 9, 4);
            infogrb = GrB_Matrix_setElement(*faces, -face_sign, 2, 5); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 3, 5); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 6, 5); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 10, 5);
            infogrb = GrB_Matrix_setElement(*faces,  face_sign, 0, 6); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 3, 6); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 7, 6); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 11, 6);

            sf[0] = (nb_tnp1 == 4) ? 2 : 1;
            sf[1] = -1;
            sf[2] = -1;
            sf[3] = 2 + *get_ith_elem_vec_int(grid->status_edge, 0);
            sf[4] = 2 + *get_ith_elem_vec_int(grid->status_edge, 1);
            sf[5] = 2 + *get_ith_elem_vec_int(grid->status_edge, 2);
            sf[6] = 2 + *get_ith_elem_vec_int(grid->status_edge, 3);
            for (k=0;k<7;k++) push_back_vec_int(&status_face, &sf[k]);

        } else {
            int n1, n2, e1, e2, e_13, e_12, e_24;
            Point2D *p2D;
            Point3D pt13D_n, pt13D_np1, pt23D_n, pt23D_np1, pt33D_n, pt33D_np1, pt43D_n, pt43D_np1;
            my_real_c lam1, lam2, lam3, lam4;

            vertex_neighbors(i_pt1, &n1, &n2, &e1, &e2);
            if (n1 == i_pt2){
                if (grid_edge_sign(grid, i_pt1, e1) > 0){ int tmp = i_pt1; i_pt1 = i_pt2; i_pt2 = tmp; }
            } else {
                if (grid_edge_sign(grid, i_pt1, e2) > 0){ int tmp = i_pt1; i_pt1 = i_pt2; i_pt2 = tmp; }
            }

            vertex_neighbors(i_pt1, &n1, &n2, &e1, &e2);
            if (n1 != i_pt2){ i_pt3 = n1; e_13 = e1; e_12 = e2; }
            else            { i_pt3 = n2; e_13 = e2; e_12 = e1; }
            vertex_neighbors(i_pt2, &n1, &n2, &e1, &e2);
            if (n1 != i_pt1){ i_pt4 = n1; e_24 = e1; }
            else            { i_pt4 = n2; e_24 = e2; }

            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt1);
            pt13D_n = (Point3D){p2D->x,p2D->y,0.0}; pt13D_np1 = (Point3D){p2D->x,p2D->y,dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt2);
            pt23D_n = (Point3D){p2D->x,p2D->y,0.0}; pt23D_np1 = (Point3D){p2D->x,p2D->y,dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt3);
            pt33D_n = (Point3D){p2D->x,p2D->y,0.0}; pt33D_np1 = (Point3D){p2D->x,p2D->y,dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt4);
            pt43D_n = (Point3D){p2D->x,p2D->y,0.0}; pt43D_np1 = (Point3D){p2D->x,p2D->y,dt};

            lam1 = level_set_tn[i_pt1]   / (level_set_tn[i_pt1]   - level_set_tn[i_pt3]);
            lam3 = level_set_tnp1[i_pt1] / (level_set_tnp1[i_pt1] - level_set_tnp1[i_pt3]);
            lam2 = level_set_tn[i_pt2]   / (level_set_tn[i_pt2]   - level_set_tn[i_pt4]);
            lam4 = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tnp1[i_pt4]);

            verts[0] = pt13D_n;
            verts[1] = pt13D_np1;
            verts[2] = pt23D_n;
            verts[3] = pt23D_np1;
            verts[4] = lerp_pt3D(lam1, pt13D_n, pt33D_n);
            verts[5] = lerp_pt3D(lam3, pt13D_np1, pt33D_np1);
            verts[6] = lerp_pt3D(lam2, pt23D_n, pt43D_n);
            verts[7] = lerp_pt3D(lam4, pt23D_np1, pt43D_np1);
            for (k=0;k<8;k++) push_back_vec_pts3D(&vertices, &verts[k]);

            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 0); infogrb = GrB_Matrix_setElement(*edges,  1, 4, 0);
            infogrb = GrB_Matrix_setElement(*edges, -1, 1, 1); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 1);
            infogrb = GrB_Matrix_setElement(*edges, -1, 2, 2); infogrb = GrB_Matrix_setElement(*edges,  1, 6, 2);
            infogrb = GrB_Matrix_setElement(*edges, -1, 3, 3); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 3);
            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 4); infogrb = GrB_Matrix_setElement(*edges,  1, 1, 4);
            infogrb = GrB_Matrix_setElement(*edges, -1, 1, 5); infogrb = GrB_Matrix_setElement(*edges,  1, 3, 5);
            infogrb = GrB_Matrix_setElement(*edges, -1, 3, 6); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 6);
            infogrb = GrB_Matrix_setElement(*edges, -1, 2, 7); infogrb = GrB_Matrix_setElement(*edges,  1, 0, 7);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 8); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 8);
            infogrb = GrB_Matrix_setElement(*edges, -1, 5, 9); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 9);
            infogrb = GrB_Matrix_setElement(*edges, -1, 7, 10); infogrb = GrB_Matrix_setElement(*edges,  1, 6, 10);
            infogrb = GrB_Matrix_setElement(*edges, -1, 6, 11); infogrb = GrB_Matrix_setElement(*edges,  1, 4, 11);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 12); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 12);

            // faces: net sign = -1 (literal negated by the final unconditional flip)
            infogrb = GrB_Matrix_setElement(*faces,  1, 4, 0); infogrb = GrB_Matrix_setElement(*faces,  1, 5, 0); infogrb = GrB_Matrix_setElement(*faces,  1, 6, 0); infogrb = GrB_Matrix_setElement(*faces,  1, 7, 0);
            infogrb = GrB_Matrix_setElement(*faces, -1, 8, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 9, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 12, 1);
            infogrb = GrB_Matrix_setElement(*faces, -1, 10, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 11, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 12, 2);
            infogrb = GrB_Matrix_setElement(*faces,  1, 0, 3); infogrb = GrB_Matrix_setElement(*faces, -1, 1, 3); infogrb = GrB_Matrix_setElement(*faces, -1, 4, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 8, 3);
            infogrb = GrB_Matrix_setElement(*faces,  1, 1, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 3, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 5, 4); infogrb = GrB_Matrix_setElement(*faces,  1, 9, 4);
            infogrb = GrB_Matrix_setElement(*faces, -1, 2, 5); infogrb = GrB_Matrix_setElement(*faces,  1, 3, 5); infogrb = GrB_Matrix_setElement(*faces, -1, 6, 5); infogrb = GrB_Matrix_setElement(*faces,  1, 10, 5);
            infogrb = GrB_Matrix_setElement(*faces, -1, 0, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 2, 6); infogrb = GrB_Matrix_setElement(*faces, -1, 7, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 11, 6);

            sf[0] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_12);
            sf[1] = -1;
            sf[2] = -1;
            sf[3] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_13);
            sf[4] = 2;
            sf[5] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_24);
            sf[6] = 1;
            for (k=0;k<7;k++) push_back_vec_int(&status_face, &sf[k]);
        }

        for (k=0;k<7;k++) infogrb = GrB_Matrix_setElement(*volumes, 1, k, 0);

    } else if (case2){
        int i_pta, i_ptb, i_ptc, i_ptd;
        Point3D verts[12];
        long int sf[12];
        int k;

        vertices    = alloc_with_capacity_vec_pts3D(12);
        infogrb = GrB_Matrix_new(edges,   GrB_INT8, 12, 20);
        infogrb = GrB_Matrix_new(faces,   GrB_INT8, 20, 12);
        infogrb = GrB_Matrix_new(volumes, GrB_INT8, 12, 2);
        status_face = alloc_with_capacity_vec_int(12);

        if (nb_tn == 1){
            // --- first: the bottom tetrahedron ---
            int n1, n2, e1, e2, e_12, e_13, e_34, e_24;
            Point2D *p2D;
            Point3D pt13D_n, pt13D_np1, pt23D_n, pt23D_np1, pt33D_n, pt33D_np1, pt43D_n, pt43D_np1;
            my_real_c lam_a, lam_b, lam_c;

            i_pta = i_pt2; i_ptb = i_pt3; i_ptc = i_pt4;

            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt1);
            pt13D_n = (Point3D){p2D->x,p2D->y,0.0}; pt13D_np1 = (Point3D){p2D->x,p2D->y,dt};

            vertex_neighbors(i_pt1, &n1, &n2, &e1, &e2);
            if (grid_edge_sign(grid, n1, e1) < 0){ i_pt2 = n1; i_pt3 = n2; e_12 = e1; e_13 = e2; }
            else                                  { i_pt3 = n1; i_pt2 = n2; e_13 = e1; e_12 = e2; }

            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt2);
            pt23D_n = (Point3D){p2D->x,p2D->y,0.0}; pt23D_np1 = (Point3D){p2D->x,p2D->y,dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt3);
            pt33D_n = (Point3D){p2D->x,p2D->y,0.0}; pt33D_np1 = (Point3D){p2D->x,p2D->y,dt};

            lam_a = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tn[i_pt2]);
            lam_b = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tn[i_pt3]);
            lam_c = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tnp1[i_pt1]);

            verts[0] = pt13D_n;
            verts[1] = lerp_pt3D(lam_a, pt13D_n, pt23D_n);
            verts[2] = lerp_pt3D(lam_b, pt13D_n, pt33D_n); // corrected: was lam_a in the Julia source
            verts[3] = lerp_pt3D(lam_c, pt13D_n, pt13D_np1);

            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 0); infogrb = GrB_Matrix_setElement(*edges,  1, 1, 0);
            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 1); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 1);
            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 2); infogrb = GrB_Matrix_setElement(*edges,  1, 3, 2);
            infogrb = GrB_Matrix_setElement(*edges, -1, 3, 3); infogrb = GrB_Matrix_setElement(*edges,  1, 1, 3);
            infogrb = GrB_Matrix_setElement(*edges, -1, 3, 4); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 4);
            infogrb = GrB_Matrix_setElement(*edges, -1, 1, 5); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 5);

            infogrb = GrB_Matrix_setElement(*faces, -1, 0, 0); infogrb = GrB_Matrix_setElement(*faces,  1, 1, 0); infogrb = GrB_Matrix_setElement(*faces, -1, 5, 0);
            infogrb = GrB_Matrix_setElement(*faces,  1, 0, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 2, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 3, 1);
            infogrb = GrB_Matrix_setElement(*faces, -1, 1, 2); infogrb = GrB_Matrix_setElement(*faces,  1, 2, 2); infogrb = GrB_Matrix_setElement(*faces,  1, 4, 2);
            infogrb = GrB_Matrix_setElement(*faces,  1, 3, 3); infogrb = GrB_Matrix_setElement(*faces, -1, 4, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 5, 3);

            sf[0] = 1;
            sf[1] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_12);
            sf[2] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_13);
            sf[3] = -1;

            // --- second, "weird" volume ---
            // We already have pt2 -> pt1 -> pt3; we only need to identify pt4, the corner opposite pt1.
            if      (i_pta != i_pt2 && i_pta != i_pt3) i_pt4 = i_pta;
            else if (i_ptb != i_pt2 && i_ptb != i_pt3) i_pt4 = i_ptb;
            else                                        i_pt4 = i_ptc;

            // The orientation is necessarily pt1->pt3->pt4->pt2.
            vertex_neighbors(i_pt4, &n1, &n2, &e1, &e2);
            if (n1 == i_pt3){ e_34 = e1; e_24 = e2; }
            else             { e_34 = e2; e_24 = e1; }

            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt4);
            pt43D_n = (Point3D){p2D->x,p2D->y,0.0}; pt43D_np1 = (Point3D){p2D->x,p2D->y,dt};

            lam_a = level_set_tn[i_pt2] / (level_set_tn[i_pt2] - level_set_tnp1[i_pt2]);
            lam_b = level_set_tn[i_pt3] / (level_set_tn[i_pt3] - level_set_tnp1[i_pt3]);
            lam_c = level_set_tn[i_pt4] / (level_set_tn[i_pt4] - level_set_tnp1[i_pt4]);

            verts[4] = pt33D_np1;
            verts[5] = lerp_pt3D(lam_b, pt33D_n, pt33D_np1);
            verts[6] = pt43D_np1;
            verts[7] = lerp_pt3D(lam_c, pt43D_n, pt43D_np1);
            verts[8] = pt23D_np1;
            verts[9] = lerp_pt3D(lam_a, pt23D_n, pt23D_np1);

            lam_a = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tnp1[i_pt1]);
            lam_b = level_set_tnp1[i_pt3] / (level_set_tnp1[i_pt3] - level_set_tnp1[i_pt1]);

            verts[10] = lerp_pt3D(lam_b, pt33D_np1, pt13D_np1);
            verts[11] = lerp_pt3D(lam_a, pt23D_np1, pt13D_np1);

            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 6);  infogrb = GrB_Matrix_setElement(*edges,  1, 5, 6);
            infogrb = GrB_Matrix_setElement(*edges, -1, 6, 7);  infogrb = GrB_Matrix_setElement(*edges,  1, 7, 7);
            infogrb = GrB_Matrix_setElement(*edges, -1, 8, 8);  infogrb = GrB_Matrix_setElement(*edges,  1, 9, 8);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 9);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 9);
            infogrb = GrB_Matrix_setElement(*edges, -1, 5, 10); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 10);
            infogrb = GrB_Matrix_setElement(*edges, -1, 6, 11); infogrb = GrB_Matrix_setElement(*edges,  1, 8, 11);
            infogrb = GrB_Matrix_setElement(*edges, -1, 7, 12); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 12);
            infogrb = GrB_Matrix_setElement(*edges, -1, 8, 13); infogrb = GrB_Matrix_setElement(*edges,  1, 11, 13);
            infogrb = GrB_Matrix_setElement(*edges, -1, 9, 14); infogrb = GrB_Matrix_setElement(*edges,  1, 11, 14);
            infogrb = GrB_Matrix_setElement(*edges, -1, 11, 15); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 15);
            infogrb = GrB_Matrix_setElement(*edges, -1, 10, 16); infogrb = GrB_Matrix_setElement(*edges,  1, 4, 16);
            infogrb = GrB_Matrix_setElement(*edges, -1, 10, 17); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 17);
            infogrb = GrB_Matrix_setElement(*edges, -1, 11, 18); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 18);
            infogrb = GrB_Matrix_setElement(*edges, -1, 11, 19); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 19);

            infogrb = GrB_Matrix_setElement(*faces, -1, 9, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 11, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 13, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 15, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 16, 4);
            infogrb = GrB_Matrix_setElement(*faces,  1, 6, 5); infogrb = GrB_Matrix_setElement(*faces,  1, 16, 5); infogrb = GrB_Matrix_setElement(*faces, -1, 17, 5);
            infogrb = GrB_Matrix_setElement(*faces, -1, 6, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 7, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 9, 6); infogrb = GrB_Matrix_setElement(*faces, -1, 10, 6);
            infogrb = GrB_Matrix_setElement(*faces, -1, 7, 7); infogrb = GrB_Matrix_setElement(*faces,  1, 8, 7); infogrb = GrB_Matrix_setElement(*faces,  1, 11, 7); infogrb = GrB_Matrix_setElement(*faces, -1, 12, 7);
            infogrb = GrB_Matrix_setElement(*faces, -1, 8, 8); infogrb = GrB_Matrix_setElement(*faces,  1, 13, 8); infogrb = GrB_Matrix_setElement(*faces, -1, 14, 8);
            infogrb = GrB_Matrix_setElement(*faces,  1, 15, 9); infogrb = GrB_Matrix_setElement(*faces,  1, 17, 9); infogrb = GrB_Matrix_setElement(*faces, -1, 18, 9);
            infogrb = GrB_Matrix_setElement(*faces,  1, 10, 10); infogrb = GrB_Matrix_setElement(*faces,  1, 18, 10); infogrb = GrB_Matrix_setElement(*faces, -1, 19, 10);
            infogrb = GrB_Matrix_setElement(*faces,  1, 12, 11); infogrb = GrB_Matrix_setElement(*faces,  1, 14, 11); infogrb = GrB_Matrix_setElement(*faces,  1, 19, 11);

            sf[4] = 2;
            sf[5] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_13);
            sf[6] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_34);
            sf[7] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_24);
            sf[8] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_12);
            sf[9] = -1; sf[10] = -1; sf[11] = -1;

        } else if (nb_tn == 3){
            // --- first: the top tetrahedron ---
            int n1, n2, e1, e2, e_24, e_34, e_13, e_12;
            Point2D *p2D;
            Point3D pt13D_n, pt13D_np1, pt23D_n, pt23D_np1, pt33D_n, pt33D_np1, pt43D_n, pt43D_np1;
            my_real_c lam_a, lam_b, lam_c;

            i_pta = i_pt1; i_ptb = i_pt2; i_ptc = i_pt3;

            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt4);
            pt43D_n = (Point3D){p2D->x,p2D->y,0.0}; pt43D_np1 = (Point3D){p2D->x,p2D->y,dt};

            vertex_neighbors(i_pt4, &n1, &n2, &e1, &e2);
            if (grid_edge_sign(grid, n1, e1) < 0){ i_pt2 = n2; i_pt3 = n1; e_24 = e2; e_34 = e1; }
            else                                  { i_pt3 = n2; i_pt2 = n1; e_34 = e2; e_24 = e1; }

            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt2);
            pt23D_n = (Point3D){p2D->x,p2D->y,0.0}; pt23D_np1 = (Point3D){p2D->x,p2D->y,dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt3);
            pt33D_n = (Point3D){p2D->x,p2D->y,0.0}; pt33D_np1 = (Point3D){p2D->x,p2D->y,dt};

            lam_a = level_set_tnp1[i_pt4] / (level_set_tnp1[i_pt4] - level_set_tnp1[i_pt2]);
            lam_b = level_set_tnp1[i_pt4] / (level_set_tnp1[i_pt4] - level_set_tnp1[i_pt3]);
            lam_c = level_set_tn[i_pt4]   / (level_set_tn[i_pt4]   - level_set_tnp1[i_pt4]);

            verts[0] = pt43D_np1;
            verts[1] = lerp_pt3D(lam_a, pt43D_np1, pt23D_np1);
            verts[2] = lerp_pt3D(lam_b, pt43D_np1, pt33D_np1); // corrected: was lam_a in the Julia source
            verts[3] = lerp_pt3D(lam_c, pt43D_n, pt43D_np1);

            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 0); infogrb = GrB_Matrix_setElement(*edges,  1, 1, 0);
            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 1); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 1);
            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 2); infogrb = GrB_Matrix_setElement(*edges,  1, 3, 2);
            infogrb = GrB_Matrix_setElement(*edges, -1, 3, 3); infogrb = GrB_Matrix_setElement(*edges,  1, 1, 3);
            infogrb = GrB_Matrix_setElement(*edges, -1, 3, 4); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 4);
            infogrb = GrB_Matrix_setElement(*edges, -1, 1, 5); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 5);

            infogrb = GrB_Matrix_setElement(*faces,  1, 0, 0); infogrb = GrB_Matrix_setElement(*faces, -1, 1, 0); infogrb = GrB_Matrix_setElement(*faces,  1, 5, 0);
            infogrb = GrB_Matrix_setElement(*faces, -1, 0, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 2, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 3, 1);
            infogrb = GrB_Matrix_setElement(*faces,  1, 1, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 2, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 4, 2);
            infogrb = GrB_Matrix_setElement(*faces, -1, 3, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 4, 3); infogrb = GrB_Matrix_setElement(*faces, -1, 5, 3);

            sf[0] = 2;
            sf[1] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_24);
            sf[2] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_34);
            sf[3] = -1;

            // --- second, "weird" volume ---
            // We already have pt2 -> pt4 -> pt3; we only need to identify pt1, the corner opposite pt4.
            if      (i_pta != i_pt2 && i_pta != i_pt3) i_pt1 = i_pta;
            else if (i_ptb != i_pt2 && i_ptb != i_pt3) i_pt1 = i_ptb;
            else                                        i_pt1 = i_ptc;

            // The orientation is necessarily pt4->pt3->pt1->pt2.
            vertex_neighbors(i_pt1, &n1, &n2, &e1, &e2);
            if (n1 == i_pt3){ e_13 = e1; e_12 = e2; }
            else             { e_13 = e2; e_12 = e1; }

            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt1);
            pt13D_n = (Point3D){p2D->x,p2D->y,0.0}; pt13D_np1 = (Point3D){p2D->x,p2D->y,dt};

            lam_a = level_set_tn[i_pt2] / (level_set_tn[i_pt2] - level_set_tnp1[i_pt2]);
            lam_b = level_set_tn[i_pt3] / (level_set_tn[i_pt3] - level_set_tnp1[i_pt3]);
            lam_c = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tnp1[i_pt1]);

            verts[4] = pt33D_n;
            verts[5] = lerp_pt3D(lam_b, pt33D_n, pt33D_np1);
            verts[6] = pt13D_n;
            verts[7] = lerp_pt3D(lam_c, pt13D_n, pt13D_np1);
            verts[8] = pt23D_n;
            verts[9] = lerp_pt3D(lam_a, pt23D_n, pt23D_np1);

            lam_a = level_set_tn[i_pt2] / (level_set_tn[i_pt2] - level_set_tn[i_pt4]);
            lam_b = level_set_tn[i_pt3] / (level_set_tn[i_pt3] - level_set_tn[i_pt4]);

            verts[10] = lerp_pt3D(lam_b, pt33D_n, pt43D_n);
            verts[11] = lerp_pt3D(lam_a, pt23D_n, pt43D_n);

            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 6);  infogrb = GrB_Matrix_setElement(*edges,  1, 5, 6);
            infogrb = GrB_Matrix_setElement(*edges, -1, 6, 7);  infogrb = GrB_Matrix_setElement(*edges,  1, 7, 7);
            infogrb = GrB_Matrix_setElement(*edges, -1, 8, 8);  infogrb = GrB_Matrix_setElement(*edges,  1, 9, 8);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 9);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 9);
            infogrb = GrB_Matrix_setElement(*edges, -1, 5, 10); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 10);
            infogrb = GrB_Matrix_setElement(*edges, -1, 6, 11); infogrb = GrB_Matrix_setElement(*edges,  1, 8, 11);
            infogrb = GrB_Matrix_setElement(*edges, -1, 7, 12); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 12);
            infogrb = GrB_Matrix_setElement(*edges, -1, 8, 13); infogrb = GrB_Matrix_setElement(*edges,  1, 11, 13);
            infogrb = GrB_Matrix_setElement(*edges, -1, 9, 14); infogrb = GrB_Matrix_setElement(*edges,  1, 11, 14);
            infogrb = GrB_Matrix_setElement(*edges, -1, 11, 15); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 15);
            infogrb = GrB_Matrix_setElement(*edges, -1, 10, 16); infogrb = GrB_Matrix_setElement(*edges,  1, 4, 16);
            infogrb = GrB_Matrix_setElement(*edges, -1, 10, 17); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 17);
            infogrb = GrB_Matrix_setElement(*edges, -1, 11, 18); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 18);
            infogrb = GrB_Matrix_setElement(*edges, -1, 11, 19); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 19);

            infogrb = GrB_Matrix_setElement(*faces, -1, 9, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 11, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 13, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 15, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 16, 4);
            infogrb = GrB_Matrix_setElement(*faces, -1, 6, 5); infogrb = GrB_Matrix_setElement(*faces, -1, 16, 5); infogrb = GrB_Matrix_setElement(*faces,  1, 17, 5);
            infogrb = GrB_Matrix_setElement(*faces,  1, 6, 6); infogrb = GrB_Matrix_setElement(*faces, -1, 7, 6); infogrb = GrB_Matrix_setElement(*faces, -1, 9, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 10, 6);
            infogrb = GrB_Matrix_setElement(*faces,  1, 7, 7); infogrb = GrB_Matrix_setElement(*faces, -1, 8, 7); infogrb = GrB_Matrix_setElement(*faces, -1, 11, 7); infogrb = GrB_Matrix_setElement(*faces,  1, 12, 7);
            infogrb = GrB_Matrix_setElement(*faces,  1, 8, 8); infogrb = GrB_Matrix_setElement(*faces, -1, 13, 8); infogrb = GrB_Matrix_setElement(*faces,  1, 14, 8);
            infogrb = GrB_Matrix_setElement(*faces, -1, 15, 9); infogrb = GrB_Matrix_setElement(*faces, -1, 17, 9); infogrb = GrB_Matrix_setElement(*faces,  1, 18, 9);
            infogrb = GrB_Matrix_setElement(*faces, -1, 10, 10); infogrb = GrB_Matrix_setElement(*faces, -1, 18, 10); infogrb = GrB_Matrix_setElement(*faces,  1, 19, 10);
            infogrb = GrB_Matrix_setElement(*faces, -1, 12, 11); infogrb = GrB_Matrix_setElement(*faces, -1, 14, 11); infogrb = GrB_Matrix_setElement(*faces, -1, 19, 11);

            sf[4] = 1;
            sf[5] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_34);
            sf[6] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_13);
            sf[7] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_12);
            sf[8] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_24);
            sf[9] = -1; sf[10] = -1; sf[11] = -1;

        } else { // nb_tn == 2
            int n1, n2, e1, e2, e_12, e_13, e_34, e_24;
            Point2D *p2D;
            Point3D pt13D_n, pt13D_np1, pt23D_n, pt23D_np1, pt33D_n, pt33D_np1, pt43D_n, pt43D_np1;
            my_real_c lam_a, lam_b, lam_c;

            if (neighs_12){
                // corner is at t^{n+1}, therefore necessarily i_pt3 or i_pt4
                if (i_pt3 == i_pt1 || i_pt3 == i_pt2){
                    i_ptc = i_pt4; i_pta = i_pt1; i_ptb = i_pt2; i_ptd = i_pt3;
                } else {
                    i_ptc = i_pt3; i_pta = i_pt1; i_ptb = i_pt2; i_ptd = i_pt4;
                }
                i_pt1 = i_ptc;

                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt1);
                pt13D_n = (Point3D){p2D->x,p2D->y,0.0}; pt13D_np1 = (Point3D){p2D->x,p2D->y,dt};

                vertex_neighbors(i_pt1, &n1, &n2, &e1, &e2);
                if (grid_edge_sign(grid, n1, e1) < 0){ i_pt2 = n1; i_pt3 = n2; e_12 = e1; e_13 = e2; }
                else                                  { i_pt3 = n1; i_pt2 = n2; e_13 = e1; e_12 = e2; }

                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt2);
                pt23D_n = (Point3D){p2D->x,p2D->y,0.0}; pt23D_np1 = (Point3D){p2D->x,p2D->y,dt};
                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt3);
                pt33D_n = (Point3D){p2D->x,p2D->y,0.0}; pt33D_np1 = (Point3D){p2D->x,p2D->y,dt};

                lam_a = level_set_tnp1[i_pt1] / (level_set_tnp1[i_pt1] - level_set_tnp1[i_pt2]);
                lam_b = level_set_tnp1[i_pt1] / (level_set_tnp1[i_pt1] - level_set_tnp1[i_pt3]);
                lam_c = level_set_tn[i_pt1]   / (level_set_tn[i_pt1]   - level_set_tnp1[i_pt1]);

                verts[0] = pt13D_np1;
                verts[1] = lerp_pt3D(lam_a, pt13D_np1, pt23D_np1);
                verts[2] = lerp_pt3D(lam_b, pt13D_np1, pt33D_np1);
                verts[3] = lerp_pt3D(lam_c, pt13D_n, pt13D_np1);

                infogrb = GrB_Matrix_setElement(*edges, -1, 0, 0); infogrb = GrB_Matrix_setElement(*edges,  1, 1, 0);
                infogrb = GrB_Matrix_setElement(*edges, -1, 0, 1); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 1);
                infogrb = GrB_Matrix_setElement(*edges, -1, 0, 2); infogrb = GrB_Matrix_setElement(*edges,  1, 3, 2);
                infogrb = GrB_Matrix_setElement(*edges, -1, 3, 3); infogrb = GrB_Matrix_setElement(*edges,  1, 1, 3);
                infogrb = GrB_Matrix_setElement(*edges, -1, 3, 4); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 4);
                infogrb = GrB_Matrix_setElement(*edges, -1, 1, 5); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 5);

                infogrb = GrB_Matrix_setElement(*faces, -1, 0, 0); infogrb = GrB_Matrix_setElement(*faces,  1, 1, 0); infogrb = GrB_Matrix_setElement(*faces, -1, 5, 0);
                infogrb = GrB_Matrix_setElement(*faces, -1, 0, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 2, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 3, 1);
                infogrb = GrB_Matrix_setElement(*faces,  1, 1, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 2, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 4, 2);
                infogrb = GrB_Matrix_setElement(*faces, -1, 3, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 4, 3); infogrb = GrB_Matrix_setElement(*faces, -1, 5, 3);

                sf[0] = 1;
                sf[1] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_12);
                sf[2] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_13);
                sf[3] = -1;

                // second volume: keep i_pt3 as the "negative" point; swap 2/3 if needed
                if (i_pt2 == i_pta || i_pt2 == i_ptb || i_pt2 == i_ptd){
                    int tmp = i_pt2; i_pt2 = i_pt3; i_pt3 = tmp;
                    p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt2);
                    pt23D_n = (Point3D){p2D->x,p2D->y,0.0}; pt23D_np1 = (Point3D){p2D->x,p2D->y,dt};
                    p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt3);
                    pt33D_n = (Point3D){p2D->x,p2D->y,0.0}; pt33D_np1 = (Point3D){p2D->x,p2D->y,dt};
                }

                if      (i_pta != i_pt2 && i_pta != i_pt3) i_pt4 = i_pta;
                else if (i_ptb != i_pt2 && i_ptb != i_pt3) i_pt4 = i_ptb;
                else                                        i_pt4 = i_ptd;

                vertex_neighbors(i_pt4, &n1, &n2, &e1, &e2);
                if (n1 == i_pt3){ e_34 = e1; e_24 = e2; }
                else             { e_34 = e2; e_24 = e1; }

                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt4);
                pt43D_n = (Point3D){p2D->x,p2D->y,0.0}; pt43D_np1 = (Point3D){p2D->x,p2D->y,dt};

                lam_a = level_set_tn[i_pt3] / (level_set_tn[i_pt3] - level_set_tn[i_pt1]);
                lam_b = level_set_tn[i_pt4] / (level_set_tn[i_pt4] - level_set_tn[i_pt2]);

                verts[4] = pt33D_n;
                verts[5] = lerp_pt3D(lam_a, pt33D_n, pt13D_n);
                verts[6] = pt43D_n;
                verts[7] = lerp_pt3D(lam_b, pt43D_n, pt23D_n);

                lam_a = level_set_tnp1[i_pt4] / (level_set_tnp1[i_pt4] - level_set_tnp1[i_pt3]);
                lam_b = level_set_tnp1[i_pt4] / (level_set_tnp1[i_pt4] - level_set_tnp1[i_pt2]);
                verts[8] = pt43D_np1;
                verts[9] = lerp_pt3D(lam_a, pt43D_np1, pt33D_np1);
                verts[10] = lerp_pt3D(lam_b, pt43D_np1, pt23D_np1);

                lam_a = level_set_tn[i_pt3] / (level_set_tn[i_pt3] - level_set_tnp1[i_pt3]);
                verts[11] = lerp_pt3D(lam_a, pt33D_n, pt33D_np1);

                infogrb = GrB_Matrix_setElement(*edges, -1, 4, 6);  infogrb = GrB_Matrix_setElement(*edges,  1, 5, 6);
                infogrb = GrB_Matrix_setElement(*edges, -1, 6, 7);  infogrb = GrB_Matrix_setElement(*edges,  1, 7, 7);
                infogrb = GrB_Matrix_setElement(*edges, -1, 8, 8);  infogrb = GrB_Matrix_setElement(*edges,  1, 9, 8);
                infogrb = GrB_Matrix_setElement(*edges, -1, 4, 9);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 9);
                infogrb = GrB_Matrix_setElement(*edges, -1, 5, 10); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 10);
                infogrb = GrB_Matrix_setElement(*edges, -1, 6, 11); infogrb = GrB_Matrix_setElement(*edges,  1, 8, 11);
                infogrb = GrB_Matrix_setElement(*edges, -1, 7, 12); infogrb = GrB_Matrix_setElement(*edges, 1, 10, 12);
                infogrb = GrB_Matrix_setElement(*edges, -1, 8, 13); infogrb = GrB_Matrix_setElement(*edges, 1, 10, 13);
                infogrb = GrB_Matrix_setElement(*edges, -1, 9, 14); infogrb = GrB_Matrix_setElement(*edges, 1, 11, 14);
                infogrb = GrB_Matrix_setElement(*edges, -1, 9, 15); infogrb = GrB_Matrix_setElement(*edges, 1, 10, 15);
                infogrb = GrB_Matrix_setElement(*edges, -1, 11, 16); infogrb = GrB_Matrix_setElement(*edges,  1, 4, 16);
                infogrb = GrB_Matrix_setElement(*edges, -1, 11, 17); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 17);
                infogrb = GrB_Matrix_setElement(*edges, -1, 10, 18); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 18);
                infogrb = GrB_Matrix_setElement(*edges, -1, 5, 19); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 19);

                infogrb = GrB_Matrix_setElement(*faces,  1, 8, 4); infogrb = GrB_Matrix_setElement(*faces,  1, 9, 4); infogrb = GrB_Matrix_setElement(*faces,  1, 11, 4); infogrb = GrB_Matrix_setElement(*faces,  1, 14, 4); infogrb = GrB_Matrix_setElement(*faces,  1, 16, 4);
                infogrb = GrB_Matrix_setElement(*faces,  1, 7, 5); infogrb = GrB_Matrix_setElement(*faces, -1, 11, 5); infogrb = GrB_Matrix_setElement(*faces,  1, 12, 5); infogrb = GrB_Matrix_setElement(*faces, -1, 13, 5);
                infogrb = GrB_Matrix_setElement(*faces,  1, 6, 6); infogrb = GrB_Matrix_setElement(*faces, -1, 7, 6); infogrb = GrB_Matrix_setElement(*faces, -1, 9, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 10, 6);
                infogrb = GrB_Matrix_setElement(*faces, -1, 8, 7); infogrb = GrB_Matrix_setElement(*faces,  1, 13, 7); infogrb = GrB_Matrix_setElement(*faces, -1, 15, 7);
                infogrb = GrB_Matrix_setElement(*faces, -1, 6, 8); infogrb = GrB_Matrix_setElement(*faces, -1, 16, 8); infogrb = GrB_Matrix_setElement(*faces,  1, 17, 8);
                infogrb = GrB_Matrix_setElement(*faces, -1, 14, 9); infogrb = GrB_Matrix_setElement(*faces, -1, 17, 9); infogrb = GrB_Matrix_setElement(*faces, -1, 19, 9);
                infogrb = GrB_Matrix_setElement(*faces,  1, 15, 10); infogrb = GrB_Matrix_setElement(*faces,  1, 18, 10); infogrb = GrB_Matrix_setElement(*faces,  1, 19, 10);
                infogrb = GrB_Matrix_setElement(*faces, -1, 10, 11); infogrb = GrB_Matrix_setElement(*faces, -1, 12, 11); infogrb = GrB_Matrix_setElement(*faces, -1, 18, 11);

                sf[4] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_34);
                sf[5] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_24);
                sf[6] = 1;
                sf[7] = 2;
                sf[8] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_13);
                sf[9] = -1; sf[10] = -1; sf[11] = -1;

            } else {
                // corner is at t^n, therefore necessarily i_pt1 or i_pt2
                if (i_pt1 == i_pt2 || i_pt1 == i_pt3){
                    i_ptc = i_pt2; i_pta = i_pt1; i_ptb = i_pt3; i_ptd = i_pt4;
                } else {
                    i_ptc = i_pt1; i_pta = i_pt2; i_ptb = i_pt3; i_ptd = i_pt4;
                }
                i_pt1 = i_ptc;

                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt1);
                pt13D_n = (Point3D){p2D->x,p2D->y,0.0}; pt13D_np1 = (Point3D){p2D->x,p2D->y,dt};

                vertex_neighbors(i_pt1, &n1, &n2, &e1, &e2);
                if (grid_edge_sign(grid, n1, e1) < 0){ i_pt2 = n2; i_pt3 = n1; e_12 = e2; e_13 = e1; }
                else                                  { i_pt3 = n2; i_pt2 = n1; e_13 = e2; e_12 = e1; }

                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt2);
                pt23D_n = (Point3D){p2D->x,p2D->y,0.0}; pt23D_np1 = (Point3D){p2D->x,p2D->y,dt};
                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt3);
                pt33D_n = (Point3D){p2D->x,p2D->y,0.0}; pt33D_np1 = (Point3D){p2D->x,p2D->y,dt};

                lam_a = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tn[i_pt2]);
                lam_b = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tn[i_pt3]);
                lam_c = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tnp1[i_pt1]);

                verts[0] = pt13D_n;
                verts[1] = lerp_pt3D(lam_a, pt13D_n, pt23D_n);
                verts[2] = lerp_pt3D(lam_b, pt13D_n, pt33D_n);
                verts[3] = lerp_pt3D(lam_c, pt13D_n, pt13D_np1);

                infogrb = GrB_Matrix_setElement(*edges, -1, 0, 0); infogrb = GrB_Matrix_setElement(*edges,  1, 1, 0);
                infogrb = GrB_Matrix_setElement(*edges, -1, 0, 1); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 1);
                infogrb = GrB_Matrix_setElement(*edges, -1, 0, 2); infogrb = GrB_Matrix_setElement(*edges,  1, 3, 2);
                infogrb = GrB_Matrix_setElement(*edges, -1, 3, 3); infogrb = GrB_Matrix_setElement(*edges,  1, 1, 3);
                infogrb = GrB_Matrix_setElement(*edges, -1, 3, 4); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 4);
                infogrb = GrB_Matrix_setElement(*edges, -1, 1, 5); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 5);

                infogrb = GrB_Matrix_setElement(*faces,  1, 0, 0); infogrb = GrB_Matrix_setElement(*faces, -1, 1, 0); infogrb = GrB_Matrix_setElement(*faces,  1, 5, 0);
                infogrb = GrB_Matrix_setElement(*faces, -1, 0, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 2, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 3, 1);
                infogrb = GrB_Matrix_setElement(*faces,  1, 1, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 2, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 4, 2);
                infogrb = GrB_Matrix_setElement(*faces, -1, 3, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 4, 3); infogrb = GrB_Matrix_setElement(*faces, -1, 5, 3);

                sf[0] = 1;
                sf[1] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_12);
                sf[2] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_13);
                sf[3] = -1;

                if (i_pt2 == i_pta || i_pt2 == i_ptb || i_pt2 == i_ptd){
                    int tmp = i_pt2; i_pt2 = i_pt3; i_pt3 = tmp;
                    p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt2);
                    pt23D_n = (Point3D){p2D->x,p2D->y,0.0}; pt23D_np1 = (Point3D){p2D->x,p2D->y,dt};
                    p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt3);
                    pt33D_n = (Point3D){p2D->x,p2D->y,0.0}; pt33D_np1 = (Point3D){p2D->x,p2D->y,dt};
                }

                if      (i_pta != i_pt2 && i_pta != i_pt3) i_pt4 = i_pta;
                else if (i_ptb != i_pt2 && i_ptb != i_pt3) i_pt4 = i_ptb;
                else                                        i_pt4 = i_ptd;

                vertex_neighbors(i_pt4, &n1, &n2, &e1, &e2);
                if (n1 == i_pt3){ e_34 = e1; e_24 = e2; }
                else             { e_34 = e2; e_24 = e1; }

                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt4);
                pt43D_n = (Point3D){p2D->x,p2D->y,0.0}; pt43D_np1 = (Point3D){p2D->x,p2D->y,dt};

                lam_a = level_set_tnp1[i_pt3] / (level_set_tnp1[i_pt3] - level_set_tnp1[i_pt1]);
                lam_b = level_set_tnp1[i_pt4] / (level_set_tnp1[i_pt4] - level_set_tnp1[i_pt2]);

                verts[4] = pt33D_np1;
                verts[5] = lerp_pt3D(lam_a, pt33D_np1, pt13D_np1);
                verts[6] = pt43D_np1;
                verts[7] = lerp_pt3D(lam_b, pt43D_np1, pt23D_np1);

                lam_a = level_set_tn[i_pt4] / (level_set_tn[i_pt4] - level_set_tn[i_pt3]);
                lam_b = level_set_tn[i_pt4] / (level_set_tn[i_pt4] - level_set_tn[i_pt2]);
                verts[8] = pt43D_n;
                verts[9] = lerp_pt3D(lam_a, pt43D_n, pt33D_n);
                verts[10] = lerp_pt3D(lam_b, pt43D_n, pt23D_n);

                lam_a = level_set_tn[i_pt3] / (level_set_tn[i_pt3] - level_set_tnp1[i_pt3]);
                verts[11] = lerp_pt3D(lam_a, pt33D_n, pt33D_np1);

                infogrb = GrB_Matrix_setElement(*edges, -1, 4, 6);  infogrb = GrB_Matrix_setElement(*edges,  1, 5, 6);
                infogrb = GrB_Matrix_setElement(*edges, -1, 6, 7);  infogrb = GrB_Matrix_setElement(*edges,  1, 7, 7);
                infogrb = GrB_Matrix_setElement(*edges, -1, 8, 8);  infogrb = GrB_Matrix_setElement(*edges,  1, 9, 8);
                infogrb = GrB_Matrix_setElement(*edges, -1, 4, 9);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 9);
                infogrb = GrB_Matrix_setElement(*edges, -1, 5, 10); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 10);
                infogrb = GrB_Matrix_setElement(*edges, -1, 6, 11); infogrb = GrB_Matrix_setElement(*edges,  1, 8, 11);
                infogrb = GrB_Matrix_setElement(*edges, -1, 7, 12); infogrb = GrB_Matrix_setElement(*edges, 1, 10, 12);
                infogrb = GrB_Matrix_setElement(*edges, -1, 8, 13); infogrb = GrB_Matrix_setElement(*edges, 1, 10, 13);
                infogrb = GrB_Matrix_setElement(*edges, -1, 9, 14); infogrb = GrB_Matrix_setElement(*edges, 1, 11, 14);
                infogrb = GrB_Matrix_setElement(*edges, -1, 9, 15); infogrb = GrB_Matrix_setElement(*edges, 1, 10, 15);
                infogrb = GrB_Matrix_setElement(*edges, -1, 11, 16); infogrb = GrB_Matrix_setElement(*edges,  1, 4, 16);
                infogrb = GrB_Matrix_setElement(*edges, -1, 11, 17); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 17);
                infogrb = GrB_Matrix_setElement(*edges, -1, 10, 18); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 18);
                infogrb = GrB_Matrix_setElement(*edges, -1, 5, 19); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 19);

                infogrb = GrB_Matrix_setElement(*faces, -1, 8, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 9, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 11, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 14, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 16, 4);
                infogrb = GrB_Matrix_setElement(*faces, -1, 7, 5); infogrb = GrB_Matrix_setElement(*faces,  1, 11, 5); infogrb = GrB_Matrix_setElement(*faces, -1, 12, 5); infogrb = GrB_Matrix_setElement(*faces,  1, 13, 5);
                infogrb = GrB_Matrix_setElement(*faces, -1, 6, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 7, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 9, 6); infogrb = GrB_Matrix_setElement(*faces, -1, 10, 6);
                infogrb = GrB_Matrix_setElement(*faces,  1, 8, 7); infogrb = GrB_Matrix_setElement(*faces, -1, 13, 7); infogrb = GrB_Matrix_setElement(*faces,  1, 15, 7);
                infogrb = GrB_Matrix_setElement(*faces,  1, 6, 8); infogrb = GrB_Matrix_setElement(*faces,  1, 16, 8); infogrb = GrB_Matrix_setElement(*faces, -1, 17, 8);
                infogrb = GrB_Matrix_setElement(*faces,  1, 14, 9); infogrb = GrB_Matrix_setElement(*faces,  1, 17, 9); infogrb = GrB_Matrix_setElement(*faces,  1, 19, 9);
                infogrb = GrB_Matrix_setElement(*faces, -1, 15, 10); infogrb = GrB_Matrix_setElement(*faces, -1, 18, 10); infogrb = GrB_Matrix_setElement(*faces, -1, 19, 10);
                infogrb = GrB_Matrix_setElement(*faces,  1, 10, 11); infogrb = GrB_Matrix_setElement(*faces,  1, 12, 11); infogrb = GrB_Matrix_setElement(*faces,  1, 18, 11);

                sf[4] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_34);
                sf[5] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_24);
                sf[6] = 2;
                sf[7] = 1;
                sf[8] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_12);
                sf[9] = -1; sf[10] = -1; sf[11] = -1;
            }
        }

        for (k = 0; k < 12; k++) push_back_vec_pts3D(&vertices, &verts[k]);
        for (k = 0; k < 12; k++) push_back_vec_int(&status_face, &sf[k]);
        // volumes[1:4,1] / volumes[4:12,2] .= 1  (row 4, 1-based, belongs to BOTH columns — preserved as in the Julia source)
        for (k = 0; k < 4;  k++) infogrb = GrB_Matrix_setElement(*volumes, 1, k, 0);
        for (k = 3; k < 12; k++) infogrb = GrB_Matrix_setElement(*volumes, 1, k, 1);

    } else if (case3){
        int i_pta, i_ptb, i_ptc, i_ptd;
        Point3D verts[12];
        long int sf[12];
        int k;

        vertices    = alloc_with_capacity_vec_pts3D(12);
        infogrb = GrB_Matrix_new(edges,   GrB_INT8, 12, 20);
        infogrb = GrB_Matrix_new(faces,   GrB_INT8, 20, 12);
        infogrb = GrB_Matrix_new(volumes, GrB_INT8, 12, 2);
        status_face = alloc_with_capacity_vec_int(12);

        if (neighs_12){
            int n1, n2, e1, e2, e_ab, e_ac, e_bd, e_cd;
            Point2D *p2D;
            Point3D ptA3D_n, ptA3D_np1, ptB3D_n, ptB3D_np1, ptC3D_n, ptC3D_np1, ptD3D_n, ptD3D_np1;
            my_real_c lam_a, lam_b;

            i_pta = i_pt1; i_ptb = i_pt2;
            vertex_neighbors(i_pt1, &n1, &n2, &e1, &e2);
            if (n1 == i_pt2) e_ab = e1;
            else             e_ab = e2;

            if (grid_edge_sign(grid, i_pt2, e_ab) < 0){
                i_pta = i_pt1; i_ptb = i_pt2;
                if (n1 != i_pt2){ i_ptc = n1; e_ac = e1; }
                else            { i_ptc = n2; e_ac = e2; }
                vertex_neighbors(i_pt2, &n1, &n2, &e1, &e2);
                if (n1 != i_pt1){ i_ptd = n1; e_bd = e1; }
                else            { i_ptd = n2; e_bd = e2; }
            } else {
                i_ptb = i_pt1; i_pta = i_pt2;
                if (n1 != i_pt2){ i_ptd = n1; e_bd = e1; }
                else            { i_ptd = n2; e_bd = e2; }
                vertex_neighbors(i_pt2, &n1, &n2, &e1, &e2);
                if (n1 != i_pt1){ i_ptc = n1; e_ac = e1; }
                else            { i_ptc = n2; e_ac = e2; }
            }

            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pta);
            ptA3D_n = (Point3D){p2D->x,p2D->y,0.0}; ptA3D_np1 = (Point3D){p2D->x,p2D->y,dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptb);
            ptB3D_n = (Point3D){p2D->x,p2D->y,0.0}; ptB3D_np1 = (Point3D){p2D->x,p2D->y,dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptc);
            ptC3D_n = (Point3D){p2D->x,p2D->y,0.0}; ptC3D_np1 = (Point3D){p2D->x,p2D->y,dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptd);
            ptD3D_n = (Point3D){p2D->x,p2D->y,0.0}; ptD3D_np1 = (Point3D){p2D->x,p2D->y,dt};

            // --- first volume ---
            lam_a = level_set_tn[i_pta] / (level_set_tn[i_pta] - level_set_tnp1[i_pta]);
            lam_b = level_set_tn[i_ptb] / (level_set_tn[i_ptb] - level_set_tnp1[i_ptb]);

            verts[0] = ptA3D_n;
            verts[1] = lerp_pt3D(lam_a, ptA3D_n, ptA3D_np1);
            verts[2] = ptB3D_n;
            verts[3] = lerp_pt3D(lam_b, ptB3D_n, ptB3D_np1);

            lam_a = level_set_tn[i_pta] / (level_set_tn[i_pta] - level_set_tn[i_ptc]);
            lam_b = level_set_tn[i_ptb] / (level_set_tn[i_ptb] - level_set_tn[i_ptd]);

            verts[4] = lerp_pt3D(lam_a, ptA3D_n, ptC3D_n);
            verts[5] = lerp_pt3D(lam_b, ptB3D_n, ptD3D_n);

            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 0); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 0);
            infogrb = GrB_Matrix_setElement(*edges, -1, 1, 1); infogrb = GrB_Matrix_setElement(*edges,  1, 3, 1);
            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 2); infogrb = GrB_Matrix_setElement(*edges,  1, 1, 2);
            infogrb = GrB_Matrix_setElement(*edges, -1, 2, 3); infogrb = GrB_Matrix_setElement(*edges,  1, 3, 3);
            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 4); infogrb = GrB_Matrix_setElement(*edges,  1, 4, 4);
            infogrb = GrB_Matrix_setElement(*edges, -1, 2, 5); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 5);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 6); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 6);
            infogrb = GrB_Matrix_setElement(*edges, -1, 1, 7); infogrb = GrB_Matrix_setElement(*edges,  1, 4, 7);
            infogrb = GrB_Matrix_setElement(*edges, -1, 3, 8); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 8);
            infogrb = GrB_Matrix_setElement(*edges, -1, 1, 9); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 9);

            // faces: net sign = -1
            infogrb = GrB_Matrix_setElement(*faces, -1, 0, 0); infogrb = GrB_Matrix_setElement(*faces,  1, 4, 0); infogrb = GrB_Matrix_setElement(*faces, -1, 5, 0); infogrb = GrB_Matrix_setElement(*faces,  1, 6, 0);
            infogrb = GrB_Matrix_setElement(*faces,  1, 0, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 1, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 2, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 3, 1);
            infogrb = GrB_Matrix_setElement(*faces,  1, 2, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 4, 2); infogrb = GrB_Matrix_setElement(*faces,  1, 7, 2);
            infogrb = GrB_Matrix_setElement(*faces, -1, 3, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 5, 3); infogrb = GrB_Matrix_setElement(*faces, -1, 8, 3);
            infogrb = GrB_Matrix_setElement(*faces,  1, 1, 4); infogrb = GrB_Matrix_setElement(*faces,  1, 8, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 9, 4);
            infogrb = GrB_Matrix_setElement(*faces, -1, 6, 5); infogrb = GrB_Matrix_setElement(*faces, -1, 7, 5); infogrb = GrB_Matrix_setElement(*faces,  1, 9, 5);

            sf[0] = 1;
            sf[1] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_ab);
            sf[2] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_ac);
            sf[3] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_bd);
            sf[4] = -1; sf[5] = -1;

            // --- second volume ---
            vertex_neighbors(i_ptc, &n1, &n2, &e1, &e2);
            if (n1 == i_ptd) e_cd = e1;
            else             e_cd = e2;

            lam_a = level_set_tn[i_ptd] / (level_set_tn[i_ptd] - level_set_tnp1[i_ptd]);
            lam_b = level_set_tn[i_ptc] / (level_set_tn[i_ptc] - level_set_tnp1[i_ptc]);

            verts[6] = ptD3D_np1;
            verts[7] = lerp_pt3D(lam_a, ptD3D_n, ptD3D_np1);
            verts[8] = ptC3D_np1;
            verts[9] = lerp_pt3D(lam_b, ptC3D_n, ptC3D_np1);

            lam_a = level_set_tnp1[i_pta] / (level_set_tnp1[i_pta] - level_set_tnp1[i_ptc]);
            lam_b = level_set_tnp1[i_ptb] / (level_set_tnp1[i_ptb] - level_set_tnp1[i_ptd]);

            verts[10] = lerp_pt3D(lam_b, ptB3D_np1, ptD3D_np1);
            verts[11] = lerp_pt3D(lam_a, ptA3D_np1, ptC3D_np1);

            infogrb = GrB_Matrix_setElement(*edges, -1, 6, 10);  infogrb = GrB_Matrix_setElement(*edges,  1, 8, 10);
            infogrb = GrB_Matrix_setElement(*edges, -1, 7, 11);  infogrb = GrB_Matrix_setElement(*edges,  1, 9, 11);
            infogrb = GrB_Matrix_setElement(*edges, -1, 6, 12);  infogrb = GrB_Matrix_setElement(*edges,  1, 7, 12);
            infogrb = GrB_Matrix_setElement(*edges, -1, 8, 13);  infogrb = GrB_Matrix_setElement(*edges,  1, 9, 13);
            infogrb = GrB_Matrix_setElement(*edges, -1, 6, 14);  infogrb = GrB_Matrix_setElement(*edges,  1, 10, 14);
            infogrb = GrB_Matrix_setElement(*edges, -1, 8, 15);  infogrb = GrB_Matrix_setElement(*edges,  1, 11, 15);
            infogrb = GrB_Matrix_setElement(*edges, -1, 10, 16); infogrb = GrB_Matrix_setElement(*edges,  1, 11, 16);
            infogrb = GrB_Matrix_setElement(*edges, -1, 7, 17);  infogrb = GrB_Matrix_setElement(*edges,  1, 10, 17);
            infogrb = GrB_Matrix_setElement(*edges, -1, 9, 18);  infogrb = GrB_Matrix_setElement(*edges,  1, 11, 18);
            infogrb = GrB_Matrix_setElement(*edges, -1, 7, 19);  infogrb = GrB_Matrix_setElement(*edges,  1, 11, 19);

            infogrb = GrB_Matrix_setElement(*faces,  1, 10, 6); infogrb = GrB_Matrix_setElement(*faces, -1, 14, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 15, 6); infogrb = GrB_Matrix_setElement(*faces, -1, 16, 6);
            infogrb = GrB_Matrix_setElement(*faces, -1, 10, 7); infogrb = GrB_Matrix_setElement(*faces,  1, 11, 7); infogrb = GrB_Matrix_setElement(*faces,  1, 12, 7); infogrb = GrB_Matrix_setElement(*faces, -1, 13, 7);
            infogrb = GrB_Matrix_setElement(*faces, -1, 12, 8); infogrb = GrB_Matrix_setElement(*faces,  1, 14, 8); infogrb = GrB_Matrix_setElement(*faces, -1, 17, 8);
            infogrb = GrB_Matrix_setElement(*faces,  1, 13, 9); infogrb = GrB_Matrix_setElement(*faces, -1, 15, 9); infogrb = GrB_Matrix_setElement(*faces,  1, 18, 9);
            infogrb = GrB_Matrix_setElement(*faces, -1, 11, 10); infogrb = GrB_Matrix_setElement(*faces, -1, 18, 10); infogrb = GrB_Matrix_setElement(*faces,  1, 19, 10);
            infogrb = GrB_Matrix_setElement(*faces,  1, 16, 11); infogrb = GrB_Matrix_setElement(*faces,  1, 17, 11); infogrb = GrB_Matrix_setElement(*faces, -1, 19, 11);

            sf[6] = 2;
            sf[7] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_cd);
            sf[8] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_bd);
            sf[9] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_ac);
            sf[10] = -1; sf[11] = -1;

        } else {
            // Here, i_pt1 == i_pt3 && i_pt2 == i_pt4
            int n1, n2, e1, e2, e_13, e_14, e_23, e_24;
            Point2D *p2D;
            Point3D pt13D_n, pt13D_np1, pt23D_n, pt23D_np1, pt33D_n, pt33D_np1, pt43D_n, pt43D_np1;
            my_real_c lam_a, lam_b;

            vertex_neighbors(i_pt1, &n1, &n2, &e1, &e2);
            if (grid_edge_sign(grid, n1, e1) < 0){ i_pt3 = n1; i_pt4 = n2; e_13 = e1; e_14 = e2; }
            else                                  { i_pt3 = n2; i_pt4 = n1; e_13 = e2; e_14 = e1; }
            vertex_neighbors(i_pt2, &n1, &n2, &e1, &e2);
            if (n1 == i_pt3){ e_23 = e1; e_24 = e2; }
            else            { e_23 = e2; e_24 = e1; }

            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt1);
            pt13D_n = (Point3D){p2D->x,p2D->y,0.0}; pt13D_np1 = (Point3D){p2D->x,p2D->y,dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt2);
            pt23D_n = (Point3D){p2D->x,p2D->y,0.0}; pt23D_np1 = (Point3D){p2D->x,p2D->y,dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt3);
            pt33D_n = (Point3D){p2D->x,p2D->y,0.0}; pt33D_np1 = (Point3D){p2D->x,p2D->y,dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt4);
            pt43D_n = (Point3D){p2D->x,p2D->y,0.0}; pt43D_np1 = (Point3D){p2D->x,p2D->y,dt};

            lam_a = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tn[i_pt3]);
            lam_b = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tn[i_pt4]);
            verts[0] = pt13D_n;
            verts[1] = lerp_pt3D(lam_a, pt13D_n, pt33D_n);
            verts[2] = lerp_pt3D(lam_b, pt13D_n, pt43D_n);

            lam_a = level_set_tnp1[i_pt1] / (level_set_tnp1[i_pt1] - level_set_tnp1[i_pt3]);
            lam_b = level_set_tnp1[i_pt1] / (level_set_tnp1[i_pt1] - level_set_tnp1[i_pt4]);
            verts[3] = pt13D_np1;
            verts[4] = lerp_pt3D(lam_a, pt13D_np1, pt33D_np1);
            verts[5] = lerp_pt3D(lam_b, pt13D_np1, pt43D_np1);

            lam_a = level_set_tn[i_pt2] / (level_set_tn[i_pt2] - level_set_tn[i_pt3]);
            lam_b = level_set_tn[i_pt2] / (level_set_tn[i_pt2] - level_set_tn[i_pt4]);
            verts[6] = pt23D_n;
            verts[7] = lerp_pt3D(lam_b, pt23D_n, pt43D_n);
            verts[8] = lerp_pt3D(lam_a, pt23D_n, pt33D_n);

            lam_a = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tnp1[i_pt3]);
            lam_b = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tnp1[i_pt4]);
            verts[9]  = pt23D_np1;
            verts[10] = lerp_pt3D(lam_b, pt23D_np1, pt43D_np1);
            verts[11] = lerp_pt3D(lam_a, pt23D_np1, pt33D_np1);

            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 0); infogrb = GrB_Matrix_setElement(*edges,  1, 1, 0);
            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 1); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 1);
            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 2); infogrb = GrB_Matrix_setElement(*edges,  1, 3, 2);
            infogrb = GrB_Matrix_setElement(*edges, -1, 3, 3); infogrb = GrB_Matrix_setElement(*edges,  1, 4, 3);
            infogrb = GrB_Matrix_setElement(*edges, -1, 3, 4); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 4);
            infogrb = GrB_Matrix_setElement(*edges, -1, 1, 5); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 5);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 6); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 6);
            infogrb = GrB_Matrix_setElement(*edges, -1, 1, 7); infogrb = GrB_Matrix_setElement(*edges,  1, 4, 7);
            infogrb = GrB_Matrix_setElement(*edges, -1, 2, 8); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 8);
            infogrb = GrB_Matrix_setElement(*edges, -1, 2, 9); infogrb = GrB_Matrix_setElement(*edges,  1, 4, 9);

            infogrb = GrB_Matrix_setElement(*edges, -1, 6, 10);  infogrb = GrB_Matrix_setElement(*edges,  1, 7, 10);
            infogrb = GrB_Matrix_setElement(*edges, -1, 6, 11);  infogrb = GrB_Matrix_setElement(*edges,  1, 8, 11);
            infogrb = GrB_Matrix_setElement(*edges, -1, 6, 12);  infogrb = GrB_Matrix_setElement(*edges,  1, 9, 12);
            infogrb = GrB_Matrix_setElement(*edges, -1, 9, 13);  infogrb = GrB_Matrix_setElement(*edges,  1, 10, 13);
            infogrb = GrB_Matrix_setElement(*edges, -1, 9, 14);  infogrb = GrB_Matrix_setElement(*edges,  1, 11, 14);
            infogrb = GrB_Matrix_setElement(*edges, -1, 7, 15);  infogrb = GrB_Matrix_setElement(*edges,  1, 8, 15);
            infogrb = GrB_Matrix_setElement(*edges, -1, 10, 16); infogrb = GrB_Matrix_setElement(*edges,  1, 11, 16);
            infogrb = GrB_Matrix_setElement(*edges, -1, 7, 17);  infogrb = GrB_Matrix_setElement(*edges,  1, 10, 17);
            infogrb = GrB_Matrix_setElement(*edges, -1, 8, 18);  infogrb = GrB_Matrix_setElement(*edges,  1, 11, 18);
            infogrb = GrB_Matrix_setElement(*edges, -1, 8, 19);  infogrb = GrB_Matrix_setElement(*edges,  1, 10, 19);

            // faces: net sign = -1
            infogrb = GrB_Matrix_setElement(*faces, -1, 0, 0); infogrb = GrB_Matrix_setElement(*faces,  1, 1, 0); infogrb = GrB_Matrix_setElement(*faces, -1, 5, 0);
            infogrb = GrB_Matrix_setElement(*faces,  1, 3, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 4, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 6, 1);
            infogrb = GrB_Matrix_setElement(*faces,  1, 0, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 2, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 3, 2); infogrb = GrB_Matrix_setElement(*faces,  1, 7, 2);
            infogrb = GrB_Matrix_setElement(*faces, -1, 1, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 2, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 4, 3); infogrb = GrB_Matrix_setElement(*faces, -1, 8, 3);
            infogrb = GrB_Matrix_setElement(*faces,  1, 5, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 7, 4); infogrb = GrB_Matrix_setElement(*faces,  1, 9, 4);
            infogrb = GrB_Matrix_setElement(*faces, -1, 6, 5); infogrb = GrB_Matrix_setElement(*faces,  1, 8, 5); infogrb = GrB_Matrix_setElement(*faces, -1, 9, 5);

            infogrb = GrB_Matrix_setElement(*faces, -1, 10, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 11, 6); infogrb = GrB_Matrix_setElement(*faces, -1, 15, 6);
            infogrb = GrB_Matrix_setElement(*faces,  1, 13, 7); infogrb = GrB_Matrix_setElement(*faces, -1, 14, 7); infogrb = GrB_Matrix_setElement(*faces,  1, 16, 7);
            infogrb = GrB_Matrix_setElement(*faces,  1, 10, 8); infogrb = GrB_Matrix_setElement(*faces, -1, 12, 8); infogrb = GrB_Matrix_setElement(*faces, -1, 13, 8); infogrb = GrB_Matrix_setElement(*faces,  1, 17, 8);
            infogrb = GrB_Matrix_setElement(*faces, -1, 11, 9); infogrb = GrB_Matrix_setElement(*faces,  1, 12, 9); infogrb = GrB_Matrix_setElement(*faces,  1, 14, 9); infogrb = GrB_Matrix_setElement(*faces, -1, 18, 9);
            infogrb = GrB_Matrix_setElement(*faces,  1, 15, 10); infogrb = GrB_Matrix_setElement(*faces, -1, 17, 10); infogrb = GrB_Matrix_setElement(*faces,  1, 19, 10);
            infogrb = GrB_Matrix_setElement(*faces, -1, 16, 11); infogrb = GrB_Matrix_setElement(*faces,  1, 18, 11); infogrb = GrB_Matrix_setElement(*faces, -1, 19, 11);

            sf[0] = 1;
            sf[1] = 2;
            sf[2] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_13);
            sf[3] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_14);
            sf[4] = -1; sf[5] = -1;
            sf[6] = 1;
            sf[7] = 2;
            sf[8] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_24);
            sf[9] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_23);
            sf[10] = -1; sf[11] = -1;
        }

        for (k = 0; k < 12; k++) push_back_vec_pts3D(&vertices, &verts[k]);
        for (k = 0; k < 12; k++) push_back_vec_int(&status_face, &sf[k]);
        // volumes[1:6,1] / volumes[6:12,2] .= 1  (row 6, 1-based, belongs to BOTH columns — preserved as in the Julia source)
        for (k = 0; k < 6;  k++) infogrb = GrB_Matrix_setElement(*volumes, 1, k, 0);
        for (k = 5; k < 12; k++) infogrb = GrB_Matrix_setElement(*volumes, 1, k, 1);

    } else if (case4) {
        Point3D verts[11];
        long int sf[12];
        int k;

        vertices    = alloc_with_capacity_vec_pts3D(11);
        infogrb = GrB_Matrix_new(edges,   GrB_INT8, 11, 21);
        infogrb = GrB_Matrix_new(faces,   GrB_INT8, 21, 12);
        infogrb = GrB_Matrix_new(volumes, GrB_INT8, 12, 1);
        status_face = alloc_with_capacity_vec_int(12);

        if (nb_tn == 1){
            int n1, n2, e1, e2;
            int i_pta, i_ptb, i_ptc, i_ptd;
            int e_12, e_13, e_23, e_34, e_14;
            int neighs_23, neighs_24;
            Point2D *p2D;
            Point3D pt13D_n, pt13D_np1, pt23D_n, pt33D_n, pt23D_np1;
            Point3D pt33D_np1, pt43D_n, pt43D_np1;
            my_real_c lam_a, lam_b, xi, eta, zeta;

            vertex_neighbors(i_pt1, &n1, &n2, &e1, &e2);
            if (grid_edge_sign(grid, n1, e1) < 0){ i_pta = n1; i_ptb = n2; e_12 = e1; e_13 = e2; }
            else                                  { i_ptb = n1; i_pta = n2; e_13 = e1; e_12 = e2; }

            lam_a = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tn[i_pta]);
            lam_b = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tn[i_ptb]);

            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt1);
            pt13D_n = (Point3D){p2D->x,p2D->y,0.0}; pt13D_np1 = (Point3D){p2D->x,p2D->y,dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pta);
            pt23D_n = (Point3D){p2D->x,p2D->y,0.0};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptb);
            pt33D_n = (Point3D){p2D->x,p2D->y,0.0};

            verts[0] = pt13D_n;
            verts[1] = lerp_pt3D(lam_a, pt13D_n, pt23D_n);
            verts[2] = lerp_pt3D(lam_b, pt13D_n, pt33D_n);
            verts[3] = pt13D_np1;

            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 0); infogrb = GrB_Matrix_setElement(*edges,  1, 1, 0);
            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 1); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 1);
            infogrb = GrB_Matrix_setElement(*edges, -1, 1, 2); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 2);
            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 3); infogrb = GrB_Matrix_setElement(*edges,  1, 3, 3);

            // faces: net sign = -1
            infogrb = GrB_Matrix_setElement(*faces, -1, 0, 0); infogrb = GrB_Matrix_setElement(*faces,  1, 1, 0); infogrb = GrB_Matrix_setElement(*faces, -1, 2, 0);

            sf[0] = 1;

            // Re-order points 2, 3, 4 so that i_pt3 is the "middle" corner point, others correctly oriented.
            neighs_23 = is_adjacent_ref_square(i_pt2, i_pt3);
            neighs_24 = is_adjacent_ref_square(i_pt2, i_pt4);
            if (!(neighs_23 && neighs_34)){ // if pt3 is not in the middle (otherwise nothing to do)
                i_pta = i_pt2; i_ptb = i_pt3; i_ptc = i_pt4;
                vertex_neighbors(i_pt3, &n1, &n2, &e1, &e2);
                if (neighs_23 && neighs_24){ // pt2 in the middle
                    if (n1 == i_pt2){
                        if (grid_edge_sign(grid, i_pt2, e1) < 0){ i_pt2 = i_ptc; i_pt3 = i_pta; i_pt4 = i_ptb; }
                        else                                     { i_pt2 = i_ptb; i_pt3 = i_pta; i_pt4 = i_ptc; }
                    } else {
                        if (grid_edge_sign(grid, i_pt2, e2) < 0){ i_pt2 = i_ptc; i_pt3 = i_pta; i_pt4 = i_ptb; }
                        else                                     { i_pt2 = i_ptb; i_pt3 = i_pta; i_pt4 = i_ptc; }
                    }
                } else { // pt4 in the middle
                    if (n1 == i_pt4){
                        if (grid_edge_sign(grid, i_pt4, e1) < 0){ i_pt2 = i_pta; i_pt3 = i_ptc; i_pt4 = i_ptb; }
                        else                                     { i_pt2 = i_ptb; i_pt3 = i_ptc; i_pt4 = i_pta; }
                    } else {
                        if (grid_edge_sign(grid, i_pt4, e2) < 0){ i_pt2 = i_pta; i_pt3 = i_ptc; i_pt4 = i_ptb; }
                        else                                     { i_pt2 = i_ptb; i_pt3 = i_ptc; i_pt4 = i_pta; }
                    }
                }
            }

            vertex_neighbors(i_pt2, &n1, &n2, &e1, &e2);
            if (n1 == i_pt3){ i_ptd = n2; e_12 = e2; e_23 = e1; }
            else            { i_ptd = n1; e_12 = e1; e_23 = e2; }
            vertex_neighbors(i_pt4, &n1, &n2, &e1, &e2);
            if (n1 == i_pt3){ e_34 = e1; e_14 = e2; }
            else            { e_14 = e1; e_34 = e2; }

            if (i_pt1 == i_pt2){
                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt3);
                pt33D_n = (Point3D){p2D->x,p2D->y,0.0}; pt33D_np1 = (Point3D){p2D->x,p2D->y,dt};
                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt4);
                pt43D_n = (Point3D){p2D->x,p2D->y,0.0}; pt43D_np1 = (Point3D){p2D->x,p2D->y,dt};

                lam_a = level_set_tn[i_pt3] / (level_set_tn[i_pt3] - level_set_tnp1[i_pt3]);
                lam_b = level_set_tn[i_pt4] / (level_set_tn[i_pt4] - level_set_tnp1[i_pt4]);

                verts[4] = pt33D_np1;
                verts[5] = lerp_pt3D(lam_a, pt33D_n, pt33D_np1);
                verts[6] = pt43D_np1;
                verts[7] = lerp_pt3D(lam_b, pt43D_n, pt43D_np1);

                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptd);
                pt23D_np1 = (Point3D){p2D->x,p2D->y,dt};
                lam_a = level_set_tnp1[i_pt1] / (level_set_tnp1[i_pt1] - level_set_tnp1[i_ptd]);
                lam_b = level_set_tnp1[i_pt4] / (level_set_tnp1[i_pt4] - level_set_tnp1[i_ptd]);

                verts[8] = lerp_pt3D(lam_a, pt13D_np1, pt23D_np1);
                verts[9] = lerp_pt3D(lam_b, pt43D_np1, pt23D_np1);

                find_0pt_Q1(level_set_tn, level_set_tnp1, &xi, &eta, &zeta);
                { Point2D pt2D_c = bilinear_pt2D(grid, xi, eta);
                  verts[10] = (Point3D){pt2D_c.x, pt2D_c.y, zeta*dt}; }

                infogrb = GrB_Matrix_setElement(*edges, -1, 3, 4);  infogrb = GrB_Matrix_setElement(*edges,  1, 4, 4);
                infogrb = GrB_Matrix_setElement(*edges, -1, 4, 5);  infogrb = GrB_Matrix_setElement(*edges,  1, 5, 5);
                infogrb = GrB_Matrix_setElement(*edges, -1, 4, 6);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 6);
                infogrb = GrB_Matrix_setElement(*edges, -1, 6, 7);  infogrb = GrB_Matrix_setElement(*edges,  1, 7, 7);
                infogrb = GrB_Matrix_setElement(*edges, -1, 3, 8);  infogrb = GrB_Matrix_setElement(*edges,  1, 8, 8);
                infogrb = GrB_Matrix_setElement(*edges, -1, 6, 9);  infogrb = GrB_Matrix_setElement(*edges,  1, 9, 9);
                infogrb = GrB_Matrix_setElement(*edges, -1, 8, 10); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 10);
                infogrb = GrB_Matrix_setElement(*edges, -1, 1, 11); infogrb = GrB_Matrix_setElement(*edges,  1, 8, 11);
                infogrb = GrB_Matrix_setElement(*edges, -1, 2, 12); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 12);
                infogrb = GrB_Matrix_setElement(*edges, -1, 5, 13); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 13);
                infogrb = GrB_Matrix_setElement(*edges, -1, 7, 14); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 14);
                infogrb = GrB_Matrix_setElement(*edges, -1, 1, 15); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 15);
                infogrb = GrB_Matrix_setElement(*edges, -1, 2, 16); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 16);
                infogrb = GrB_Matrix_setElement(*edges, -1, 5, 17); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 17);
                infogrb = GrB_Matrix_setElement(*edges, -1, 7, 18); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 18);
                infogrb = GrB_Matrix_setElement(*edges, -1, 8, 19); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 19);
                infogrb = GrB_Matrix_setElement(*edges, -1, 9, 20); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 20);

                infogrb = GrB_Matrix_setElement(*faces, -1, 4, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 6, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 8, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 9, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 10, 1);
                infogrb = GrB_Matrix_setElement(*faces,  1, 0, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 3, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 8, 2); infogrb = GrB_Matrix_setElement(*faces,  1, 11, 2);
                infogrb = GrB_Matrix_setElement(*faces, -1, 1, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 3, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 4, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 5, 3); infogrb = GrB_Matrix_setElement(*faces, -1, 12, 3);
                infogrb = GrB_Matrix_setElement(*faces, -1, 5, 4); infogrb = GrB_Matrix_setElement(*faces,  1, 6, 4); infogrb = GrB_Matrix_setElement(*faces,  1, 7, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 13, 4);
                infogrb = GrB_Matrix_setElement(*faces, -1, 7, 5); infogrb = GrB_Matrix_setElement(*faces,  1, 9, 5); infogrb = GrB_Matrix_setElement(*faces, -1, 14, 5);
                infogrb = GrB_Matrix_setElement(*faces,  1, 2, 6); infogrb = GrB_Matrix_setElement(*faces, -1, 15, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 16, 6);
                infogrb = GrB_Matrix_setElement(*faces,  1, 12, 7); infogrb = GrB_Matrix_setElement(*faces, -1, 16, 7); infogrb = GrB_Matrix_setElement(*faces,  1, 17, 7);
                infogrb = GrB_Matrix_setElement(*faces,  1, 13, 8); infogrb = GrB_Matrix_setElement(*faces, -1, 17, 8); infogrb = GrB_Matrix_setElement(*faces,  1, 18, 8);
                infogrb = GrB_Matrix_setElement(*faces,  1, 14, 9); infogrb = GrB_Matrix_setElement(*faces, -1, 18, 9); infogrb = GrB_Matrix_setElement(*faces,  1, 20, 9);
                infogrb = GrB_Matrix_setElement(*faces, -1, 10, 10); infogrb = GrB_Matrix_setElement(*faces,  1, 19, 10); infogrb = GrB_Matrix_setElement(*faces, -1, 20, 10);
                infogrb = GrB_Matrix_setElement(*faces, -1, 11, 11); infogrb = GrB_Matrix_setElement(*faces,  1, 15, 11); infogrb = GrB_Matrix_setElement(*faces, -1, 19, 11);

                sf[1] = 2;
                sf[2] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_12);
                sf[3] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_23);
                sf[4] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_34);
                sf[5] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_14);
                for (k = 6; k < 12; k++) sf[k] = -1;

            } else if (i_pt1 == i_pt3){
                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt2);
                pt23D_n = (Point3D){p2D->x,p2D->y,0.0}; pt23D_np1 = (Point3D){p2D->x,p2D->y,dt};
                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt4);
                pt43D_n = (Point3D){p2D->x,p2D->y,0.0}; pt43D_np1 = (Point3D){p2D->x,p2D->y,dt};

                lam_a = level_set_tn[i_pt2] / (level_set_tn[i_pt2] - level_set_tnp1[i_pt2]);
                lam_b = level_set_tn[i_pt4] / (level_set_tn[i_pt4] - level_set_tnp1[i_pt4]);

                verts[4] = pt23D_np1;
                verts[5] = lerp_pt3D(lam_a, pt23D_n, pt23D_np1);
                verts[6] = pt43D_np1;
                verts[7] = lerp_pt3D(lam_b, pt43D_n, pt43D_np1);

                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptd);
                pt33D_np1 = (Point3D){p2D->x,p2D->y,dt};
                lam_a = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tnp1[i_ptd]);
                lam_b = level_set_tnp1[i_pt4] / (level_set_tnp1[i_pt4] - level_set_tnp1[i_ptd]);

                verts[8] = lerp_pt3D(lam_a, pt23D_np1, pt33D_np1);
                verts[9] = lerp_pt3D(lam_b, pt43D_np1, pt33D_np1);

                find_0pt_Q1(level_set_tn, level_set_tnp1, &xi, &eta, &zeta);
                { Point2D pt2D_c = bilinear_pt2D(grid, xi, eta);
                  verts[10] = (Point3D){pt2D_c.x, pt2D_c.y, zeta*dt}; }

                infogrb = GrB_Matrix_setElement(*edges, -1, 3, 4);  infogrb = GrB_Matrix_setElement(*edges,  1, 4, 4);
                infogrb = GrB_Matrix_setElement(*edges, -1, 4, 5);  infogrb = GrB_Matrix_setElement(*edges,  1, 5, 5);
                infogrb = GrB_Matrix_setElement(*edges, -1, 3, 6);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 6);
                infogrb = GrB_Matrix_setElement(*edges, -1, 6, 7);  infogrb = GrB_Matrix_setElement(*edges,  1, 7, 7);
                infogrb = GrB_Matrix_setElement(*edges, -1, 4, 8);  infogrb = GrB_Matrix_setElement(*edges,  1, 8, 8);
                infogrb = GrB_Matrix_setElement(*edges, -1, 6, 9);  infogrb = GrB_Matrix_setElement(*edges,  1, 9, 9);
                infogrb = GrB_Matrix_setElement(*edges, -1, 8, 10); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 10);
                infogrb = GrB_Matrix_setElement(*edges, -1, 1, 11); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 11);
                infogrb = GrB_Matrix_setElement(*edges, -1, 2, 12); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 12);
                infogrb = GrB_Matrix_setElement(*edges, -1, 5, 13); infogrb = GrB_Matrix_setElement(*edges,  1, 8, 13);
                infogrb = GrB_Matrix_setElement(*edges, -1, 7, 14); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 14);
                infogrb = GrB_Matrix_setElement(*edges, -1, 1, 15); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 15);
                infogrb = GrB_Matrix_setElement(*edges, -1, 2, 16); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 16);
                infogrb = GrB_Matrix_setElement(*edges, -1, 5, 17); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 17);
                infogrb = GrB_Matrix_setElement(*edges, -1, 7, 18); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 18);
                infogrb = GrB_Matrix_setElement(*edges, -1, 8, 19); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 19);
                infogrb = GrB_Matrix_setElement(*edges, -1, 9, 20); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 20);

                infogrb = GrB_Matrix_setElement(*faces,  1, 4, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 6, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 8, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 9, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 10, 1);
                infogrb = GrB_Matrix_setElement(*faces,  1, 0, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 3, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 4, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 5, 2); infogrb = GrB_Matrix_setElement(*faces,  1, 11, 2);
                infogrb = GrB_Matrix_setElement(*faces, -1, 1, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 3, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 6, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 7, 3); infogrb = GrB_Matrix_setElement(*faces, -1, 12, 3);
                infogrb = GrB_Matrix_setElement(*faces, -1, 7, 4); infogrb = GrB_Matrix_setElement(*faces,  1, 9, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 14, 4);
                infogrb = GrB_Matrix_setElement(*faces,  1, 5, 5); infogrb = GrB_Matrix_setElement(*faces, -1, 8, 5); infogrb = GrB_Matrix_setElement(*faces,  1, 13, 5);
                infogrb = GrB_Matrix_setElement(*faces,  1, 2, 6); infogrb = GrB_Matrix_setElement(*faces, -1, 15, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 16, 6);
                infogrb = GrB_Matrix_setElement(*faces,  1, 12, 7); infogrb = GrB_Matrix_setElement(*faces, -1, 16, 7); infogrb = GrB_Matrix_setElement(*faces,  1, 18, 7);
                infogrb = GrB_Matrix_setElement(*faces, -1, 10, 8); infogrb = GrB_Matrix_setElement(*faces,  1, 18, 8); infogrb = GrB_Matrix_setElement(*faces, -1, 20, 8);
                infogrb = GrB_Matrix_setElement(*faces,  1, 14, 9); infogrb = GrB_Matrix_setElement(*faces, -1, 17, 9); infogrb = GrB_Matrix_setElement(*faces,  1, 20, 9);
                infogrb = GrB_Matrix_setElement(*faces, -1, 11, 10); infogrb = GrB_Matrix_setElement(*faces,  1, 17, 10); infogrb = GrB_Matrix_setElement(*faces, -1, 19, 10);
                infogrb = GrB_Matrix_setElement(*faces, -1, 13, 11); infogrb = GrB_Matrix_setElement(*faces,  1, 15, 11); infogrb = GrB_Matrix_setElement(*faces,  1, 19, 11);

                sf[1] = 2;
                sf[2] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_23);
                sf[3] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_34);
                sf[4] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_14);
                sf[5] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_12);
                for (k = 6; k < 12; k++) sf[k] = -1;

            } else if (i_pt1 == i_pt4){
                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt2);
                pt23D_n = (Point3D){p2D->x,p2D->y,0.0}; pt23D_np1 = (Point3D){p2D->x,p2D->y,dt};
                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt3);
                pt33D_n = (Point3D){p2D->x,p2D->y,0.0}; pt33D_np1 = (Point3D){p2D->x,p2D->y,dt};

                lam_a = level_set_tn[i_pt3] / (level_set_tn[i_pt3] - level_set_tnp1[i_pt3]);
                lam_b = level_set_tn[i_pt2] / (level_set_tn[i_pt2] - level_set_tnp1[i_pt2]);

                verts[4] = pt23D_np1;
                verts[5] = lerp_pt3D(lam_b, pt23D_n, pt23D_np1);
                verts[6] = pt33D_np1;
                verts[7] = lerp_pt3D(lam_a, pt33D_n, pt33D_np1);

                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptd);
                pt43D_np1 = (Point3D){p2D->x,p2D->y,dt};
                lam_a = level_set_tnp1[i_pt1] / (level_set_tnp1[i_pt1] - level_set_tnp1[i_ptd]);
                lam_b = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tnp1[i_ptd]);

                verts[8] = lerp_pt3D(lam_a, pt13D_np1, pt43D_np1);
                verts[9] = lerp_pt3D(lam_b, pt23D_np1, pt43D_np1);

                find_0pt_Q1(level_set_tn, level_set_tnp1, &xi, &eta, &zeta);
                { Point2D pt2D_c = bilinear_pt2D(grid, xi, eta);
                  verts[10] = (Point3D){pt2D_c.x, pt2D_c.y, zeta*dt}; }

                infogrb = GrB_Matrix_setElement(*edges, -1, 4, 4);  infogrb = GrB_Matrix_setElement(*edges,  1, 9, 4);
                infogrb = GrB_Matrix_setElement(*edges, -1, 4, 5);  infogrb = GrB_Matrix_setElement(*edges,  1, 5, 5);
                infogrb = GrB_Matrix_setElement(*edges, -1, 4, 6);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 6);
                infogrb = GrB_Matrix_setElement(*edges, -1, 6, 7);  infogrb = GrB_Matrix_setElement(*edges,  1, 7, 7);
                infogrb = GrB_Matrix_setElement(*edges, -1, 3, 8);  infogrb = GrB_Matrix_setElement(*edges,  1, 8, 8);
                infogrb = GrB_Matrix_setElement(*edges, -1, 3, 9);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 9);
                infogrb = GrB_Matrix_setElement(*edges, -1, 8, 10); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 10);
                infogrb = GrB_Matrix_setElement(*edges, -1, 1, 11); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 11);
                infogrb = GrB_Matrix_setElement(*edges, -1, 2, 12); infogrb = GrB_Matrix_setElement(*edges,  1, 8, 12);
                infogrb = GrB_Matrix_setElement(*edges, -1, 5, 13); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 13);
                infogrb = GrB_Matrix_setElement(*edges, -1, 5, 14); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 14);
                infogrb = GrB_Matrix_setElement(*edges, -1, 1, 15); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 15);
                infogrb = GrB_Matrix_setElement(*edges, -1, 2, 16); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 16);
                infogrb = GrB_Matrix_setElement(*edges, -1, 5, 17); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 17);
                infogrb = GrB_Matrix_setElement(*edges, -1, 7, 18); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 18);
                infogrb = GrB_Matrix_setElement(*edges, -1, 8, 19); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 19);
                infogrb = GrB_Matrix_setElement(*edges, -1, 9, 20); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 20);

                infogrb = GrB_Matrix_setElement(*faces,  1, 4, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 6, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 8, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 9, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 10, 1);
                infogrb = GrB_Matrix_setElement(*faces, -1, 5, 2); infogrb = GrB_Matrix_setElement(*faces,  1, 6, 2); infogrb = GrB_Matrix_setElement(*faces,  1, 7, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 13, 2);
                infogrb = GrB_Matrix_setElement(*faces,  1, 0, 3); infogrb = GrB_Matrix_setElement(*faces, -1, 3, 3); infogrb = GrB_Matrix_setElement(*faces, -1, 7, 3); infogrb = GrB_Matrix_setElement(*faces, -1, 9, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 11, 3);
                infogrb = GrB_Matrix_setElement(*faces, -1, 1, 4); infogrb = GrB_Matrix_setElement(*faces,  1, 3, 4); infogrb = GrB_Matrix_setElement(*faces,  1, 8, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 12, 4);
                infogrb = GrB_Matrix_setElement(*faces, -1, 4, 5); infogrb = GrB_Matrix_setElement(*faces,  1, 5, 5); infogrb = GrB_Matrix_setElement(*faces,  1, 14, 5);
                infogrb = GrB_Matrix_setElement(*faces,  1, 2, 6); infogrb = GrB_Matrix_setElement(*faces, -1, 15, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 16, 6);
                infogrb = GrB_Matrix_setElement(*faces,  1, 12, 7); infogrb = GrB_Matrix_setElement(*faces, -1, 16, 7); infogrb = GrB_Matrix_setElement(*faces,  1, 19, 7);
                infogrb = GrB_Matrix_setElement(*faces, -1, 10, 8); infogrb = GrB_Matrix_setElement(*faces,  1, 19, 8); infogrb = GrB_Matrix_setElement(*faces, -1, 20, 8);
                infogrb = GrB_Matrix_setElement(*faces, -1, 14, 9); infogrb = GrB_Matrix_setElement(*faces,  1, 18, 9); infogrb = GrB_Matrix_setElement(*faces, -1, 20, 9);
                infogrb = GrB_Matrix_setElement(*faces,  1, 13, 10); infogrb = GrB_Matrix_setElement(*faces, -1, 17, 10); infogrb = GrB_Matrix_setElement(*faces,  1, 18, 10);
                infogrb = GrB_Matrix_setElement(*faces, -1, 11, 11); infogrb = GrB_Matrix_setElement(*faces,  1, 15, 11); infogrb = GrB_Matrix_setElement(*faces, -1, 17, 11);

                sf[1] = 2;
                sf[2] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_23);
                sf[3] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_34);
                sf[4] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_14);
                sf[5] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_12);
                for (k = 6; k < 12; k++) sf[k] = -1;

            } else {
                printf("Error: case 4 while it should be case 2.\n");
                return;
            }

        } else if (nb_tn == 3){
            int n1, n2, e1, e2;
            int i_pta, i_ptb, i_ptc, i_ptd;
            int e_12, e_13, e_23, e_34, e_14;
            int neighs_23, neighs_13;
            Point2D *p2D;
            Point3D pt13D_n, pt13D_np1, pt23D_n, pt23D_np1, pt33D_n, pt33D_np1, pt43D_n, pt43D_np1;
            my_real_c lam_a, lam_b, xi, eta, zeta;

            vertex_neighbors(i_pt4, &n1, &n2, &e1, &e2);
            if (grid_edge_sign(grid, n1, e1) < 0){ i_pta = n1; i_ptb = n2; e_12 = e1; e_13 = e2; }
            else                                  { i_ptb = n1; i_pta = n2; e_13 = e1; e_12 = e2; }

            lam_a = level_set_tnp1[i_pt4] / (level_set_tnp1[i_pt4] - level_set_tnp1[i_pta]);
            lam_b = level_set_tnp1[i_pt4] / (level_set_tnp1[i_pt4] - level_set_tnp1[i_ptb]);

            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt4);
            pt43D_n = (Point3D){p2D->x,p2D->y,0.0}; pt43D_np1 = (Point3D){p2D->x,p2D->y,dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pta);
            pt23D_n = (Point3D){p2D->x,p2D->y,0.0};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptb);
            pt33D_n = (Point3D){p2D->x,p2D->y,0.0};
            // (note: only the t^{n+1} coordinates of pta/ptb are actually needed below)
            pt23D_np1 = (Point3D){pt23D_n.x,pt23D_n.y,dt};
            pt33D_np1 = (Point3D){pt33D_n.x,pt33D_n.y,dt};

            verts[0] = pt43D_np1;
            verts[1] = lerp_pt3D(lam_a, pt43D_np1, pt23D_np1);
            verts[2] = lerp_pt3D(lam_b, pt43D_np1, pt33D_np1);
            verts[3] = pt43D_n;

            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 0); infogrb = GrB_Matrix_setElement(*edges,  1, 1, 0);
            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 1); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 1);
            infogrb = GrB_Matrix_setElement(*edges, -1, 1, 2); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 2);
            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 3); infogrb = GrB_Matrix_setElement(*edges,  1, 3, 3);

            // faces: net sign = -1
            infogrb = GrB_Matrix_setElement(*faces,  1, 0, 0); infogrb = GrB_Matrix_setElement(*faces, -1, 1, 0); infogrb = GrB_Matrix_setElement(*faces,  1, 2, 0);

            sf[0] = 2;

            // Re-order points 1, 2, 3 so that i_pt2 is the "middle" corner point, others correctly oriented.
            neighs_23 = is_adjacent_ref_square(i_pt2, i_pt3);
            neighs_13 = is_adjacent_ref_square(i_pt1, i_pt3);
            if (!(neighs_12 && neighs_23)){ // if pt2 is not in the middle (otherwise nothing to do)
                i_pta = i_pt1; i_ptb = i_pt2; i_ptc = i_pt3;
                vertex_neighbors(i_pt2, &n1, &n2, &e1, &e2);
                if (neighs_12 && neighs_13){ // pt1 in the middle
                    if (n1 == i_pt1){
                        if (grid_edge_sign(grid, i_pt1, e1) < 0){ i_pt1 = i_ptc; i_pt2 = i_pta; i_pt3 = i_ptb; }
                        else                                     { i_pt1 = i_ptb; i_pt2 = i_pta; i_pt3 = i_ptc; }
                    } else {
                        if (grid_edge_sign(grid, i_pt1, e2) < 0){ i_pt1 = i_ptc; i_pt2 = i_pta; i_pt3 = i_ptb; }
                        else                                     { i_pt1 = i_ptb; i_pt2 = i_pta; i_pt3 = i_ptc; }
                    }
                } else { // pt3 in the middle
                    if (n1 == i_pt3){
                        if (grid_edge_sign(grid, i_pt3, e1) < 0){ i_pt1 = i_pta; i_pt2 = i_ptc; i_pt3 = i_ptb; }
                        else                                     { i_pt1 = i_ptb; i_pt2 = i_ptc; i_pt3 = i_pta; }
                    } else {
                        if (grid_edge_sign(grid, i_pt3, e2) < 0){ i_pt1 = i_pta; i_pt2 = i_ptc; i_pt3 = i_ptb; }
                        else                                     { i_pt1 = i_ptb; i_pt2 = i_ptc; i_pt3 = i_pta; }
                    }
                }
            }

            vertex_neighbors(i_pt1, &n1, &n2, &e1, &e2);
            if (n1 == i_pt2){ i_ptd = n2; e_14 = e2; e_12 = e1; }
            else            { i_ptd = n1; e_14 = e1; e_12 = e2; }
            vertex_neighbors(i_pt3, &n1, &n2, &e1, &e2);
            if (n1 == i_pt2){ e_23 = e1; e_34 = e2; }
            else            { e_34 = e1; e_23 = e2; }

            if (i_pt4 == i_pt1){
                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt2);
                pt23D_n = (Point3D){p2D->x,p2D->y,0.0}; pt23D_np1 = (Point3D){p2D->x,p2D->y,dt};
                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt3);
                pt33D_n = (Point3D){p2D->x,p2D->y,0.0}; pt33D_np1 = (Point3D){p2D->x,p2D->y,dt};

                lam_a = level_set_tn[i_pt2] / (level_set_tn[i_pt2] - level_set_tnp1[i_pt2]);
                lam_b = level_set_tn[i_pt3] / (level_set_tn[i_pt3] - level_set_tnp1[i_pt3]);

                verts[4] = pt23D_n;
                verts[5] = lerp_pt3D(lam_a, pt23D_n, pt23D_np1);
                verts[6] = pt33D_n;
                verts[7] = lerp_pt3D(lam_b, pt33D_n, pt33D_np1);

                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptd);
                pt13D_n = (Point3D){p2D->x,p2D->y,0.0};
                lam_a = level_set_tn[i_pt4] / (level_set_tn[i_pt4] - level_set_tn[i_ptd]);
                lam_b = level_set_tn[i_pt3] / (level_set_tn[i_pt3] - level_set_tn[i_ptd]);

                verts[8] = lerp_pt3D(lam_a, pt43D_n, pt13D_n);
                verts[9] = lerp_pt3D(lam_b, pt33D_n, pt13D_n);

                find_0pt_Q1(level_set_tn, level_set_tnp1, &xi, &eta, &zeta);
                { Point2D pt2D_c = bilinear_pt2D(grid, xi, eta);
                  verts[10] = (Point3D){pt2D_c.x, pt2D_c.y, zeta*dt}; }

                infogrb = GrB_Matrix_setElement(*edges, -1, 3, 4);  infogrb = GrB_Matrix_setElement(*edges,  1, 4, 4);
                infogrb = GrB_Matrix_setElement(*edges, -1, 4, 5);  infogrb = GrB_Matrix_setElement(*edges,  1, 5, 5);
                infogrb = GrB_Matrix_setElement(*edges, -1, 4, 6);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 6);
                infogrb = GrB_Matrix_setElement(*edges, -1, 6, 7);  infogrb = GrB_Matrix_setElement(*edges,  1, 7, 7);
                infogrb = GrB_Matrix_setElement(*edges, -1, 3, 8);  infogrb = GrB_Matrix_setElement(*edges,  1, 8, 8);
                infogrb = GrB_Matrix_setElement(*edges, -1, 6, 9);  infogrb = GrB_Matrix_setElement(*edges,  1, 9, 9);
                infogrb = GrB_Matrix_setElement(*edges, -1, 8, 10); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 10);
                infogrb = GrB_Matrix_setElement(*edges, -1, 1, 11); infogrb = GrB_Matrix_setElement(*edges,  1, 8, 11);
                infogrb = GrB_Matrix_setElement(*edges, -1, 2, 12); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 12);
                infogrb = GrB_Matrix_setElement(*edges, -1, 5, 13); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 13);
                infogrb = GrB_Matrix_setElement(*edges, -1, 7, 14); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 14);
                infogrb = GrB_Matrix_setElement(*edges, -1, 1, 15); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 15);
                infogrb = GrB_Matrix_setElement(*edges, -1, 2, 16); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 16);
                infogrb = GrB_Matrix_setElement(*edges, -1, 5, 17); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 17);
                infogrb = GrB_Matrix_setElement(*edges, -1, 7, 18); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 18);
                infogrb = GrB_Matrix_setElement(*edges, -1, 8, 19); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 19);
                infogrb = GrB_Matrix_setElement(*edges, -1, 9, 20); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 20);

                infogrb = GrB_Matrix_setElement(*faces, -1, 4, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 6, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 8, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 9, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 10, 1);
                infogrb = GrB_Matrix_setElement(*faces, -1, 0, 2); infogrb = GrB_Matrix_setElement(*faces,  1, 3, 2); infogrb = GrB_Matrix_setElement(*faces,  1, 8, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 11, 2);
                infogrb = GrB_Matrix_setElement(*faces,  1, 1, 3); infogrb = GrB_Matrix_setElement(*faces, -1, 3, 3); infogrb = GrB_Matrix_setElement(*faces, -1, 4, 3); infogrb = GrB_Matrix_setElement(*faces, -1, 5, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 12, 3);
                infogrb = GrB_Matrix_setElement(*faces,  1, 5, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 6, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 7, 4); infogrb = GrB_Matrix_setElement(*faces,  1, 13, 4);
                infogrb = GrB_Matrix_setElement(*faces,  1, 7, 5); infogrb = GrB_Matrix_setElement(*faces, -1, 9, 5); infogrb = GrB_Matrix_setElement(*faces,  1, 14, 5);
                infogrb = GrB_Matrix_setElement(*faces, -1, 2, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 15, 6); infogrb = GrB_Matrix_setElement(*faces, -1, 16, 6);
                infogrb = GrB_Matrix_setElement(*faces, -1, 12, 7); infogrb = GrB_Matrix_setElement(*faces,  1, 16, 7); infogrb = GrB_Matrix_setElement(*faces, -1, 17, 7);
                infogrb = GrB_Matrix_setElement(*faces, -1, 13, 8); infogrb = GrB_Matrix_setElement(*faces,  1, 17, 8); infogrb = GrB_Matrix_setElement(*faces, -1, 18, 8);
                infogrb = GrB_Matrix_setElement(*faces, -1, 14, 9); infogrb = GrB_Matrix_setElement(*faces,  1, 18, 9); infogrb = GrB_Matrix_setElement(*faces, -1, 20, 9);
                infogrb = GrB_Matrix_setElement(*faces,  1, 10, 10); infogrb = GrB_Matrix_setElement(*faces, -1, 19, 10); infogrb = GrB_Matrix_setElement(*faces,  1, 20, 10);
                infogrb = GrB_Matrix_setElement(*faces,  1, 11, 11); infogrb = GrB_Matrix_setElement(*faces, -1, 15, 11); infogrb = GrB_Matrix_setElement(*faces,  1, 19, 11);

                sf[1] = 1;
                sf[2] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_14);
                sf[3] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_12);
                sf[4] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_23);
                sf[5] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_34);
                for (k = 6; k < 12; k++) sf[k] = -1;

            } else if (i_pt4 == i_pt2){
                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt1);
                pt13D_n = (Point3D){p2D->x,p2D->y,0.0}; pt13D_np1 = (Point3D){p2D->x,p2D->y,dt};
                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt3);
                pt33D_n = (Point3D){p2D->x,p2D->y,0.0}; pt33D_np1 = (Point3D){p2D->x,p2D->y,dt};

                lam_a = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tnp1[i_pt1]);
                lam_b = level_set_tn[i_pt3] / (level_set_tn[i_pt3] - level_set_tnp1[i_pt3]);

                verts[4] = pt13D_n;
                verts[5] = lerp_pt3D(lam_a, pt13D_n, pt13D_np1);
                verts[6] = pt33D_n;
                verts[7] = lerp_pt3D(lam_b, pt33D_n, pt33D_np1);

                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptd);
                pt23D_n = (Point3D){p2D->x,p2D->y,0.0};
                lam_a = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tn[i_ptd]);
                lam_b = level_set_tn[i_pt3] / (level_set_tn[i_pt3] - level_set_tn[i_ptd]);

                verts[8] = lerp_pt3D(lam_a, pt13D_n, pt23D_n);
                verts[9] = lerp_pt3D(lam_b, pt33D_n, pt23D_n);

                find_0pt_Q1(level_set_tn, level_set_tnp1, &xi, &eta, &zeta);
                { Point2D pt2D_c = bilinear_pt2D(grid, xi, eta);
                  verts[10] = (Point3D){pt2D_c.x, pt2D_c.y, zeta*dt}; }

                infogrb = GrB_Matrix_setElement(*edges, -1, 3, 4);  infogrb = GrB_Matrix_setElement(*edges,  1, 4, 4);
                infogrb = GrB_Matrix_setElement(*edges, -1, 4, 5);  infogrb = GrB_Matrix_setElement(*edges,  1, 5, 5);
                infogrb = GrB_Matrix_setElement(*edges, -1, 3, 6);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 6);
                infogrb = GrB_Matrix_setElement(*edges, -1, 6, 7);  infogrb = GrB_Matrix_setElement(*edges,  1, 7, 7);
                infogrb = GrB_Matrix_setElement(*edges, -1, 4, 8);  infogrb = GrB_Matrix_setElement(*edges,  1, 8, 8);
                infogrb = GrB_Matrix_setElement(*edges, -1, 6, 9);  infogrb = GrB_Matrix_setElement(*edges,  1, 9, 9);
                infogrb = GrB_Matrix_setElement(*edges, -1, 8, 10); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 10);
                infogrb = GrB_Matrix_setElement(*edges, -1, 1, 11); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 11);
                infogrb = GrB_Matrix_setElement(*edges, -1, 2, 12); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 12);
                infogrb = GrB_Matrix_setElement(*edges, -1, 5, 13); infogrb = GrB_Matrix_setElement(*edges,  1, 8, 13);
                infogrb = GrB_Matrix_setElement(*edges, -1, 7, 14); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 14);
                infogrb = GrB_Matrix_setElement(*edges, -1, 1, 15); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 15);
                infogrb = GrB_Matrix_setElement(*edges, -1, 2, 16); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 16);
                infogrb = GrB_Matrix_setElement(*edges, -1, 5, 17); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 17);
                infogrb = GrB_Matrix_setElement(*edges, -1, 7, 18); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 18);
                infogrb = GrB_Matrix_setElement(*edges, -1, 8, 19); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 19);
                infogrb = GrB_Matrix_setElement(*edges, -1, 9, 20); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 20);

                infogrb = GrB_Matrix_setElement(*faces,  1, 4, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 6, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 8, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 9, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 10, 1);
                infogrb = GrB_Matrix_setElement(*faces,  1, 0, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 3, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 4, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 5, 2); infogrb = GrB_Matrix_setElement(*faces,  1, 11, 2);
                infogrb = GrB_Matrix_setElement(*faces, -1, 1, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 3, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 6, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 7, 3); infogrb = GrB_Matrix_setElement(*faces, -1, 12, 3);
                infogrb = GrB_Matrix_setElement(*faces, -1, 7, 4); infogrb = GrB_Matrix_setElement(*faces,  1, 9, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 14, 4);
                infogrb = GrB_Matrix_setElement(*faces,  1, 5, 5); infogrb = GrB_Matrix_setElement(*faces, -1, 8, 5); infogrb = GrB_Matrix_setElement(*faces,  1, 13, 5);
                infogrb = GrB_Matrix_setElement(*faces,  1, 2, 6); infogrb = GrB_Matrix_setElement(*faces, -1, 15, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 16, 6);
                infogrb = GrB_Matrix_setElement(*faces,  1, 12, 7); infogrb = GrB_Matrix_setElement(*faces, -1, 16, 7); infogrb = GrB_Matrix_setElement(*faces,  1, 18, 7);
                infogrb = GrB_Matrix_setElement(*faces, -1, 10, 8); infogrb = GrB_Matrix_setElement(*faces,  1, 18, 8); infogrb = GrB_Matrix_setElement(*faces, -1, 20, 8);
                infogrb = GrB_Matrix_setElement(*faces,  1, 14, 9); infogrb = GrB_Matrix_setElement(*faces, -1, 17, 9); infogrb = GrB_Matrix_setElement(*faces,  1, 20, 9);
                infogrb = GrB_Matrix_setElement(*faces, -1, 11, 10); infogrb = GrB_Matrix_setElement(*faces,  1, 17, 10); infogrb = GrB_Matrix_setElement(*faces, -1, 19, 10);
                infogrb = GrB_Matrix_setElement(*faces, -1, 13, 11); infogrb = GrB_Matrix_setElement(*faces,  1, 15, 11); infogrb = GrB_Matrix_setElement(*faces,  1, 19, 11);

                sf[1] = 1;
                sf[2] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_12);
                sf[3] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_23);
                sf[4] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_34);
                sf[5] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_14);
                for (k = 6; k < 12; k++) sf[k] = -1;

            } else if (i_pt4 == i_pt3){
                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt1);
                pt13D_n = (Point3D){p2D->x,p2D->y,0.0}; pt13D_np1 = (Point3D){p2D->x,p2D->y,dt};
                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt2);
                pt23D_n = (Point3D){p2D->x,p2D->y,0.0}; pt23D_np1 = (Point3D){p2D->x,p2D->y,dt};

                lam_a = level_set_tn[i_pt2] / (level_set_tn[i_pt2] - level_set_tnp1[i_pt2]);
                lam_b = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tnp1[i_pt1]);

                verts[4] = pt13D_n;
                verts[5] = lerp_pt3D(lam_b, pt13D_n, pt13D_np1);
                verts[6] = pt23D_n;
                verts[7] = lerp_pt3D(lam_a, pt23D_n, pt23D_np1);

                p2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptd);
                pt33D_n = (Point3D){p2D->x,p2D->y,0.0};
                lam_a = level_set_tn[i_pt4] / (level_set_tn[i_pt4] - level_set_tn[i_ptd]);
                lam_b = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tn[i_ptd]);

                verts[8] = lerp_pt3D(lam_a, pt43D_n, pt33D_n);
                verts[9] = lerp_pt3D(lam_b, pt13D_n, pt33D_n);

                find_0pt_Q1(level_set_tn, level_set_tnp1, &xi, &eta, &zeta);
                { Point2D pt2D_c = bilinear_pt2D(grid, xi, eta);
                  verts[10] = (Point3D){pt2D_c.x, pt2D_c.y, zeta*dt}; }

                infogrb = GrB_Matrix_setElement(*edges, -1, 4, 4);  infogrb = GrB_Matrix_setElement(*edges,  1, 9, 4);
                infogrb = GrB_Matrix_setElement(*edges, -1, 4, 5);  infogrb = GrB_Matrix_setElement(*edges,  1, 5, 5);
                infogrb = GrB_Matrix_setElement(*edges, -1, 4, 6);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 6);
                infogrb = GrB_Matrix_setElement(*edges, -1, 6, 7);  infogrb = GrB_Matrix_setElement(*edges,  1, 7, 7);
                infogrb = GrB_Matrix_setElement(*edges, -1, 3, 8);  infogrb = GrB_Matrix_setElement(*edges,  1, 8, 8);
                infogrb = GrB_Matrix_setElement(*edges, -1, 3, 9);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 9);
                infogrb = GrB_Matrix_setElement(*edges, -1, 8, 10); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 10);
                infogrb = GrB_Matrix_setElement(*edges, -1, 1, 11); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 11);
                infogrb = GrB_Matrix_setElement(*edges, -1, 2, 12); infogrb = GrB_Matrix_setElement(*edges,  1, 8, 12);
                infogrb = GrB_Matrix_setElement(*edges, -1, 5, 13); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 13);
                infogrb = GrB_Matrix_setElement(*edges, -1, 5, 14); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 14);
                infogrb = GrB_Matrix_setElement(*edges, -1, 1, 15); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 15);
                infogrb = GrB_Matrix_setElement(*edges, -1, 2, 16); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 16);
                infogrb = GrB_Matrix_setElement(*edges, -1, 5, 17); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 17);
                infogrb = GrB_Matrix_setElement(*edges, -1, 7, 18); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 18);
                infogrb = GrB_Matrix_setElement(*edges, -1, 8, 19); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 19);
                infogrb = GrB_Matrix_setElement(*edges, -1, 9, 20); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 20);

                infogrb = GrB_Matrix_setElement(*faces, -1, 4, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 6, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 8, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 9, 1); infogrb = GrB_Matrix_setElement(*faces,  1, 10, 1);
                infogrb = GrB_Matrix_setElement(*faces,  1, 5, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 6, 2); infogrb = GrB_Matrix_setElement(*faces, -1, 7, 2); infogrb = GrB_Matrix_setElement(*faces,  1, 13, 2);
                infogrb = GrB_Matrix_setElement(*faces, -1, 0, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 3, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 7, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 9, 3); infogrb = GrB_Matrix_setElement(*faces, -1, 11, 3);
                infogrb = GrB_Matrix_setElement(*faces,  1, 1, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 3, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 8, 4); infogrb = GrB_Matrix_setElement(*faces,  1, 12, 4);
                infogrb = GrB_Matrix_setElement(*faces,  1, 4, 5); infogrb = GrB_Matrix_setElement(*faces, -1, 5, 5); infogrb = GrB_Matrix_setElement(*faces, -1, 14, 5);
                infogrb = GrB_Matrix_setElement(*faces, -1, 2, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 15, 6); infogrb = GrB_Matrix_setElement(*faces, -1, 16, 6);
                infogrb = GrB_Matrix_setElement(*faces, -1, 12, 7); infogrb = GrB_Matrix_setElement(*faces,  1, 16, 7); infogrb = GrB_Matrix_setElement(*faces, -1, 19, 7);
                infogrb = GrB_Matrix_setElement(*faces,  1, 10, 8); infogrb = GrB_Matrix_setElement(*faces, -1, 19, 8); infogrb = GrB_Matrix_setElement(*faces,  1, 20, 8);
                infogrb = GrB_Matrix_setElement(*faces,  1, 14, 9); infogrb = GrB_Matrix_setElement(*faces, -1, 18, 9); infogrb = GrB_Matrix_setElement(*faces,  1, 20, 9);
                infogrb = GrB_Matrix_setElement(*faces, -1, 13, 10); infogrb = GrB_Matrix_setElement(*faces,  1, 17, 10); infogrb = GrB_Matrix_setElement(*faces, -1, 18, 10);
                infogrb = GrB_Matrix_setElement(*faces,  1, 11, 11); infogrb = GrB_Matrix_setElement(*faces, -1, 15, 11); infogrb = GrB_Matrix_setElement(*faces,  1, 17, 11);

                sf[1] = 1;
                sf[2] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_12);
                sf[3] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_23);
                sf[4] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_34);
                sf[5] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_14);
                for (k = 6; k < 12; k++) sf[k] = -1;

            } else {
                printf("Error: case 4 while it should be case 2.\n");
                return;
            }

        } else { // nb_tn == 2
            // We must identify the points: we will have i_ptb == i_ptc,
            // i_pta being the second point at t^n, i_ptd being the second point at t^{n+1}
            int n1, n2, e1, e2;
            int i_pta, i_ptb, i_ptc, i_ptd;
            int e_13, e_12, e_24, e_34;
            Point2D *p2D;
            Point3D pt13D_n, pt13D_np1, pt23D_n, pt23D_np1, pt33D_n, pt33D_np1, pt43D_n, pt43D_np1;
            my_real_c lam1, xi, eta, zeta;
            int face_sign;

            if (i_pt1 == i_pt3){
                i_pta = i_pt2; i_ptb = i_pt1; i_ptc = i_pt3; i_ptd = i_pt4;
            } else if (i_pt1 == i_pt4){
                i_pta = i_pt2; i_ptb = i_pt1; i_ptc = i_pt4; i_ptd = i_pt3;
            } else if (i_pt2 == i_pt3){
                i_pta = i_pt1; i_ptb = i_pt2; i_ptc = i_pt3; i_ptd = i_pt4;
            } else if (i_pt2 == i_pt4){
                i_pta = i_pt1; i_ptb = i_pt2; i_ptc = i_pt4; i_ptd = i_pt3;
            } else {
                printf("Error: a case has been forgotten.\n");
                return;
            }

            vertex_neighbors(i_pta, &n1, &n2, &e1, &e2);
            if (n1 != i_ptb){ i_pt1 = n1; e_13 = e1; e_12 = e2; }
            else            { i_pt1 = n2; e_13 = e2; e_12 = e1; }
            vertex_neighbors(i_ptb, &n1, &n2, &e1, &e2);
            if (n1 != i_pta){ i_pt2 = n1; e_24 = e1; }
            else            { i_pt2 = n2; e_24 = e2; }
            vertex_neighbors(i_pt1, &n1, &n2, &e1, &e2);
            if (n1 == i_pt2) e_34 = e1;
            else             e_34 = e2;

            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pta);
            pt13D_n = (Point3D){p2D->x,p2D->y,0.0}; pt13D_np1 = (Point3D){p2D->x,p2D->y,dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_ptb);
            pt23D_n = (Point3D){p2D->x,p2D->y,0.0}; pt23D_np1 = (Point3D){p2D->x,p2D->y,dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt1);
            pt33D_n = (Point3D){p2D->x,p2D->y,0.0}; pt33D_np1 = (Point3D){p2D->x,p2D->y,dt};
            p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt2);
            pt43D_n = (Point3D){p2D->x,p2D->y,0.0}; pt43D_np1 = (Point3D){p2D->x,p2D->y,dt};

            lam1 = level_set_tn[i_pta] / (level_set_tn[i_pta] - level_set_tn[i_pt1]);
            verts[0] = lerp_pt3D(lam1, pt13D_n, pt33D_n);
            verts[1] = pt13D_n;

            lam1 = level_set_tn[i_ptb] / (level_set_tn[i_ptb] - level_set_tn[i_pt2]);
            verts[2] = pt23D_n;
            verts[3] = lerp_pt3D(lam1, pt23D_n, pt43D_n);

            lam1 = level_set_tn[i_pta] / (level_set_tn[i_pta] - level_set_tnp1[i_pta]);
            verts[4] = lerp_pt3D(lam1, pt13D_n, pt13D_np1);

            lam1 = level_set_tnp1[i_ptb] / (level_set_tnp1[i_ptb] - level_set_tnp1[i_pta]);
            verts[5] = lerp_pt3D(lam1, pt13D_np1, pt23D_np1);
            verts[6] = pt23D_np1;

            lam1 = level_set_tnp1[i_pt2] / (level_set_tnp1[i_pt2] - level_set_tnp1[i_pt1]);
            verts[7] = pt43D_np1;
            verts[8] = lerp_pt3D(lam1, pt43D_np1, pt33D_np1);

            lam1 = level_set_tn[i_pt2] / (level_set_tn[i_pt2] - level_set_tnp1[i_pt2]);
            verts[9] = lerp_pt3D(lam1, pt43D_n, pt43D_np1);

            find_0pt_Q1(level_set_tn, level_set_tnp1, &xi, &eta, &zeta);
            { Point2D pt2D_c = bilinear_pt2D(grid, xi, eta);
              verts[10] = (Point3D){pt2D_c.x, pt2D_c.y, zeta*dt}; }

            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 0);  infogrb = GrB_Matrix_setElement(*edges,  1, 1, 0);
            infogrb = GrB_Matrix_setElement(*edges, -1, 1, 1);  infogrb = GrB_Matrix_setElement(*edges,  1, 2, 1);
            infogrb = GrB_Matrix_setElement(*edges, -1, 1, 2);  infogrb = GrB_Matrix_setElement(*edges,  1, 4, 2);
            infogrb = GrB_Matrix_setElement(*edges, -1, 2, 3);  infogrb = GrB_Matrix_setElement(*edges,  1, 3, 3);
            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 4);  infogrb = GrB_Matrix_setElement(*edges,  1, 3, 4);
            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 5);  infogrb = GrB_Matrix_setElement(*edges,  1, 4, 5);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 6);  infogrb = GrB_Matrix_setElement(*edges,  1, 5, 6);
            infogrb = GrB_Matrix_setElement(*edges, -1, 2, 7);  infogrb = GrB_Matrix_setElement(*edges,  1, 6, 7);
            infogrb = GrB_Matrix_setElement(*edges, -1, 3, 8);  infogrb = GrB_Matrix_setElement(*edges,  1, 9, 8);
            infogrb = GrB_Matrix_setElement(*edges, -1, 7, 9);  infogrb = GrB_Matrix_setElement(*edges,  1, 9, 9);
            infogrb = GrB_Matrix_setElement(*edges, -1, 8, 10); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 10);
            infogrb = GrB_Matrix_setElement(*edges, -1, 5, 11); infogrb = GrB_Matrix_setElement(*edges,  1, 6, 11);
            infogrb = GrB_Matrix_setElement(*edges, -1, 6, 12); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 12);
            infogrb = GrB_Matrix_setElement(*edges, -1, 7, 13); infogrb = GrB_Matrix_setElement(*edges,  1, 8, 13);
            infogrb = GrB_Matrix_setElement(*edges, -1, 5, 14); infogrb = GrB_Matrix_setElement(*edges,  1, 8, 14);
            infogrb = GrB_Matrix_setElement(*edges, -1, 0, 15); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 15);
            infogrb = GrB_Matrix_setElement(*edges, -1, 3, 16); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 16);
            infogrb = GrB_Matrix_setElement(*edges, -1, 4, 17); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 17);
            infogrb = GrB_Matrix_setElement(*edges, -1, 5, 18); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 18);
            infogrb = GrB_Matrix_setElement(*edges, -1, 8, 19); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 19);
            infogrb = GrB_Matrix_setElement(*edges, -1, 9, 20); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 20);

            face_sign = ((i_pta == i_ptb + 1) || (i_pta == 0 && i_ptb == 3)) ? 1 : -1;

            infogrb = GrB_Matrix_setElement(*faces, -face_sign, 0, 0); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 1, 0); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 3, 0); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 4, 0);
            infogrb = GrB_Matrix_setElement(*faces,  face_sign, 0, 1); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 2, 1); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 5, 1);
            infogrb = GrB_Matrix_setElement(*faces, -face_sign, 2, 2); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 1, 2); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 6, 2); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 7, 2); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 11, 2);
            infogrb = GrB_Matrix_setElement(*faces,  face_sign, 3, 3); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 7, 3); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 8, 3); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 9, 3); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 12, 3);
            infogrb = GrB_Matrix_setElement(*faces,  face_sign, 9, 4); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 10, 4); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 13, 4);
            infogrb = GrB_Matrix_setElement(*faces,  face_sign, 11, 5); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 12, 5); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 13, 5); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 14, 5);
            infogrb = GrB_Matrix_setElement(*faces, -face_sign, 4, 6); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 15, 6); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 16, 6);
            infogrb = GrB_Matrix_setElement(*faces, -face_sign, 8, 7); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 16, 7); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 20, 7);
            infogrb = GrB_Matrix_setElement(*faces,  face_sign, 10, 8); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 19, 8); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 20, 8);
            infogrb = GrB_Matrix_setElement(*faces,  face_sign, 14, 9); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 18, 9); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 19, 9);
            infogrb = GrB_Matrix_setElement(*faces,  face_sign, 6, 10); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 17, 10); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 18, 10);
            infogrb = GrB_Matrix_setElement(*faces,  face_sign, 5, 11); infogrb = GrB_Matrix_setElement(*faces, -face_sign, 15, 11); infogrb = GrB_Matrix_setElement(*faces,  face_sign, 17, 11);

            sf[0] = 1;
            sf[1] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_13);
            sf[2] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_12);
            sf[3] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_24);
            sf[4] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_34);
            sf[5] = 2;
            for (k = 6; k < 12; k++) sf[k] = -1;
        }

        for (k = 0; k < 11; k++) push_back_vec_pts3D(&vertices, &verts[k]);
        for (k = 0; k < 12; k++) push_back_vec_int(&status_face, &sf[k]);
        for (k = 0; k < 12; k++) infogrb = GrB_Matrix_setElement(*volumes, 1, k, 0);

    } else if (case5){
        int n1, n2, e1, e2;
        int e_13, e_14, e_23, e_24;
        Point2D *p2D;
        Point3D pt13D_n, pt13D_np1, pt23D_n, pt23D_np1, pt33D_n, pt33D_np1, pt43D_n, pt43D_np1;
        my_real_c lam_a, lam_b, lam_c;
        Point3D verts[16];
        long int sf[16];
        int k;

        vertices    = alloc_with_capacity_vec_pts3D(16);
        infogrb = GrB_Matrix_new(edges,   GrB_INT8, 16, 24);
        infogrb = GrB_Matrix_new(faces,   GrB_INT8, 24, 16);
        infogrb = GrB_Matrix_new(volumes, GrB_INT8, 16, 4);
        status_face = alloc_with_capacity_vec_int(16);

        // Keep pt1 and pt2 as they are; make sure pt3 and pt4 are correctly ordered.
        vertex_neighbors(i_pt1, &n1, &n2, &e1, &e2);
        if (n1 == i_pt3){
            if (grid_edge_sign(grid, i_pt1, e1) < 0){ int tmp = i_pt3; i_pt3 = i_pt4; i_pt4 = tmp; e_13 = e2; e_14 = e1; }
            else                                     { e_13 = e1; e_14 = e2; }
        } else {
            if (grid_edge_sign(grid, i_pt1, e2) < 0){ int tmp = i_pt3; i_pt3 = i_pt4; i_pt4 = tmp; e_13 = e1; e_14 = e2; }
            else                                     { e_13 = e2; e_14 = e1; }
        }
        vertex_neighbors(i_pt2, &n1, &n2, &e1, &e2);
        if (n1 == i_pt3){ e_23 = e1; e_24 = e2; }
        else            { e_23 = e2; e_24 = e1; }

        p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt1);
        pt13D_n = (Point3D){p2D->x,p2D->y,0.0}; pt13D_np1 = (Point3D){p2D->x,p2D->y,dt};
        p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt2);
        pt23D_n = (Point3D){p2D->x,p2D->y,0.0}; pt23D_np1 = (Point3D){p2D->x,p2D->y,dt};
        p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt3);
        pt33D_n = (Point3D){p2D->x,p2D->y,0.0}; pt33D_np1 = (Point3D){p2D->x,p2D->y,dt};
        p2D = get_ith_elem_vec_pts2D(grid->vertices, i_pt4);
        pt43D_n = (Point3D){p2D->x,p2D->y,0.0}; pt43D_np1 = (Point3D){p2D->x,p2D->y,dt};

        // --- corner tetrahedron at pt1 ---
        lam_a = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tn[i_pt3]);
        lam_b = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tn[i_pt4]);
        lam_c = level_set_tn[i_pt1] / (level_set_tn[i_pt1] - level_set_tnp1[i_pt1]);

        verts[0] = pt13D_n;
        verts[1] = lerp_pt3D(lam_a, pt13D_n, pt33D_n);
        verts[2] = lerp_pt3D(lam_b, pt13D_n, pt43D_n);
        verts[3] = lerp_pt3D(lam_c, pt13D_n, pt13D_np1);

        infogrb = GrB_Matrix_setElement(*edges, -1, 0, 0); infogrb = GrB_Matrix_setElement(*edges,  1, 1, 0);
        infogrb = GrB_Matrix_setElement(*edges, -1, 0, 1); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 1);
        infogrb = GrB_Matrix_setElement(*edges, -1, 0, 2); infogrb = GrB_Matrix_setElement(*edges,  1, 3, 2);
        infogrb = GrB_Matrix_setElement(*edges, -1, 3, 3); infogrb = GrB_Matrix_setElement(*edges,  1, 1, 3);
        infogrb = GrB_Matrix_setElement(*edges, -1, 3, 4); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 4);
        infogrb = GrB_Matrix_setElement(*edges, -1, 1, 5); infogrb = GrB_Matrix_setElement(*edges,  1, 2, 5);

        infogrb = GrB_Matrix_setElement(*faces, -1, 0, 0); infogrb = GrB_Matrix_setElement(*faces,  1, 1, 0); infogrb = GrB_Matrix_setElement(*faces, -1, 5, 0);
        infogrb = GrB_Matrix_setElement(*faces,  1, 0, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 2, 1); infogrb = GrB_Matrix_setElement(*faces, -1, 3, 1);
        infogrb = GrB_Matrix_setElement(*faces, -1, 1, 2); infogrb = GrB_Matrix_setElement(*faces,  1, 2, 2); infogrb = GrB_Matrix_setElement(*faces,  1, 4, 2);
        infogrb = GrB_Matrix_setElement(*faces,  1, 3, 3); infogrb = GrB_Matrix_setElement(*faces, -1, 4, 3); infogrb = GrB_Matrix_setElement(*faces,  1, 5, 3);

        sf[0] = 1;
        sf[1] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_13);
        sf[2] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_14);
        sf[3] = -1;

        // --- corner tetrahedron at pt2 ---
        lam_a = level_set_tn[i_pt2] / (level_set_tn[i_pt2] - level_set_tn[i_pt3]);
        lam_b = level_set_tn[i_pt2] / (level_set_tn[i_pt2] - level_set_tn[i_pt4]);
        lam_c = level_set_tn[i_pt2] / (level_set_tn[i_pt2] - level_set_tnp1[i_pt2]);

        verts[4] = pt23D_n;
        verts[5] = lerp_pt3D(lam_b, pt23D_n, pt43D_n);
        verts[6] = lerp_pt3D(lam_a, pt23D_n, pt33D_n);
        verts[7] = lerp_pt3D(lam_c, pt23D_n, pt23D_np1);

        infogrb = GrB_Matrix_setElement(*edges, -1, 4, 6); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 6);
        infogrb = GrB_Matrix_setElement(*edges, -1, 4, 7); infogrb = GrB_Matrix_setElement(*edges,  1, 6, 7);
        infogrb = GrB_Matrix_setElement(*edges, -1, 4, 8); infogrb = GrB_Matrix_setElement(*edges,  1, 7, 8);
        infogrb = GrB_Matrix_setElement(*edges, -1, 7, 9); infogrb = GrB_Matrix_setElement(*edges,  1, 5, 9);
        infogrb = GrB_Matrix_setElement(*edges, -1, 7, 10); infogrb = GrB_Matrix_setElement(*edges,  1, 6, 10);
        infogrb = GrB_Matrix_setElement(*edges, -1, 5, 11); infogrb = GrB_Matrix_setElement(*edges,  1, 6, 11);

        infogrb = GrB_Matrix_setElement(*faces, -1, 6, 4); infogrb = GrB_Matrix_setElement(*faces,  1, 7, 4); infogrb = GrB_Matrix_setElement(*faces, -1, 11, 4);
        infogrb = GrB_Matrix_setElement(*faces,  1, 6, 5); infogrb = GrB_Matrix_setElement(*faces, -1, 8, 5); infogrb = GrB_Matrix_setElement(*faces, -1, 9, 5);
        infogrb = GrB_Matrix_setElement(*faces, -1, 7, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 8, 6); infogrb = GrB_Matrix_setElement(*faces,  1, 10, 6);
        infogrb = GrB_Matrix_setElement(*faces,  1, 9, 7); infogrb = GrB_Matrix_setElement(*faces, -1, 10, 7); infogrb = GrB_Matrix_setElement(*faces,  1, 11, 7);

        sf[4] = 1;
        sf[5] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_24);
        sf[6] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_23);
        sf[7] = -1;

        // --- volume at pt3 (t^{n+1} side) ---
        lam_a = level_set_tnp1[i_pt3] / (level_set_tnp1[i_pt3] - level_set_tnp1[i_pt2]);
        lam_b = level_set_tnp1[i_pt3] / (level_set_tnp1[i_pt3] - level_set_tnp1[i_pt1]);
        lam_c = level_set_tn[i_pt3]   / (level_set_tn[i_pt3]   - level_set_tnp1[i_pt3]);

        verts[8]  = pt33D_np1;
        verts[9]  = lerp_pt3D(lam_a, pt33D_np1, pt23D_np1);
        verts[10] = lerp_pt3D(lam_b, pt33D_np1, pt13D_np1);
        verts[11] = lerp_pt3D(lam_c, pt33D_n, pt33D_np1);

        infogrb = GrB_Matrix_setElement(*edges, -1, 8, 12); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 12);
        infogrb = GrB_Matrix_setElement(*edges, -1, 8, 13); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 13);
        infogrb = GrB_Matrix_setElement(*edges, -1, 8, 14); infogrb = GrB_Matrix_setElement(*edges,  1, 11, 14);
        infogrb = GrB_Matrix_setElement(*edges, -1, 11, 15); infogrb = GrB_Matrix_setElement(*edges,  1, 9, 15);
        infogrb = GrB_Matrix_setElement(*edges, -1, 11, 16); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 16);
        infogrb = GrB_Matrix_setElement(*edges, -1, 9, 17); infogrb = GrB_Matrix_setElement(*edges,  1, 10, 17);

        infogrb = GrB_Matrix_setElement(*faces,  1, 12, 8); infogrb = GrB_Matrix_setElement(*faces, -1, 13, 8); infogrb = GrB_Matrix_setElement(*faces,  1, 17, 8);
        infogrb = GrB_Matrix_setElement(*faces, -1, 12, 9); infogrb = GrB_Matrix_setElement(*faces,  1, 14, 9); infogrb = GrB_Matrix_setElement(*faces,  1, 15, 9);
        infogrb = GrB_Matrix_setElement(*faces,  1, 13, 10); infogrb = GrB_Matrix_setElement(*faces, -1, 14, 10); infogrb = GrB_Matrix_setElement(*faces, -1, 16, 10);
        infogrb = GrB_Matrix_setElement(*faces, -1, 15, 11); infogrb = GrB_Matrix_setElement(*faces,  1, 16, 11); infogrb = GrB_Matrix_setElement(*faces, -1, 17, 11);

        sf[8] = 2;
        sf[9] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_23);
        sf[10] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_13);
        sf[11] = -1;

        lam_a = level_set_tnp1[i_pt4] / (level_set_tnp1[i_pt4] - level_set_tnp1[i_pt2]);
        lam_b = level_set_tnp1[i_pt4] / (level_set_tnp1[i_pt4] - level_set_tnp1[i_pt1]);
        lam_c = level_set_tn[i_pt4]   / (level_set_tn[i_pt4]   - level_set_tnp1[i_pt4]);

        verts[12] = pt43D_np1;
        verts[13] = lerp_pt3D(lam_b, pt43D_np1, pt13D_np1);
        verts[14] = lerp_pt3D(lam_a, pt43D_np1, pt23D_np1);
        verts[15] = lerp_pt3D(lam_c, pt43D_n, pt43D_np1);

        infogrb = GrB_Matrix_setElement(*edges, -1, 12, 18); infogrb = GrB_Matrix_setElement(*edges,  1, 13, 18);
        infogrb = GrB_Matrix_setElement(*edges, -1, 12, 19); infogrb = GrB_Matrix_setElement(*edges,  1, 14, 19);
        infogrb = GrB_Matrix_setElement(*edges, -1, 12, 20); infogrb = GrB_Matrix_setElement(*edges,  1, 15, 20);
        infogrb = GrB_Matrix_setElement(*edges, -1, 15, 21); infogrb = GrB_Matrix_setElement(*edges,  1, 13, 21);
        infogrb = GrB_Matrix_setElement(*edges, -1, 15, 22); infogrb = GrB_Matrix_setElement(*edges,  1, 14, 22);
        infogrb = GrB_Matrix_setElement(*edges, -1, 13, 23); infogrb = GrB_Matrix_setElement(*edges,  1, 14, 23);

        infogrb = GrB_Matrix_setElement(*faces,  1, 18, 12); infogrb = GrB_Matrix_setElement(*faces, -1, 19, 12); infogrb = GrB_Matrix_setElement(*faces,  1, 23, 12);
        infogrb = GrB_Matrix_setElement(*faces, -1, 18, 13); infogrb = GrB_Matrix_setElement(*faces,  1, 20, 13); infogrb = GrB_Matrix_setElement(*faces,  1, 21, 13);
        infogrb = GrB_Matrix_setElement(*faces,  1, 19, 14); infogrb = GrB_Matrix_setElement(*faces, -1, 20, 14); infogrb = GrB_Matrix_setElement(*faces, -1, 22, 14);
        infogrb = GrB_Matrix_setElement(*faces, -1, 21, 15); infogrb = GrB_Matrix_setElement(*faces,  1, 22, 15); infogrb = GrB_Matrix_setElement(*faces, -1, 23, 15);

        sf[12] = 2;
        sf[13] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_14);
        sf[14] = 2 + *get_ith_elem_vec_int(grid->status_edge, e_24);
        sf[15] = -1;

        for (k = 0; k < 16; k++) push_back_vec_pts3D(&vertices, &verts[k]);
        for (k = 0; k < 16; k++) push_back_vec_int(&status_face, &sf[k]);
        for (k = 0;  k < 4;  k++) infogrb = GrB_Matrix_setElement(*volumes, 1, k, 0);
        for (k = 3;  k < 8;  k++) infogrb = GrB_Matrix_setElement(*volumes, 1, k, 1);
        for (k = 7;  k < 12; k++) infogrb = GrB_Matrix_setElement(*volumes, 1, k, 2);
        for (k = 11; k < 16; k++) infogrb = GrB_Matrix_setElement(*volumes, 1, k, 3);

    } else {
        printf("Error: that should not happen - have we forgotten some cases?\n");
        return;
    }

    built = new_Polyhedron3D_vefvs(vertices, edges, faces, volumes, status_face);
    copy_Polyhedron3D(built, clipped3D);

    dealloc_Polyhedron3D(built); free(built);
    dealloc_vec_pts3D(vertices); free(vertices);
    dealloc_vec_int(status_face); free(status_face);
    GrB_free(edges);   free(edges);
    GrB_free(faces);   free(faces);
    GrB_free(volumes); free(volumes);
}

/// @brief Compute the effective areas in `grid` when occupied by `clipped3D` represeting a moving polygon.
/// @details clipped3D is supposed to be clipped by grid at input.
/// @param grid [IN] non-moving polygon, clipper.
/// @param clipped3D [IN] moving polygon.
/// @param dt [IN] time-step
/// @param nb_pts [IN] Number of points in 2D polygon used to build clipped3D
/// @param id_cell [IN] Cell identificator of grid
/// @param lambdas_arr [OUT] Array of effective area for each edge of `grid`. Allocated inside the function.
/// @param big_lambda_n [OUT] First cell: effective area at time t^n. Second cell: area of grid - effective area. Allocated inside the function.
/// @param big_lambda_np1 [OUT] First cell: effective area at time t^n+`dt`. Second cell: area of grid - effective area. Allocated inside the function.
/// @param normals_ptr [OUT] list of normal vectors of the faces of `clipped3D` inside `grid`.
/// @param edge_indices [OUT] list of edge indices clipped inside grid.
/// @param is_narrowband [OUT] True if the intersection of `grid` and `clipped` is not empty at time t^n or t^n+`dt`, false otherwise.
void compute_lambdas2D_noclip(const Polygon2D* grid, const Polyhedron3D *clipped3D, const my_real_c dt, \
                        Array_double **lambdas_arr, Vector_double** big_lambda_n, Vector_double** big_lambda_np1, \
                        Vector_points3D **normals_ptr, Vector_int64 **edge_indices, bool *is_narrowband){
    const unsigned int nb_regions = 2;
    my_real_c *val = (my_real_c*)malloc(sizeof(my_real_c));
    my_real_c nm;
    Point3D *area, *occupied;
    Vector_points3D *occupied_area = NULL;
    //Polygon2D *mini_clipped;
    Polyhedron3D *cell3D = NULL;
    //long *sfj;
    Point3D *pt3D, *local_l, zero_pt;
    bool local_narrowband;
    //Vector_points2D *vec_move_grid;
    my_real_c *vec_move_gridx = NULL, *vec_move_gridy = NULL;
    Vector_points3D *surfaces = NULL;
    GrB_Index nb_edge, i, j, k;
    //GrB_Index nb_clipped_faces;
    Vector_points3D *lambdas3D = NULL; //TODO : Change this when nb_regions>2
    Array_points3D *local_lambdas = NULL;
    
    GrB_Matrix_ncols(&nb_edge, *(grid->edges));
    *lambdas_arr = alloc_with_capacity_arr_double(nb_edge, nb_regions); //All set to 0
    local_lambdas = alloc_with_capacity_arr_pts3D(nb_edge + 2, nb_regions-1); //All set to 0
    occupied_area = alloc_with_capacity_vec_pts3D(nb_edge + 2);

    *big_lambda_n = alloc_with_capacity_vec_double(nb_regions);
    *big_lambda_np1 = alloc_with_capacity_vec_double(nb_regions);
    *val = 0.;
    for (i=0; i<nb_regions; i++){
        set_ith_elem_vec_double(*big_lambda_n, i, val);
        set_ith_elem_vec_double(*big_lambda_np1, i, val);
    }
    
    zero_pt = (Point3D){0., 0., 0.};
    for (i=0; i<nb_edge + 2; i++){
        set_ith_elem_vec_pts3D(occupied_area, i, &zero_pt);
    }

    *is_narrowband = false;

    if ((clipped3D) && (clipped3D->vertices->size>2)){
        vec_move_gridx = calloc(grid->vertices->size, sizeof(my_real_c));
        vec_move_gridy = calloc(grid->vertices->size, sizeof(my_real_c));
        cell3D = build_space2D_time_cell(grid, vec_move_gridx, vec_move_gridy, grid->vertices->size, dt, false, NULL);
        surfaces = points3D_from_matrix(surfaces_poly3D(cell3D)); //Could be changed: cell3D is not really needed, this is just the length of each 
        k = 0; //should be k = clipped->status_edge[i], or another variable to indicate what region covers face nb i.
        
        compute_lambdas2D_time_clipped(nb_edge + 2, clipped3D, &lambdas3D, normals_ptr, edge_indices, &local_narrowband);

        *is_narrowband |= local_narrowband;
        
        //λe += λ[2:end]
        for (j=0; j<lambdas3D->size; j++){
            pt3D = get_ith_elem_vec_pts3D(lambdas3D, j); 
            local_l = get_ijth_elem_arr_pts3D(local_lambdas, j, k);
            local_l->x += pt3D->x;
            local_l->y += pt3D->y;
            local_l->t += pt3D->t;

            local_l = get_ith_elem_vec_pts3D(occupied_area, j);
            local_l->x += pt3D->x;
            local_l->y += pt3D->y;
            local_l->t += pt3D->t;
        }
        

        i = 0;
        area = get_ith_elem_vec_pts3D(surfaces, i);
        occupied = get_ith_elem_vec_pts3D(occupied_area, i);
        *val = fmax(0., norm_pt3D(*area) - norm_pt3D(*occupied));

        set_ith_elem_vec_double(*big_lambda_n, 0, val);
        for(k = 1; k<nb_regions; k++){
            nm = norm_pt3D(*get_ijth_elem_arr_pts3D(local_lambdas, i, k-1));
            set_ith_elem_vec_double(*big_lambda_n, k, &nm);
        }
        

        i = 1;
        area = get_ith_elem_vec_pts3D(surfaces, i);
        occupied = get_ith_elem_vec_pts3D(occupied_area, i);
        *val = fmax(0., norm_pt3D(*area) - norm_pt3D(*occupied));
        set_ith_elem_vec_double(*big_lambda_np1, 0, val);
        for(k = 1; k<nb_regions; k++){
            nm = norm_pt3D(*get_ijth_elem_arr_pts3D(local_lambdas, i, k-1));
            set_ith_elem_vec_double(*big_lambda_np1, k, &nm);
        }

        for (i=2; i<nb_edge + 2; i++){
            area = get_ith_elem_vec_pts3D(surfaces, i);
            occupied = get_ith_elem_vec_pts3D(occupied_area, i);
            *val = fmax(0., norm_pt3D(*area) - norm_pt3D(*occupied)) / dt;
            set_ijth_elem_arr_double(*lambdas_arr, i-2, 0, val);
            for(k = 1; k<nb_regions; k++){
                nm = norm_pt3D(*get_ijth_elem_arr_pts3D(local_lambdas, i, k-1)) / dt;
                //if (nm > 0.13){
                //    printf("Warning: large effective area on face %lu of grid cell, value = %lf\n", i-2, nm);
                //    print_vec_pt3D(*lambdas3D);
                //}
                set_ijth_elem_arr_double(*lambdas_arr, i-2, k, &nm); 
            }
        }
    } else {
        i = 0;
        *val = norm_pt3D(*get_ith_elem_vec_pts3D(surfaces, i));
        set_ith_elem_vec_double(*big_lambda_n, 0, val);
        nm = 0.;
        for(k = 1; k<nb_regions; k++){
            set_ith_elem_vec_double(*big_lambda_n, k, &nm);
        }

        i = 1;
        *val = norm_pt3D(*get_ith_elem_vec_pts3D(surfaces, i));
        set_ith_elem_vec_double(*big_lambda_np1, 0, val);
        for(k = 1; k<nb_regions; k++){
            set_ith_elem_vec_double(*big_lambda_np1, k, &nm);
        }

        for (i=2; i<nb_edge + 2; i++){
            *val = norm_pt3D(*get_ith_elem_vec_pts3D(surfaces, i)) / dt;
            set_ijth_elem_arr_double(*lambdas_arr, i-2, 0, val);

            for(k = 1; k<nb_regions; k++){
                set_ijth_elem_arr_double(*lambdas_arr, i-2, k, &nm); 
            }
        }
        if (!(*normals_ptr)){
            *normals_ptr = alloc_empty_vec_pts3D();
        }
    }

    dealloc_Polyhedron3D(cell3D); free(cell3D);
    dealloc_vec_pts3D(surfaces); free(surfaces);
    dealloc_vec_pts3D(occupied_area); free(occupied_area);
    dealloc_vec_pts3D(lambdas3D); free(lambdas3D);
    dealloc_arr_pts3D(local_lambdas); free(local_lambdas);
    if (vec_move_gridx) free(vec_move_gridx);
    if (vec_move_gridy) free(vec_move_gridy);
    if (val) free(val);
}
