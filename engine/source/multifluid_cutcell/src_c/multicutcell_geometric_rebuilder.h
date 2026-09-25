#ifndef MULTICUTCELL_GEOMETRIC_REBUILDER_H
#define MULTICUTCELL_GEOMETRIC_REBUILDER_H

#include "my_real_c.inc"
#include "Polygon2D.h"
#include "Polyhedron3D.h"
#include "array_double.h"
#include "vector_double.h"
#include "vector_Points.h"
#include "vector_int.h"

void rebuild_polyhedron_level_set(const Polygon2D* grid, Polyhedron3D* clipped3D, my_real_c dt, \
                                const my_real_c* level_set_tn, const my_real_c* level_set_tnp1, \
                                long long int nb_tn, long long int nb_tnp1, long long int *is_reversed);

void compute_lambdas2D_noclip(const Polygon2D* grid, const Polyhedron3D *clipped3D, const my_real_c dt, \
                        Array_double **lambdas_arr, Vector_double** big_lambda_n, Vector_double** big_lambda_np1, \
                        Vector_points3D **normals_ptr, Vector_int64 **edge_indices, bool *is_narrowband);

#endif