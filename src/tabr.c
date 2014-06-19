
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "tabr.h"

#define NO_INDEX (-1)

#define REALLOC_FACTOR 2
#define N_MAX_WAITING_DEFAULT 10
#define FAVOR_EDGES_WITH_ROOTS 0
#define CONVEXITY_THRESHOLD 0.01

/* Nomenclature: NUMBER refers to indices running from 0 to 2 (edges
 * in triangle) or 0 to 1 (vertices in edges).  INDEX refers to
 * indices into the vertex, edge, or triangle tables, and run from 0
 * to n_vertices, n_edges, or n_triangles.
 */

static inline double
hypot2_sqr(double x, double y)
{
  return x*x + y*y;
}

static inline double
hypot3_sqr(double x, double y, double z)
{
  return x*x + y*y + z*z;
}

static double
dist2_sqr(const double u[2], const double v[2])
{
  double dx = u[0] - v[0];
  double dy = u[1] - v[1];

  return hypot2_sqr(dx, dy);
}

static double
transdot2(const double origin[2], 
	  const double a[2], const double b[2], 
	  int normalize)
{
  double u[2], v[2];

  u[0] = a[0] - origin[0];
  u[1] = a[1] - origin[1];
  v[0] = b[0] - origin[0];
  v[1] = b[1] - origin[1];

  double norm = 1.0;

  if (normalize)
    {
      double u2 = hypot2_sqr(u[0], u[1]);
      double v2 = hypot2_sqr(v[0], v[1]);
      norm = sqrt(u2 * v2);
    }

  return (u[0]*v[0] + u[1]*v[1]) / norm;
}

static double
transcross2(const double origin[2], 
	    const double a[2], const double b[2], 
	    int normalize)
{
  double u[2], v[2];

  u[0] = a[0] - origin[0];
  u[1] = a[1] - origin[1];
  v[0] = b[0] - origin[0];
  v[1] = b[1] - origin[1];

  double norm = 1.0;

  if (normalize)
    {
      double u2 = hypot2_sqr(u[0], u[1]);
      double v2 = hypot2_sqr(v[0], v[1]);
      norm = sqrt(u2 * v2);
    }

  return (u[0]*v[1] - u[1]*v[0]) / norm;
}

static void
cross3(const double v0[3],
       const double v1[3],
       double cross[3])
{
  cross[0] = v0[1] * v1[2] - v0[2] * v1[1];
  cross[1] = v0[2] * v1[0] - v0[0] * v1[2];
  cross[2] = v0[0] * v1[1] - v0[1] * v1[0];
}

static double
dot3(const double v0[3], const double v1[3])
{
  return v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2];
}

static void
normalize3(double v[3])
{
  double s = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

  if (s == 0.0) return;

  v[0] /= s;
  v[1] /= s;
  v[2] /= s;
}

static void
transcross3(const double origin[3],
	    const double v0[3],
	    const double v1[3],
	    double result[3])
{
  double u0[3], u1[3];
  int i;

  for (i = 0; i < 3; ++i)
    {
      u0[i] = v0[i] - origin[i];
      u1[i] = v1[i] - origin[i];
    }

  cross3(u0, u1, result);
}

static double
triangle3_area(const double v0[3],
	       const double v1[3],
	       const double v2[3])
{
  double u[3], w[3], uxw[3];
  int i;
  
  for (i = 0; i < 3; ++i)
    {
      u[i] = v1[i] - v0[i];
      w[i] = v2[i] - v0[i];
    }

  cross3(u, w, uxw);

  return 0.5 * sqrt(hypot3_sqr(uxw[0], uxw[1], uxw[2]));
}

/* Calculate the cosine of the angles between the normals of two
 * triangles (v0, v2, v1) and (v0, v1, v3). */
double
cos_angle_triangle_normals(const double v0[3], const double v1[3],
			   const double v2[3], const double v3[3])
{
  double normal_left[3];
  double normal_right[3];

  transcross3(v0, v1, v2, normal_left);
  transcross3(v0, v3, v1, normal_right);
  normalize3(normal_left);
  normalize3(normal_right);

  return dot3(normal_left, normal_right);
}


static inline int
next_number(int n)
{
  return n == 2 ? 0 : n + 1;
}

static inline int
prev_number(int n)
{
  return n == 0 ? 2 : n - 1;
}

static inline int
has_daughter_edge(tabr_t *self, int edge_index)
{
  if (self->edge[edge_index].daughter_edge_index[0] != -1)
    {
      assert(self->edge[edge_index].daughter_edge_index[1] != -1);

      return 1;
    }
  return 0;
}

static inline void
get_edge_in_triangle(tabr_t *self, 
		     int triangle_index, int edge_number,
		     int *edge_index, int *edge_orientation)
{
  *edge_index = self->triangle[triangle_index].edge_index[edge_number];
  *edge_orientation = self->triangle[triangle_index].edge_orientation[edge_number];
}

static inline int
get_counter_clockwise_vertex_in_edge(tabr_t *self, int edge_index, int edge_orientation)
{
  return self->edge[edge_index].vertex_index[edge_orientation == 1 ? 1 : 0];
}

static inline int
get_clockwise_vertex_in_edge(tabr_t *self, int edge_index, int edge_orientation)
{
  return self->edge[edge_index].vertex_index[edge_orientation == 1 ? 0 : 1];
}

/* Given triangle, fetch its vertices in CCW order */
static int
get_triangle_vertex_indices(tabr_t *self, int triangle_index, int vertex[3])
{
  assert(triangle_index >= 0);
  assert(triangle_index < self->n_triangles);

  int e0 = self->triangle[triangle_index].edge_index[0];
  int e1 = self->triangle[triangle_index].edge_index[1];
  int e2 = self->triangle[triangle_index].edge_index[2];

  int o0 = self->triangle[triangle_index].edge_orientation[0];
  int o1 = self->triangle[triangle_index].edge_orientation[1];
  int o2 = self->triangle[triangle_index].edge_orientation[2];

  vertex[0] = self->edge[e0].vertex_index[o0 == 1 ? 0 : 1];
  vertex[1] = self->edge[e1].vertex_index[o1 == 1 ? 0 : 1];
  vertex[2] = self->edge[e2].vertex_index[o2 == 1 ? 0 : 1];

  assert(vertex[0] == self->edge[e2].vertex_index[o2 == 1 ? 1 : 0]);
  assert(vertex[1] == self->edge[e0].vertex_index[o0 == 1 ? 1 : 0]);
  assert(vertex[2] == self->edge[e1].vertex_index[o1 == 1 ? 1 : 0]);

  return 0;
}

/* Given triangle and edge indices, find the number of the edge within
 * the triangle where the edge appears (if it appears). */
static int
find_edge_in_triangle(tabr_t *self, int triangle_index, int edge_index)
{
  int edge_number;
  for (edge_number = 0; edge_number < 3; ++edge_number)
    {
      if (self->triangle[triangle_index].edge_index[edge_number] == edge_index)
	{
	  return edge_number;
	}
    }

  return -1;
}

/* Given edge, find its neighbourhood of triangles, edges, and
 * vertices. If v0 and v1 are vertices of this edge, then v2 is the
 * vertex on the left of edge (v0, v1) (may not exist), and v3 is the
 * vertex on the right of (v0, v1) (if exists). Edge 0 is (v0, v2),
 * edge 1 is (v2, v1), edge 2 is (v0, v3), and edge 3 is (v3, v1).
 */
static void
get_edge_neighbourhood(tabr_t *self, int edge_index, 
		       int triangle[2],
		       int edge[4],
		       int edge_orientation[4],
		       int vertex[4])
{
  int edge_number;

  vertex[0] = self->edge[edge_index].vertex_index[0];
  vertex[1] = self->edge[edge_index].vertex_index[1];

  triangle[0] = self->edge[edge_index].triangle_index[0];
  triangle[1] = self->edge[edge_index].triangle_index[1];

  int left_right;
  for (left_right = 0; left_right < 2; ++left_right)
    {
      int fst_edgenum = 2*left_right + 0;
      int snd_edgenum = 2*left_right + 1;
      int mid_vertnum = 2 + left_right;

      if (triangle[left_right] != NO_INDEX)
	{
	  edge_number = find_edge_in_triangle(self, triangle[left_right], 
					      edge_index);
	  
	  assert(edge_number != NO_INDEX);
	  
	  int this_edge_index, this_edge_orientation;
	  
	  get_edge_in_triangle(self, triangle[left_right], edge_number,
			       &this_edge_index, &this_edge_orientation);

	  assert(this_edge_index == edge_index);

	  if (left_right == 0)
	    {
	      get_edge_in_triangle(self, 
				   triangle[left_right], prev_number(edge_number), 
				   &edge[fst_edgenum], 
				   &edge_orientation[fst_edgenum]);
	  
	      get_edge_in_triangle(self, 
				   triangle[left_right], next_number(edge_number), 
				   &edge[snd_edgenum],
				   &edge_orientation[snd_edgenum]);
	      
	      vertex[mid_vertnum] 
		= get_clockwise_vertex_in_edge(self, 
					       edge[fst_edgenum], 
					       edge_orientation[fst_edgenum]);
	  
	      assert(vertex[mid_vertnum] 
		     == get_counter_clockwise_vertex_in_edge(self, 
							     edge[snd_edgenum], 
							     edge_orientation[snd_edgenum]));
	    }
	  else
	    {
	      get_edge_in_triangle(self, 
				   triangle[left_right], next_number(edge_number), 
				   &edge[fst_edgenum], 
				   &edge_orientation[fst_edgenum]);
	  
	      get_edge_in_triangle(self, 
				   triangle[left_right], prev_number(edge_number), 
				   &edge[snd_edgenum], 
				   &edge_orientation[snd_edgenum]);


	      vertex[mid_vertnum] 
		= get_counter_clockwise_vertex_in_edge(self, 
						       edge[fst_edgenum], 
						       edge_orientation[fst_edgenum]);
	      
	      assert(vertex[mid_vertnum] 
		     == get_clockwise_vertex_in_edge(self, 
						     edge[snd_edgenum], 
						     edge_orientation[snd_edgenum]));
	    }

	}
      else
	{
	  edge[snd_edgenum] = edge[fst_edgenum] = NO_INDEX;
	  vertex[mid_vertnum] = NO_INDEX;
	}
    }
}

static int
more_vertices_ok(tabr_t *self)
{
  int new_n_max_vertices = REALLOC_FACTOR * self->n_max_vertices;
  tabr_vertex_t *new_vertex;

  new_vertex = realloc(self->vertex, new_n_max_vertices * sizeof(*new_vertex));
  if (new_vertex == NULL) return 0;

  int i;
  for (i = self->n_max_vertices; i < new_n_max_vertices; ++i)
    {
      new_vertex[i].position = NULL;
    }

  for (i = self->n_max_vertices; i < new_n_max_vertices; ++i)
    {
      new_vertex[i].position = malloc((2 + self->n_z) * sizeof(double));
      if (new_vertex[i].position == NULL) return 0;
    }

  self->n_max_vertices = new_n_max_vertices;
  self->vertex         = new_vertex;

  return 1;
}

static int
more_edges_ok(tabr_t *self)
{
  int new_n_max_edges = REALLOC_FACTOR * self->n_max_edges;
  tabr_edge_t *new_edge;

  new_edge = realloc(self->edge, new_n_max_edges * sizeof(*new_edge));
  if (new_edge == NULL) return 0;

  self->n_max_edges = new_n_max_edges;
  self->edge        = new_edge;

  return 1;
}

static int
more_triangles_ok(tabr_t *self)
{
  int new_n_max_triangles = REALLOC_FACTOR * self->n_max_triangles;
  tabr_triangle_t *new_triangle;

  new_triangle = realloc(self->triangle, new_n_max_triangles * sizeof(*new_triangle));
  if (new_triangle == NULL) return 0;

  self->n_max_triangles = new_n_max_triangles;
  self->triangle        = new_triangle;

  return 1;
}

static int
mod_vertex_ok(tabr_t *self, int vertex_index, double x, double y)
{
  self->vertex[vertex_index].position[0] = x;
  self->vertex[vertex_index].position[1] = y;
  self->vertex[vertex_index].has_data    = 0;

  return 1;
}

static int
add_vertex_ok(tabr_t *self, double x, double y)
{
  if (self->n_vertices == self->n_max_vertices)
    {
      if (!more_vertices_ok(self)) return 0;
    }

  return mod_vertex_ok(self, self->n_vertices++, x, y);
}

static int
mod_edge_ok(tabr_t *self, int edge_index, 
	    int vi1, int vi2, int ti1, int ti2,
	    int dei1, int dei2, int pei)
{
  self->edge[edge_index].vertex_index[0] = vi1;
  self->edge[edge_index].vertex_index[1] = vi2;

  self->edge[edge_index].triangle_index[0] = ti1;
  self->edge[edge_index].triangle_index[1] = ti2;

  self->edge[edge_index].daughter_edge_index[0] = dei1;
  self->edge[edge_index].daughter_edge_index[1] = dei2;

  self->edge[edge_index].parent_edge_index = pei;

  return 1;
}

static int
add_edge_ok(tabr_t *self, int vi1, int vi2, int ti1, int ti2, int pei)
{
  if (self->n_edges == self->n_max_edges)
    {
      if (!more_edges_ok(self)) return 0;
    }

  return mod_edge_ok(self, self->n_edges++, 
		     vi1, vi2, ti1, ti2, -1, -1, pei);
}

static int
mod_triangle_ok(tabr_t *self, int triangle_index,
		int ei1, int ei2, int ei3,
		int oi1, int oi2, int oi3)
{
  self->triangle[triangle_index].edge_index[0] = ei1;
  self->triangle[triangle_index].edge_index[1] = ei2;
  self->triangle[triangle_index].edge_index[2] = ei3;

  self->triangle[triangle_index].edge_orientation[0] = oi1;
  self->triangle[triangle_index].edge_orientation[1] = oi2;
  self->triangle[triangle_index].edge_orientation[2] = oi3;

  return 1;
}

static int
add_triangle_ok(tabr_t *self, 
		int ei1, int ei2, int ei3,
		int oi1, int oi2, int oi3)
{
  if (self->n_triangles == self->n_max_triangles)
    {
      if (!more_triangles_ok(self)) return 0;
    }

  return mod_triangle_ok(self, self->n_triangles++,
			 ei1, ei2, ei3,
			 oi1, oi2, oi3);
}

tabr_t *
tabr_alloc(double x_min, double x_max, int x_n,
	   double y_min, double y_max, int y_n,
	   int n_rec, int n_z)
{
  int i, j;

  if (x_n <= 1 || y_n <= 1 || n_rec <= 0) return NULL;
  if (x_min >= x_max || y_min >= y_max) return NULL;
  if (n_z < 1) return NULL;

  int n_triangles_per_row = 2 * (x_n - 1);
  int n_edges_per_row     = 3 * (x_n - 1) + 1;
  int n_rows              = y_n - 1;
  int n_columns           = x_n - 1;

  int n_vertices  = x_n * y_n;
  int n_edges     = n_rows * n_edges_per_row + n_columns;
  int n_triangles = n_rows * n_triangles_per_row;

  int overhead_multiplier = 3;
  int n_max_vertices  = overhead_multiplier * n_vertices;
  int n_max_edges     = overhead_multiplier * n_edges;
  int n_max_triangles = overhead_multiplier * n_triangles;

  tabr_t *self;
  self = malloc(sizeof(*self));
  if (self == NULL) return NULL;

  self->n_z = n_z;

  self->vertex   = NULL;
  self->edge     = NULL;
  self->triangle = NULL;

  self->vertex = malloc(n_max_vertices * sizeof(*self->vertex));
  if (self->vertex == NULL) goto fail;

  for (i = 0; i < n_max_vertices; ++i) self->vertex[i].position = NULL;
  for (i = 0; i < n_max_vertices; ++i)
    {
      self->vertex[i].position = malloc((2 + n_z) * sizeof(double));
      if (self->vertex[i].position == NULL) goto fail;
    }

  self->edge = malloc(n_max_edges * sizeof(*self->edge));
  if (self->edge == NULL) goto fail;

  self->triangle = malloc(n_max_triangles * sizeof(*self->triangle));
  if (self->triangle == NULL) goto fail;


  self->n_vertices  = 0;
  self->n_edges     = 0;
  self->n_triangles = 0;

  self->n_vertices_pending = 0;
  self->n_waiting          = 0;
  self->n_max_waiting      = N_MAX_WAITING_DEFAULT;
  self->n_computed         = 0;

  self->n_max_vertices  = n_max_vertices;
  self->n_max_edges     = n_max_edges;
  self->n_max_triangles = n_max_triangles;

  double delta_x = (x_max - x_min) / (x_n - 1);
  double delta_y = (y_max - y_min) / (y_n - 1);

  self->delta_x_max = delta_x;
  self->delta_y_max = delta_y;

  self->split_root_edges = 0;
  self->split_big_small_triangle_pairs = 0;
  self->split_big_small_ratio_limit = 1.0 / 4.0;
  self->split_angled_triangles = 1;
  self->split_angled_cos_angle_limit = cos(2 * M_PI * 0.1);
  self->split_propagate_two_cuts = 1;

  self->n_mesh_smoothing_passes = 0;

  self->n_max_recursion = n_rec;
  self->n_recursion     = 0;

  self->x_min = x_min;
  self->x_max = x_max;
  self->y_min = y_min;
  self->y_max = y_max;

  self->z_min = 0.0;
  self->z_max = 0.0;

  self->z_sum     = 0.0;
  self->z_sqr_sum = 0.0;

  double x, y;
  for (j = 0; j < y_n; ++j)
    {
      y = y_min + j * delta_y;

      for (i = 0; i < x_n; ++i)
	{
	  x = x_min + i * delta_x;
	  if (!add_vertex_ok(self, x, y)) goto fail;
	}
    }

  assert(self->n_vertices == n_vertices);

  for (j = 0; j < n_rows; ++j)
    {
      for (i = 0; i < n_columns; ++i)
	{
	  int this_vertex     = i + j * x_n;
	  int right_vertex    = this_vertex + 1;
	  int up_vertex       = this_vertex + x_n;
	  int up_right_vertex = this_vertex + x_n + 1;

	  int left_edge     = self->n_edges;
	  int bottom_edge   = self->n_edges + 1;
	  int diagonal_edge = self->n_edges + 2;
	  int right_edge    = self->n_edges + 3;
	  int top_edge      = j == n_rows - 1 ? n_edges - n_columns + i : self->n_edges + n_edges_per_row + 1;

	  int lo_triangle     = self->n_triangles;
	  int up_triangle     = self->n_triangles + 1;
	  int bottom_triangle = j == 0             ? -1 : self->n_triangles - n_triangles_per_row + 1;
	  int top_triangle    = j == n_rows - 1    ? -1 : self->n_triangles + n_triangles_per_row;
	  int left_triangle   = i == 0             ? -1 : self->n_triangles - 1;
	  int right_triangle  = i == n_columns - 1 ? -1 : self->n_triangles + 2;

	  if (!add_edge_ok(self,    up_vertex,  this_vertex, lo_triangle,   left_triangle, -1)) goto fail;
	  if (!add_edge_ok(self,  this_vertex, right_vertex, lo_triangle, bottom_triangle, -1)) goto fail;
	  if (!add_edge_ok(self, right_vertex,    up_vertex, lo_triangle,     up_triangle, -1)) goto fail;
	  if (i == n_columns - 1)
	    {
	      if (!add_edge_ok(self, up_right_vertex, right_vertex, right_triangle, up_triangle, -1)) goto fail;
	    }
	  
	  if (!add_triangle_ok(self,  left_edge, bottom_edge, diagonal_edge, +1, +1, +1)) goto fail;
	  if (!add_triangle_ok(self, right_edge,    top_edge, diagonal_edge, -1, -1, -1)) goto fail;
	}
    }

  for (i = 0; i < n_columns; ++i)
    {
      int this_vertex  = i + (y_n - 1) * x_n;
      int right_vertex = this_vertex + 1;
      int bottom_triangle = self->n_triangles - n_triangles_per_row + 2 * i + 1;
      
      if (!add_edge_ok(self, this_vertex, right_vertex, -1, bottom_triangle, -1)) goto fail;
    }

  assert(self->n_edges == n_edges);
  assert(self->n_triangles == n_triangles);

  self->state = TABR_RUNNING;

  return self;

 fail:

  tabr_free(self);
  return NULL;
}

void
tabr_set_propagate_splits(tabr_t *self, int propagate_splits)
{
  self->split_propagate_two_cuts = propagate_splits;
}

void
tabr_set_refine_roots(tabr_t *self, int on_off)
{
  self->split_root_edges = on_off;
}

void
tabr_set_max_angle(tabr_t *self, double angle)
{
  if (angle >= M_PI)
    {
      self->split_angled_triangles = 0;
    }
  else
    {
      self->split_angled_triangles = 1;
      self->split_angled_cos_angle_limit = cos(angle);
    }
}

void
tabr_set_mesh_smoothing_passes(tabr_t *self, int n_mesh_smoothing_passes)
{
  self->n_mesh_smoothing_passes = n_mesh_smoothing_passes;
}


void
tabr_free(tabr_t *self)
{
  if (self->vertex)
    {
      int i;
      for (i = 0; i < self->n_max_vertices; ++i)
	{
	  if (self->vertex[i].position) free(self->vertex[i].position);
	}
      free(self->vertex);
    }
  if (self->edge) free(self->edge);
  if (self->triangle) free(self->triangle);
  free(self);
}

int
tabr_success(tabr_t *self)
{
  return self->state == TABR_SUCCESS;
}



/* Calculate the roundness of a triangle, here the ratio of the
 * shortest and longest edges (squared, to avoid unnecessary
 * sqrts). */
static double
roundness(const double v1[2], const double v2[2], const double v3[2], int *short_edge)
{
  double d12, d23, d31, dmin, dmax;

  d12 = hypot2_sqr(v1[0] - v2[0], v1[1] - v2[1]);
  d23 = hypot2_sqr(v2[0] - v3[0], v2[1] - v3[1]);
  d31 = hypot2_sqr(v3[0] - v1[0], v3[1] - v1[1]);

  *short_edge = 0;
  dmax = dmin = d12;
  if (d23 < dmin) { dmin = d23; *short_edge = 1; }
  else if (d23 > dmax) dmax = d23;

  if (d31 < dmin) { dmin = d31; *short_edge = 2; }
  else if (d31 > dmax) dmax = d31;

  if (dmax == 0.0) return 1.0;

  return dmin / dmax;
}

int
linsol(int n, double *a, double *b)
{
  int i, j, k;
  int swap[n];

  for (i = 0; i < n; ++i) swap[i] = i;

  for (i = 0; i < n; ++i)
    {
      int j_max;
      double a_max;

      j_max = i;
      a_max = fabs(a[n*i + i]);
      
      for (j = i + 1; j < n; ++j)
	{
	  const double a_this = fabs(a[n*i + j]);
	  if (a_this > a_max) { a_max = a_this; j_max = j; }
	}
      
      swap[i] = j_max;
      swap[j_max] = i;
      
      for (k = 0; k < n; ++k)
	{
	  const double a_temp = a[n*k + i];
	  a[n*k + i] = a[n*k + j_max];
	  a[n*k + j_max] = a_temp;
	}

      const double a_ii = a[n*i + i];
      if (a_ii == 0.0) return 1;

      b[i] /= a_ii;

      for (j = i + 1; j < n; ++j) a[n*i + j] /= a_ii;

      for (j = i + 1; j < n; ++j)
	{
	  const double g = a[n*j + i];

	  b[j] -= g * b[i];

	  for (k = i + 1; k < n; ++k)
	    {
	      a[n*j + k] -= g * a[n*i + k];
	    }
	}
    }

  for (i = n - 1; i >= 0; --i)
    {
      for (j = 0; j < i; ++j)
	{
	  b[j] -= b[i] * a[n*j + i];
	}
    }


  for (i = n - 1; i >= 0; --i)
    {
      if (swap[i] != i)
	{
	  int swap_i = swap[i];
	  double b_temp = b[i];
	  b[i] = b[swap_i];
	  b[swap_i] = b_temp;
	}
    }

  return 0;
}

/* Given three*/
static double
bilinear_root(const double r0[3], const double r1[3],
	      const double r2[3], const double r3[3])
{
  assert(r0[2] * r3[2] <= 0.0);

  double x_scale = 1.0;
  double y_scale = 1.0;
  double z_scale = 1.0;

  double x1 = (r1[0] - r0[0]) / x_scale;
  double x2 = (r2[0] - r0[0]) / x_scale;
  double x3 = (r3[0] - r0[0]) / x_scale;

  double y1 = (r1[1] - r0[1]) / y_scale;
  double y2 = (r2[1] - r0[1]) / y_scale;
  double y3 = (r3[1] - r0[1]) / y_scale;

  double z0 = r0[2] / z_scale;
  double z1 = r1[2] / z_scale;
  double z2 = r2[2] / z_scale;
  double z3 = r3[2] / z_scale;

  double a[9], b[3];

  a[0] = x1;  a[1] = y1;  a[2] = x1 * y1;
  a[3] = x2;  a[4] = y2;  a[5] = x2 * y2;
  a[6] = x3;  a[7] = y3;  a[8] = x3 * y3;

  b[0] = z1 - z0;
  b[1] = z2 - z0;
  b[2] = z3 - z0;

  int status;
  status = linsol(3, a, b);

  double A = z0;
  double B = b[0];
  double C = b[1];
  double D = b[2];

  /* A + B x + C y + D x y = z */

  /* Let (x, y) = (s * x3, s * y3) and z = 0:
   *
   *   A + B * x3 * s + C * y3 * s + D * x3 * y3 * s**2 = 0
   *
   *   A + (B * x3 + C * y3) * s + D * x3 * y3 * s**2 = 0
   *
   **/

  double alpha = D * x3 * y3;
  double beta  = B * x3 + C * y3;
  double gamma = A;

  double disc = beta * beta - 4.0 * alpha * gamma;

  if (disc < 0.0) return 0.5;
  // assert(disc >= 0.0);

  double s = beta < 0.0 ? -1 : 1;
  double Q = beta + s * sqrt(disc);

  double p0 = -Q / (2.0 * alpha);
  double p1 = -2.0 * gamma / Q;

  if (0.0 < p0 && p0 < 1.0) return p0;
  if (0.0 < p1 && p1 < 1.0) return p1;

  return 0.5;
}

static double
bilin(tabr_t *self, int edge_index)
{
  int t1 = self->edge[edge_index].triangle_index[0];
  int t2 = self->edge[edge_index].triangle_index[1];
  int v0 = self->edge[edge_index].vertex_index[0];
  int v1 = self->edge[edge_index].vertex_index[1];

  int en1, en2, nen1, nen2;
  
  for (en1 = 0; en1 < 3; ++en1) if (self->triangle[t1].edge_index[en1] == edge_index) break;
  assert(en1 != 3);
  
  for (en2 = 0; en2 < 3; ++en2) if (self->triangle[t2].edge_index[en2] == edge_index) break;
  assert(en2 != 3);
  
  nen1 = en1 == 2 ? 0 : en1 + 1;
  nen2 = en2 == 2 ? 0 : en2 + 1;
  
  int nei1 = self->triangle[t1].edge_index[nen1];
  int nei2 = self->triangle[t2].edge_index[nen2];
  int neo1 = self->triangle[t1].edge_orientation[nen1];
  int neo2 = self->triangle[t2].edge_orientation[nen2];
  
  int vleft = self->edge[nei1].vertex_index[neo1 == 1 ? 1 : 0];
  int vrite = self->edge[nei2].vertex_index[neo2 == 1 ? 1 : 0];
  
  assert(v1 == self->edge[nei1].vertex_index[neo1 == 1 ? 0 : 1]);
  assert(v0 == self->edge[nei2].vertex_index[neo2 == 1 ? 0 : 1]);
  
  return bilinear_root(self->vertex[v0].position,
		       self->vertex[vleft].position,
		       self->vertex[vrite].position,
		       self->vertex[v1].position);
}

static double
illinois(tabr_t *self, int edge_index)
{
  int v0 = self->edge[edge_index].vertex_index[0];
  int v1 = self->edge[edge_index].vertex_index[1];

  double z0 = self->vertex[v0].position[2];
  double z1 = self->vertex[v1].position[2];

  assert(z0 * z1 < 0.0);
  
  int side = 0;
  double gamma = 1.0;

  while (self->edge[edge_index].parent_edge_index != -1)
    {
      int parent_index = self->edge[edge_index].parent_edge_index;
      if (self->edge[parent_index].daughter_edge_index[0] == edge_index)
	{
	  if (side > 0) break;
	  if (side <= -1) gamma *= 0.5;
	  side--;
	}
      else
	{
	  assert(self->edge[parent_index].daughter_edge_index[1] == edge_index);
	  if (side < 0) break;
	  if (side >= 1) gamma *= 0.5;
	  side++;
	}

      edge_index = parent_index;
    }

  double w0, w1;
  if (side < 0)
    {
      w0 = gamma * z0; w1 = z1;
    }
  else if (side > 0)
    {
      w0 = z0; w1 = gamma * z1;
    }
  else
    {
      w0 = z0; w1 = z1;
    }

  assert(w0 * w1 < 0.0);

  return -w0 / (w1 - w0);
}

/* Return true if edge (v0, v1) should be split. Also given vertices
 * v2 and v3 to the left and right of (v0, v1). */
static int
test_edge(tabr_t *self, 
	  const double v0[3], const double v1[3],
	  const double v2[3], const double v3[3],
	  int splits[5])
{
  if (self->split_root_edges)
    {
      if (v0[2] * v1[2] < 0.0) splits[0] = 1;

      if (v2)
	{
	  if (v0[2] * v2[2] < 0.0 || v1[2] * v2[2] < 0.0) splits[0] = 1;
	}
      
      if (v3)
	{
	  if (v0[2] * v3[2] < 0.0 || v1[2] * v3[2] < 0.0) splits[0] = 1;
	}

      if (v2 && v3)
	{
	  if (v2[2] * v3[2] < 0.0) splits[0] = 1;
	}
    }

  if (v2 == NULL || v3 == NULL) return 0;

  if (self->split_angled_triangles)
    {
      if (cos_angle_triangle_normals(v0, v1, v2, v3) < self->split_angled_cos_angle_limit)
	{
	  splits[0] = 1;
	  splits[1] = splits[2] = 1;
	  splits[3] = splits[4] = 1;
	}
    }

  if (self->split_big_small_triangle_pairs)
    {
      double area_left, area_right, area_min_max_ratio;

      area_left  = triangle3_area(v0, v1, v2);
      area_right = triangle3_area(v0, v1, v3);

      if (area_left > area_right)
	{
	  area_min_max_ratio = area_right / area_left;

	  if (area_min_max_ratio < self->split_big_small_ratio_limit)
	    {
	      splits[1] = splits[2] = 1;
	    }
	}
      else
	{
	  area_min_max_ratio = area_left / area_right;

	  if (area_min_max_ratio < self->split_big_small_ratio_limit)
	    {
	      splits[3] = splits[4] = 1;
	    }
	}

    }

  return 0;
}


/* /\* Return true if edge should be split *\/ */
/* static int */
/* test_edge(tabr_t *self, int edge_index, double *split_position) */
/* { */
/*   int v0 = self->edge[edge_index].vertex_index[0]; */
/*   int v1 = self->edge[edge_index].vertex_index[1]; */

/*   if (self->vertex[v0].has_data == 0 || self->vertex[v1].has_data == 0) return 0; */

/*   /\*   if (fabs(self->vertex[v0].position[0] - self->vertex[v1].position[0]) < self->delta_x_min && *\/ */
/*   /\*       fabs(self->vertex[v0].position[1] - self->vertex[v1].position[1]) < self->delta_y_min) return 0; *\/ */
  
/*   /\* Per edge tests here *\/ */

/*   /\* XXX FOR NOW, JUST DO A ROOT TEST *\/ */

/*   if (self->vertex[v0].position[2] * self->vertex[v1].position[2] < 0.0) */
/*     { */
/*       double bisection = 0.5; */
/*       double root_estimate_illinois; */
/*       double root_estimate_bilinear; */

/*       /\* Illinois root estimate *\/ */

/*       root_estimate_illinois = illinois(self, edge_index); */
/*       if (root_estimate_illinois <= 0.0 || root_estimate_illinois >= 1.0) root_estimate_illinois = 0.5; */

/*       /\* Bilinear interpolation *\/ */

/*       if (self->edge[edge_index].triangle_index[0] != -1 && */
/* 	  self->edge[edge_index].triangle_index[1] != -1) */
/* 	{ */
/* 	  root_estimate_bilinear = bilin(self, edge_index); */
/* 	} */
/*       else */
/* 	{ */
/* 	  root_estimate_bilinear = 0.5; */
/* 	} */

/*       double weight = ((double) self->n_recursion) / (self->n_max_recursion - 1); */
/*       weight *= weight; */

/*       if (self->n_recursion > self->n_max_recursion - N_ROOT_REFINEMENTS) weight = 1.0; */
/*       else weight = 0.0; */

/*       *split_position = (1.0 - weight) * bisection + weight * root_estimate_illinois; */

      
/* 	/\* */
/*       if (self->n_recursion < self->n_max_recursion - 3) */
/*       // if (self->n_recursion % 3 != 2) */
/*       // if (0) */
/*       	{ */
/*       	  *split_position = 0.5; */
/*       	} */
/*       else */
/*       	{ */
/* 	} */
/* 	*\/ */

/*       return 1; */
/*     } */

/*   return 0;   */
/* } */

/* Try to split an edge. relpos sets the relative position of the
 * split, 0 < relpos < 1, zero meaning the first vertex, one the
 * second. */
static int
split_edge_ok(tabr_t *self, int edge_index, double relpos)
{
  assert(self->edge[edge_index].daughter_edge_index[0] == -1);
  assert(self->edge[edge_index].daughter_edge_index[1] == -1);
  assert(0.0 < relpos && relpos < 1.0);

  int vi1 = self->edge[edge_index].vertex_index[0];
  int vi2 = self->edge[edge_index].vertex_index[1];

  double w1 = relpos;
  double w0 = 1.0 - relpos;
  
  double new_x = w0 * self->vertex[vi1].position[0] + w1 * self->vertex[vi2].position[0];
  double new_y = w0 * self->vertex[vi1].position[1] + w1 * self->vertex[vi2].position[1];

  int new_vertex = self->n_vertices;

  int dei1 = self->n_edges;
  int dei2 = self->n_edges + 1;

  if (!add_vertex_ok(self, new_x, new_y)) return 0;

  if (!add_edge_ok(self, vi1, new_vertex, 
		   self->edge[edge_index].triangle_index[0],
		   self->edge[edge_index].triangle_index[1],
		   edge_index)) return 0;
  if (!add_edge_ok(self, new_vertex, vi2,
		   self->edge[edge_index].triangle_index[0],
		   self->edge[edge_index].triangle_index[1],
		   edge_index)) return 0;

  self->edge[edge_index].daughter_edge_index[0] = dei1;
  self->edge[edge_index].daughter_edge_index[1] = dei2;

  return 1;
}

static int
compute_job(tabr_t *self, tabr_job_t *job)
{
  assert(self->n_vertices_pending < self->n_vertices);

  job->type   = TABR_JOB_COMPUTE;
  job->x      = self->vertex[self->n_vertices_pending].position[0];
  job->y      = self->vertex[self->n_vertices_pending].position[1];
  job->result = TABR_JOB_ERROR;
  job->handle = self->n_vertices_pending;
  
  self->n_vertices_pending++;

  return 1;
}

static int
wait_job(tabr_t *self, tabr_job_t *job)
{
  job->type   = TABR_JOB_WAIT;
  job->result = TABR_JOB_SUCCESS;

  self->n_waiting++;

  return 1;
}

/* Split triangle triangle_index over all edges. */
int
split_all_ok(tabr_t *self, int triangle_index)
{
  /*
   *
   *
   *           /\
   *  e0   e00/t0\e21   e2
   *         /____\
   *        /\ c  /\
   *    e01/t1\  /t2\e20
   *      /____\/____\
   *       e10    e11
   *           
   *           e1
   */
  int edge_0 = self->triangle[triangle_index].edge_index[0];
  int edge_1 = self->triangle[triangle_index].edge_index[1];
  int edge_2 = self->triangle[triangle_index].edge_index[2];
  
  int or_0 = self->triangle[triangle_index].edge_orientation[0];
  int or_1 = self->triangle[triangle_index].edge_orientation[1];
  int or_2 = self->triangle[triangle_index].edge_orientation[2];
  
  /* Outer edges; These inherit the orientations of the parent
   * edges, and so the orientation flags or_n can be used to
   * determine the directions of these. */
  
  int edge_00 = self->edge[edge_0].daughter_edge_index[or_0 == 1 ? 0 : 1];
  int edge_01 = self->edge[edge_0].daughter_edge_index[or_0 == 1 ? 1 : 0];
  
  int edge_10 = self->edge[edge_1].daughter_edge_index[or_1 == 1 ? 0 : 1];
  int edge_11 = self->edge[edge_1].daughter_edge_index[or_1 == 1 ? 1 : 0];
  
  int edge_20 = self->edge[edge_2].daughter_edge_index[or_2 == 1 ? 0 : 1];
  int edge_21 = self->edge[edge_2].daughter_edge_index[or_2 == 1 ? 1 : 0];

  assert(edge_00 != -1 && edge_01 != -1);
  assert(edge_10 != -1 && edge_11 != -1);
  assert(edge_20 != -1 && edge_21 != -1);

  
  /* Tip vertices */
  
  int vertex_t0 = self->edge[edge_00].vertex_index[or_0 == 1 ? 0 : 1];
  int vertex_t1 = self->edge[edge_10].vertex_index[or_1 == 1 ? 0 : 1];
  int vertex_t2 = self->edge[edge_20].vertex_index[or_2 == 1 ? 0 : 1];

  assert(vertex_t0 == self->edge[edge_0].vertex_index[or_0 == 1 ? 0 : 1]);
  assert(vertex_t1 == self->edge[edge_1].vertex_index[or_1 == 1 ? 0 : 1]);
  assert(vertex_t2 == self->edge[edge_2].vertex_index[or_2 == 1 ? 0 : 1]);

  assert(vertex_t0 == self->edge[edge_21].vertex_index[or_2 == 1 ? 1 : 0]);
  assert(vertex_t1 == self->edge[edge_01].vertex_index[or_0 == 1 ? 1 : 0]);
  assert(vertex_t2 == self->edge[edge_11].vertex_index[or_1 == 1 ? 1 : 0]);
  
  /* Center vertices */
  
  int vertex_c0 = self->edge[edge_00].vertex_index[or_0 == 1 ? 1 : 0];
  int vertex_c1 = self->edge[edge_10].vertex_index[or_1 == 1 ? 1 : 0];
  int vertex_c2 = self->edge[edge_20].vertex_index[or_2 == 1 ? 1 : 0];

  assert(vertex_c0 == self->edge[edge_01].vertex_index[or_0 == 1 ? 0 : 1]);
  assert(vertex_c1 == self->edge[edge_11].vertex_index[or_1 == 1 ? 0 : 1]);
  assert(vertex_c2 == self->edge[edge_21].vertex_index[or_2 == 1 ? 0 : 1]);
  
  /* Outer, neighbouring triangles */
  
  int triangle_00 = self->edge[edge_00].triangle_index[or_0 == 1 ? 1 : 0];
  int triangle_01 = self->edge[edge_01].triangle_index[or_0 == 1 ? 1 : 0];
  int triangle_10 = self->edge[edge_10].triangle_index[or_1 == 1 ? 1 : 0];
  int triangle_11 = self->edge[edge_11].triangle_index[or_1 == 1 ? 1 : 0];
  int triangle_20 = self->edge[edge_20].triangle_index[or_2 == 1 ? 1 : 0];
  int triangle_21 = self->edge[edge_21].triangle_index[or_2 == 1 ? 1 : 0];
  
  assert(triangle_index == self->edge[edge_00].triangle_index[or_0 == 1 ? 0 : 1]);
  assert(triangle_index == self->edge[edge_01].triangle_index[or_0 == 1 ? 0 : 1]);
  assert(triangle_index == self->edge[edge_10].triangle_index[or_1 == 1 ? 0 : 1]);
  assert(triangle_index == self->edge[edge_11].triangle_index[or_1 == 1 ? 0 : 1]);
  assert(triangle_index == self->edge[edge_20].triangle_index[or_2 == 1 ? 0 : 1]);
  assert(triangle_index == self->edge[edge_21].triangle_index[or_2 == 1 ? 0 : 1]);

  /* New triangles; to be allocated */
  
  int triangle_t0 = self->n_triangles;
  int triangle_t1 = self->n_triangles + 1;
  int triangle_t2 = self->n_triangles + 2;
  int triangle_c  = triangle_index;
  
  /* Inner edges; to be allocated */
  
  int edge_i0 = self->n_edges;
  int edge_i1 = self->n_edges + 1;
  int edge_i2 = self->n_edges + 2;
  
  if (!add_edge_ok(self, vertex_c0, vertex_c2, triangle_t0, triangle_c, -1)) return 0;
  if (!add_edge_ok(self, vertex_c1, vertex_c0, triangle_t1, triangle_c, -1)) return 0;
  if (!add_edge_ok(self, vertex_c2, vertex_c1, triangle_t2, triangle_c, -1)) return 0;
  
  if (!add_triangle_ok(self, edge_00, edge_i0, edge_21, or_0,    1, or_2)) return 0;
  if (!add_triangle_ok(self, edge_01, edge_10, edge_i1, or_0, or_1,    1)) return 0;
  if (!add_triangle_ok(self, edge_11, edge_20, edge_i2, or_1, or_2,    1)) return 0;
  if (!mod_triangle_ok(self, triangle_c, 
		       edge_i0, edge_i1, edge_i2, -1, -1, -1)) return 0;
		       
  
  /* Update edge data */
  
  self->edge[edge_00].triangle_index[or_0 == 1 ? 0 : 1] = triangle_t0;
  self->edge[edge_01].triangle_index[or_0 == 1 ? 0 : 1] = triangle_t1;
  self->edge[edge_10].triangle_index[or_1 == 1 ? 0 : 1] = triangle_t1;
  self->edge[edge_11].triangle_index[or_1 == 1 ? 0 : 1] = triangle_t2;
  self->edge[edge_20].triangle_index[or_2 == 1 ? 0 : 1] = triangle_t2;
  self->edge[edge_21].triangle_index[or_2 == 1 ? 0 : 1] = triangle_t0;

  return 1;
}

/* Split two edges. A new edge needs to be created. */
static int
split_two_ok(tabr_t *self, int triangle_index, int intact_edge_number)
{
  /* Visualize the triangle upright with the intact edge as the
   * base. Need to determine new triangle configuration: Either lower
   * left to upper right cut (LLUR or LU) or upper left to lower right
   * cut (ULLR or UL).
   *
   **/

  /* Intact edge index and orientation */
  int iei = self->triangle[triangle_index].edge_index[intact_edge_number];
  int ieo = self->triangle[triangle_index].edge_orientation[intact_edge_number];

  /* Left/right edge number/index/orientation (new edges below the cut) */
  int len = intact_edge_number == 0 ? 2 : intact_edge_number - 1;
  int ren = intact_edge_number == 2 ? 0 : intact_edge_number + 1;

  int lei_parent = self->triangle[triangle_index].edge_index[len];
  int rei_parent = self->triangle[triangle_index].edge_index[ren];
  int leo = self->triangle[triangle_index].edge_orientation[len];
  int reo = self->triangle[triangle_index].edge_orientation[ren];

  assert(self->edge[lei_parent].daughter_edge_index[0] != -1);
  assert(self->edge[rei_parent].daughter_edge_index[0] != -1);

  int lei = self->edge[lei_parent].daughter_edge_index[leo == 1 ? 1 : 0];
  int rei = self->edge[rei_parent].daughter_edge_index[reo == 1 ? 0 : 1];

  /* Also get left and right edges of the "tip" triangle. */
  int tip_lei = self->edge[lei_parent].daughter_edge_index[leo == 1 ? 0 : 1];
  int tip_rei = self->edge[rei_parent].daughter_edge_index[reo == 1 ? 1 : 0];

  /* Left and right edge splitting vertices. */
  int lvi = self->edge[lei].vertex_index[leo == 1 ? 0 : 1];
  int rvi = self->edge[rei].vertex_index[reo == 1 ? 1 : 0];

  /* Intact edge left and right vertices */
  int ilvi = self->edge[iei].vertex_index[ieo == 1 ? 0 : 1];
  int irvi = self->edge[iei].vertex_index[ieo == 1 ? 1 : 0];

  assert(ilvi == self->edge[lei].vertex_index[leo == 1 ? 1 : 0]);
  assert(irvi == self->edge[rei].vertex_index[reo == 1 ? 0 : 1]);
	 
  /* Tip vertex */
  int tip_vi = self->edge[tip_lei].vertex_index[leo == 1 ? 0 : 1];

  assert(tip_vi == self->edge[tip_rei].vertex_index[reo == 1 ? 1 : 0]);

  /* New triangle indices. Make this triangle the tip triangle, and
   * allocate two more. First comes the one that shares an edge with
   * the tip. */

  int tip_ti = triangle_index;
  int mid_ti = self->n_triangles;
  int bot_ti = self->n_triangles + 1;
  
  /* New edges, one cutting the tip triangle from the bottom,
   * and second cutting the bottom part into the middle and
   * bottom triangles.  */
  int tip_edge = self->n_edges;
  int mid_edge = self->n_edges + 1;

  if (!add_edge_ok(self, lvi, rvi, tip_ti, mid_ti, -1)) return 0;
 
  /* Next figure out how to do the bottom cut. Aim to make as round
   * triangles as possible. */

  int do_lower_left_to_upper_right 
    = dist2_sqr(self->vertex[ilvi].position,
	    self->vertex[rvi].position)
    < dist2_sqr(self->vertex[irvi].position,
	    self->vertex[lvi].position);
  
  if (do_lower_left_to_upper_right)
    {
      /* Cut from LOWER LEFT to UPPER RIGHT */

      if (!add_edge_ok(self, ilvi, rvi, mid_ti, bot_ti, -1)) return 0;

      if (!add_triangle_ok(self, mid_edge, tip_edge, lei, 1, -1, leo)) return 0;
      if (!add_triangle_ok(self, mid_edge, iei, rei, -1, ieo, reo)) return 0;

      /* Update edges with new triangle indices */
      self->edge[tip_lei].triangle_index[leo == 1 ? 0 : 1] = tip_ti;
      self->edge[  lei  ].triangle_index[leo == 1 ? 0 : 1] = mid_ti;
      self->edge[  iei  ].triangle_index[ieo == 1 ? 0 : 1] = bot_ti;
      self->edge[  rei  ].triangle_index[reo == 1 ? 0 : 1] = bot_ti;
      self->edge[tip_rei].triangle_index[reo == 1 ? 0 : 1] = tip_ti;
    }
  else
    {
      /* Cut from UPPER LEFT to LOWER RIGHT */

      if (!add_edge_ok(self, lvi, irvi, mid_ti, bot_ti, -1)) return 0;

      if (!add_triangle_ok(self, mid_edge, rei, tip_edge, 1, reo, -1)) return 0;
      if (!add_triangle_ok(self, mid_edge, lei, iei, -1, leo, ieo)) return 0;

      self->edge[tip_lei].triangle_index[leo == 1 ? 0 : 1] = tip_ti;
      self->edge[  lei  ].triangle_index[leo == 1 ? 0 : 1] = bot_ti;
      self->edge[  iei  ].triangle_index[ieo == 1 ? 0 : 1] = bot_ti;
      self->edge[  rei  ].triangle_index[reo == 1 ? 0 : 1] = mid_ti;
      self->edge[tip_rei].triangle_index[reo == 1 ? 0 : 1] = tip_ti;
    }

  /* Modify this triangle, making it the tip */
  if (!mod_triangle_ok(self, tip_ti, tip_edge, tip_rei, tip_lei, 1, reo, leo)) return 0;


  /* Done */

  return 1;
}

/* Split triangle in two, along line from edge broken_edge_number to
 * the opposing vertex. */
static int
split_one_ok(tabr_t *self, int triangle_index, int broken_edge_number)
{
  int next_edge_number = broken_edge_number == 2 ? 0 : broken_edge_number + 1;
  int prev_edge_number = broken_edge_number == 0 ? 2 : broken_edge_number - 1;
  
  int split_edge_index       = self->triangle[triangle_index].edge_index[broken_edge_number];
  int split_edge_orientation = self->triangle[triangle_index].edge_orientation[broken_edge_number];
  int next_edge_index        = self->triangle[triangle_index].edge_index[next_edge_number];
  int next_edge_orientation  = self->triangle[triangle_index].edge_orientation[next_edge_number];
  int prev_edge_index        = self->triangle[triangle_index].edge_index[prev_edge_number];
  int prev_edge_orientation  = self->triangle[triangle_index].edge_orientation[prev_edge_number];

  int split_left_edge_index  = self->edge[split_edge_index].daughter_edge_index[split_edge_orientation == 1 ? 0 : 1];
  int split_right_edge_index = self->edge[split_edge_index].daughter_edge_index[split_edge_orientation == 1 ? 1 : 0];

  assert(split_left_edge_index != -1);
  assert(split_right_edge_index != -1);
  
  int left_vertex_index     = self->edge[split_edge_index].vertex_index[split_edge_orientation == 1 ? 0 : 1];
  int right_vertex_index    = self->edge[split_edge_index].vertex_index[split_edge_orientation == 1 ? 1 : 0];

  assert(left_vertex_index  == self->edge[prev_edge_index].vertex_index[prev_edge_orientation == 1 ? 1 : 0]);
  assert(right_vertex_index == self->edge[next_edge_index].vertex_index[next_edge_orientation == 1 ? 0 : 1]);
  
  int center_vertex_index   = self->edge[split_left_edge_index].vertex_index[split_edge_orientation == 1 ? 1 : 0];
  int opposing_vertex_index = self->edge[next_edge_index].vertex_index[next_edge_orientation == 1 ? 1 : 0];

  assert(center_vertex_index   == self->edge[split_right_edge_index].vertex_index[split_edge_orientation == 1 ? 0 : 1]);
  assert(opposing_vertex_index == self->edge[prev_edge_index].vertex_index[prev_edge_orientation == 1 ? 0 : 1]);
  
  int left_outer_triangle   = self->edge[prev_edge_index].triangle_index[prev_edge_orientation == 1 ? 1 : 0];
  int right_outer_triangle  = self->edge[next_edge_index].triangle_index[next_edge_orientation == 1 ? 1 : 0];
  int left_bottom_triangle  = self->edge[ split_left_edge_index].triangle_index[split_edge_orientation == 1 ? 1 : 0];
  int right_bottom_triangle = self->edge[split_right_edge_index].triangle_index[split_edge_orientation == 1 ? 1 : 0];

  assert(triangle_index == self->edge[prev_edge_index].triangle_index[prev_edge_orientation == 1 ? 0 : 1]);
  assert(triangle_index == self->edge[next_edge_index].triangle_index[next_edge_orientation == 1 ? 0 : 1]);
  assert(triangle_index == self->edge[ split_left_edge_index].triangle_index[split_edge_orientation == 1 ? 0 : 1]);
  assert(triangle_index == self->edge[split_right_edge_index].triangle_index[split_edge_orientation == 1 ? 0 : 1]);
  
  int splitting_edge = self->n_edges;
  int left_triangle  = self->n_triangles;
  int right_triangle = triangle_index;
  
  if (!add_edge_ok(self, 
		   center_vertex_index, opposing_vertex_index, 
		   left_triangle, right_triangle, -1)) return 0;
  
  if (!add_triangle_ok(self, prev_edge_index, split_left_edge_index, splitting_edge, prev_edge_orientation, split_edge_orientation, 1)) return 0;
  if (!mod_triangle_ok(self, right_triangle, 
		       split_right_edge_index, next_edge_index, splitting_edge,
		       split_edge_orientation, next_edge_orientation, -1)) return 0;

  self->edge[ split_left_edge_index].triangle_index[split_edge_orientation == 1 ? 0 : 1] = left_triangle;
  self->edge[split_right_edge_index].triangle_index[split_edge_orientation == 1 ? 0 : 1] = right_triangle;
  self->edge[prev_edge_index].triangle_index[prev_edge_orientation == 1 ? 0 : 1] = left_triangle;
  self->edge[next_edge_index].triangle_index[next_edge_orientation == 1 ? 0 : 1] = right_triangle;

  return 1;
}


static void
compute_normalized_positions(tabr_t *self, 
			     int *vertex, int n_normalize_vertices,
			     double normalized_position[][3])
{
  int i, j;
  double scale[3];

  /* TODO: Option to choose normalization method. */

  double z_mean = self->z_sum / self->n_vertices;
  double z_sqr_mean = self->z_sqr_sum / self->n_vertices;
  double z_var = z_sqr_mean - z_mean*z_mean;

  if (z_var < 0.0) z_var = 0.0;

  scale[0] = self->x_max - self->x_min;
  scale[1] = self->y_max - self->y_min;
  scale[2] = sqrt(z_var); // self->z_max - self->z_min;

  if (scale[2] == 0.0) scale[2] = 1.0;

  for (i = 0; i < n_normalize_vertices; ++i)
    {
      if (vertex[i] == NO_INDEX) continue;
      for (j = 0; j < 3; ++j)
	{
	  normalized_position[i][j] = self->vertex[vertex[i]].position[j] / scale[j];
	}
    }
}

static int
test_and_split_edge_ok(tabr_t *self, int edge_index)
{
  if (has_daughter_edge(self, edge_index)) return 1;

  int triangle[2];
  int edge[4];
  int edge_orientation[4];
  int vertex[4];
  int splits[5] = {0};

  get_edge_neighbourhood(self, edge_index, 
			 triangle, edge, edge_orientation, vertex);

  double normalized_position[4][3];

  compute_normalized_positions(self, vertex, 4, normalized_position);

  test_edge(self, 
	    normalized_position[0], 
	    normalized_position[1],
	    vertex[2] != NO_INDEX ? normalized_position[2] : NULL, 
	    vertex[3] != NO_INDEX ? normalized_position[3] : NULL,
	    splits);

  if (splits[0]) { if (!split_edge_ok(self, edge_index, 0.5)) return 0; }
  if (splits[1] && !has_daughter_edge(self, edge[0])) { if (!split_edge_ok(self, edge[0], 0.5)) return 0; }
  if (splits[2] && !has_daughter_edge(self, edge[1])) { if (!split_edge_ok(self, edge[1], 0.5)) return 0; }
  if (splits[3] && !has_daughter_edge(self, edge[2])) { if (!split_edge_ok(self, edge[2], 0.5)) return 0; }
  if (splits[4] && !has_daughter_edge(self, edge[3])) { if (!split_edge_ok(self, edge[3], 0.5)) return 0; }

  return 1;
}

static int
refine_triangle_propagate_splits_ok(tabr_t *self, int triangle_index)
{
  int edge_number;
  int n_broken_edges = 0;
  int intact_edge_number;
  
  for (edge_number = 0; edge_number < 3; ++edge_number)
    {
      int edge_index = self->triangle[triangle_index].edge_index[edge_number];
      
      if (self->edge[edge_index].daughter_edge_index[0] != -1)
	{
	  n_broken_edges++;
	}
      else
	{
	  intact_edge_number = edge_number;
	}
    }

  if (n_broken_edges != 2) return 1;

  /* If two edges have been split, split the third one as well. Then
   * move on to the triangle on the other side of the split edge, and
   * see if that has now 2 broken edges. Proceed recursively. */

  int intact_edge_index       = self->triangle[triangle_index].edge_index[intact_edge_number];
  int intact_edge_orientation = self->triangle[triangle_index].edge_orientation[intact_edge_number];

  if (!split_edge_ok(self, intact_edge_index, 0.5)) return 0;
  
  int neighbour_triangle_index = self->edge[intact_edge_index].triangle_index[intact_edge_orientation == 1 ? 1 : 0];

  assert(self->edge[intact_edge_index].triangle_index[intact_edge_orientation == 1 ? 0 : 1] == triangle_index);

  if (neighbour_triangle_index == -1) return 1;

  /* Tail recursive call. */
  return refine_triangle_propagate_splits_ok(self, neighbour_triangle_index);
}


static int
refine_triangle_create_triangles_ok(tabr_t *self, int triangle_index)
{
  int edge_number;
  int n_broken_edges = 0;
  int broken_edge_number;
  int intact_edge_number;
  
  for (edge_number = 0; edge_number < 3; ++edge_number)
    {
      int edge_index = self->triangle[triangle_index].edge_index[edge_number];
      
      if (self->edge[edge_index].daughter_edge_index[0] != -1) 
	{
	  broken_edge_number = edge_number;
	  n_broken_edges++;
	}
      else
	{
	  intact_edge_number = edge_number;
	}
    }

  if (n_broken_edges == 2)
    {
      if (!split_two_ok(self, triangle_index, intact_edge_number)) return 0;
    }
  else if (n_broken_edges == 3)
    {
      if (!split_all_ok(self, triangle_index)) return 0;
    }
  else if (n_broken_edges == 1)
    {
      if (!split_one_ok(self, triangle_index, broken_edge_number)) return 0;
    }

  return 1;
}

/* Given four vertices in counter clockwise order, return true if the
 * resulting polygon is strictly convex (convex and no collinear
 * edges). Some tolerance is allowed by the threshold parameter.  If
 * threshold < 0.0, some concavity is allowed, whereas a positive
 * threshold requires more strict convexity. 
 */
static int
is_strictly_convex(const double v0[2], const double v1[2],
		   const double v2[2], const double v3[2],
		   double threshold)
{
  if (transcross2(v0, v1, v3, 1) <= threshold) return 0;
  if (transcross2(v1, v2, v0, 1) <= threshold) return 0;
  if (transcross2(v2, v3, v1, 1) <= threshold) return 0;
  if (transcross2(v3, v0, v2, 1) <= threshold) return 0;
  return 1;
}

static int
test_and_rotate_edge_ok(tabr_t *self, int edge_index)
{
  if (has_daughter_edge(self, edge_index)) return 1;

  int triangle[2];
  int edge[4];
  int edge_orientation[4];
  int vertex[4];

  get_edge_neighbourhood(self, edge_index, triangle, edge, edge_orientation, vertex);

  if (triangle[0] == NO_INDEX || triangle[1] == NO_INDEX) return 1;

  /* Choose edge so that angle between the triangles is as small as
   * possible.*/

  double normalized_position[4][3];
  compute_normalized_positions(self, vertex, 4, normalized_position);


  double cos_angle_current 
    = cos_angle_triangle_normals(normalized_position[0],
				 normalized_position[1],
				 normalized_position[2],
				 normalized_position[3]);
  double cos_angle_rotated
    = cos_angle_triangle_normals(normalized_position[3],
				 normalized_position[2],
				 normalized_position[0],
				 normalized_position[1]);

  if (cos_angle_current > cos_angle_rotated) return 1;
      
  /* But do not let edge become very long */

  double current_length_sqr, rotated_length_sqr;
  
  current_length_sqr = dist2_sqr(normalized_position[0],
				 normalized_position[1]);
  rotated_length_sqr = dist2_sqr(normalized_position[3],
				 normalized_position[2]);
  
  if (rotated_length_sqr > 5.0 * current_length_sqr) return 1;
  
  /* The neighbourhood must be convex and the rotated edge must not
     coincide with a pre-existing edge */

  if (!is_strictly_convex(self->vertex[vertex[0]].position, 
			  self->vertex[vertex[3]].position,
			  self->vertex[vertex[1]].position, 
			  self->vertex[vertex[2]].position,
			  CONVEXITY_THRESHOLD)) return 1;

  /* Edge is rotated counter-clockwise */
  
  /* Update outer perimeter edge triangle indices */
  self->edge[edge[1]].triangle_index[edge_orientation[1] == 1 ? 0 : 1] = triangle[1];
  self->edge[edge[2]].triangle_index[edge_orientation[2] == 1 ? 0 : 1] = triangle[0];

  /* Update this edge */
  if (!mod_edge_ok(self, edge_index,
		   vertex[3], vertex[2], triangle[0], triangle[1], -1, -1, -1)) return 0;
  
  /* Update triangles */
  if (!mod_triangle_ok(self, triangle[0],
		       edge[2], edge_index, edge[0], edge_orientation[2], 1, edge_orientation[0])) return 0;

  if (!mod_triangle_ok(self, triangle[1],
		       edge[3], edge[1], edge_index, edge_orientation[3], edge_orientation[1], -1)) return 0;

  return 1;
}

int
tabr_get_job(tabr_t *self,
	     tabr_job_t *job)
{
  /* Quick exit if bad/finished state. */
  if (self->state) return 0;

  /* If some data not yet computed, give some of that. Proceed when
   * all pending data is computed. */
  if (self->n_vertices > self->n_vertices_pending) return compute_job(self, job);
  if (self->n_vertices > self->n_computed) return wait_job(self, job);

  /* All data is ready.  */

  int triangle_index;
  int edge_index;


  /* Tidy the mesh. */

  int repeat;
  for (repeat = 0; repeat < self->n_mesh_smoothing_passes; ++repeat)
    {
      for (edge_index = 0; edge_index < self->n_edges; ++edge_index)
	{
	  if (!test_and_rotate_edge_ok(self, edge_index)) goto fail;
	}
    }

  /* Reached maximum recursion? */
  if (self->n_max_recursion == self->n_recursion)
    {
      if (self->n_vertices == self->n_computed) self->state = TABR_SUCCESS;
      return 0;
    }

  /* Mesh could be finer. Refine mesh. */

  int n_triangles_start = self->n_triangles;
  int n_edges_start = self->n_edges;

  /* First pass: Test and split edges. */
  for (edge_index = 0; edge_index < n_edges_start; ++edge_index)
    {
      if (!test_and_split_edge_ok(self, edge_index)) goto fail;
    }


  /* Second pass: Propagate splits. If two edges are split, then the
   * third one is split as well. The propagate split function is
   * recursive. If a new edge is split, the function calls itself with
   * the split-edge neighbour triangle as the argument. */
  if (self->split_propagate_two_cuts)
    {
      for (triangle_index = 0; triangle_index < n_triangles_start; triangle_index++)
	{
	  if (!refine_triangle_propagate_splits_ok(self, triangle_index)) goto fail;
	}
    }

  /* Third pass: Create new triangles. */
  for (triangle_index = 0; triangle_index < n_triangles_start; triangle_index++)
    {
      if (!refine_triangle_create_triangles_ok(self, triangle_index)) goto fail;
    }

  self->n_recursion++;

  /* Mesh refinement done.*/

  return tabr_get_job(self, job);

 fail:
  
  self->state = TABR_ERROR;
  return 0;
}

int
tabr_job_wait(tabr_job_t *job)
{
  return job->type == TABR_JOB_WAIT;
}

void
tabr_job_get_data(tabr_job_t *job,
		  double *x, double *y)
{
  if (job->type != TABR_JOB_COMPUTE)
    {
      job->result = TABR_JOB_ERROR;
      return;
    }

  *x = job->x;
  *y = job->y;
  
}


int
tabr_job_finish_wait(tabr_t *self,
		     tabr_job_t *job)
{
  if (job->type != TABR_JOB_WAIT)
    {
      self->state = TABR_ERROR;
      return 1;
    }

  assert(self->n_waiting > 0);

  if (self->n_waiting < 1)
    {
      self->state = TABR_ERROR;
      return 1;
    }

  self->n_waiting--;  

  return 0;
}
		     

int
tabr_job_finish(tabr_t *self,
		tabr_job_t *job,
		const double *z)
{
  assert(job->type == TABR_JOB_COMPUTE);

  if (job->type != TABR_JOB_COMPUTE || z == NULL)
    {
      self->state = TABR_ERROR;
      return 1;
    }

  int i_vertex = job->handle;
  
  assert(self->vertex[i_vertex].position[0] == job->x);
  assert(self->vertex[i_vertex].position[1] == job->y);
  
  memcpy(self->vertex[i_vertex].position + 2, z, sizeof(double) * self->n_z);
  self->vertex[i_vertex].has_data    = 1;

  if (self->n_computed == 0)
    {
      self->z_min = self->z_max = z[0];
    }
  else
    {
      if (z[0] > self->z_max) self->z_max = z[0];
      else if (z[0] < self->z_min) self->z_min = z[0];
    }

  self->z_sum     += z[0];
  self->z_sqr_sum += z[0] * z[0];
  
  self->n_computed++;

  return (self->n_computed == self->n_vertices) && self->n_waiting;
}

int
tabr_sanity_check_ok(tabr_t *self)
{
  /* Things that can be checked:
   *
   * 1. Edges form closed path
   * 2. Edges have correct triangle indices
   * 
   */
  int i;

  for (i = 0; i < self->n_triangles; ++i)
    {
      /* 1. Closed path test: First edge ends where second begins etc.  */
      int e1 = self->triangle[i].edge_index[0];
      int e2 = self->triangle[i].edge_index[1];
      int e3 = self->triangle[i].edge_index[2];
      int o1 = self->triangle[i].edge_orientation[0];
      int o2 = self->triangle[i].edge_orientation[1];
      int o3 = self->triangle[i].edge_orientation[2];

      assert(/* Edge 1 ends where edge 2 begins */ self->edge[e1].vertex_index[o1 == 1 ? 1 : 0] == self->edge[e2].vertex_index[o2 == 1 ? 0 : 1]);
      assert(/* Edge 2 ends where edge 3 begins */ self->edge[e2].vertex_index[o2 == 1 ? 1 : 0] == self->edge[e3].vertex_index[o3 == 1 ? 0 : 1]);
      assert(/* Edge 3 ends where edge 1 begins */ self->edge[e3].vertex_index[o3 == 1 ? 1 : 0] == self->edge[e1].vertex_index[o1 == 1 ? 0 : 1]);

      /* 2. Edges have correct triangle indices: Edges of triangle i should have links back to the same triangle i. */

      assert(/* Edge 1 has link to this triangle */ self->edge[e1].triangle_index[0] == i || self->edge[e1].triangle_index[1] == i);
      assert(/* Edge 2 has link to this triangle */ self->edge[e2].triangle_index[0] == i || self->edge[e2].triangle_index[1] == i);
      assert(/* Edge 3 has link to this triangle */ self->edge[e3].triangle_index[0] == i || self->edge[e3].triangle_index[1] == i);
    }

  return 1;
}

void
tabr_dump_vertex_data(tabr_t *self, FILE *file,
		      const char *sep, const char *lf)
{
  int i;
  int j;

  for (i = 0; i < self->n_vertices; ++i)
    {
      fprintf(file, "%24.16le%s%24.16le", 
	      self->vertex[i].position[0], sep,
	      self->vertex[i].position[1]);

      for (j = 0; j < self->n_z; ++j)
	{
	  fprintf(file, "%s%24.16le", 
		  sep, self->vertex[i].position[2 + j]);
	}

      fprintf(file, "%s", lf);
    }
}

/*
void 
tabr_dump_edge_data(tabr_t *self, FILE *file,
		      const char *sep, const char *lf)
{
  int i;

  for (i = 0; i < self->n_edges; ++i)
    {
      if (has_daughter_edge(self, i)) continue;

      fprintf(file, "%6d%s%6d%s",
	      self->edge[i].vertex_index[
    }
}
*/

void
tabr_dump_triangle_data(tabr_t *self, FILE *file,
			const char *sep, const char *lf)
{
  int i;
  int vertex[3];

  for (i = 0; i < self->n_triangles; ++i)
    {
      get_triangle_vertex_indices(self, i, vertex);
      
      fprintf(file, "%6d%s%6d%s%6d%s",
	      vertex[0], sep,
	      vertex[1], sep,
	      vertex[2], lf);
    }
}
 
