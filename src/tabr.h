
#ifndef TABR_H
#define TABR_H

#include <stdio.h>

/* General idea: Evaluate z = f(x, y) on a mesh of (x, y) values,
 * where the (x, y) points are chosen adaptively: Points should be
 * dense where e.g. f(x, y) changes rapidly, passes through a critical
 * value, etc.
 *
 * Here, z can be a vector, though the default mesh refinement
 * functions only refer to the first value in z.
 *
 *
 * MULTITHREADING
 * 
 * The adaptive sampling method should allow for several threads to
 * work on the problem concurrently: A worker thread asks (while
 * holding say a mutex lock) for new (x, y) pair. These are supplied
 * by the mesh sampling object, along with a handle identifying the
 * request. These values are internally marked as pending for
 * data. The worker thread then evaluates z = f(x, y) and (again
 * locking a mutex) provides the output along with the handle to the
 * adaptive sampler object. This continues until no new (x, y) pairs
 * can be supplied (everything is done, or an error occurred).
 *
 * Under the hood, a x,y data request may trigger a mesh refinement.
 * This occurs if all z values have been computed and the mesh has not
 * been refined upto requested fineness. The thread making the data
 * request alone computes the refinement (it should hold a mutex lock
 * and block requests from other threads).
 *
 * Mesh refinement may be blocked: If all x,y values have been passed
 * to worker threads, but not all z values are ready yet, the
 * refinement cannot take place. In this case, instead of x,y data, a
 * WAIT request is returned. The waiting thread should ideally use a
 * condition variable to wait for a signal to try again. The return
 * value of the function used to pass ready z values to the refinement
 * object indicates if waiting threads should be signaled to wake up.
 *
 * If POSIX threads were enabled at compile time, tabr_threadman.h has
 * functions for automating the thread/data management. Also see
 * examples directory.
 *
 *
 * SINGLE-THREADED CODE
 *
 * All that single-threaded code needs to do is to set up a simple
 * loop that pulls (x,y) values, computes z = f(x,y), and passes the
 * result to the mesh refinement structure. See examples directory for
 * basic usage.
 *
 *
 */

/* 
 * TODO: 
 *
 * 1. Custom edge/triangle split tests
 * 2. Better failed evaluation handling
 * 3. Accessor functions for the computed data. Now just 
 *    data dump functions.
 * 4. Control for the normalization of data entering split
 *    tests. Also: Current normalization code has a lot of redundant
 *    computation. Not pretty.
 * 5. Parametric surfaces: (x, y, z, ...) = f(u, v).
 * 5. There are some hard-coded parameter values, e.g. in the edge
 *    rotation function (test_and_rotate_edge_ok, triangle roundness 
 *    and polygon convexity tests).
 * 6. Better documentation
 * 7. Speed optimization.
 *
 */

/* Mesh vertex structure. */
typedef struct {
  double *position; /*< Vertex x, y, and z coordinates, z may be
  		        vector, position data is then x, y, z[0],
  		        z[1], ..., z[n-1] array where n is the number
  		        of elements, specified in the creation of the
  		        mesh refinement object. */
  int has_data;
} tabr_vertex_t;

typedef struct {
  int vertex_index[2];         /*< Indices to the vertex table,
				   describing the start and end points
				   of the edge. */
  int triangle_index[2];       /*< Indices to the triangles that use
				   this edge. When looking from vertex
				   0 to 1, triangle 0 is on the left
				   and 1 on the right. Set to -1 if no
				   triangles on that side of the
				   edge. This data is used to navigate
				   triangles and edges, finding
				   neighbours, etc. */
  int daughter_edge_index[2];  /*< When an edge is split, two daughter
				   edges are formed. This array
				   contains the edge table indices of
				   those edges. Set to -1 if no
				   daughters. */
  int parent_edge_index;       /*< Index of the parent edge. If this
				   edge was created by splitting a
				   pre-existing edge, this index links
				   back to that edge; -1 if
				   otherwise. This may be removed in
				   the future. */
} tabr_edge_t;

typedef struct {
  int edge_index[3];           /*< Indices to the edge table. Edges
				   are to be listed in
				   counter-clockwise order. */
  int edge_orientation[3];     /*< Is 1 if vertices of the edge are
				   listed in counter-clocwise order;
				   -1 if otherwise. */
} tabr_triangle_t;


/* Adaptive mesh refinement algorithm state. */
typedef enum {
  TABR_RUNNING = 0,
  TABR_SUCCESS,
  TABR_ERROR
} tabr_state_t;

/* Adaptive mesh struct, holding the mesh itself plus all the
 * variables controlling the refinement process. */
typedef struct {
  int n_z;              /*< Length of the z ( = f(x, t) ) data vector. */
  int n_vertices;       /*< Number of vertices currently in the mesh */
  int n_edges;          /*< Number of edges currently in the mesh */
  int n_triangles;      /*< Number of triangles currently in the mesh */

  int n_max_vertices;   /*< Max. number of vertices the vertex table
  			    can accomodate. Increased as needed. */ 
  int n_max_edges;      /*< Max. number of edges the edge table can
			    accomodate. Increased as needed. */
  int n_max_triangles;  /*< Max. number of triangles the triangle
			    table can accomodate. Increased as
			    needed. */

  tabr_vertex_t   *vertex;    /*< Vertex table, holding the mesh vertices. */
  tabr_edge_t     *edge;      /*< Edge table, holding the mesh edges. */
  tabr_triangle_t *triangle;  /*< Triangle table, holding the mesh triangles. */

  int n_vertices_pending;     /*< Number of vertices scheduled to be
				  computed, but which have not yet
				  been completed. */
  
  double delta_y_max;
  double delta_x_max;
  double delta_y_min;
  double delta_x_min;

  int split_root_edges;
  int split_big_small_triangle_pairs;
  int split_angled_triangles;
  int split_propagate_two_cuts;
  double split_big_small_ratio_limit;
  double split_angled_cos_angle_limit;

  int n_mesh_smoothing_passes;

  int n_max_recursion;      /*< Maximum number of mesh refinement
			        recursions. */
  int n_recursion;          /*< Number of refinement recursions
			        done. */

  double x_min;
  double x_max;
  double y_min;
  double y_max;
  double z_min;
  double z_max;

  double z_sum;
  double z_sqr_sum;

  int n_max_waiting;
  int n_waiting;           /*< Number of WAIT jobs returned but not
			       finished. */
  int n_computed;          /*< Number of z = f(x,y) evaluations
			       done. */

  tabr_state_t state;

} tabr_t;


/* Allocate a new mesh refinement object. Takes as arguments the x and
 * y ranges, the number of initial points in each direction, and the
 * number of mesh refinement recursions to do. Lastly, n_z specifies
 * the dimension of the result z = f(x, y).  After allocation, one may
 * optionally fine-tune the algorithm using the
 * tabr_set_... functions. After that, the tabr_t object is ready. */
tabr_t *
tabr_alloc(double x_min, double x_max, int x_n,
	   double y_min, double y_max, int y_n,
	   int n_rec, int n_z);


/* Set root searching on or off. If on_off = 1, the algorithm will
 * split all edges that cross a root, allowing the mesh become finer
 * around zeros of the function. */
void
tabr_set_refine_roots(tabr_t *self, int on_off);


/* Set the maximum angle between adjacent triangles. Setting angle >
 * pi turns off angle testing. */
void
tabr_set_max_angle(tabr_t *self, double angle);


/* Set the number of smoothing passes the mesh refinement algorithm
 * does. After the mesh is completely evaluated, and before the mesh
 * is refined, the algorithm attempts to rotate edges so that the
 * resulting mesh is as smooth as possible. Setting smoothing passes
 * to zero skips this; more than one smoothing pass may improve the
 * grid. 
 *
 * This may have unpredictable effects, and may also lead to
 * artifacts.
 */
void
tabr_set_mesh_smoothing_passes(tabr_t *self, int n_mesh_smoothing_passes);


/* Set the behaviour for splitting edges of triangles with two cut
 * edges. Setting propagate_splits = 1 causes the algorithm to force
 * splitting of the third edge of a triangle with two split edges.  If
 * propagate_splits = 0, a new edge is created from one of the corners
 * to the one of the midpoints of the split edges. */
void
tabr_set_propagate_splits(tabr_t *self, int propagate_splits);


/* Destroy the mesh refinement structure. */
void
tabr_free(tabr_t *self);


void
tabr_dump_vertex_data(tabr_t *self, FILE *file,
		      const char *sep, const char *lf);

void
tabr_dump_triangle_data(tabr_t *self, FILE *file,
			const char *sep, const char *lf);


/* Return true, if the computation has successfully finished. False
 * return value indicates either unfinished computation, or an
 * error. */
int
tabr_success(tabr_t *self);

/* Return true if mesh is internally consistent. This function is for
 * debugging purposes, and simply contains asserts that check out if
 * all is well. */
int
tabr_sanity_check_ok(tabr_t *self);


typedef enum {
  TABR_JOB_COMPUTE,
  TABR_JOB_WAIT
} tabr_job_type_t;

typedef enum {
  TABR_JOB_SUCCESS = 0,
  TABR_JOB_ERROR
} tabr_job_result_t;


typedef struct {
  tabr_job_type_t type;
  tabr_job_result_t result;
  double x;
  double y;
  int handle;
} tabr_job_t;


/* Try to get new x,y data. Returns false if all computation is
 * done. Otherwise, a new task is given in job. In single threaded
 * applications, job is always a new x, y pair, values of which can be
 * extracted with tabr_job_get_data. In multithreaded code, the job
 * may be a WAIT, meaning that missing data is blocking further
 * progress and the thread must wait. A wait job is tested using
 * tabr_job_wait. The wait may be a fixed time, but preferrably a
 * condition variable. Multithreaded code should check the return
 * value of tabr_job_finish, which will indicate whether waiting
 * threads can be signalled to wake up. 
 *
 * All jobs need to be finished either by calling tabr_job_finish or
 * tabr_job_finish_wait.
 *
 */
int
tabr_get_job(tabr_t *self,
	     tabr_job_t *job);


/* Return true if job is a wait job, meaning that the calling thread
 * must wait until more data has arrived. If only a single thread is
 * making calls to tabr_get_job, a job is never a wait. */
int
tabr_job_wait(tabr_job_t *job);

void
tabr_job_get_data(tabr_job_t *job,
		  double *x, double *y);

/* Finish a job. Updates the tabr object with data from the job
 * object, return true if threads waiting for data should be signalled
 * to wake up. Single-threaded applications can of course neglect
 * this. */
int
tabr_job_finish(tabr_t *self,
		tabr_job_t *job,
		const double *z);

/* Finish a wait job. Updates the tabr object so that the number of
 * waiting threads is current. As tabr_finish_job, will return true if
 * threads need to be signaled, though this only occurs in a case of
 * unrecoverable errors (something has gone wrong, and all other
 * waiting threads must wake up and exit). */
int
tabr_job_finish_wait(tabr_t *self,
		     tabr_job_t *job);



#endif
