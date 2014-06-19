

#ifndef TABR_THREADMAN_H
#define TABR_THREADMAN_H

#include <stdio.h>
#include <pthread.h>

#include "tabr.h"

/* Routines for automatic thread management using POSIX threads.
 *
 * The user merely needs to supply a function that evaluates z = f(x,
 * y, g), where g is some global work data shared by all threads, OR a
 * function z = f(x, y, l, g), where l is per thread work data,
 * allocated by a user supplied function.
 *
 * Alternatively, the user can also create his own threads, and use
 * the thread manager functions to fetch input x,y data from the mesh
 * refinement object. The manager takes care of mutex locking and
 * waiting on condition variables.
 *
 */

/* Most basic user supplied f : z = f(x, y) function type. Takes as
 * input the x and y values, and evaluates the result to z. The
 * function may depend on additional parameters which are passed
 * through the global_pars parameter (user makes sure that its
 * contents is not modified); this may be NULL if no additional
 * parameters are needed.  
 */
typedef int (*tabr_threadman_func_t)(double x, double y, double *z, 
				     const void *global_pars);


/* Improved user supplied f : z = f(x, y). In addition to depending on
 * some global parameters, it is supposed that the function needs some
 * local workspace in order to evaluate. Since multiple function
 * evaulations are occurring concurrently, each thread needs its own
 * workspace.  This is passed via the local_work parameter. The
 * function types tabr_threadman_alloc_local_t and
 * tabr_threadman_free_local_t allocate and free these local
 * workspaces, and are used by the autorun functions described below. 
 */
typedef int (*tabr_threadman_func_local_t)(double x, double y, double *z, void *local_work, const void *global_work);
typedef void* (*tabr_threadman_alloc_local_t)(const void *global_pars);
typedef void  (*tabr_threadman_free_local_t)(void *local_work);

struct tabr_threadman;
typedef struct tabr_threadman tabr_threadman_t;

typedef struct {
  int thread_number;
  tabr_threadman_t *manager;
  tabr_t *tabr;
  tabr_threadman_func_t func;
  tabr_threadman_func_local_t func_environ;
  const void *pars;
  void *environ;
  double *z;
} tabr_threadman_worker_pars_t;

struct tabr_threadman {
  pthread_mutex_t mutex;
  pthread_cond_t wait_condvar;

  FILE *logstream;

  int create_threads;
  int n_threads;
  tabr_job_t *job;

  int autorun_init_ok;
  pthread_t *thread;
  tabr_threadman_worker_pars_t *pars;
};

#define TABR_THREADMAN_CREATE_THREADS 1
#define TABR_THREADMAN_DO_NOT_CREATE_THREADS 0

/* Allocate a thread management object. Here n_threads is the number
 * of threads to use, and create_threads = 1, if the manager takes
 * care of creating the threads; = 0 if user starts the threads
 * herself. */
tabr_threadman_t *
tabr_threadman_alloc(int n_threads, int create_threads);


/* Deallocate the thread manager. */
void
tabr_threadman_free(tabr_threadman_t *self);


/* Turn on data streaming. If the log stream is set, the threadman
 * manager will send the computed x,y,z data to the specified stream
 * in plain ASCII format. This can be useful for e.g. monitoring the
 * progress of the computation. By default, streaming is off; this can
 * also be done by giving NULL as the stream parameter. This function
 * should be called prior to starting worker threads or calling the
 * autorun functions. 
 *
 * Output is written when the threadmanagers mutex is locked; it is up
 * to the user to make sure the file is not in use by other threads.
 **/
void
tabr_threadman_set_logstream(tabr_threadman_t *self, FILE *stream);


/* Get new data using the thread manager object. Takes as input
 * arguments the manager and mesh refinement objects, plus the number
 * of the calling thread (used to identify the data requests, 0 <=
 * thread_number < n_threads, specified at tabr_threadman_alloc).  On
 * output, if return value was non-zero, x and y contain the new input
 * x, y values. If return value is zero, computation is done (either
 * due completion or error). The thread manager takes care of locking
 * mutexes, and waiting on condition variables. 
 */
int
tabr_threadman_get_data(tabr_threadman_t *manager, 
			tabr_t *tabr,
			int thread_number,
			double *x, double *y);


/* Pass computed data to the mesh refinement object, using the thread
 * manager. As inputs, the thread manager and mesh refinement objects,
 * the number of the calling thread, a flag indicating error (1) or
 * success (0), and the computed value z. 
 */
void
tabr_threadman_finish(tabr_threadman_t *manager,
		      tabr_t *tabr,
		      int thread_number,
		      int error,
		      const double *z);


/* Do a full automated mesh refinement run with given function z =
 * f(x, y). The function f is given as
 *
 * int func(double x, double y, double *z, const void *pars),
 *
 * where x, y are the inputs and z is the value to be computed, and
 * pars are some global parameters for the function (may be NULL).  It
 * is user's responsibility to make sure that what is pointed to by
 * pars is not modified during the run. 
 *
 * If the function needs to allocate per thread workspace (e.g. extra
 * storage that can be modified by func), then one should use the
 * function tabr_threadman_autorun_locals instead. 
 *
 * Return zero for failure.
 */
int
tabr_threadman_autorun(tabr_threadman_t *manager,
		       tabr_t *tabr,
		       tabr_threadman_func_t func,
		       const void *pars);


/* Do a full automated mesh refinement run with given function z = f(x, y).
 * The function f is given as
 *
 * int func(double x, double y, double *z, void *local_work, const void *global_pars),
 *
 * where x, y are the inputs and z is the value to be computed;
 * local_work is per thread workspace that can be modified during the
 * functions evaluation, global_pars are global parameters (user must
 * make sure global_pars are not modified during execution). Return
 * zero for success, non-zero for failure.
 *
 * The local_work structures are allocated and destroyed automatically
 * with the functions alloc_local and free_local.
 *
 * Return zero for failure.
 */
int
tabr_threadman_autorun_locals(tabr_threadman_t *manager,
			      tabr_t *tabr,
			      tabr_threadman_func_local_t func_local,
			      tabr_threadman_alloc_local_t alloc_local,
			      tabr_threadman_free_local_t free_local,
			      const void *globvar);
#endif
