
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tabr.h"
#include "tabr_threadman.h"


int
fun(double x, double y, double *z, const void *params)
{
  double alpha = *(const double *)params;
  double envelope;

  envelope = exp(-alpha * (x * x + y * y));
  *z = envelope * exp(cos(5 * M_PI * x) * sin(2 * M_PI * y));

  return 0;
}


int
main(int argc, char **argv)
{
  int ok, status;
  int n_threads;
  tabr_t *tabr;
  tabr_job_t job;
  tabr_threadman_t *threadman;

  double alpha;

  /* Initialization */

  alpha = 4.0;

  n_threads = 3;

  tabr = tabr_alloc(-1.0, +1.0, 10,
		    -1.0, +1.0, 10,
		    6, 1);
  if (tabr == NULL)
    {
      fprintf(stderr, "Failed creating mesh refinement object\n");
      goto fail;
    }

  threadman = tabr_threadman_alloc(n_threads, TABR_THREADMAN_CREATE_THREADS);
  if (threadman == NULL)
    {
      fprintf(stderr, "Failed creating thread manager object\n");
      goto fail;
    }


  tabr_set_max_angle(tabr, M_PI / 20.0);
  tabr_threadman_set_logstream(threadman, stdout);


  /* Main evaluation loop */

  ok = tabr_threadman_autorun(threadman, tabr, fun, &alpha);
  if (!ok)
    {
      fprintf(stderr, "Adaptive mesh evaluation of the function failed\n");
      goto fail;
    }


  /* Clean-up and exit. */
  
  tabr_threadman_free(threadman);
  tabr_free(tabr);

  return 0;

 fail:

  printf("fail.\n");
  return -1;
}


