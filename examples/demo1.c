
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tabr.h"


double
fun(double x, double y)
{
  const double envelope = exp(-4.0 * (x * x + y * y));
  return envelope * exp(cos(5 * M_PI * x) * sin(2 * M_PI * y));
}


int
main(int argc, char **argv)
{
  tabr_t *tabr;
  tabr_job_t job;


  /* Initialization */

  tabr = tabr_alloc(-1.0, +1.0, 10,
		    -1.0, +1.0, 10,
		    6, 1);
  if (tabr == NULL)
    {
      fprintf(stderr, "Failed creating mesh refinement object\n");
      goto fail;
    }

  tabr_set_max_angle(tabr, M_PI / 20.0);


  /* Main evaluation loop */

  while (tabr_get_job(tabr, &job))
    {
      double x, y, z;
      tabr_job_get_data(&job, &x, &y);

      z = fun(x, y);

      printf(" %25.15lf %25.15lf %25.15lf\n", x, y, z);

      tabr_job_finish(tabr, &job, &z);
    }


  /* Clean-up and exit. */
  
  tabr_free(tabr);

  return 0;

 fail:

  printf("fail.\n");
  return -1;
}


