
#include <stdlib.h>

#include "tabr_threadman.h"

tabr_threadman_t *
tabr_threadman_alloc(int n_threads, int create_threads)
{
  tabr_threadman_t *self;

  self = malloc(sizeof(*self));
  if (self == NULL) return NULL;

  self->logstream = NULL;

  self->create_threads = create_threads;
  self->n_threads = n_threads;
  self->job = malloc(n_threads * sizeof(*self->job));
  if (self->job == NULL)
    {
      free(self);
      return NULL;
    }
  
  if (create_threads)
    {
      self->thread = malloc(n_threads * sizeof(*self->thread));
      if (self->thread == NULL)
	{
	  free(self->job);
	  free(self);
	  return NULL;
	}

      self->pars = malloc(n_threads * sizeof(*self->pars));
      if (self->pars == NULL)
	{
	  free(self->thread);
	  free(self->job);
	  free(self);
	  return NULL;
	}
    }
  else
    {
      self->thread = NULL;
      self->pars   = NULL;
    }

  pthread_mutex_init(&self->mutex, NULL);
  pthread_cond_init(&self->wait_condvar, NULL);

  return self;
}

void
tabr_threadman_free(tabr_threadman_t *self)
{
  pthread_mutex_destroy(&self->mutex);
  pthread_cond_destroy(&self->wait_condvar);

  if (self->logstream) fprintf(self->logstream, "\n\n\n");

  if (self->thread) free(self->thread);
  if (self->pars)   free(self->pars);

  free(self->job);
  free(self);
}

void
tabr_threadman_set_logstream(tabr_threadman_t *self, FILE *stream)
{
  self->logstream = stream;
}


int
tabr_threadman_get_data(tabr_threadman_t *manager, 
			tabr_t *tabr,
			int thread_number,
			double *x, double *y)
{
  pthread_mutex_lock(&manager->mutex);

  int got_job;
  tabr_job_t *job = &manager->job[thread_number];

  while ((got_job = tabr_get_job(tabr, job)) && tabr_job_wait(job))
    {
      pthread_cond_wait(&manager->wait_condvar, &manager->mutex);
      tabr_job_finish_wait(tabr, job);
    }

  pthread_mutex_unlock(&manager->mutex);

  if (got_job)
    {
      tabr_job_get_data(job, x, y);
    }

  return got_job;
}


void
tabr_threadman_finish(tabr_threadman_t *manager,
		      tabr_t *tabr,
		      int thread_number,
		      int error_code,
		      const double *z)
{
  pthread_mutex_lock(&manager->mutex);

  tabr_job_t *job = &manager->job[thread_number];

  if (manager->logstream)
    {
      fprintf(manager->logstream, " %25.15le %25.15le", job->x, job->y);
      int i;
      for (i = 0; i < tabr->n_z; ++i)
	{
	  fprintf(manager->logstream, " %25.15le", z[i]);
	}
      fprintf(manager->logstream, "\n");
    }

  int wakeup = tabr_job_finish(tabr, job, z);
  
  if (wakeup)
    {
      pthread_cond_broadcast(&manager->wait_condvar);
    }

  pthread_mutex_unlock(&manager->mutex);
}

static void *
worker(void *ptr)
{
  int status;
  tabr_threadman_worker_pars_t *pars = ptr;
  double x, y;

  /* Master thread holds the mutex lock while starting workers. All
   * workers must wait until the master has finished the
   * initialization. */
  pthread_mutex_lock(&pars->manager->mutex);
  status = !pars->manager->autorun_init_ok;
  pthread_mutex_unlock(&pars->manager->mutex);

  /* If autorun init failed when starting later threads, just exit. */
  if (status)
    {
      pthread_exit(NULL);
    }

  while (tabr_threadman_get_data(pars->manager, 
				 pars->tabr,
				 pars->thread_number,
				 &x, &y))
    {
      if (pars->func)
	{
	  status = pars->func(x, y, pars->z, pars->pars);
	}
      else
	{
	  status = pars->func_environ(x, y, pars->z, pars->environ, pars->pars);
	}

      tabr_threadman_finish(pars->manager, 
			    pars->tabr, 
			    pars->thread_number, 
			    status, 
			    pars->z);
    }

  pthread_exit(NULL);
}


static int 
autorun(tabr_threadman_t *manager,
	tabr_t *tabr,
	tabr_threadman_func_t func,
	tabr_threadman_func_local_t func_environ,
	tabr_threadman_alloc_local_t alloc_environ,
	tabr_threadman_free_local_t free_environ,
	const void *globvar)
{
  if (!manager->create_threads) return 0;

  manager->autorun_init_ok = 0;

  int status;
  pthread_attr_t thread_attr;
  int n_z = tabr->n_z;

  pthread_attr_init(&thread_attr);
  pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);


  /* The spawned worker threads will start by trying to lock the
   * mutex. By holding the mutex, this main thread can block the
   * workers until initialization is finished. If the initialization
   * fails for whatever reasons, the workers then quit gracefully by
   * checking the state of the thread manager. */

  pthread_mutex_lock(&manager->mutex);

  int i_thread;
  for (i_thread = 0; i_thread < manager->n_threads; ++i_thread) 
    {
      manager->pars[i_thread].environ = NULL;
      manager->pars[i_thread].z       = NULL;
    }

  for (i_thread = 0; i_thread < manager->n_threads; ++i_thread)
    {
      manager->pars[i_thread].thread_number = i_thread;
      manager->pars[i_thread].manager = manager;
      manager->pars[i_thread].tabr = tabr;
      manager->pars[i_thread].pars = globvar;
      manager->pars[i_thread].func = func;
      manager->pars[i_thread].z    = malloc(sizeof(double) * n_z);

      if (manager->pars[i_thread].z == NULL)
	{
	  manager->autorun_init_ok = 0;
	  break;
	}
      
      if (func_environ && alloc_environ)
	{
	  manager->pars[i_thread].func_environ = func_environ;
	  manager->pars[i_thread].environ = alloc_environ(globvar);

	  if (manager->pars[i_thread].environ == NULL)
	    {
	      manager->autorun_init_ok = 0;
	      break;
	    }
	}

      status = pthread_create(&manager->thread[i_thread], 
			      &thread_attr, 
			      worker, (void *) &manager->pars[i_thread]);
      if (status)
	{
	  manager->autorun_init_ok = 0;
	  break;
	}
    }
  
  int n_threads_started = i_thread;

  if (n_threads_started == manager->n_threads) manager->autorun_init_ok = 1;
  
  pthread_mutex_unlock(&manager->mutex);

  for (i_thread = 0; i_thread < n_threads_started; ++i_thread)
    {
      void *retval;
      pthread_join(manager->thread[i_thread], &retval);
    }
  
  for (i_thread = 0; i_thread < manager->n_threads; ++i_thread) 
    {
      if (manager->pars[i_thread].z)
	{
	  free(manager->pars[i_thread].z);
	}
      
      if (manager->pars[i_thread].environ)
	{
	  free_environ(manager->pars[i_thread].environ);
	}
    }
  
  
  if (!manager->autorun_init_ok) return 0;
  
  return tabr_success(tabr);
  
}
	
int
tabr_threadman_autorun(tabr_threadman_t *manager,
		       tabr_t *tabr,
		       tabr_threadman_func_t func,
		       const void *pars)
{
  return autorun(manager, tabr, 
		 func, 
		 NULL, NULL, NULL, 
		 pars);
}
	

int
tabr_threadman_autorun_locals(tabr_threadman_t *manager,
			      tabr_t *tabr,
			      tabr_threadman_func_local_t func_environ,
			      tabr_threadman_alloc_local_t alloc_environ,
			      tabr_threadman_free_local_t free_environ,
			      const void *globvar)
{
  return autorun(manager, tabr, 
		 NULL, 
		 func_environ, alloc_environ, free_environ, 
		 globvar);
}




