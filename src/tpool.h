#ifndef CAMIL_TPOOL_H
#define CAMIL_TPOOL_H


#ifdef __cplusplus
extern "C" {
#endif

#include "logger.h"
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <stddef.h>


struct tpool_job {
	void (*fn)(void *arg);
	void *arg;
};

struct tpool {
	pthread_t *threads;       // worker threads, NULL in inline mode
	int nthreads;             // number of workers, 0 means inline mode
	struct tpool_job *queue;  // circular buffer of pending jobs
	size_t capacity;          // number of slots in queue
	size_t head;              // index of the next job to pop
	size_t count;             // number of pending jobs
	size_t running;           // number of jobs currently executing
	int shutdown;             // set by tpool_destroy() to release workers
	pthread_mutex_t mutex;    // guards every field above
	pthread_cond_t not_empty; // signalled when a job is queued
	pthread_cond_t not_full;  // signalled when a queue slot frees up
	pthread_cond_t idle;      // signalled when queue drains and jobs finish
};

// starts a pool with nthreads workers and a queue of capacity jobs
int tpool_init(struct tpool *pool, int nthreads, size_t capacity);

// queues fn(arg) for execution, blocking while the queue is full
void tpool_submit(struct tpool *pool, void (*fn)(void *), void *arg);

// blocks until every submitted job has finished
void tpool_wait(struct tpool *pool);

// waits for outstanding jobs, joins the workers and releases all resources
void tpool_destroy(struct tpool *pool);

// number of worker threads the pool actually runs, 0 for inline mode
int tpool_nthreads(const struct tpool *pool);

#ifdef __cplusplus
}
#endif

#endif
