#include "tpool.h"


// body of every worker thread: pop a job, run it outside the lock, repeat until the pool is shut down and the queue is empty
static void *tpool_worker(void *arg) {
	struct tpool *pool = (struct tpool *)arg;

	for (;;) {
		struct tpool_job job;

		pthread_mutex_lock(&pool->mutex);
		while (pool->count == 0 && !pool->shutdown) {
			pthread_cond_wait(&pool->not_empty, &pool->mutex);
		}
		if (pool->count == 0 && pool->shutdown) {
			pthread_mutex_unlock(&pool->mutex);
			break;
		}

		job = pool->queue[pool->head];
		pool->head = (pool->head + 1) % pool->capacity;
		pool->count--;
		pool->running++;
		pthread_cond_signal(&pool->not_full);
		pthread_mutex_unlock(&pool->mutex);

		job.fn(job.arg);

		pthread_mutex_lock(&pool->mutex);
		pool->running--;
		if (pool->count == 0 && pool->running == 0) {
			pthread_cond_broadcast(&pool->idle);
		}
		pthread_mutex_unlock(&pool->mutex);
	}

	return NULL;
}

int tpool_init(struct tpool *pool, int nthreads, size_t capacity) {
	int i;

	memset(pool, 0, sizeof(*pool));

	if (capacity < 1) {
		capacity = 1;
	}

	pool->capacity = capacity;
	pool->queue = (struct tpool_job *)malloc(capacity * sizeof(struct tpool_job));
	if (pool->queue == NULL) {
		log_error("thread pool: out of memory while allocating the job queue");
		return -1;
	}

	pthread_mutex_init(&pool->mutex, NULL);
	pthread_cond_init(&pool->not_empty, NULL);
	pthread_cond_init(&pool->not_full, NULL);
	pthread_cond_init(&pool->idle, NULL);

	// Inline mode: no workers, jobs run on the submitting thread.
	if (nthreads <= 1) {
		pool->nthreads = 0;
		return 0;
	}

	pool->threads = (pthread_t *)malloc((size_t)nthreads * sizeof(pthread_t));
	if (pool->threads == NULL) {
		log_error("thread pool: out of memory while allocating thread handles");
		free(pool->queue);
		pthread_mutex_destroy(&pool->mutex);
		pthread_cond_destroy(&pool->not_empty);
		pthread_cond_destroy(&pool->not_full);
		pthread_cond_destroy(&pool->idle);
		memset(pool, 0, sizeof(*pool));
		return -1;
	}

	for (i = 0; i < nthreads; i++) {
		if (pthread_create(&pool->threads[i], NULL, tpool_worker, pool) != 0) {
			// Keep the workers that did start; a smaller pool is still usable.
			log_warn("thread pool: could only start %d of %d threads", i, nthreads);
			break;
		}
	}
	pool->nthreads = i;

	// not a single worker could be started, fall back to inline execution
	if (pool->nthreads == 0) {
		free(pool->threads);
		pool->threads = NULL;
	}

	return 0;
}

void tpool_submit(struct tpool *pool, void (*fn)(void *), void *arg) {
	size_t tail;

	if (pool->nthreads == 0) {
		fn(arg);
		return;
	}

	pthread_mutex_lock(&pool->mutex);
	while (pool->count == pool->capacity) {
		pthread_cond_wait(&pool->not_full, &pool->mutex);
	}

	tail = (pool->head + pool->count) % pool->capacity;
	pool->queue[tail].fn = fn;
	pool->queue[tail].arg = arg;
	pool->count++;

	pthread_cond_signal(&pool->not_empty);
	pthread_mutex_unlock(&pool->mutex);
}

void tpool_wait(struct tpool *pool) {
	if (pool->nthreads == 0) {
		return;
	}

	pthread_mutex_lock(&pool->mutex);
	while (pool->count > 0 || pool->running > 0) {
		pthread_cond_wait(&pool->idle, &pool->mutex);
	}
	pthread_mutex_unlock(&pool->mutex);
}

void tpool_destroy(struct tpool *pool) {
	int i;

	if (pool->queue == NULL) {
		return;
	}

	tpool_wait(pool);

	pthread_mutex_lock(&pool->mutex);
	pool->shutdown = 1;
	pthread_cond_broadcast(&pool->not_empty);
	pthread_mutex_unlock(&pool->mutex);

	for (i = 0; i < pool->nthreads; i++) {
		pthread_join(pool->threads[i], NULL);
	}

	free(pool->threads);
	free(pool->queue);
	pthread_mutex_destroy(&pool->mutex);
	pthread_cond_destroy(&pool->not_empty);
	pthread_cond_destroy(&pool->not_full);
	pthread_cond_destroy(&pool->idle);
	memset(pool, 0, sizeof(*pool));
}

int tpool_nthreads(const struct tpool *pool) {
	return pool->nthreads;
}
