#include "logger.h"


// guards the output stream so that concurrently logging worker threads do not interleave partial lines
static pthread_mutex_t log_mutex = PTHREAD_MUTEX_INITIALIZER;

// set to nonzero by log_set_verbose(). Without it only warnings and errors
// are printed; progress (INFO) and detail (DEBUG) need --verbose.
static int log_verbose = 0;

// wall clock reading taken the first time log_elapsed() runs; all reported timestamps are relative to it
static double log_start = -1.0;

// returns the current time in seconds
static double log_now(void) {
	struct timespec ts;
#if defined(_POSIX_C_SOURCE) && _POSIX_C_SOURCE >= 199309L
	if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0) {
		return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
	}
#endif
	if (timespec_get(&ts, TIME_UTC) == TIME_UTC) {
		return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
	}
	return 0.0;
}

double log_elapsed(void) {
	double now = log_now();
	if (log_start < 0.0) {
		log_start = now;
	}
	return now - log_start;
}

void log_set_verbose(int verbose) {
	log_verbose = verbose;
}

int log_is_verbose(void) {
	return log_verbose;
}

// shared body of the four public entry points
static void log_emit(const char *level, const char *fmt, va_list ap) {
	double elapsed = log_elapsed();

	pthread_mutex_lock(&log_mutex);
	fprintf(stderr, "[%8.3f] [%s] ", elapsed, level);
	vfprintf(stderr, fmt, ap);
	fputc('\n', stderr);
	fflush(stderr);
	pthread_mutex_unlock(&log_mutex);
}

void log_info(const char *fmt, ...) {
	va_list ap;
	if (!log_verbose) {
		return;
	}
	va_start(ap, fmt);
	log_emit("INFO", fmt, ap);
	va_end(ap);
}

void log_warn(const char *fmt, ...) {
	va_list ap;
	va_start(ap, fmt);
	log_emit("WARN", fmt, ap);
	va_end(ap);
}

void log_error(const char *fmt, ...) {
	va_list ap;
	va_start(ap, fmt);
	log_emit("ERROR", fmt, ap);
	va_end(ap);
}

void log_debug(const char *fmt, ...) {
	va_list ap;
	if (!log_verbose) {
		return;
	}
	va_start(ap, fmt);
	log_emit("DEBUG", fmt, ap);
	va_end(ap);
}
