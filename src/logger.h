#ifndef CAMIL_LOGGER_H
#define CAMIL_LOGGER_H

#define _POSIX_C_SOURCE 200809L // this is must to prevent a bug

#ifdef __cplusplus
extern "C" {
#endif

#include <pthread.h>
#include <stdarg.h>
#include <stdio.h>
#include <time.h>

#if defined(__GNUC__)
#define CAMIL_PRINTF(fmt_index, first_arg) \
	__attribute__((format(printf, fmt_index, first_arg)))
#else
#define CAMIL_PRINTF(fmt_index, first_arg)
#endif


// enables (nonzero) or disables (zero) debug level messages
void log_set_verbose(int verbose);

// returns nonzero when debug level messages are enabled
int log_is_verbose(void);

// emits a progress message. Intended for milestones a user cares about
void log_info(const char *fmt, ...) CAMIL_PRINTF(1, 2);

// emits a recoverable problem, for example an unreadable input file that is skipped rather than aborting the whole run
void log_warn(const char *fmt, ...) CAMIL_PRINTF(1, 2);

// emits a fatal or near fatal problem
void log_error(const char *fmt, ...) CAMIL_PRINTF(1, 2);

// emits a message that is only shown in verbose mode
void log_debug(const char *fmt, ...) CAMIL_PRINTF(1, 2);

// seconds elapsed since the first call to this function (which the logger performs during the first message), using a monotonic clock when available
double log_elapsed(void);

#ifdef __cplusplus
}
#endif

#endif
