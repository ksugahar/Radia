/*
 * HACApK Logging Macros for Radia Integration
 *
 * This header provides logging macros that can be configured to
 * output to console (for standalone use) or integrate with Radia's
 * logging system.
 *
 * Copyright (c) 2025 Radia Project
 * License: MIT
 */

#ifndef HACAPK_LOG_H_INCLUDED
#define HACAPK_LOG_H_INCLUDED

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Log levels
 */
#define HACAPK_LOG_SILENT  0
#define HACAPK_LOG_ERROR   1
#define HACAPK_LOG_WARN    2
#define HACAPK_LOG_INFO    3
#define HACAPK_LOG_DEBUG   4

/*
 * Global log level (can be set at runtime)
 * Default: HACAPK_LOG_INFO for normal operation
 */
#ifndef HACAPK_LOG_LEVEL
#define HACAPK_LOG_LEVEL HACAPK_LOG_INFO
#endif

/*
 * Log output macros
 *
 * In Radia integration mode, these could be redirected to radTSend
 * For now, they use simple console output
 */

#ifdef HACAPK_RADIA_LOGGING
/* Radia-style logging (silent by default, errors to stderr) */
#define HACAPK_LOG_ERROR_MSG(fmt, ...) \
    fprintf(stderr, "[HACApK Error] " fmt "\n", ##__VA_ARGS__)

#define HACAPK_LOG_WARN_MSG(fmt, ...) \
    do { /* Warnings suppressed in Radia mode */ } while(0)

#define HACAPK_LOG_INFO_MSG(fmt, ...) \
    do { /* Info suppressed in Radia mode */ } while(0)

#define HACAPK_LOG_DEBUG_MSG(fmt, ...) \
    do { /* Debug suppressed in Radia mode */ } while(0)

#else
/* Standalone logging with level control */
#define HACAPK_LOG_ERROR_MSG(fmt, ...) \
    fprintf(stderr, "[HACApK Error] " fmt "\n", ##__VA_ARGS__)

#define HACAPK_LOG_WARN_MSG(fmt, ...) \
    do { if (HACAPK_LOG_LEVEL >= HACAPK_LOG_WARN) \
        fprintf(stderr, "[HACApK Warn] " fmt "\n", ##__VA_ARGS__); } while(0)

#define HACAPK_LOG_INFO_MSG(fmt, ...) \
    do { if (HACAPK_LOG_LEVEL >= HACAPK_LOG_INFO) \
        printf("[HACApK] " fmt "\n", ##__VA_ARGS__); } while(0)

#define HACAPK_LOG_DEBUG_MSG(fmt, ...) \
    do { if (HACAPK_LOG_LEVEL >= HACAPK_LOG_DEBUG) \
        printf("[HACApK Debug] " fmt "\n", ##__VA_ARGS__); } while(0)

#endif /* HACAPK_RADIA_LOGGING */

/*
 * Conditional logging based on print_level parameter
 * (matches original HACApK behavior)
 */
#define HACAPK_LOG_IF(print_level, level, fmt, ...) \
    do { if ((print_level) >= (level)) HACAPK_LOG_INFO_MSG(fmt, ##__VA_ARGS__); } while(0)

#define HACAPK_LOG_IF_DEBUG(print_level, fmt, ...) \
    do { if ((print_level) >= 2) HACAPK_LOG_DEBUG_MSG(fmt, ##__VA_ARGS__); } while(0)

#ifdef __cplusplus
}
#endif

#endif /* HACAPK_LOG_H_INCLUDED */
