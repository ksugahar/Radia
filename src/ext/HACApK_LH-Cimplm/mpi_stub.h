/*
 * MPI Stub Header for Single-Process HACApK
 *
 * This header provides minimal MPI function stubs to allow HACApK
 * to compile and run without actual MPI support.
 *
 * For Radia integration, we run on a single process only.
 *
 * Copyright (c) 2025 Radia Project
 * License: MIT
 */

#ifndef MPI_STUB_H_INCLUDED
#define MPI_STUB_H_INCLUDED

#include <stddef.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

/* MPI Types */
typedef int MPI_Comm;
typedef int MPI_Datatype;
typedef int MPI_Op;
typedef int MPI_Status;
typedef int MPI_Request;

/* MPI Constants */
#define MPI_COMM_WORLD 0
#define MPI_SUCCESS 0
#define MPI_DOUBLE 0
#define MPI_INT 1
#define MPI_CHAR 2
#define MPI_BYTE 3
#define MPI_SUM 0
#define MPI_MAX 1
#define MPI_MIN 2
#define MPI_UNDEFINED (-1)
#define MPI_STATUS_IGNORE ((MPI_Status*)0)

/* Convert Fortran communicator to C - just return the same value */
static inline MPI_Comm MPI_Comm_f2c(int fortran_comm) {
    (void)fortran_comm;
    return MPI_COMM_WORLD;
}

/* Get rank - always return 0 for single process */
static inline int MPI_Comm_rank(MPI_Comm comm, int *rank) {
    (void)comm;
    *rank = 0;
    return MPI_SUCCESS;
}

/* Get size - always return 1 for single process */
static inline int MPI_Comm_size(MPI_Comm comm, int *size) {
    (void)comm;
    *size = 1;
    return MPI_SUCCESS;
}

/* Initialize MPI - no-op */
static inline int MPI_Init(int *argc, char ***argv) {
    (void)argc;
    (void)argv;
    return MPI_SUCCESS;
}

/* Finalize MPI - no-op */
static inline int MPI_Finalize(void) {
    return MPI_SUCCESS;
}

/* Barrier - no-op for single process */
static inline int MPI_Barrier(MPI_Comm comm) {
    (void)comm;
    return MPI_SUCCESS;
}

/* Allreduce - just copy for single process */
static inline int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                                MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
    (void)datatype;
    (void)op;
    (void)comm;

    size_t elem_size;
    switch (datatype) {
        case 0: elem_size = sizeof(double); break;  /* MPI_DOUBLE */
        case 1: elem_size = sizeof(int); break;     /* MPI_INT */
        default: elem_size = 1; break;
    }

    if (sendbuf != recvbuf) {
        memcpy(recvbuf, sendbuf, count * elem_size);
    }
    return MPI_SUCCESS;
}

/* Reduce - just copy to root for single process */
static inline int MPI_Reduce(const void *sendbuf, void *recvbuf, int count,
                             MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
    (void)root;
    return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
}

/* Bcast - no-op for single process */
static inline int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype,
                            int root, MPI_Comm comm) {
    (void)buffer;
    (void)count;
    (void)datatype;
    (void)root;
    (void)comm;
    return MPI_SUCCESS;
}

/* Scatter - just copy for single process */
static inline int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                              void *recvbuf, int recvcount, MPI_Datatype recvtype,
                              int root, MPI_Comm comm) {
    (void)sendcount;
    (void)sendtype;
    (void)recvtype;
    (void)root;
    (void)comm;

    size_t elem_size;
    switch (recvtype) {
        case 0: elem_size = sizeof(double); break;
        case 1: elem_size = sizeof(int); break;
        default: elem_size = 1; break;
    }

    if (sendbuf != recvbuf) {
        memcpy(recvbuf, sendbuf, recvcount * elem_size);
    }
    return MPI_SUCCESS;
}

/* Gather - just copy for single process */
static inline int MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                             void *recvbuf, int recvcount, MPI_Datatype recvtype,
                             int root, MPI_Comm comm) {
    (void)sendcount;
    (void)sendtype;
    (void)recvtype;
    (void)root;
    (void)comm;

    size_t elem_size;
    switch (sendtype) {
        case 0: elem_size = sizeof(double); break;
        case 1: elem_size = sizeof(int); break;
        default: elem_size = 1; break;
    }

    if (sendbuf != recvbuf) {
        memcpy(recvbuf, sendbuf, recvcount * elem_size);
    }
    return MPI_SUCCESS;
}

/* Allgather - just copy for single process */
static inline int MPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                void *recvbuf, int recvcount, MPI_Datatype recvtype,
                                MPI_Comm comm) {
    return MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, 0, comm);
}

/* Allgatherv - just copy for single process */
static inline int MPI_Allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                 void *recvbuf, const int *recvcounts, const int *displs,
                                 MPI_Datatype recvtype, MPI_Comm comm) {
    (void)recvcounts;
    (void)displs;
    return MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, sendcount, recvtype, 0, comm);
}

/* Comm_split - return same communicator */
static inline int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm) {
    (void)comm;
    (void)color;
    (void)key;
    *newcomm = MPI_COMM_WORLD;
    return MPI_SUCCESS;
}

/* Comm_free - no-op */
static inline int MPI_Comm_free(MPI_Comm *comm) {
    (void)comm;
    return MPI_SUCCESS;
}

/* Wtime - return wall clock time */
static inline double MPI_Wtime(void) {
#ifdef _WIN32
    #include <windows.h>
    LARGE_INTEGER freq, count;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&count);
    return (double)count.QuadPart / (double)freq.QuadPart;
#else
    #include <sys/time.h>
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
#endif
}

/* Abort - exit with error */
static inline int MPI_Abort(MPI_Comm comm, int errorcode) {
    (void)comm;
    exit(errorcode);
    return MPI_SUCCESS;  /* Never reached */
}

/* Irecv, Isend - no-op for single process */
static inline int MPI_Irecv(void *buf, int count, MPI_Datatype datatype,
                            int source, int tag, MPI_Comm comm, MPI_Request *request) {
    (void)buf; (void)count; (void)datatype; (void)source; (void)tag; (void)comm;
    *request = 0;
    return MPI_SUCCESS;
}

static inline int MPI_Isend(const void *buf, int count, MPI_Datatype datatype,
                            int dest, int tag, MPI_Comm comm, MPI_Request *request) {
    (void)buf; (void)count; (void)datatype; (void)dest; (void)tag; (void)comm;
    *request = 0;
    return MPI_SUCCESS;
}

static inline int MPI_Waitall(int count, MPI_Request *array_of_requests,
                              MPI_Status *array_of_statuses) {
    (void)count; (void)array_of_requests; (void)array_of_statuses;
    return MPI_SUCCESS;
}

static inline int MPI_Wait(MPI_Request *request, MPI_Status *status) {
    (void)request; (void)status;
    return MPI_SUCCESS;
}

#ifdef __cplusplus
}
#endif

#endif /* MPI_STUB_H_INCLUDED */
