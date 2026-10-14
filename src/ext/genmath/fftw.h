/* fftw.h - MKL DFTI implementation replacing FFTW2 */
/*
 * This file provides FFTW2-compatible API using Intel MKL DFTI.
 * Original FFTW2 code was GPL licensed from MIT.
 * This implementation uses Intel MKL for FFT operations.
 *
 * Modified for Radia project - 2026
 */

#ifndef FFTW_H
#define FFTW_H

#include <stdlib.h>
#include <stdio.h>
#include <mkl_dfti.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Define for using single precision */
#define FFTW_ENABLE_FLOAT 1

/* our real numbers */
#ifdef FFTW_ENABLE_FLOAT
typedef float fftw_real;
#define MKL_FFT_PRECISION DFTI_SINGLE
#else
typedef double fftw_real;
#define MKL_FFT_PRECISION DFTI_DOUBLE
#endif

/*********************************************
 * Complex numbers and operations
 *********************************************/
typedef struct {
	fftw_real re, im;
} fftw_complex;
#define c_re(c)  ((c).re)
#define c_im(c)  ((c).im)

typedef enum {
	FFTW_FORWARD = -1, FFTW_BACKWARD = 1
} fftw_direction;

/* backward compatibility with FFTW-1.3 */
typedef fftw_complex FFTW_COMPLEX;
typedef fftw_real FFTW_REAL;

/*********************************************
 * Success or failure status
 *********************************************/
typedef enum {
	FFTW_SUCCESS = 0, FFTW_FAILURE = -1
} fftw_status;

/*****************************
 *        Plans (MKL DFTI descriptors)
 *****************************/

/* 1D FFT plan structure */
struct fftw_plan_struct {
	DFTI_DESCRIPTOR_HANDLE handle;
	int n;
	fftw_direction dir;
	int flags;
	int is_inplace;
};
typedef struct fftw_plan_struct* fftw_plan;

/* N-dimensional FFT plan structure */
struct fftwnd_plan_struct {
	DFTI_DESCRIPTOR_HANDLE handle;
	int rank;
	int* n;
	fftw_direction dir;
	int flags;
	int is_inplace;
};
typedef struct fftwnd_plan_struct* fftwnd_plan;

/* flags for the planner */
#define FFTW_ESTIMATE (0)
#define FFTW_MEASURE  (1)
#define FFTW_OUT_OF_PLACE (0)
#define FFTW_IN_PLACE (8)
#define FFTW_USE_WISDOM (16)
#define FFTW_THREADSAFE (128)
#define FFTWND_FORCE_BUFFERED (256)
#define FFTW_NO_VECTOR_RECURSE (512)

/*****************************
 *  1D FFT Functions
 *****************************/

/* Create 1D FFT plan */
static inline fftw_plan fftw_create_plan(int n, fftw_direction dir, int flags)
{
	fftw_plan plan = (fftw_plan)malloc(sizeof(struct fftw_plan_struct));
	if (!plan) return NULL;

	plan->n = n;
	plan->dir = dir;
	plan->flags = flags;
	plan->is_inplace = (flags & FFTW_IN_PLACE) ? 1 : 0;

	MKL_LONG status = DftiCreateDescriptor(&plan->handle, MKL_FFT_PRECISION, DFTI_COMPLEX, 1, (MKL_LONG)n);
	if (status != DFTI_NO_ERROR) {
		free(plan);
		return NULL;
	}

	/* Set placement */
	if (plan->is_inplace) {
		DftiSetValue(plan->handle, DFTI_PLACEMENT, DFTI_INPLACE);
	} else {
		DftiSetValue(plan->handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	}

	/* MKL uses different scaling convention
	 * FFTW: forward has no scaling, backward has 1/n scaling
	 * MKL default: no scaling for either direction
	 * We need to handle this in the transform function
	 */

	status = DftiCommitDescriptor(plan->handle);
	if (status != DFTI_NO_ERROR) {
		DftiFreeDescriptor(&plan->handle);
		free(plan);
		return NULL;
	}

	return plan;
}

/* Destroy 1D FFT plan */
static inline void fftw_destroy_plan(fftw_plan plan)
{
	if (plan) {
		if (plan->handle) {
			DftiFreeDescriptor(&plan->handle);
		}
		free(plan);
	}
}

/* Execute 1D FFT
 * howmany: number of transforms to perform
 * in: input data (complex interleaved)
 * istride: stride between elements in a single transform
 * idist: distance between first elements of consecutive transforms
 * out: output data
 * ostride: output stride
 * odist: output distance
 */
static inline void fftw(fftw_plan plan, int howmany,
                        fftw_complex* in, int istride, int idist,
                        fftw_complex* out, int ostride, int odist)
{
	if (!plan || !in) return;

	/* For batch transforms with strides, we need to process one by one
	 * or reconfigure the descriptor */
	for (int i = 0; i < howmany; i++) {
		fftw_complex* in_ptr = in + i * idist;
		fftw_complex* out_ptr = out + i * odist;

		if (plan->dir == FFTW_FORWARD) {
			if (plan->is_inplace || in_ptr == out_ptr) {
				DftiComputeForward(plan->handle, in_ptr);
			} else {
				DftiComputeForward(plan->handle, in_ptr, out_ptr);
			}
		} else {
			if (plan->is_inplace || in_ptr == out_ptr) {
				DftiComputeBackward(plan->handle, in_ptr);
			} else {
				DftiComputeBackward(plan->handle, in_ptr, out_ptr);
			}
		}
	}
}

/* Execute single 1D FFT */
static inline void fftw_one(fftw_plan plan, fftw_complex* in, fftw_complex* out)
{
	fftw(plan, 1, in, 1, 0, out, 1, 0);
}

/*****************************
 *  2D FFT Functions
 *****************************/

/* Create 2D FFT plan */
static inline fftwnd_plan fftw2d_create_plan(int nx, int ny, fftw_direction dir, int flags)
{
	fftwnd_plan plan = (fftwnd_plan)malloc(sizeof(struct fftwnd_plan_struct));
	if (!plan) return NULL;

	plan->rank = 2;
	plan->n = (int*)malloc(2 * sizeof(int));
	if (!plan->n) {
		free(plan);
		return NULL;
	}
	plan->n[0] = nx;  /* Row dimension (number of rows) */
	plan->n[1] = ny;  /* Column dimension (number of columns) */
	plan->dir = dir;
	plan->flags = flags;
	plan->is_inplace = (flags & FFTW_IN_PLACE) ? 1 : 0;

	MKL_LONG dims[2] = { (MKL_LONG)nx, (MKL_LONG)ny };
	MKL_LONG status = DftiCreateDescriptor(&plan->handle, MKL_FFT_PRECISION, DFTI_COMPLEX, 2, dims);
	if (status != DFTI_NO_ERROR) {
		free(plan->n);
		free(plan);
		return NULL;
	}

	/* Set placement */
	if (plan->is_inplace) {
		DftiSetValue(plan->handle, DFTI_PLACEMENT, DFTI_INPLACE);
	} else {
		DftiSetValue(plan->handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	}

	status = DftiCommitDescriptor(plan->handle);
	if (status != DFTI_NO_ERROR) {
		DftiFreeDescriptor(&plan->handle);
		free(plan->n);
		free(plan);
		return NULL;
	}

	return plan;
}

/* Create 3D FFT plan */
static inline fftwnd_plan fftw3d_create_plan(int nx, int ny, int nz, fftw_direction dir, int flags)
{
	fftwnd_plan plan = (fftwnd_plan)malloc(sizeof(struct fftwnd_plan_struct));
	if (!plan) return NULL;

	plan->rank = 3;
	plan->n = (int*)malloc(3 * sizeof(int));
	if (!plan->n) {
		free(plan);
		return NULL;
	}
	plan->n[0] = nx;
	plan->n[1] = ny;
	plan->n[2] = nz;
	plan->dir = dir;
	plan->flags = flags;
	plan->is_inplace = (flags & FFTW_IN_PLACE) ? 1 : 0;

	MKL_LONG dims[3] = { (MKL_LONG)nx, (MKL_LONG)ny, (MKL_LONG)nz };
	MKL_LONG status = DftiCreateDescriptor(&plan->handle, MKL_FFT_PRECISION, DFTI_COMPLEX, 3, dims);
	if (status != DFTI_NO_ERROR) {
		free(plan->n);
		free(plan);
		return NULL;
	}

	if (plan->is_inplace) {
		DftiSetValue(plan->handle, DFTI_PLACEMENT, DFTI_INPLACE);
	} else {
		DftiSetValue(plan->handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	}

	status = DftiCommitDescriptor(plan->handle);
	if (status != DFTI_NO_ERROR) {
		DftiFreeDescriptor(&plan->handle);
		free(plan->n);
		free(plan);
		return NULL;
	}

	return plan;
}

/* Create N-dimensional FFT plan */
static inline fftwnd_plan fftwnd_create_plan(int rank, const int* n, fftw_direction dir, int flags)
{
	fftwnd_plan plan = (fftwnd_plan)malloc(sizeof(struct fftwnd_plan_struct));
	if (!plan) return NULL;

	plan->rank = rank;
	plan->n = (int*)malloc(rank * sizeof(int));
	if (!plan->n) {
		free(plan);
		return NULL;
	}

	MKL_LONG* dims = (MKL_LONG*)malloc(rank * sizeof(MKL_LONG));
	if (!dims) {
		free(plan->n);
		free(plan);
		return NULL;
	}

	for (int i = 0; i < rank; i++) {
		plan->n[i] = n[i];
		dims[i] = (MKL_LONG)n[i];
	}

	plan->dir = dir;
	plan->flags = flags;
	plan->is_inplace = (flags & FFTW_IN_PLACE) ? 1 : 0;

	MKL_LONG status = DftiCreateDescriptor(&plan->handle, MKL_FFT_PRECISION, DFTI_COMPLEX, (MKL_LONG)rank, dims);
	free(dims);

	if (status != DFTI_NO_ERROR) {
		free(plan->n);
		free(plan);
		return NULL;
	}

	if (plan->is_inplace) {
		DftiSetValue(plan->handle, DFTI_PLACEMENT, DFTI_INPLACE);
	} else {
		DftiSetValue(plan->handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	}

	status = DftiCommitDescriptor(plan->handle);
	if (status != DFTI_NO_ERROR) {
		DftiFreeDescriptor(&plan->handle);
		free(plan->n);
		free(plan);
		return NULL;
	}

	return plan;
}

/* Destroy N-dimensional FFT plan */
static inline void fftwnd_destroy_plan(fftwnd_plan plan)
{
	if (plan) {
		if (plan->handle) {
			DftiFreeDescriptor(&plan->handle);
		}
		if (plan->n) {
			free(plan->n);
		}
		free(plan);
	}
}

/* Execute N-dimensional FFT */
static inline void fftwnd(fftwnd_plan plan, int howmany,
                          fftw_complex* in, int istride, int idist,
                          fftw_complex* out, int ostride, int odist)
{
	if (!plan || !in) return;

	for (int i = 0; i < howmany; i++) {
		fftw_complex* in_ptr = in + i * idist;
		fftw_complex* out_ptr = out + i * odist;

		if (plan->dir == FFTW_FORWARD) {
			if (plan->is_inplace || in_ptr == out_ptr) {
				DftiComputeForward(plan->handle, in_ptr);
			} else {
				DftiComputeForward(plan->handle, in_ptr, out_ptr);
			}
		} else {
			if (plan->is_inplace || in_ptr == out_ptr) {
				DftiComputeBackward(plan->handle, in_ptr);
			} else {
				DftiComputeBackward(plan->handle, in_ptr, out_ptr);
			}
		}
	}
}

/* Execute single N-dimensional FFT */
static inline void fftwnd_one(fftwnd_plan p, fftw_complex* in, fftw_complex* out)
{
	fftwnd(p, 1, in, 1, 0, out, 1, 0);
}

/* Printing functions (stub implementations) */
static inline void fftw_print_plan(fftw_plan plan) { (void)plan; }
static inline void fftwnd_print_plan(fftwnd_plan plan) { (void)plan; }
static inline void fftw_fprint_plan(FILE* f, fftw_plan plan) { (void)f; (void)plan; }
static inline void fftwnd_fprint_plan(FILE* f, fftwnd_plan plan) { (void)f; (void)plan; }

/* Memory functions (use standard malloc/free) */
static inline void* fftw_malloc(size_t n) { return malloc(n); }
static inline void fftw_free(void* p) { free(p); }

/* Size query */
static inline size_t fftw_sizeof_fftw_real(void) { return sizeof(fftw_real); }

/* Die function */
static inline void fftw_die(const char* s)
{
	fprintf(stderr, "fftw: %s\n", s);
	exit(1);
}

/* Wisdom functions (not supported with MKL - stubs) */
static inline void fftw_forget_wisdom(void) {}
static inline void fftw_export_wisdom(void (*emitter)(char c, void*), void* data) { (void)emitter; (void)data; }
static inline fftw_status fftw_import_wisdom(int (*g)(void*), void* data) { (void)g; (void)data; return FFTW_SUCCESS; }
static inline void fftw_export_wisdom_to_file(FILE* f) { (void)f; }
static inline fftw_status fftw_import_wisdom_from_file(FILE* f) { (void)f; return FFTW_SUCCESS; }
static inline char* fftw_export_wisdom_to_string(void) { return NULL; }
static inline fftw_status fftw_import_wisdom_from_string(const char* s) { (void)s; return FFTW_SUCCESS; }

#define FFTW_HAS_WISDOM
#define FFTW_HAS_PLAN_SPECIFIC
#define FFTW_HAS_FPRINT_PLAN
#define FFTWND_HAS_PRINT_PLAN

#ifdef __cplusplus
}
#endif

#endif /* FFTW_H */
