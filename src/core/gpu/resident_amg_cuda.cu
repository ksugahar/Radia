// SPDX-License-Identifier: Apache-2.0
// No FE or TaskManager ownership: uploaded algebraic Galerkin hierarchy only.
#include "resident_amg_host.hpp"
#include <algorithm>
#include <cmath>
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <memory>
#include <stdexcept>
#include <sstream>
namespace {
void ck(cudaError_t e) {
  if (e != cudaSuccess)
    throw std::runtime_error(cudaGetErrorString(e));
}
void cb(cublasStatus_t e) {
  if (e != CUBLAS_STATUS_SUCCESS)
    throw std::runtime_error("cuBLAS AMG failure");
}
template <class T> struct Buffer {
  T *p = nullptr;
  size_t n;
  int device;
  explicit Buffer(size_t s) : n(s) {
    ck(cudaGetDevice(&device));
    ck(cudaMalloc(&p, std::max<size_t>(s, 1) * sizeof(T)));
  }
  ~Buffer() {
    if (p) {
      int old = -1;
      cudaGetDevice(&old);
      cudaSetDevice(device);
      cudaFree(p);
      if (old >= 0 && old != device)
        cudaSetDevice(old);
    }
  }
  void upload(const std::vector<T> &v) {
    if (v.size() != n)
      throw std::runtime_error("upload size");
    ck(cudaMemcpy(p, v.data(), n * sizeof(T), cudaMemcpyHostToDevice));
  }
  Buffer(const Buffer &) = delete;
};
thread_local cudaStream_t workStream = nullptr;
constexpr int block = 256;
__global__ void mult(int n, const int *rp, const int *ci, const double *a,
                     const double *x, double *y) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) {
    double sum = 0;
    for (int k = rp[i]; k < rp[i + 1]; ++k)
      sum += a[k] * x[ci[k]];
    y[i] = sum;
  }
}
// A warp shares long CSR rows; short transfer rows retain the scalar kernel.
// Every participating lane reaches the shuffle, including the final partial
// group of rows, so inactive rows never leave a partial warp reduction.
__global__ void multWarp(int n, const int *rp, const int *ci, const double *a,
                         const double *x, double *y) {
  const size_t thread = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
  const int row = int(thread / 32), lane = int(thread % 32);
  double sum = 0;
  if (row < n)
    for (size_t k = size_t(rp[row]) + lane; k < size_t(rp[row + 1]); k += 32)
      sum += a[k] * x[ci[k]];
  for (int offset = 16; offset; offset /= 2)
    sum += __shfl_down_sync(0xffffffffu, sum, offset);
  if (row < n && lane == 0)
    y[row] = sum;
}
__global__ void residual(int n, const double *b, const double *ax, double *r) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n)
    r[i] = b[i] - ax[i];
}
__global__ void smooth(int n, const double *b, const double *ax,
                       const double *d, double *x) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n)
    x[i] += d[i] * (b[i] - ax[i]);
}
__global__ void add(int n, const double *y, double *x) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n)
    x[i] += y[i];
}
__global__ void dense(int n, const double *a, const double *b, double *x) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) {
    double s = 0;
    for (int j = 0; j < n; ++j)
      s += a[size_t(i) * n + j] * b[j];
    x[i] = s;
  }
}
__global__ void update(int n, double alpha, const double *p, const double *ap,
                       double *x, double *r) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) {
    x[i] += alpha * p[i];
    r[i] -= alpha * ap[i];
  }
}
__global__ void direction(int n, double beta, const double *z, double *p) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n)
    p[i] = z[i] + beta * p[i];
}
struct IterationState {
  double rho, alpha, beta, bnorm, tolerance;
  int iteration, limit, active, error;
};
__global__ void graphAlpha(IterationState *s, const double *v) {
  if (!s->active)
    return;
  if (!(v[0] > 0) || !isfinite(v[0]) || !(s->rho > 0)) {
    s->error = 1;
    s->active = 0;
    return;
  }
  s->alpha = s->rho / v[0];
  if (!isfinite(s->alpha)) {
    s->error = 2;
    s->active = 0;
  }
}
__global__ void graphUpdate(int n, const IterationState *s, const double *p,
                            const double *ap, double *x, double *r) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n && s->active) {
    x[i] += s->alpha * p[i];
    r[i] -= s->alpha * ap[i];
  }
}
__global__ void graphConverge(IterationState *s, const double *v) {
  if (!s->active)
    return;
  ++s->iteration;
  if (!isfinite(v[1])) {
    s->error = 3;
    s->active = 0;
  } else if (v[1] / s->bnorm <= s->tolerance)
    s->active = 0;
  else if (s->iteration >= s->limit) {
    s->error = 4;
    s->active = 0;
  }
}
__global__ void graphBeta(IterationState *s, const double *v) {
  if (!s->active)
    return;
  if (!(v[2] > 0) || !isfinite(v[2])) {
    s->error = 1;
    s->active = 0;
    return;
  }
  s->beta = v[2] / s->rho;
  s->rho = v[2];
  if (!isfinite(s->beta)) {
    s->error = 2;
    s->active = 0;
  }
}
__global__ void graphDirection(int n, const IterationState *s, const double *z,
                               double *p) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n && s->active)
    p[i] = z[i] + s->beta * p[i];
}
struct CSR {
  int n, m;
  Buffer<int> row, col;
  Buffer<double> val;
  explicit CSR(const AmgCSR &a)
      : n(a.n), m(a.m), row(a.row.size()), col(a.col.size()),
        val(a.val.size()) {
    row.upload(a.row);
    col.upload(a.col);
    val.upload(a.val);
  }
  void apply(const double *x, double *y) {
    if (n && col.n / size_t(n) >= 16)
      multWarp<<<(n + block / 32 - 1) / (block / 32), block, 0, workStream>>>(
          n, row.p, col.p, val.p, x, y);
    else
      mult<<<(n + block - 1) / block, block, 0, workStream>>>(
          n, row.p, col.p, val.p, x, y);
  }
  size_t bytes() const { return row.n * 4 + col.n * 4 + val.n * 8; }
};
struct Level {
  CSR a;
  std::unique_ptr<CSR> p, r;
  Buffer<double> d, b, x, tmp;
  explicit Level(const AmgLevel &s)
      : a(s.a), d(s.smoother.size()), b(s.a.n), x(s.a.n), tmp(s.a.n) {
    d.upload(s.smoother);
    if (s.p.n) {
      p = std::make_unique<CSR>(s.p);
      r = std::make_unique<CSR>(s.r);
    }
  }
  size_t bytes() const {
    return a.bytes() + (p ? p->bytes() + r->bytes() : 0) +
           (d.n + b.n + x.n + tmp.n) * 8;
  }
};
struct Context {
  std::vector<std::unique_ptr<Level>> levels;
  Buffer<double> inverse, b, x, r, p, z, ap, original_b, solution;
  Buffer<double> scalars{3};
  Buffer<IterationState> state{1};
  cudaGraph_t batchGraph = nullptr;
  cudaGraphExec_t batchExec = nullptr;
  int n, device;
  unsigned solves = 0;
  cublasHandle_t blas = nullptr;
  cudaStream_t stream = nullptr;
  cudaGraph_t graph = nullptr;
  cudaGraphExec_t executable = nullptr;
  Context(const std::vector<AmgLevel> &s, const std::vector<double> &inv)
      : inverse(inv.size()), b(s.front().a.n), x(s.front().a.n),
        r(s.front().a.n), p(s.front().a.n), z(s.front().a.n), ap(s.front().a.n),
        original_b(s.front().a.n), solution(s.front().a.n), n(s.front().a.n) {
    ck(cudaGetDevice(&device));
    inverse.upload(inv);
    for (auto &a : s)
      levels.emplace_back(std::make_unique<Level>(a));
    ck(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    workStream = stream;
    try {
      cb(cublasCreate(&blas));
      cb(cublasSetStream(blas, stream));
    } catch (...) {
      if (blas) cublasDestroy(blas);
      cudaStreamDestroy(stream);
      workStream = nullptr;
      throw;
    }
  }
  ~Context() {
    int old = -1;
    cudaGetDevice(&old);
    cudaSetDevice(device);
    if (batchExec)
      cudaGraphExecDestroy(batchExec);
    if (batchGraph)
      cudaGraphDestroy(batchGraph);
    if (executable)
      cudaGraphExecDestroy(executable);
    if (graph)
      cudaGraphDestroy(graph);
    if (blas)
      cublasDestroy(blas);
    if (stream)
      cudaStreamDestroy(stream);
    levels.clear();
    if (old >= 0 && old != device)
      cudaSetDevice(old);
  }
  void selected() {
    workStream = stream;
    int d;
    ck(cudaGetDevice(&d));
    if (d != device)
      throw std::runtime_error("AMG context CUDA device mismatch");
  }
  void cycle(size_t i) {
    auto &l = *levels[i];
    int k = (l.a.n + block - 1) / block;
    if (i + 1 == levels.size()) {
      dense<<<k, block, 0, workStream>>>(l.a.n, inverse.p, l.b.p, l.x.p);
      return;
    }
    ck(cudaMemsetAsync(l.x.p, 0, l.x.n * 8, stream));
    for (int j = 0; j < 2; ++j) {
      l.a.apply(l.x.p, l.tmp.p);
      smooth<<<k, block, 0, workStream>>>(l.a.n, l.b.p, l.tmp.p, l.d.p, l.x.p);
    }
    l.a.apply(l.x.p, l.tmp.p);
    residual<<<k, block, 0, workStream>>>(l.a.n, l.b.p, l.tmp.p, l.tmp.p);
    l.r->apply(l.tmp.p, levels[i + 1]->b.p);
    cycle(i + 1);
    l.p->apply(levels[i + 1]->x.p, l.tmp.p);
    add<<<k, block, 0, workStream>>>(l.a.n, l.tmp.p, l.x.p);
    for (int j = 0; j < 2; ++j) {
      l.a.apply(l.x.p, l.tmp.p);
      smooth<<<k, block, 0, workStream>>>(l.a.n, l.b.p, l.tmp.p, l.d.p, l.x.p);
    }
  }
  void precondition() { ck(cudaGraphLaunch(executable, stream)); }
  void capture() {
    ck(cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal));
    try {
      ck(cudaMemcpyAsync(levels[0]->b.p, r.p, n * 8, cudaMemcpyDeviceToDevice,
                         stream));
      cycle(0);
      ck(cudaMemcpyAsync(z.p, levels[0]->x.p, n * 8, cudaMemcpyDeviceToDevice,
                         stream));
      ck(cudaStreamEndCapture(stream, &graph));
      ck(cudaGraphInstantiate(&executable, graph, nullptr, nullptr, 0));
    } catch (...) {
      cudaGraph_t aborted = nullptr;
      cudaStreamEndCapture(stream, &aborted);
      if (aborted)
        cudaGraphDestroy(aborted);
      throw;
    }
  }
  void inlinePrecondition() {
    ck(cudaMemcpyAsync(levels[0]->b.p, r.p, n * 8, cudaMemcpyDeviceToDevice,
                       stream));
    cycle(0);
    ck(cudaMemcpyAsync(z.p, levels[0]->x.p, n * 8, cudaMemcpyDeviceToDevice,
                       stream));
  }
  void captureBatch() {
    cb(cublasSetPointerMode(blas, CUBLAS_POINTER_MODE_DEVICE));
    ck(cudaMemsetAsync(r.p, 0, n * 8, stream));
    ck(cudaMemsetAsync(p.p, 0, n * 8, stream));
    cb(cublasDdot(blas, n, r.p, 1, r.p, 1, scalars.p));
    cb(cublasDnrm2(blas, n, r.p, 1, scalars.p + 1));
    ck(cudaStreamSynchronize(stream));
    ck(cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal));
    try {
      for (int i = 0; i < 8; ++i) {
        levels[0]->a.apply(p.p, ap.p);
        cb(cublasDdot(blas, n, p.p, 1, ap.p, 1, scalars.p));
        graphAlpha<<<1, 1, 0, stream>>>(state.p, scalars.p);
        graphUpdate<<<(n + block - 1) / block, block, 0, stream>>>(
            n, state.p, p.p, ap.p, x.p, r.p);
        cb(cublasDnrm2(blas, n, r.p, 1, scalars.p + 1));
        graphConverge<<<1, 1, 0, stream>>>(state.p, scalars.p);
        inlinePrecondition();
        cb(cublasDdot(blas, n, r.p, 1, z.p, 1, scalars.p + 2));
        graphBeta<<<1, 1, 0, stream>>>(state.p, scalars.p);
        graphDirection<<<(n + block - 1) / block, block, 0, stream>>>(
            n, state.p, z.p, p.p);
      }
      ck(cudaStreamEndCapture(stream, &batchGraph));
      ck(cudaGraphInstantiate(&batchExec, batchGraph, nullptr, nullptr, 0));
    } catch (...) {
      cudaGraph_t aborted = nullptr;
      cudaStreamEndCapture(stream, &aborted);
      if (aborted)
        cudaGraphDestroy(aborted);
      throw;
    }
    cb(cublasSetPointerMode(blas, CUBLAS_POINTER_MODE_HOST));
  }
  double dot(const double *a, const double *b) {
    double v;
    cb(cublasDdot(blas, n, a, 1, b, 1, &v));
    if (!std::isfinite(v))
      throw std::runtime_error("AMG nonfinite dot");
    return v;
  }
  double norm(const double *a) {
    double v;
    cb(cublasDnrm2(blas, n, a, 1, &v));
    if (!std::isfinite(v))
      throw std::runtime_error("AMG nonfinite norm");
    return v;
  }
  size_t bytes() const {
    size_t s =
        (inverse.n + b.n + x.n + r.n + p.n + z.n + ap.n + original_b.n +
         solution.n + scalars.n) * 8 +
        sizeof(IterationState);
    for (auto &l : levels)
      s += l->bytes();
    return s;
  }
};
} // namespace
void *amg_create(const std::vector<AmgLevel> &s,
                 const std::vector<double> &inv) {
  int device;
  ck(cudaGetDevice(&device));
  cudaDeviceProp properties{};
  ck(cudaGetDeviceProperties(&properties, device));
  if (properties.major < 7 || (properties.major == 7 && properties.minor < 5))
    throw std::runtime_error("Resident AMG requires compute capability 7.5 or newer");
  auto c = std::make_unique<Context>(s, inv);
  c->capture();
  c->captureBatch();
  ck(cudaDeviceSynchronize());
  return c.release();
}
AmgStats amg_solve(void *handle, const double *rhs, double *out, double tol,
                   int limit) {
  auto &c = *static_cast<Context *>(handle);
  c.selected();
  ck(cudaMemcpyAsync(c.b.p, rhs, c.n * 8, cudaMemcpyHostToDevice, c.stream));
  ck(cudaMemcpyAsync(c.original_b.p, c.b.p, c.n * 8, cudaMemcpyDeviceToDevice,
                     c.stream));
  double bn = c.norm(c.b.p);
  int iter = 0;
  unsigned batches = 0, refinements = 0;
  if (bn == 0 && std::any_of(rhs, rhs + c.n, [](double v) { return v != 0; }))
    throw std::runtime_error("RHS underflow");
  auto pcg = [&]() {
    ck(cudaMemcpyAsync(c.r.p, c.b.p, c.n * 8, cudaMemcpyDeviceToDevice, c.stream));
    ck(cudaMemsetAsync(c.x.p, 0, c.n * 8, c.stream));
    double current_bn = c.norm(c.b.p);
    if (current_bn == 0)
      return;
    int remaining = limit - iter;
    if (remaining < 1)
      throw std::runtime_error("AMG refinement total iteration limit");
    c.precondition();
    ck(cudaMemcpyAsync(c.p.p, c.z.p, c.n * 8, cudaMemcpyDeviceToDevice,
                       c.stream));
    double rho = c.dot(c.r.p, c.z.p);
    if (!(rho > 0))
      throw std::runtime_error("AMG nonpositive preconditioner");
    IterationState status{rho, 0, 0, current_bn, tol, 0, remaining, 1, 0};
    ck(cudaMemcpyAsync(c.state.p, &status, sizeof(status),
                       cudaMemcpyHostToDevice, c.stream));
    for (int i = 0; i < remaining; i += 8) {
      ck(cudaGraphLaunch(c.batchExec, c.stream));
      ++batches;
      ck(cudaMemcpyAsync(&status, c.state.p, sizeof(status),
                         cudaMemcpyDeviceToHost, c.stream));
      ck(cudaStreamSynchronize(c.stream));
      if (!status.active)
        break;
    }
    if (status.error)
      throw std::runtime_error("AMG graph PCG breakdown or iteration limit");
    iter += status.iteration;
  };
  pcg();
  ck(cudaMemcpyAsync(c.solution.p, c.x.p, c.n * 8, cudaMemcpyDeviceToDevice, c.stream));
  auto original_residual = [&]() {
    c.levels[0]->a.apply(c.solution.p, c.ap.p);
    residual<<<(c.n + block - 1) / block, block, 0, workStream>>>(
        c.n, c.original_b.p, c.ap.p, c.r.p);
    ck(cudaGetLastError());
    return bn > 0 ? c.norm(c.r.p) / bn : c.norm(c.r.p);
  };
  double rr = original_residual(), initial_rr = rr;
  const double gate = std::min(std::max(10 * tol, 1e-12), 1e-6);
  // A recursively updated Krylov residual can drift from b-A*x. Restart
  // against that original residual on the SAME resident operator, using the
  // SAME requested PCG tolerance. The separate original-system certificate
  // is min(max(10*tol, 1e-12), 1e-6), recorded in the returned statistics.
  while (std::isfinite(rr) && rr > gate && refinements < 3) {
    ck(cudaMemcpyAsync(c.b.p, c.r.p, c.n * 8, cudaMemcpyDeviceToDevice, c.stream));
    pcg();
    add<<<(c.n + block - 1) / block, block, 0, workStream>>>(c.n, c.x.p, c.solution.p);
    ck(cudaGetLastError());
    ++refinements;
    rr = original_residual();
  }
  if (iter > limit || !std::isfinite(rr) || rr > gate) {
    std::ostringstream message;
    message.precision(17);
    message << "AMG original residual gate: residual=" << rr << ", gate=" << gate
            << ", initial=" << initial_rr << ", rhs_norm=" << bn
            << ", refinements=" << refinements << ", iterations=" << iter;
    throw std::runtime_error(message.str());
  }
  ck(cudaMemcpyAsync(out, c.solution.p, c.n * 8, cudaMemcpyDeviceToHost, c.stream));
  ck(cudaStreamSynchronize(c.stream));
  ++c.solves;
  return {iter, rr, c.bytes(), c.solves, refinements, initial_rr, batches, gate};
}
void amg_release(void *p) { delete static_cast<Context *>(p); }
