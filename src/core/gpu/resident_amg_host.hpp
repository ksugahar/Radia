// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <cstddef>
#include <vector>
struct AmgCSR {
  int n = 0, m = 0;
  std::vector<int> row, col;
  std::vector<double> val;
};
struct AmgLevel {
  AmgCSR a, p, r;
  std::vector<double> smoother;
};
struct AmgStats {
  int iterations;
  double residual;
  size_t bytes;
  unsigned solves;
  unsigned residual_refinements;
  double initial_original_residual;
  unsigned graph_batches;
  double certificate_gate;
};
void *amg_create(const std::vector<AmgLevel> &, const std::vector<double> &);
AmgStats amg_solve(void *, const double *, double *, double, int);
void amg_release(void *);
