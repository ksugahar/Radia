/* SPDX-License-Identifier: MPL-2.0 */
#pragma once
#include <cstddef>
#include <vector>
// Plain checked CSR boundary; NGSolve headers never enter the CUDA parser.
struct AMSSetupCSR { int rows,cols,nnz; const int* rowptr; const int* columns; const double* values; };
struct AMSSetupProduct { std::vector<int> rowptr,columns; std::vector<double> values; };
struct AMSSetupGPUStats { size_t products=0,upload_bytes=0,download_bytes=0,peak_owned_bytes=0; double seconds=0; };
void* AMSSetupGPUCreate();
void AMSSetupGPUDestroy(void*) noexcept;
AMSSetupProduct AMSSetupGPURestrict(void*,const void* identity,const AMSSetupCSR&,const AMSSetupCSR&,const AMSSetupCSR&);
AMSSetupGPUStats AMSSetupGPUStatistics(void*);
