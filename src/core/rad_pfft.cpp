/*
 * rad_pfft.cpp
 *
 * Precorrected FFT (pFFT) acceleration for integral equation solver
 * Implementation using Intel MKL DFTI
 *
 * Part of Radia project
 */

#include "rad_pfft.h"
#include "rad_constants.h"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace radia {

radTPfft::radTPfft()
    : nx_(0), ny_(0), nz_(0)
    , nx2_(0), ny2_(0), nz2_(0)
    , totalGridSize_(0)
    , gridSpacing_(0)
    , interpOrder_(2)
    , numPanels_(0)
    , fftHandle_(nullptr)
    , fftInitialized_(false)
    , initialized_(false)
    , kernelSetup_(false)
{
    gridOrigin_ = TVector3d(0, 0, 0);
}

radTPfft::~radTPfft() {
    Cleanup();
}

void radTPfft::Cleanup() {
    DestroyMKLFFT();

    gridKernel_.clear();
    gridWork_.clear();
    projMatrix_.clear();
    interpMatrix_.clear();
    directPairs_.clear();
    directValues_.clear();
    gridValues_.clear();

    initialized_ = false;
    kernelSetup_ = false;
    numPanels_ = 0;
}

void radTPfft::InitializeMKLFFT() {
    if (fftInitialized_) {
        DestroyMKLFFT();
    }

    MKL_LONG dims[3] = { static_cast<MKL_LONG>(nz2_),
                         static_cast<MKL_LONG>(ny2_),
                         static_cast<MKL_LONG>(nx2_) };

    MKL_LONG status = DftiCreateDescriptor(&fftHandle_, DFTI_DOUBLE,
                                           DFTI_COMPLEX, 3, dims);
    if (status != DFTI_NO_ERROR) {
        throw std::runtime_error("DftiCreateDescriptor failed");
    }

    // In-place transform
    status = DftiSetValue(fftHandle_, DFTI_PLACEMENT, DFTI_INPLACE);
    if (status != DFTI_NO_ERROR) {
        DftiFreeDescriptor(&fftHandle_);
        throw std::runtime_error("DftiSetValue DFTI_PLACEMENT failed");
    }

    // Backward scaling: 1/N normalization
    double scale = 1.0 / static_cast<double>(totalGridSize_);
    status = DftiSetValue(fftHandle_, DFTI_BACKWARD_SCALE, scale);
    if (status != DFTI_NO_ERROR) {
        DftiFreeDescriptor(&fftHandle_);
        throw std::runtime_error("DftiSetValue DFTI_BACKWARD_SCALE failed");
    }

    status = DftiCommit(fftHandle_);
    if (status != DFTI_NO_ERROR) {
        DftiFreeDescriptor(&fftHandle_);
        throw std::runtime_error("DftiCommit failed");
    }

    fftInitialized_ = true;
}

void radTPfft::DestroyMKLFFT() {
    if (fftInitialized_ && fftHandle_ != nullptr) {
        DftiFreeDescriptor(&fftHandle_);
        fftHandle_ = nullptr;
        fftInitialized_ = false;
    }
}

void radTPfft::FFT3D_Forward(std::vector<Complex>& data) {
    if (!fftInitialized_) {
        throw std::runtime_error("FFT not initialized");
    }

    MKL_LONG status = DftiComputeForward(fftHandle_, data.data());
    if (status != DFTI_NO_ERROR) {
        throw std::runtime_error("DftiComputeForward failed");
    }
}

void radTPfft::FFT3D_Backward(std::vector<Complex>& data) {
    if (!fftInitialized_) {
        throw std::runtime_error("FFT not initialized");
    }

    MKL_LONG status = DftiComputeBackward(fftHandle_, data.data());
    if (status != DFTI_NO_ERROR) {
        throw std::runtime_error("DftiComputeBackward failed");
    }
}

void radTPfft::Initialize(const std::vector<TVector3d>& panelCenters,
                          const std::vector<double>& panelAreas,
                          double gridSpacing,
                          int interpOrder) {
    Cleanup();

    if (panelCenters.empty()) {
        throw std::invalid_argument("panelCenters cannot be empty");
    }

    numPanels_ = static_cast<int>(panelCenters.size());
    interpOrder_ = interpOrder;

    // Compute bounding box
    TVector3d minPt = panelCenters[0];
    TVector3d maxPt = panelCenters[0];

    for (const auto& p : panelCenters) {
        minPt.x = std::min(minPt.x, p.x);
        minPt.y = std::min(minPt.y, p.y);
        minPt.z = std::min(minPt.z, p.z);
        maxPt.x = std::max(maxPt.x, p.x);
        maxPt.y = std::max(maxPt.y, p.y);
        maxPt.z = std::max(maxPt.z, p.z);
    }

    // Auto-compute grid spacing if not specified
    if (gridSpacing <= 0) {
        // Use average panel spacing
        double avgArea = 0;
        if (!panelAreas.empty()) {
            avgArea = std::accumulate(panelAreas.begin(), panelAreas.end(), 0.0)
                      / panelAreas.size();
            gridSpacing = std::sqrt(avgArea);
        } else {
            // Estimate from bounding box
            TVector3d extent;
            extent.x = maxPt.x - minPt.x;
            extent.y = maxPt.y - minPt.y;
            extent.z = maxPt.z - minPt.z;
            double maxExtent = std::max({extent.x, extent.y, extent.z});
            gridSpacing = maxExtent / std::cbrt(static_cast<double>(numPanels_));
        }
    }
    gridSpacing_ = gridSpacing;

    // Add margin for interpolation stencil
    double margin = (interpOrder + 1) * gridSpacing;
    minPt.x -= margin;
    minPt.y -= margin;
    minPt.z -= margin;
    maxPt.x += margin;
    maxPt.y += margin;
    maxPt.z += margin;

    gridOrigin_ = minPt;

    // Compute grid dimensions
    nx_ = static_cast<int>(std::ceil((maxPt.x - minPt.x) / gridSpacing)) + 1;
    ny_ = static_cast<int>(std::ceil((maxPt.y - minPt.y) / gridSpacing)) + 1;
    nz_ = static_cast<int>(std::ceil((maxPt.z - minPt.z) / gridSpacing)) + 1;

    // Extended grid for Toeplitz -> Circulant embedding
    nx2_ = 2 * nx_;
    ny2_ = 2 * ny_;
    nz2_ = 2 * nz_;
    totalGridSize_ = nx2_ * ny2_ * nz2_;

    // Allocate grid buffers
    gridKernel_.resize(totalGridSize_, Complex(0, 0));
    gridWork_.resize(totalGridSize_, Complex(0, 0));

    // Initialize MKL FFT
    InitializeMKLFFT();

    // Build projection and interpolation matrices
    BuildProjectionMatrix(panelCenters, panelAreas);
    BuildInterpolationMatrix(panelCenters);

    // Identify near-field pairs
    IdentifyNearFieldPairs(panelCenters, 3.0 * gridSpacing);

    initialized_ = true;
}

void radTPfft::BuildProjectionMatrix(const std::vector<TVector3d>& centers,
                                      const std::vector<double>& areas) {
    // Build sparse projection matrix P
    // P projects panel values to grid points using interpolation weights

    projMatrix_.nRows = numPanels_;
    projMatrix_.nCols = totalGridSize_;
    projMatrix_.rowPtr.resize(numPanels_ + 1);
    projMatrix_.rowPtr[0] = 0;

    int stencilSize = 1;
    for (int d = 0; d < 3; ++d) {
        stencilSize *= (2 * interpOrder_);
    }

    projMatrix_.colIdx.reserve(numPanels_ * stencilSize);
    projMatrix_.values.reserve(numPanels_ * stencilSize);

    for (int p = 0; p < numPanels_; ++p) {
        double gx, gy, gz;
        WorldToGrid(centers[p], gx, gy, gz);

        int ix0 = static_cast<int>(std::floor(gx)) - interpOrder_ + 1;
        int iy0 = static_cast<int>(std::floor(gy)) - interpOrder_ + 1;
        int iz0 = static_cast<int>(std::floor(gz)) - interpOrder_ + 1;

        // Panel weight (area)
        double panelWeight = (p < static_cast<int>(areas.size())) ? areas[p] : 1.0;

        for (int dz = 0; dz < 2 * interpOrder_; ++dz) {
            for (int dy = 0; dy < 2 * interpOrder_; ++dy) {
                for (int dx = 0; dx < 2 * interpOrder_; ++dx) {
                    int ix = ix0 + dx;
                    int iy = iy0 + dy;
                    int iz = iz0 + dz;

                    // Wrap to grid with periodic boundary
                    ix = ((ix % nx2_) + nx2_) % nx2_;
                    iy = ((iy % ny2_) + ny2_) % ny2_;
                    iz = ((iz % nz2_) + nz2_) % nz2_;

                    // Linear interpolation weight
                    double wx = 1.0 - std::abs(gx - (ix0 + dx));
                    double wy = 1.0 - std::abs(gy - (iy0 + dy));
                    double wz = 1.0 - std::abs(gz - (iz0 + dz));

                    wx = std::max(0.0, wx);
                    wy = std::max(0.0, wy);
                    wz = std::max(0.0, wz);

                    double w = wx * wy * wz * panelWeight;

                    if (w > 1e-14) {
                        int gridIdx = GridIndex(ix, iy, iz);
                        projMatrix_.colIdx.push_back(gridIdx);
                        projMatrix_.values.push_back(w);
                    }
                }
            }
        }
        projMatrix_.rowPtr[p + 1] = static_cast<int>(projMatrix_.colIdx.size());
    }
}

void radTPfft::BuildInterpolationMatrix(const std::vector<TVector3d>& centers) {
    // Interpolation matrix: grid -> panels
    // For now, use same structure as projection (transpose-like)

    interpMatrix_.nRows = numPanels_;
    interpMatrix_.nCols = totalGridSize_;
    interpMatrix_.rowPtr.resize(numPanels_ + 1);
    interpMatrix_.rowPtr[0] = 0;

    int stencilSize = 1;
    for (int d = 0; d < 3; ++d) {
        stencilSize *= (2 * interpOrder_);
    }

    interpMatrix_.colIdx.reserve(numPanels_ * stencilSize);
    interpMatrix_.values.reserve(numPanels_ * stencilSize);

    for (int p = 0; p < numPanels_; ++p) {
        double gx, gy, gz;
        WorldToGrid(centers[p], gx, gy, gz);

        int ix0 = static_cast<int>(std::floor(gx)) - interpOrder_ + 1;
        int iy0 = static_cast<int>(std::floor(gy)) - interpOrder_ + 1;
        int iz0 = static_cast<int>(std::floor(gz)) - interpOrder_ + 1;

        for (int dz = 0; dz < 2 * interpOrder_; ++dz) {
            for (int dy = 0; dy < 2 * interpOrder_; ++dy) {
                for (int dx = 0; dx < 2 * interpOrder_; ++dx) {
                    int ix = ix0 + dx;
                    int iy = iy0 + dy;
                    int iz = iz0 + dz;

                    ix = ((ix % nx2_) + nx2_) % nx2_;
                    iy = ((iy % ny2_) + ny2_) % ny2_;
                    iz = ((iz % nz2_) + nz2_) % nz2_;

                    double wx = 1.0 - std::abs(gx - (ix0 + dx));
                    double wy = 1.0 - std::abs(gy - (iy0 + dy));
                    double wz = 1.0 - std::abs(gz - (iz0 + dz));

                    wx = std::max(0.0, wx);
                    wy = std::max(0.0, wy);
                    wz = std::max(0.0, wz);

                    double w = wx * wy * wz;

                    if (w > 1e-14) {
                        int gridIdx = GridIndex(ix, iy, iz);
                        interpMatrix_.colIdx.push_back(gridIdx);
                        interpMatrix_.values.push_back(w);
                    }
                }
            }
        }
        interpMatrix_.rowPtr[p + 1] = static_cast<int>(interpMatrix_.colIdx.size());
    }
}

void radTPfft::IdentifyNearFieldPairs(const std::vector<TVector3d>& centers,
                                       double threshold) {
    // Find panel pairs within threshold distance
    // These require direct computation (near-field correction)

    double thresh2 = threshold * threshold;

    directPairs_.clear();

    for (int i = 0; i < numPanels_; ++i) {
        for (int j = 0; j < numPanels_; ++j) {
            double dx = centers[i].x - centers[j].x;
            double dy = centers[i].y - centers[j].y;
            double dz = centers[i].z - centers[j].z;
            double r2 = dx*dx + dy*dy + dz*dz;

            if (r2 < thresh2) {
                directPairs_.push_back({i, j});
            }
        }
    }

    directValues_.resize(directPairs_.size(), Complex(0, 0));
    gridValues_.resize(directPairs_.size(), Complex(0, 0));
}

void radTPfft::SetupKernel(const std::function<Complex(double)>& greenFunc) {
    if (!initialized_) {
        throw std::runtime_error("pFFT not initialized");
    }

    // Setup kernel on extended grid (Toeplitz -> Circulant embedding)
    std::fill(gridKernel_.begin(), gridKernel_.end(), Complex(0, 0));

    for (int iz = 0; iz < nz2_; ++iz) {
        int dz = (iz < nz_) ? iz : iz - nz2_;
        for (int iy = 0; iy < ny2_; ++iy) {
            int dy = (iy < ny_) ? iy : iy - ny2_;
            for (int ix = 0; ix < nx2_; ++ix) {
                int dx = (ix < nx_) ? ix : ix - nx2_;

                double r = gridSpacing_ * std::sqrt(
                    static_cast<double>(dx*dx + dy*dy + dz*dz));

                int idx = GridIndex(ix, iy, iz);

                if (r > 1e-15) {
                    gridKernel_[idx] = greenFunc(r);
                } else {
                    // Self-term: handle separately or set to 0
                    gridKernel_[idx] = Complex(0, 0);
                }
            }
        }
    }

    // FFT the kernel
    FFT3D_Forward(gridKernel_);

    kernelSetup_ = true;
}

void radTPfft::SetupKernelMQS() {
    // Magneto-quasi-static: G(r) = 1/(4*pi*r)
    SetupKernel([](double r) -> Complex {
        return Complex(INV_FOUR_PI / r, 0);
    });
}

void radTPfft::SetupKernelFullWave(double frequency, double eps, double mu) {
    // Full-wave: G(r) = exp(-jkr) / (4*pi*r)
    double omega = 2.0 * PI * frequency;
    double k = omega * std::sqrt(eps * mu);

    SetupKernel([k](double r) -> Complex {
        Complex jkr(0, k * r);
        return std::exp(-jkr) * INV_FOUR_PI / r;
    });
}

void radTPfft::SetupDirectCorrection(const std::function<Complex(double)>& greenFunc,
                                      const std::vector<TVector3d>& panelCenters,
                                      double thresholdFactor) {
    if (!initialized_) {
        throw std::runtime_error("pFFT not initialized");
    }

    // Compute direct values and grid contribution for near-field pairs
    for (size_t k = 0; k < directPairs_.size(); ++k) {
        int i = directPairs_[k].first;
        int j = directPairs_[k].second;

        double dx = panelCenters[i].x - panelCenters[j].x;
        double dy = panelCenters[i].y - panelCenters[j].y;
        double dz = panelCenters[i].z - panelCenters[j].z;
        double r = std::sqrt(dx*dx + dy*dy + dz*dz);

        if (r > 1e-15) {
            directValues_[k] = greenFunc(r);
        } else {
            // Self-term: needs special handling
            directValues_[k] = Complex(0, 0);
        }

        // TODO: Compute and subtract grid contribution for correction
        gridValues_[k] = Complex(0, 0);
    }
}

void radTPfft::MatVec(const std::vector<Complex>& x, std::vector<Complex>& y) {
    if (!initialized_ || !kernelSetup_) {
        throw std::runtime_error("pFFT not properly initialized");
    }

    if (static_cast<int>(x.size()) != numPanels_) {
        throw std::invalid_argument("Input vector size mismatch");
    }

    y.resize(numPanels_);
    std::fill(y.begin(), y.end(), Complex(0, 0));

    // Step 1: Project panels to grid
    std::fill(gridWork_.begin(), gridWork_.end(), Complex(0, 0));

    for (int i = 0; i < numPanels_; ++i) {
        for (int k = projMatrix_.rowPtr[i]; k < projMatrix_.rowPtr[i + 1]; ++k) {
            int gridIdx = projMatrix_.colIdx[k];
            double w = projMatrix_.values[k];
            gridWork_[gridIdx] += w * x[i];
        }
    }

    // Step 2: Forward FFT
    FFT3D_Forward(gridWork_);

    // Step 3: Multiply by kernel in frequency domain
    for (int i = 0; i < totalGridSize_; ++i) {
        gridWork_[i] *= gridKernel_[i];
    }

    // Step 4: Backward FFT
    FFT3D_Backward(gridWork_);

    // Step 5: Interpolate grid to panels (far-field contribution)
    for (int i = 0; i < numPanels_; ++i) {
        for (int k = interpMatrix_.rowPtr[i]; k < interpMatrix_.rowPtr[i + 1]; ++k) {
            int gridIdx = interpMatrix_.colIdx[k];
            double w = interpMatrix_.values[k];
            y[i] += w * gridWork_[gridIdx];
        }
    }

    // Step 6: Add direct (near-field) correction
    for (size_t k = 0; k < directPairs_.size(); ++k) {
        int i = directPairs_[k].first;
        int j = directPairs_[k].second;

        // Direct value minus grid contribution (precorrection)
        Complex correction = directValues_[k] - gridValues_[k];
        y[i] += correction * x[j];
    }
}

} // namespace radia
