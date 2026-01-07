/*
 * rad_pfft.h
 *
 * Precorrected FFT (pFFT) acceleration for integral equation solver
 * Based on FastImp methodology, using Intel MKL DFTI for FFT
 *
 * References:
 * [1] Z. Zhu, B. Song, J.K. White, "Algorithms in FastImp", IEEE TCAD, 2005
 *
 * Part of Radia project
 */

#ifndef RAD_PFFT_H
#define RAD_PFFT_H

#include <mkl.h>
#include <complex>
#include <vector>
#include <functional>
#include <stdexcept>
#include "gmvect.h"

namespace radia {

// Complex type alias
using Complex = std::complex<double>;

/**
 * @brief Sparse matrix in CSR (Compressed Sparse Row) format
 */
struct SparseMatrixCSR {
    std::vector<int> rowPtr;      // Row pointers (size: nRows + 1)
    std::vector<int> colIdx;      // Column indices
    std::vector<double> values;   // Non-zero values
    int nRows = 0;
    int nCols = 0;

    void clear() {
        rowPtr.clear();
        colIdx.clear();
        values.clear();
        nRows = nCols = 0;
    }
};

/**
 * @brief pFFT acceleration engine for integral equation matrix-vector products
 *
 * Implements the precorrected FFT algorithm:
 * y = A*x = (P^T * F^{-1} * diag(G_hat) * F * P) * x + D * x
 *
 * where:
 *   P: projection matrix (panels -> grid)
 *   F: 3D FFT
 *   G_hat: FFT of Green's function on grid
 *   D: direct (near-field) correction matrix
 */
class radTPfft {
public:
    radTPfft();
    ~radTPfft();

    // Non-copyable
    radTPfft(const radTPfft&) = delete;
    radTPfft& operator=(const radTPfft&) = delete;

    /**
     * @brief Initialize pFFT grid and matrices
     * @param panelCenters Centers of surface panels
     * @param panelAreas Areas of panels (for weighting)
     * @param gridSpacing Grid spacing (auto-computed if <= 0)
     * @param interpOrder Interpolation order (default: 2)
     */
    void Initialize(const std::vector<TVector3d>& panelCenters,
                    const std::vector<double>& panelAreas,
                    double gridSpacing = -1.0,
                    int interpOrder = 2);

    /**
     * @brief Setup kernel (Green's function) on grid
     * @param greenFunc Function: r -> G(r), where r is distance
     */
    void SetupKernel(const std::function<Complex(double)>& greenFunc);

    /**
     * @brief Setup kernel for specific frequency
     * @param frequency Frequency in Hz
     * @param eps Permittivity
     * @param mu Permeability
     */
    void SetupKernelFullWave(double frequency, double eps, double mu);

    /**
     * @brief Setup kernel for MQS (quasi-static)
     */
    void SetupKernelMQS();

    /**
     * @brief Compute matrix-vector product y = A * x
     * @param x Input vector (size: numPanels)
     * @param y Output vector (size: numPanels)
     */
    void MatVec(const std::vector<Complex>& x, std::vector<Complex>& y);

    /**
     * @brief Setup direct (near-field) correction matrix
     * @param greenFunc Green's function for direct computation
     * @param panelCenters Panel centers
     * @param threshold Distance threshold for near-field (multiple of grid spacing)
     */
    void SetupDirectCorrection(const std::function<Complex(double)>& greenFunc,
                               const std::vector<TVector3d>& panelCenters,
                               double thresholdFactor = 3.0);

    /**
     * @brief Get number of panels
     */
    int NumPanels() const { return numPanels_; }

    /**
     * @brief Get grid dimensions
     */
    void GetGridDims(int& nx, int& ny, int& nz) const {
        nx = nx_; ny = ny_; nz = nz_;
    }

    /**
     * @brief Check if initialized
     */
    bool IsInitialized() const { return initialized_; }

    /**
     * @brief Release all resources
     */
    void Cleanup();

private:
    // Grid dimensions
    int nx_, ny_, nz_;           // Original grid size
    int nx2_, ny2_, nz2_;        // Extended grid (2x for Toeplitz->Circulant)
    int totalGridSize_;          // nx2_ * ny2_ * nz2_

    // Grid origin and spacing
    TVector3d gridOrigin_;
    double gridSpacing_;
    int interpOrder_;

    // Number of panels
    int numPanels_;

    // MKL FFT handle
    DFTI_DESCRIPTOR_HANDLE fftHandle_;
    bool fftInitialized_;

    // Grid data
    std::vector<Complex> gridKernel_;     // FFT of kernel (frequency domain)
    std::vector<Complex> gridWork_;       // Working buffer

    // Projection matrix (panels -> grid)
    SparseMatrixCSR projMatrix_;

    // Interpolation matrix (grid -> panels)
    SparseMatrixCSR interpMatrix_;

    // Direct correction (near-field pairs)
    std::vector<std::pair<int, int>> directPairs_;
    std::vector<Complex> directValues_;
    std::vector<Complex> gridValues_;     // Grid contribution to subtract

    // State
    bool initialized_;
    bool kernelSetup_;

    // Private methods
    void InitializeMKLFFT();
    void DestroyMKLFFT();

    void FFT3D_Forward(std::vector<Complex>& data);
    void FFT3D_Backward(std::vector<Complex>& data);

    void BuildProjectionMatrix(const std::vector<TVector3d>& centers,
                               const std::vector<double>& areas);
    void BuildInterpolationMatrix(const std::vector<TVector3d>& centers);

    void IdentifyNearFieldPairs(const std::vector<TVector3d>& centers,
                                double threshold);

    int GridIndex(int ix, int iy, int iz) const {
        return ix + nx2_ * (iy + ny2_ * iz);
    }

    void GridCoords(int idx, int& ix, int& iy, int& iz) const {
        ix = idx % nx2_;
        int tmp = idx / nx2_;
        iy = tmp % ny2_;
        iz = tmp / ny2_;
    }

    // Convert world coordinates to grid indices
    void WorldToGrid(const TVector3d& point, double& gx, double& gy, double& gz) const {
        gx = (point.x - gridOrigin_.x) / gridSpacing_;
        gy = (point.y - gridOrigin_.y) / gridSpacing_;
        gz = (point.z - gridOrigin_.z) / gridSpacing_;
    }
};

} // namespace radia

#endif // RAD_PFFT_H
