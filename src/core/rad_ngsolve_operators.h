#ifndef RAD_NGSOLVE_OPERATORS_H
#define RAD_NGSOLVE_OPERATORS_H

#include "rad_hacapk_hdiv.h"

#include <basematrix.hpp>
#include <vvector.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace radia::ngsolve_bridge {

class HDivDemagMatrix final : public ngla::BaseMatrix {
public:
    explicit HDivDemagMatrix(std::shared_ptr<RadHACApKChargeGram> gram)
        : gram_(std::move(gram)), ndof_(gram_ ? gram_->ConfiguredNFace() : 0) {
        if (!gram_ || !gram_->HasConfiguredChargeMap())
            throw std::runtime_error(
                "HDivDemagMatrix: charge map is not configured");
    }

    void Mult(const ngla::BaseVector& x, ngla::BaseVector& y) const override {
        auto xv = x.FV<double>();
        auto yv = y.FV<double>();
        if (xv.Size() != ndof_ || yv.Size() != ndof_)
            throw std::runtime_error(
                "HDivDemagMatrix.Mult: vector size mismatch");
        gram_->ApplyConfiguredDemag(xv.Data(), yv.Data(), true);
    }

    void MultAdd(double scale, const ngla::BaseVector& x,
                 ngla::BaseVector& y) const override {
        auto xv = x.FV<double>();
        auto yv = y.FV<double>();
        if (xv.Size() != ndof_ || yv.Size() != ndof_)
            throw std::runtime_error(
                "HDivDemagMatrix.MultAdd: vector size mismatch");
        gram_->ApplyConfiguredDemagAdd(scale, xv.Data(), yv.Data(), true);
    }

    void MultTransAdd(double scale, const ngla::BaseVector& x,
                      ngla::BaseVector& y) const override {
        MultAdd(scale, x, y);
    }

    int VHeight() const override { return ndof_; }
    int VWidth() const override { return ndof_; }
    ngla::AutoVector CreateRowVector() const override {
        return std::make_unique<ngla::VVector<double>>(ndof_);
    }
    ngla::AutoVector CreateColVector() const override {
        return std::make_unique<ngla::VVector<double>>(ndof_);
    }

private:
    std::shared_ptr<RadHACApKChargeGram> gram_;
    int ndof_ = 0;
};

class ComplexDiagonalInverseMatrix final : public ngla::BaseMatrix {
public:
    explicit ComplexDiagonalInverseMatrix(
        std::vector<std::complex<double>> inverse_diagonal)
        : inverse_diagonal_(std::move(inverse_diagonal)) {
        if (inverse_diagonal_.empty())
            throw std::invalid_argument("diagonal inverse must not be empty");
    }

    void Mult(const ngla::BaseVector& x, ngla::BaseVector& y) const override {
        auto xv = x.FVComplex();
        auto yv = y.FVComplex();
        if (xv.Size() != static_cast<int>(inverse_diagonal_.size()) ||
            yv.Size() != static_cast<int>(inverse_diagonal_.size()))
            throw std::runtime_error(
                "ComplexDiagonalInverse.Mult: vector size mismatch");
        for (int i = 0; i < xv.Size(); ++i)
            yv[i] = inverse_diagonal_[static_cast<std::size_t>(i)] * xv[i];
    }

    int VHeight() const override {
        return static_cast<int>(inverse_diagonal_.size());
    }
    int VWidth() const override {
        return static_cast<int>(inverse_diagonal_.size());
    }
    bool IsComplex() const override { return true; }
    ngla::AutoVector CreateRowVector() const override {
        return std::make_unique<ngla::VVector<ngcore::Complex>>(VHeight());
    }
    ngla::AutoVector CreateColVector() const override {
        return std::make_unique<ngla::VVector<ngcore::Complex>>(VWidth());
    }

private:
    std::vector<std::complex<double>> inverse_diagonal_;
};

class ProjectedBaseMatrix final : public ngla::BaseMatrix {
public:
    ProjectedBaseMatrix(std::shared_ptr<ngla::BaseMatrix> parent,
                        std::vector<std::complex<double>> projection,
                        int parent_size, int reduced_size)
        : parent_(std::move(parent)), projection_(std::move(projection)),
          parent_size_(parent_size), reduced_size_(reduced_size) {
        if (!parent_ || parent_size_ < 1 || reduced_size_ < 1 ||
            parent_->VHeight() != parent_size_ ||
            parent_->VWidth() != parent_size_ ||
            projection_.size() !=
                static_cast<std::size_t>(parent_size_) * reduced_size_)
            throw std::invalid_argument("ProjectedBaseMatrix shape mismatch");
    }

    void Mult(const ngla::BaseVector& x, ngla::BaseVector& y) const override {
        auto xv = x.FVComplex();
        auto yv = y.FVComplex();
        if (xv.Size() != reduced_size_ || yv.Size() != reduced_size_)
            throw std::runtime_error(
                "ProjectedBaseMatrix.Mult: vector size mismatch");
        std::vector<std::complex<double>> parent_x(
            static_cast<std::size_t>(parent_size_), 0.0);
        for (int row = 0; row < parent_size_; ++row) {
            const auto* projection_row =
                &projection_[static_cast<std::size_t>(row) * reduced_size_];
            for (int col = 0; col < reduced_size_; ++col)
                parent_x[static_cast<std::size_t>(row)] +=
                    projection_row[col] * xv[col];
        }
        const auto parent_y = ApplyParent(parent_x);
        for (int col = 0; col < reduced_size_; ++col) {
            std::complex<double> value = 0.0;
            for (int row = 0; row < parent_size_; ++row)
                value += std::conj(projection_[
                    static_cast<std::size_t>(row) * reduced_size_ + col]) *
                    parent_y[static_cast<std::size_t>(row)];
            yv[col] = value;
        }
    }

    int VHeight() const override { return reduced_size_; }
    int VWidth() const override { return reduced_size_; }
    bool IsComplex() const override { return true; }
    ngla::AutoVector CreateRowVector() const override {
        return std::make_unique<ngla::VVector<ngcore::Complex>>(reduced_size_);
    }
    ngla::AutoVector CreateColVector() const override {
        return std::make_unique<ngla::VVector<ngcore::Complex>>(reduced_size_);
    }

private:
    std::vector<std::complex<double>> ApplyParent(
        const std::vector<std::complex<double>>& input) const {
        std::vector<std::complex<double>> output(
            static_cast<std::size_t>(parent_size_));
        if (parent_->IsComplex()) {
            auto parent_x = parent_->CreateColVector();
            auto parent_y = parent_->CreateRowVector();
            auto xv = parent_x.FVComplex();
            auto yv = parent_y.FVComplex();
            for (int i = 0; i < parent_size_; ++i) xv[i] = input[i];
            parent_->Mult(*parent_x, *parent_y);
            for (int i = 0; i < parent_size_; ++i) output[i] = yv[i];
            return output;
        }
        auto parent_x = parent_->CreateColVector();
        auto parent_y = parent_->CreateRowVector();
        auto xv = parent_x.FV<double>();
        auto yv = parent_y.FV<double>();
        for (int i = 0; i < parent_size_; ++i) xv[i] = input[i].real();
        parent_->Mult(*parent_x, *parent_y);
        for (int i = 0; i < parent_size_; ++i) output[i] = yv[i];
        for (int i = 0; i < parent_size_; ++i) xv[i] = input[i].imag();
        parent_->Mult(*parent_x, *parent_y);
        for (int i = 0; i < parent_size_; ++i)
            output[i] += std::complex<double>(0.0, yv[i]);
        return output;
    }

    std::shared_ptr<ngla::BaseMatrix> parent_;
    std::vector<std::complex<double>> projection_;
    int parent_size_ = 0;
    int reduced_size_ = 0;
};

class ReducedBlockMatrix final : public ngla::BaseMatrix {
public:
    struct Term {
        std::shared_ptr<ngla::BaseMatrix> matrix;
        int start = 0;
        int stop = 0;
        std::complex<double> scale = 1.0;
    };

    ReducedBlockMatrix(std::vector<std::complex<double>> dense, int size,
                       std::vector<Term> terms)
        : dense_(std::move(dense)), size_(size), terms_(std::move(terms)) {
        if (size_ < 1 ||
            dense_.size() != static_cast<std::size_t>(size_) * size_)
            throw std::invalid_argument(
                "ReducedBlockMatrix dense correction must be a non-empty square matrix");
        for (const auto& term : terms_) {
            if (!term.matrix || term.start < 0 || term.stop <= term.start ||
                term.stop > size_ ||
                term.matrix->VHeight() != term.stop - term.start ||
                term.matrix->VWidth() != term.stop - term.start)
                throw std::invalid_argument(
                    "ReducedBlockMatrix embedded term shape mismatch");
            if (!std::isfinite(term.scale.real()) ||
                !std::isfinite(term.scale.imag()))
                throw std::invalid_argument(
                    "ReducedBlockMatrix term scale must be finite");
        }
    }

    void Mult(const ngla::BaseVector& x, ngla::BaseVector& y) const override {
        auto xv = x.FVComplex();
        auto yv = y.FVComplex();
        if (xv.Size() != size_ || yv.Size() != size_)
            throw std::runtime_error(
                "ReducedBlockMatrix.Mult: vector size mismatch");
        for (int row = 0; row < size_; ++row) {
            std::complex<double> value = 0.0;
            const auto* dense_row =
                &dense_[static_cast<std::size_t>(row) * size_];
            for (int col = 0; col < size_; ++col)
                value += dense_row[col] * xv[col];
            yv[row] = value;
        }
        for (const auto& term : terms_) ApplyTerm(term, xv, yv);
    }

    void MultAdd(double scale, const ngla::BaseVector& x,
                 ngla::BaseVector& y) const override {
        auto temporary = CreateRowVector();
        Mult(x, *temporary);
        y += scale * *temporary;
    }

    int VHeight() const override { return size_; }
    int VWidth() const override { return size_; }
    bool IsComplex() const override { return true; }
    ngla::AutoVector CreateRowVector() const override {
        return std::make_unique<ngla::VVector<ngcore::Complex>>(size_);
    }
    ngla::AutoVector CreateColVector() const override {
        return std::make_unique<ngla::VVector<ngcore::Complex>>(size_);
    }

    std::shared_ptr<ComplexDiagonalInverseMatrix> DiagonalPreconditioner(
        double relative_floor = 1.0e-14) const {
        if (!std::isfinite(relative_floor) || relative_floor <= 0.0)
            throw std::invalid_argument("relative_floor must be positive");
        double scale = 0.0;
        for (int i = 0; i < size_; ++i)
            scale = std::max(
                scale, std::abs(dense_[static_cast<std::size_t>(i) * size_ + i]));
        if (scale == 0.0)
            throw std::runtime_error(
                "reduced block diagonal is zero; provide a physical resistance/material block");
        const double floor = relative_floor * scale;
        std::vector<std::complex<double>> inverse(
            static_cast<std::size_t>(size_));
        for (int i = 0; i < size_; ++i) {
            const auto value =
                dense_[static_cast<std::size_t>(i) * size_ + i];
            if (std::abs(value) <= floor)
                throw std::runtime_error(
                    "reduced block diagonal is singular; provide a physical resistance/material block");
            inverse[static_cast<std::size_t>(i)] = 1.0 / value;
        }
        return std::make_shared<ComplexDiagonalInverseMatrix>(
            std::move(inverse));
    }

    int TermCount() const { return static_cast<int>(terms_.size()); }

private:
    static void ApplyTerm(const Term& term,
                          ngbla::FlatVector<ngcore::Complex> x,
                          ngbla::FlatVector<ngcore::Complex> y) {
        const int size = term.stop - term.start;
        if (term.matrix->IsComplex()) {
            auto term_x = term.matrix->CreateColVector();
            auto term_y = term.matrix->CreateRowVector();
            auto xv = term_x.FVComplex();
            auto yv = term_y.FVComplex();
            for (int i = 0; i < size; ++i) xv[i] = x[term.start + i];
            term.matrix->Mult(*term_x, *term_y);
            for (int i = 0; i < size; ++i)
                y[term.start + i] += term.scale * yv[i];
            return;
        }

        auto term_x = term.matrix->CreateColVector();
        auto term_y = term.matrix->CreateRowVector();
        auto xv = term_x.FV<double>();
        auto yv = term_y.FV<double>();
        std::vector<double> real_part(static_cast<std::size_t>(size));
        std::vector<double> imag_part(static_cast<std::size_t>(size));
        for (int i = 0; i < size; ++i)
            xv[i] = x[term.start + i].real();
        term.matrix->Mult(*term_x, *term_y);
        for (int i = 0; i < size; ++i) real_part[i] = yv[i];
        for (int i = 0; i < size; ++i)
            xv[i] = x[term.start + i].imag();
        term.matrix->Mult(*term_x, *term_y);
        for (int i = 0; i < size; ++i) imag_part[i] = yv[i];
        for (int i = 0; i < size; ++i)
            y[term.start + i] += term.scale *
                std::complex<double>(real_part[i], imag_part[i]);
    }

    std::vector<std::complex<double>> dense_;
    int size_ = 0;
    std::vector<Term> terms_;
};

} // namespace radia::ngsolve_bridge

#endif
