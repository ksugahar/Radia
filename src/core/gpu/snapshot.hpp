/* SPDX-License-Identifier: MPL-2.0 */
#pragma once
#include <climits>
#include <cstddef>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace radia_sparse {
static_assert(sizeof(int) == 4, "Sparse CUDA requires int32 host indices");
// Owned host snapshots and borrowed device views share a neutral algebraic
// boundary. No NGSolve, Python, MATLAB or CUDA headers enter this interface.
struct Node {
    enum Kind { null, number, logical, text, record, list, integers, reals } kind = null;
    double scalar = 0;
    bool boolean = false;
    std::string word;
    std::map<std::string, Node> fields;
    std::vector<Node> children;
    std::vector<int> int_values;
    std::vector<double> real_values;
    const void* device = nullptr;
    size_t count = 0;
    const Node& at(const char* name) const {
        if (kind != record || !fields.count(name))
            throw std::invalid_argument(std::string("missing snapshot field: ") + name);
        return fields.at(name);
    }
    double real() const {
        if (kind != number || !std::isfinite(scalar))
            throw std::invalid_argument("finite real snapshot scalar required");
        return scalar;
    }
    int integer() const {
        double value = real();
        if (value < 0 || value > INT_MAX || value != std::floor(value))
            throw std::invalid_argument("nonnegative int32 snapshot scalar required");
        return int(value);
    }
    bool flag() const {
        if (kind != logical) throw std::invalid_argument("logical snapshot flag required");
        return boolean;
    }
    const std::string& string() const {
        if (kind != text) throw std::invalid_argument("snapshot text required");
        return word;
    }
    const std::vector<Node>& items() const {
        if (kind != list) throw std::invalid_argument("snapshot list required");
        return children;
    }
    template<class V> const V* data(Kind expected, size_t size) const {
        if (kind != expected || count != size || (size && !device))
            throw std::invalid_argument("snapshot device array type/length mismatch");
        return static_cast<const V*>(device);
    }
};
struct Stats { int iterations = 0; double residual = 0; };
void* create(const Node&);
Stats solve(void*, const double*, size_t, double*, double, int);
void apply(void*, const double*, size_t, double*);
void release(void*) noexcept;
}
