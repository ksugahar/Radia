// Host-only runtime ownership. No Python/MEX/CUDA dependency or physics.
#pragma once
#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <utility>

namespace radia_ih_runtime {
// Policy owns typed input validation, canonical equality and device-poison query.
// All language objects must be copied before entering this runtime.
template<class Domain, class Controls, class Policy>
class Owner final {
public:
    using Snapshot = decltype(std::declval<const Domain&>().snapshot());
    struct Envelope {
        std::uint64_t generation, sample_id, epoch;
        bool poisoned;
        Snapshot snapshot;
    };
    explicit Owner(std::unique_ptr<Domain> domain) : domain_(std::move(domain)) {
        static_assert(std::is_trivially_copyable<Controls>::value,
                      "Controls must be normalized owned scalars");
        static_assert(noexcept(Policy::equal(std::declval<const Controls&>(),
                                            std::declval<const Controls&>())),
                      "Control equality must not throw");
        static_assert(noexcept(Policy::poisoned(std::declval<const Domain&>())),
                      "Poison query must not throw");
        static_assert(noexcept(std::declval<const Domain&>().snapshot()),
                      "Host snapshot capture must not throw");
        static_assert(std::is_nothrow_copy_assignable<Controls>::value,
                      "Controls assignment must not throw after commit");
        static_assert(std::is_nothrow_move_assignable<Snapshot>::value,
                      "Snapshot adoption must not throw after commit");
        static_assert(std::is_nothrow_copy_constructible<Snapshot>::value,
                      "Snapshot envelope capture must not throw");
        static_assert(noexcept(std::declval<Domain&>().close()),
                      "Domain close must not throw");
        if (!domain_) throw std::invalid_argument("Null IH domain");
        committed_ = domain_->snapshot();
        if (!committed_) throw std::invalid_argument("Null initial IH snapshot");
    }
    Owner(const Owner&) = delete;
    Owner& operator=(const Owner&) = delete;
    ~Owner() noexcept { close(); }

    Envelope advance(std::uint64_t generation, std::uint64_t sample,
                     const Controls& controls) {
        Policy::validate(controls); // No language callback, allocation or borrowed data.
        std::lock_guard<std::mutex> lock(mutex_);
        require_live(generation);
        require_usable();
        if (sample == sample_ && has_controls_) {
            if (!Policy::equal(controls, controls_))
                throw std::invalid_argument("IH sample controls changed");
            return envelope();
        }
        if (sample_ == maximum || epoch_ == maximum)
            throw std::overflow_error("IH sample/epoch overflow");
        if (sample != sample_ + 1)
            throw std::invalid_argument("IH sample must be the next committed sample");
        // Domain guarantees a nonnull immutable snapshot and noexcept commit.
        // Nothing after its successful return allocates or calls language code.
        auto next = domain_->advance(controls);
        committed_ = std::move(next);
        controls_ = controls;
        has_controls_ = true;
        sample_ = sample;
        ++epoch_;
        return envelope();
    }
    Envelope reset(std::uint64_t generation) {
        std::lock_guard<std::mutex> lock(mutex_);
        require_live(generation);
        require_usable();
        if (generation_ == maximum || epoch_ == maximum)
            throw std::overflow_error("IH generation/epoch overflow");
        auto next = domain_->reset();
        committed_ = std::move(next);
        ++generation_;
        ++epoch_;
        sample_ = 0;
        has_controls_ = false;
        return envelope();
    }
    Envelope outputs(std::uint64_t generation) const {
        std::lock_guard<std::mutex> lock(mutex_);
        require_live(generation);
        return envelope(); // Poisoned owners retain diagnostics, never compute.
    }
    void close() noexcept {
        std::lock_guard<std::mutex> lock(mutex_);
        if (!closed_) {
            domain_->close();
            closed_ = true;
        }
    }
private:
    static constexpr auto maximum = std::numeric_limits<std::uint64_t>::max();
    void require_live(std::uint64_t generation) const {
        if (closed_) throw std::invalid_argument("Closed IH owner");
        if (generation != generation_) throw std::invalid_argument("Stale IH generation");
    }
    void require_usable() const {
        if (Policy::poisoned(*domain_)) throw std::runtime_error("Poisoned IH owner");
    }
    Envelope envelope() const noexcept {
        return {generation_, sample_, epoch_, Policy::poisoned(*domain_), committed_};
    }
    std::unique_ptr<Domain> domain_;
    mutable std::mutex mutex_;
    Snapshot committed_;
    Controls controls_{};
    std::uint64_t generation_{0}, sample_{0}, epoch_{0};
    bool has_controls_{false}, closed_{false};
};

// Namespace tags MUST differ across runtime modules and handle types. Tokens are
// never reused within this registry lifetime; overflow rejects, never wraps.
// A binding owns exactly one registry per namespace and never transfers tokens.
template<class RuntimeOwner, std::uint32_t NamespaceTag>
class Registry final {
    static_assert(NamespaceTag != 0, "Checked IH handles require a namespace tag");
public:
    std::uint64_t insert(std::shared_ptr<RuntimeOwner> owner) {
        if (!owner) throw std::invalid_argument("Null IH owner");
        std::lock_guard<std::mutex> lock(mutex_);
        if (next_ == 0) throw std::overflow_error("IH handle namespace exhausted");
        const auto token = (std::uint64_t(NamespaceTag) << 32) | next_;
        entries_.emplace(token, std::move(owner)); // Failure does not consume token.
        ++next_;
        return token;
    }
    std::shared_ptr<RuntimeOwner> lookup(std::uint64_t token) const {
        std::lock_guard<std::mutex> lock(mutex_);
        return checked(token);
    }
    void destroy(std::uint64_t token) {
        std::shared_ptr<RuntimeOwner> owner;
        {
            std::lock_guard<std::mutex> lock(mutex_);
            owner = checked(token);
            entries_.erase(token);
        }
        owner->close(); // Never hold global registry lock across domain work.
    }
private:
    std::shared_ptr<RuntimeOwner> checked(std::uint64_t token) const {
        if ((token >> 32) != NamespaceTag)
            throw std::invalid_argument("Wrong IH handle runtime/type");
        const auto item = entries_.find(token);
        if (item == entries_.end()) throw std::invalid_argument("Stale IH handle");
        return item->second;
    }
    mutable std::mutex mutex_;
    std::unordered_map<std::uint64_t, std::shared_ptr<RuntimeOwner>> entries_;
    std::uint32_t next_{1};
};
} // namespace radia_ih_runtime
