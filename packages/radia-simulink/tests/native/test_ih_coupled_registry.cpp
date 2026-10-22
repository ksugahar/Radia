// Standalone host contract test; fake physics is NOT numerical/GPU acceptance.
#include "../../native/ih_coupled_registry.hpp"
#include <atomic>
#include <cmath>
#include <iostream>
#include <thread>

struct Controls { double current, ambient; int mode; };
struct Snapshot { unsigned long long revision; double value; };
struct FakeDomain {
    std::shared_ptr<const Snapshot> accepted = std::make_shared<Snapshot>(Snapshot{0, 0.});
    bool fail{}, poison{}, closed{};
    unsigned calls{}, resets{}, closes{};
    std::shared_ptr<const Snapshot> snapshot() const noexcept { return accepted; }
    std::shared_ptr<const Snapshot> advance(const Controls& c) {
        ++calls;
        if (fail || poison) throw std::runtime_error("injected trial failure");
        std::shared_ptr<const Snapshot> next = std::make_shared<Snapshot>(Snapshot{accepted->revision + 1, c.current});
        accepted.swap(next);
        return accepted;
    }
    std::shared_ptr<const Snapshot> reset() {
        ++resets;
        if (fail) throw std::runtime_error("injected reset failure");
        std::shared_ptr<const Snapshot> next = std::make_shared<Snapshot>(Snapshot{accepted->revision + 1, 0.});
        accepted.swap(next);
        return accepted;
    }
    void close() noexcept { if (!closed) { closed = true; ++closes; } }
};
struct Policy {
    static void validate(const Controls& c) {
        if (!std::isfinite(c.current) || !std::isfinite(c.ambient) || c.ambient <= 0 ||
            (c.mode == 0 ? c.current <= 0 : c.mode != 1 || c.current != 0))
            throw std::invalid_argument("invalid controls");
    }
    static bool equal(const Controls& a, const Controls& b) noexcept {
        return a.current == b.current && a.ambient == b.ambient && a.mode == b.mode;
    }
    static bool poisoned(const FakeDomain& d) noexcept { return d.poison; }
};
using Owner = radia_ih_runtime::Owner<FakeDomain, Controls, Policy>;
void check(bool value) { if (!value) throw std::runtime_error("contract failed"); }
template<class F> void rejects(F&& f) {
    bool rejected = false;
    try { f(); } catch (const std::exception&) { rejected = true; }
    check(rejected);
}
int main() {
    auto d = std::make_unique<FakeDomain>(); auto* observed = d.get();
    auto owner = std::make_shared<Owner>(std::move(d));
    radia_ih_runtime::Registry<Owner, 0x49485059> python;
    radia_ih_runtime::Registry<Owner, 0x49484d58> mex;
    auto token = python.insert(owner);
    check(python.lookup(token) == owner);
    rejects([&] { mex.lookup(token); });
    auto initial = owner->outputs(0);
    check(initial.sample_id == 0 && initial.epoch == 0 && initial.snapshot->revision == 0);
    rejects([&] { owner->advance(0, 0, {1., 300., 0}); });
    rejects([&] { owner->advance(0, 2, {1., 300., 0}); });
    rejects([&] { owner->advance(0, 1, {-1., 300., 0}); });
    check(observed->calls == 0);
    auto first = owner->advance(0, 1, {2., 300., 0});
    check(first.epoch == 1 && first.snapshot->revision == 1);
    auto replay = owner->advance(0, 1, {2., 300., 0});
    check(replay.snapshot == first.snapshot && observed->calls == 1);
    rejects([&] { owner->advance(0, 1, {3., 300., 0}); });
    rejects([&] { owner->advance(0, 1, {2., 301., 0}); });
    // Throwing in the language publisher occurs after runtime commit.
    rejects([&] { owner->advance(0, 2, {3., 300., 0}); throw std::bad_alloc(); });
    auto retrieved = owner->advance(0, 2, {3., 300., 0});
    check(observed->calls == 2 && retrieved.sample_id == 2 && retrieved.epoch == 2);
    observed->fail = true;
    rejects([&] { owner->advance(0, 3, {4., 300., 0}); });
    rejects([&] { owner->reset(0); });
    auto unchanged = owner->outputs(0);
    check(unchanged.snapshot == retrieved.snapshot && unchanged.epoch == 2);
    observed->fail = false;
    auto reset = owner->reset(0);
    check(reset.generation == 1 && reset.sample_id == 0 && reset.epoch == 3 &&
          reset.snapshot->revision == 3);
    rejects([&] { owner->outputs(0); });
    // Concurrent identical calls commit exactly once under the owner mutex.
    std::atomic<unsigned> ready{0}; std::atomic<bool> go{false};
    auto run = [&] { ++ready; while (!go.load()) std::this_thread::yield();
                    owner->advance(1, 1, {5., 300., 0}); };
    std::thread a(run), b(run);
    while (ready != 2) std::this_thread::yield();
    go = true; a.join(); b.join();
    auto concurrent = owner->outputs(1);
    check(concurrent.snapshot->revision == 4 && concurrent.epoch == 4);
    observed->poison = true;
    auto poisoned = owner->outputs(1);
    check(poisoned.poisoned && poisoned.snapshot == concurrent.snapshot);
    rejects([&] { owner->advance(1, 2, {6., 300., 0}); });
    rejects([&] { owner->reset(1); });
    python.destroy(token);
    check(observed->closed && observed->closes == 1);
    owner->close(); check(observed->closes == 1);
    rejects([&] { python.lookup(token); });
    rejects([&] { python.destroy(token); });
    rejects([&] { owner->outputs(1); });
    check(first.snapshot->value == 2. && concurrent.snapshot->value == 5.);
    auto isolated = std::make_shared<Owner>(std::make_unique<FakeDomain>());
    auto token2 = python.insert(isolated); check(token2 != token);
    check(isolated->advance(0, 1, {0., 300., 1}).sample_id == 1);
    python.destroy(token2);
    std::cout << "IH registry host contracts PASS (fake domain; no GPU/physics)\n";
}
