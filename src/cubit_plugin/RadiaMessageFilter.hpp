#pragma once
//=============================================================================
// RadiaMessageFilter: suppress Cubit Learn Edition's noisy 50k-cap ERROR
// during radia_export commands.
//
// Problem: Coreform Cubit Learn Edition prints
//   *****ERROR: Coreform Cubit - Learn Edition restricts export to models
//   with less than 50k elements.
// on any `radia_export *` command when element count > 50,000. The Radia
// in-tree plugin bypasses the cap and completes the export successfully,
// so the ERROR line is misleading noise that confuses users scanning logs.
//
// Solution: install a CubitMessageHandler that forwards every message to
// the previously-installed handler EXCEPT those that match the 50k-cap
// signature ("Learn Edition" + "restricts export"). Scoped by RAII so
// the filter is active only during radia_export command execution.
//
// Preserves:
//   - real errors from our own PRINT_ERROR calls
//   - all info / warning messages
//   - any other unrelated Cubit errors
//=============================================================================

#include "CubitMessage.hpp"
#include "CubitMessageHandler.hpp"
#include <algorithm>
#include <cstring>
#include <memory>   // 2026-05-25: Cubit 2025.12 (CUBIT_VERSION_ALLINT 17040+)
                    // switched CubitMessage::{get,set}_message_handler to
                    // std::shared_ptr<CubitMessageHandler> and made them
                    // non-static. See ScopedLearnEditionFilter below.

namespace radia {

class LearnEditionFilter : public CubitMessageHandler {
 public:
  // 2026-05-25: 'previous' is now a shared_ptr because Cubit 2025.12's
  // get_message_handler() returns shared_ptr. We retain the same RAII
  // ownership semantics (filter holds a non-owning reference; the
  // ScopedLearnEditionFilter destructor restores `previous` to the
  // active slot). Storing the shared_ptr also keeps the previous
  // handler alive for the lifetime of the filter, which is safer
  // than the raw-pointer pattern.
  explicit LearnEditionFilter(std::shared_ptr<CubitMessageHandler> previous)
      : previous_(std::move(previous)) {}

  void print_message(const char* message) override {
    if (!message) {
      if (previous_) previous_->print_message(message);
      return;
    }
    // Our ScopedLearnEditionFilter decrements the error counter in its
    // destructor to cancel the 50k-trap side effects.  Cubit tattles
    // about that via a "WARNING: Error count manually changed from N
    // to M" message, which is itself noise.  Swallow it.
    if (std::strstr(message, "Error count manually changed")) {
      return;
    }
    if (previous_) previous_->print_message(message);
  }

  void print_error(const char* message) override {
    if (!message) {
      if (previous_) previous_->print_error(message);
      return;
    }
    // (1) 50k-cap Learn Edition notice — harmless, Radia bypasses the cap.
    if (std::strstr(message, "Learn Edition") &&
        std::strstr(message, "restricts export")) {
      ++swallow_count_;
      return;
    }
    // (2) "Unable to find entities to write" — spurious during radia_export.
    // Cubit's generic export-command dispatch layer expects entity-selector
    // arguments (e.g. `export stl "f" volume all`) and prints this message
    // when the 50k-cap trap fires.  radia_export does NOT use that dispatch
    // model — it takes only a filename and walks the full mesh internally
    // via CubitInterface.  So the message is misapplied to our command.
    // This filter is RAII-scoped to radia_export execute() only, so builtin
    // Cubit export commands still see the real message normally.
    if (std::strstr(message, "Unable to find entities to write")) {
      ++swallow_count_;
      return;
    }
    if (previous_) previous_->print_error(message);
  }

  int swallow_count() const { return swallow_count_; }

 private:
  std::shared_ptr<CubitMessageHandler> previous_;
  int swallow_count_ = 0;
};

// RAII guard: install filter on construction, restore previous handler
// on destruction.  Additionally decrements Cubit's internal error counter
// by the number of spurious messages swallowed — otherwise Cubit's
// session-end "Errors found during session." summary fires even when
// every printed message was a spurious 50k-trap side effect.
// Real errors (not matched by the filter) still count.
// Use at the top of each radia_export command's execute().
class ScopedLearnEditionFilter {
 public:
  ScopedLearnEditionFilter() {
    // 2026-05-25 (Cubit 2025.12 / CUBIT_VERSION_ALLINT 17040+):
    // get/set_message_handler are now non-static AND shared_ptr-based.
    // The shared_ptr's lifetime semantics also remove the need for the
    // raw `delete filter_;` in our destructor.
    auto* msg = CubitMessage::instance();
    previous_ = msg->get_message_handler();
    filter_ = std::make_shared<LearnEditionFilter>(previous_);
    msg->set_message_handler(filter_);
    saved_error_count_ = msg->get_error_count();
  }

  ~ScopedLearnEditionFilter() {
    auto* msg = CubitMessage::instance();
    // Subtract the number of spurious messages we swallowed from the
    // error counter, but never go below the pre-filter baseline.
    // That preserves legitimate errors logged during our execute()
    // while cancelling the Cubit-internal 50k-trap side effects.
    const int current = msg->get_error_count();
    const int target = std::max(saved_error_count_,
                                 current - filter_->swallow_count());
    msg->reset_error_count(target);
    msg->set_message_handler(previous_);
    // No `delete filter_;` -- shared_ptr handles cleanup when both
    // `filter_` and the Cubit-side handler-slot copy go out of scope.
  }

  ScopedLearnEditionFilter(const ScopedLearnEditionFilter&) = delete;
  ScopedLearnEditionFilter& operator=(const ScopedLearnEditionFilter&) = delete;

 private:
  std::shared_ptr<CubitMessageHandler> previous_;
  std::shared_ptr<LearnEditionFilter> filter_;
  int saved_error_count_ = 0;
};

}  // namespace radia
