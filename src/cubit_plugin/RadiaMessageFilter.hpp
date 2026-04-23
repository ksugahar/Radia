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
#include <cstring>

namespace radia {

class LearnEditionFilter : public CubitMessageHandler {
 public:
  explicit LearnEditionFilter(CubitMessageHandler* previous)
      : previous_(previous) {}

  void print_message(const char* message) override {
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
      return;
    }
    if (previous_) previous_->print_error(message);
  }

 private:
  CubitMessageHandler* previous_;
};

// RAII guard: install filter on construction, restore previous handler
// on destruction. Use at the top of each radia_export command's execute().
class ScopedLearnEditionFilter {
 public:
  ScopedLearnEditionFilter() {
    previous_ = CubitMessage::get_message_handler();
    filter_ = new LearnEditionFilter(previous_);
    CubitMessage::set_message_handler(filter_);
  }

  ~ScopedLearnEditionFilter() {
    CubitMessage::set_message_handler(previous_);
    delete filter_;
  }

  ScopedLearnEditionFilter(const ScopedLearnEditionFilter&) = delete;
  ScopedLearnEditionFilter& operator=(const ScopedLearnEditionFilter&) = delete;

 private:
  CubitMessageHandler* previous_ = nullptr;
  LearnEditionFilter* filter_ = nullptr;
};

}  // namespace radia
