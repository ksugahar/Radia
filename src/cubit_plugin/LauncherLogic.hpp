// Pure-logic functions used by the Radia-NGSolve launcher dialog
// (RadiaComp.cpp, .ccl) AND by the headless verification command
// (VerifyLauncherCommand.cpp, .ccm).
//
// Single source of truth — keeping these in a shared header is what
// makes the headless test (`radia_export verify_launcher`) actually
// catch bugs in the GUI launcher path (2026-04-14 hardening: see
// memory/feedback_gui_test_with_deploy.md).

#pragma once

#include "CubitInterface.hpp"
#include <set>
#include <string>
#include <vector>

namespace RadiaLauncherLogic {

// ----------------------------------------------------------------
// get_model_labels — list block + sideset names, lower-cased,
// deduplicated. Used to validate dialog "[OK] / [MISSING]" rows.
//
// IMPORTANT: sidesets MUST use get_exodus_entity_name (not
// get_entity_name — that returns "" for sidesets and made the dialog
// falsely flag source/sink as MISSING on 2026-04-14 / 100号機).
// ----------------------------------------------------------------
inline std::vector<std::string> get_model_labels()
{
  std::set<std::string> names;

  std::vector<int> bids = CubitInterface::parse_cubit_list("block", "all");
  for (int bid : bids) {
    std::string n = CubitInterface::get_block_name(bid);
    if (!n.empty()) {
      for (auto &c : n) c = (char)tolower((unsigned char)c);
      names.insert(n);
    }
  }

  std::vector<int> sids = CubitInterface::parse_cubit_list("sideset", "all");
  for (int sid : sids) {
    std::string n = CubitInterface::get_exodus_entity_name("sideset", sid);
    if (!n.empty()) {
      for (auto &c : n) c = (char)tolower((unsigned char)c);
      names.insert(n);
    }
  }

  return std::vector<std::string>(names.begin(), names.end());
}

// ----------------------------------------------------------------
// default_vol_path — derive a usable .vol output path from a .jou
// path. Returns "" if the input is unusable (empty, "/", or paths
// that would produce "/.vol" / ".vol"). Caller should fall back to
// a user-home default in that case.
// ----------------------------------------------------------------
inline std::string default_vol_path(const std::string& jouPath)
{
  if (jouPath.empty() || jouPath == "/" || jouPath == "\\")
    return "";

  std::string p = jouPath;
  // Normalize separators.
  for (auto &c : p) if (c == '\\') c = '/';

  // Strip trailing .jou (case-insensitive).
  if (p.size() >= 4) {
    std::string suf = p.substr(p.size() - 4);
    for (auto &c : suf) c = (char)tolower((unsigned char)c);
    if (suf == ".jou")
      p = p.substr(0, p.size() - 4);
  }

  // Final basename empty? Reject.
  auto slash = p.find_last_of('/');
  std::string base = (slash == std::string::npos) ? p : p.substr(slash + 1);
  if (base.empty())
    return "";

  return p + ".vol";
}

// ----------------------------------------------------------------
// validate_labels — every name in *need* must be in *have* (case
// insensitive). Returns the missing names (empty = OK).
// ----------------------------------------------------------------
inline std::vector<std::string> validate_labels(
    const std::vector<std::string>& need,
    const std::vector<std::string>& have)
{
  std::set<std::string> haveSet;
  for (auto &n : have) {
    std::string lc = n;
    for (auto &c : lc) c = (char)tolower((unsigned char)c);
    haveSet.insert(lc);
  }
  std::vector<std::string> missing;
  for (auto &n : need) {
    std::string lc = n;
    for (auto &c : lc) c = (char)tolower((unsigned char)c);
    if (!haveSet.count(lc))
      missing.push_back(n);
  }
  return missing;
}

}  // namespace RadiaLauncherLogic
