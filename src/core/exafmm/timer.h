#ifndef timer_h
#define timer_h
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <atomic>
#include <chrono>

// Use C++ chrono for cross-platform timing (no sys/time.h dependency)
// This avoids Windows timeval conflicts with winsock.h

namespace exafmm_t {
  static const int stringLength = 20;           //!< Length of formatted string
  static const int decimal = 7;                 //!< Decimal precision
  static const int wait = 100;                  //!< Waiting time between output of different ranks
  static const int dividerLength = stringLength + decimal + 9;  // length of output section divider

  // Use chrono time_point instead of timeval
  using clock_type = std::chrono::high_resolution_clock;
  using time_point = clock_type::time_point;

  // C++17 inline variables to prevent multiple definitions
  inline std::atomic<long long> flop{0};
  inline time_point current_time;
  inline std::map<std::string, time_point> timer;

  inline void print(std::string s) {
    // if (!VERBOSE | (MPIRANK != 0)) return;
    s += " ";
    std::cout << "--- " << std::setw(stringLength) << std::left
              << std::setfill('-') << s << std::setw(decimal+1) << "-"
              << std::setfill(' ') << std::endl;
  }

  template<typename T>
  inline void print(std::string s, T v, bool fixed=true) {
    std::cout << std::setw(stringLength) << std::left << s << " : ";
    if(fixed)
      std::cout << std::setprecision(decimal) << std::fixed << std::scientific;
    else
      std::cout << std::setprecision(1) << std::scientific;
    std::cout << v << std::endl;
  }

  inline void print_divider(std::string s) {
    s.insert(0, " ");
    s.append(" ");
    int halfLength = (dividerLength - static_cast<int>(s.length())) / 2;
    std::cout << std::string(halfLength, '-') << s
              << std::string(dividerLength-halfLength-static_cast<int>(s.length()), '-') << std::endl;
  }

  inline void add_flop(long long n) {
    flop.fetch_add(n, std::memory_order_relaxed);
  }

  inline void start(std::string event) {
    current_time = clock_type::now();
    timer[event] = current_time;
  }

  inline double stop(std::string event, bool verbose=true) {
    current_time = clock_type::now();
    std::chrono::duration<double> elapsed = current_time - timer[event];
    double eventTime = elapsed.count();
    if (verbose)
      print(event, eventTime);
    return eventTime;
  }
}
#endif
