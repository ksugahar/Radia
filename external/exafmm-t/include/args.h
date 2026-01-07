#ifndef args_h
#define args_h
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <omp.h>

// Note: getopt.h is not available on Windows MSVC
// Radia doesn't use command-line argument parsing, so we provide a simplified Args class

namespace exafmm_t {
  // Simplified Args class without getopt dependency
  // Used only for default parameter initialization in Radia

  class Args {
  public:
    int ncrit;
    const char * distribution;
    double k;
    int maxlevel;
    int numBodies;
    int P;
    int threads;

    // Default constructor with reasonable defaults
    Args(int argc=0, char ** argv=nullptr) :
      ncrit(64),
      distribution("cube"),
      k(20),
      maxlevel(5),
      numBodies(1000000),
      P(4),
      threads(omp_get_max_threads()) {
      // No command-line parsing - use defaults
      (void)argc;
      (void)argv;
    }

    void print(int stringLength=20) {
      std::cout << std::setw(stringLength) << std::fixed << std::left
                << "ncrit" << " : " << ncrit << std::endl
                << std::setw(stringLength)
                << "distribution" << " : " << distribution << std::endl
                << std::setw(stringLength)
                << "numBodies" << " : " << numBodies << std::endl
                << std::setw(stringLength)
                << "P" << " : " << this->P << std::endl
                << std::setw(stringLength)
                << "maxlevel" << " : " << maxlevel << std::endl
                << std::setw(stringLength)
                << "threads" << " : " << threads << std::endl
                << std::setw(stringLength)
                << "wavenumber" << " : " << k << std::endl;
    }
  };
}
#endif
