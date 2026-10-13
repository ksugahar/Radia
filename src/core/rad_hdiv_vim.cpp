/* rad_hdiv_vim.cpp -- analytic charge-potential kernels for BDM1 HDiv-VIM. */
#include "rad_hdiv_vim.h"
#include <cmath>
#include <array>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <map>
#include <vector>

#ifndef RADIA_FIELD_SERIAL_ONLY
#include "rad_parallel.h"
#endif

#define RAD_FIELD_NAMESPACE rad_hdiv
#define RAD_FIELD_KERNEL
#define RAD_FIELD_POLY_CAP 18
#include "rad_hdiv_field_math.inc"
#undef RAD_FIELD_POLY_CAP
#undef RAD_FIELD_KERNEL
#undef RAD_FIELD_NAMESPACE
