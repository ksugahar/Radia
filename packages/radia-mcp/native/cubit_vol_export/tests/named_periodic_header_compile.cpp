// Guards: the actual production matcher compiles against the Netgen SDK
// typed-index headers before linking or running a native contract.
#include "NamedPeriodicIdentification.hpp"
#include <type_traits>

static_assert(std::is_same<decltype(radia_cubit_export::periodic_point_index(1)),
                           netgen::PointIndex>::value, "Native point index type required");
static_assert(std::is_same<decltype(radia_cubit_export::periodic_face_index(1)),
                           netgen::FaceRegionIndex>::value, "Native region index type required");

auto named_periodic_header_compile_entry =
    &radia_cubit_export::apply_named_periodic_identification;
