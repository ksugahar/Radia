// Guards: the actual production matcher compiles against the full Netgen
// typed-index headers before linking or running a native contract.
#include "NamedPeriodicIdentification.hpp"
#include <type_traits>

static_assert(std::is_same<decltype(cubit_mesh_export::periodic_point_index(1)),
                           netgen::PointIndex>::value, "Native point index type required");
static_assert(std::is_same<decltype(cubit_mesh_export::periodic_face_index(1)),
                           netgen::FaceRegionIndex>::value, "Native region index type required");

auto named_periodic_header_compile_entry =
    &cubit_mesh_export::apply_named_periodic_identification;
