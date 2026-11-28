#pragma once
// Netgen meshing headers for translation units that also include the Cubit SDK.
// Coreform's CubitDefines.h defines NT on Windows, and Netgen 6.2.2608 uses NT
// as a member name (DynStrideLayout::NT), so the macro is hidden around Netgen.
#pragma push_macro("NT")
#undef NT
#include <meshing.hpp>
#pragma pop_macro("NT")
