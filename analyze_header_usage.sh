#!/bin/bash

echo "=== Header File Usage Analysis ==="
echo ""

declare -A files=(
	["radappl.h"]="rad_application.h"
	["radauxst.h"]="rad_auxiliary_structures.h"
	["radcnvrg.h"]="rad_convergence.h"
	["radexcep.h"]="rad_exception.h"
	["radg.h"]="rad_geometry_base.h"
	["radg3da1.h"]="rad_geometry_3d_aux.h"
	["radhandl.h"]="rad_handle.h"
	["radmtra1.h"]="rad_material_aux.h"
	["radopnam.h"]="rad_operation_names.h"
	["radstlon.h"]="rad_string_long.h"
	["radtrans.h"]="rad_transform_def.h"
	["radyield.h"]="rad_yield.h"
)

for old in "${!files[@]}"; do
	new="${files[$old]}"
	count=$(grep -r "$old" src/ --include="*.cpp" --include="*.h" 2>/dev/null | wc -l)
	printf "%-20s -> %-35s : %3d references\n" "$old" "$new" "$count"
done | sort

echo ""
echo "Total files to rename: ${#files[@]}"
