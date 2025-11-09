#!/bin/bash

# Phase 1: Rename header files to new naming convention
# This script:
# 1. Renames header files using git mv
# 2. Updates all #include directives
# 3. Shows what was changed

set -e  # Exit on error

echo "=== Phase 1: Header File Renaming ===="
echo ""

# Define mappings (old -> new)
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

echo "Step 1: Renaming files..."
echo ""

for old in "${!files[@]}"; do
	new="${files[$old]}"

	if [ -f "src/core/$old" ]; then
		echo "  Renaming: $old -> $new"
		git mv "src/core/$old" "src/core/$new"
	else
		echo "  WARNING: File not found: src/core/$old"
	fi
done

echo ""
echo "Step 2: Updating #include directives..."
echo ""

# Update all include statements in source files
for old in "${!files[@]}"; do
	new="${files[$old]}"

	# Count occurrences
	count=$(grep -r "#include \"$old\"" src/ --include="*.cpp" --include="*.h" 2>/dev/null | wc -l)

	if [ $count -gt 0 ]; then
		echo "  Updating $count references: $old -> $new"

		# Use find + sed for cross-platform compatibility
		find src/ \( -name "*.cpp" -o -name "*.h" \) -type f -exec sed -i "s|#include \"$old\"|#include \"$new\"|g" {} +
	fi
done

echo ""
echo "=== Renaming Complete ==="
echo ""
echo "Next steps:"
echo "1. Review changes: git status"
echo "2. Build: cmake --build build --config Release"
echo "3. Test: pytest tests/"
echo "4. Commit: git commit -m 'Rename header files to new naming convention'"
