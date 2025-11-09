#!/bin/bash

# Final batch: Rename remaining files without rad_ prefix
# This completes the naming convention unification

set -e  # Exit on error

echo "=== Final Batch: Remaining File Renaming ===="
echo ""

# Define mappings (old -> new)
declare -A files=(
	["radfield_compute.cpp"]="rad_field_compute.cpp"
	["radgeom_types.h"]="rad_geom_types.h"
	["radhmat_field.cpp"]="rad_hmat_field.cpp"
	["radhmat_update.cpp"]="rad_hmat_update.cpp"
	["radintrc_hmat.cpp"]="rad_intrc_hmat.cpp"
	["radintrc_hmat.h"]="rad_intrc_hmat.h"
	["radobj_creation.cpp"]="rad_obj_creation.cpp"
	["radobj_subdivision.cpp"]="rad_obj_subdivision.cpp"
	["radpoly_analytical.cpp"]="rad_poly_analytical.cpp"
	["radpoly_analytical.h"]="rad_poly_analytical.h"
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

	# Only process header files for include updates
	if [[ "$old" == *.h ]]; then
		# Count occurrences
		count=$(grep -r "#include \"$old\"" src/ --include="*.cpp" --include="*.h" 2>/dev/null | wc -l)

		if [ $count -gt 0 ]; then
			echo "  Updating $count references: $old -> $new"

			# Use find + sed for cross-platform compatibility
			find src/ \( -name "*.cpp" -o -name "*.h" \) -type f -exec sed -i "s|#include \"$old\"|#include \"$new\"|g" {} +
		fi
	fi
done

echo ""
echo "=== Renaming Complete ==="
echo ""
echo "Next steps:"
echo "1. Review changes: git status"
echo "2. Update CMakeLists.txt for all renamed .cpp files"
echo "3. Build: cmake --build build --config Release"
echo "4. Test: pytest tests/"
echo "5. Commit: git commit -m 'Final batch: Complete rad_ prefix unification'"
