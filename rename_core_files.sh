#!/bin/bash

# Phase 2: Rename core implementation files to new naming convention
# This script:
# 1. Renames core .cpp and .h files using git mv
# 2. Updates all #include directives
# 3. Updates CMakeLists.txt if needed
# 4. Shows what was changed

set -e  # Exit on error

echo "=== Phase 2: Core Implementation File Renaming ===="
echo ""

# Define mappings (old -> new)
# Note: radinter.cpp is skipped as it's already renamed to rad_interaction.cpp
declare -A files=(
	["radcast.cpp"]="rad_type_cast.cpp"
	["radcast.h"]="rad_type_cast.h"
	["radg3d.cpp"]="rad_geometry_3d.cpp"
	["radg3d.h"]="rad_geometry_3d.h"
	["radgroup.cpp"]="rad_group.cpp"
	["radgroup.h"]="rad_group.h"
	["radmater.cpp"]="rad_material_def.cpp"
	["radmater.h"]="rad_material_def.h"
	["radmaterial.cpp"]="rad_material_impl.cpp"
	["radsend.cpp"]="rad_serialization.cpp"
	["radsend.h"]="rad_serialization.h"
	["radtransform.cpp"]="rad_transform_impl.cpp"
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
echo "2. Build: cmake --build build --config Release"
echo "3. Test: pytest tests/"
echo "4. Commit: git commit -m 'Phase 2: Rename core implementation files'"
