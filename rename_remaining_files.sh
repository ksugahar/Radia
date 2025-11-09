#!/bin/bash

# Phase 3: Rename remaining legacy files
# This script:
# 1. Renames final .cpp and .h files using git mv
# 2. Updates all #include directives
# 3. Updates CMakeLists.txt if needed

set -e  # Exit on error

echo "=== Phase 3: Final File Renaming ===="
echo ""

# Define mappings (old -> new)
declare -A files=(
	["radhmat.cpp"]="rad_hmatrix.cpp"
	["radhmat.h"]="rad_hmatrix.h"
	["radinter.cpp"]="rad_c_interface.cpp"
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
echo "2. Update CMakeLists.txt for radhmat.cpp -> rad_hmatrix.cpp and radinter.cpp -> rad_c_interface.cpp"
echo "3. Build: cmake --build build --config Release"
echo "4. Test: pytest tests/"
echo "5. Commit: git commit -m 'Phase 3: Rename final legacy files'"
