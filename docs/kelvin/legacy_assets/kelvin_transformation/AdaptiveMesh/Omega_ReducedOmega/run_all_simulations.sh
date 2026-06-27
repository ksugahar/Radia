#!/bin/bash
# Bash script to run all 18 simulation files (Sphere_3D + Cylinder_3D)
# Each model: order=2,3,4 x 3 methods (Refine_all_elements, Refine_with_zz_estimator, metric_based)

BASE_DIR="S:/NGSolve/NGSolve/2025_12_14_Kelvin変換/AdaptiveMesh/Omega_ReducedOmega"

echo "============================================================"
echo "Starting All Simulations (Sphere_3D + Cylinder_3D)"
echo "Start Time: $(date)"
echo "============================================================"
echo ""

# Array of simulation directories and files
declare -a SIMULATIONS=(
    # Sphere_3D - order=2
    "Sphere_3D/order=2/Refine_all_elements|Sphere_3D_adaptive_with_Kelvin.py"
    "Sphere_3D/order=2/Refine_with_zz_estimator|Sphere_3D_adaptive_with_Kelvin.py"
    "Sphere_3D/order=2/metric_based|Sphere_3D_metric_with_Kelvin.py"
    # Sphere_3D - order=3
    "Sphere_3D/order=3/Refine_all_elements|Sphere_3D_adaptive_with_Kelvin.py"
    "Sphere_3D/order=3/Refine_with_zz_estimator|Sphere_3D_adaptive_with_Kelvin.py"
    "Sphere_3D/order=3/metric_based|Sphere_3D_metric_with_Kelvin.py"
    # Sphere_3D - order=4
    "Sphere_3D/order=4/Refine_all_elements|Sphere_3D_adaptive_with_Kelvin.py"
    "Sphere_3D/order=4/Refine_with_zz_estimator|Sphere_3D_adaptive_with_Kelvin.py"
    "Sphere_3D/order=4/metric_based|Sphere_3D_metric_with_Kelvin.py"
    # Cylinder_3D - order=2
    "Cylinder_3D/order=2/Refine_all_elements|Cylinder_3D_adaptive_with_Kelvin.py"
    "Cylinder_3D/order=2/Refine_with_zz_estimator|Cylinder_3D_adaptive_with_Kelvin.py"
    "Cylinder_3D/order=2/metric_based|Cylinder_3D_metric_with_Kelvin.py"
    # Cylinder_3D - order=3
    "Cylinder_3D/order=3/Refine_all_elements|Cylinder_3D_adaptive_with_Kelvin.py"
    "Cylinder_3D/order=3/Refine_with_zz_estimator|Cylinder_3D_adaptive_with_Kelvin.py"
    "Cylinder_3D/order=3/metric_based|Cylinder_3D_metric_with_Kelvin.py"
    # Cylinder_3D - order=4
    "Cylinder_3D/order=4/Refine_all_elements|Cylinder_3D_adaptive_with_Kelvin.py"
    "Cylinder_3D/order=4/Refine_with_zz_estimator|Cylinder_3D_adaptive_with_Kelvin.py"
    "Cylinder_3D/order=4/metric_based|Cylinder_3D_metric_with_Kelvin.py"
)

TOTAL=${#SIMULATIONS[@]}
COMPLETED=0
FAILED=0
FAILED_NAMES=""

echo "Total simulations to run: $TOTAL"
echo ""

for SIM in "${SIMULATIONS[@]}"; do
    DIR_PART="${SIM%%|*}"
    FILE_PART="${SIM##*|}"

    COMPLETED=$((COMPLETED + 1))
    FULL_DIR="$BASE_DIR/$DIR_PART"

    echo "------------------------------------------------------------"
    echo "[$COMPLETED/$TOTAL] $DIR_PART"
    echo "File: $FILE_PART"
    echo "Started: $(date '+%H:%M:%S')"
    echo ""

    # Change to directory and run
    cd "$FULL_DIR"
    START_TIME=$(date +%s)

    python "$FILE_PART"
    EXIT_CODE=$?

    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))

    if [ $EXIT_CODE -eq 0 ]; then
        echo ""
        echo "  SUCCESS"
    else
        echo ""
        echo "  FAILED (exit code: $EXIT_CODE)"
        FAILED=$((FAILED + 1))
        FAILED_NAMES="$FAILED_NAMES  - $DIR_PART\n"
    fi

    echo "  Duration: ${DURATION}s"
    REMAINING=$((TOTAL - COMPLETED))
    echo "  Progress: $COMPLETED/$TOTAL completed, $REMAINING remaining"
    echo ""
done

echo "============================================================"
echo "All Simulations Completed"
echo "============================================================"
echo ""
echo "Summary:"
echo "  Total simulations: $TOTAL"
echo "  Successful: $((TOTAL - FAILED))"
echo "  Failed: $FAILED"

if [ $FAILED -gt 0 ]; then
    echo ""
    echo "Failed simulations:"
    echo -e "$FAILED_NAMES"
fi

echo ""
echo "End Time: $(date)"
echo ""
