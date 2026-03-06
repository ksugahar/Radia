// Test to check sizeof various Radia structures
// Compile and run to identify which structure is ~150 bytes

#include <iostream>
#include "../src/core/rad_geometry_3d.h"
#include "../src/core/radvect.h"

int main() {
    std::cout << "=== Sizeof Radia Structures ===" << std::endl;
    std::cout << std::endl;

    std::cout << "TVector3d: " << sizeof(TVector3d) << " bytes" << std::endl;
    std::cout << "radTField: " << sizeof(radTField) << " bytes" << std::endl;
    std::cout << "radTFieldKey: " << sizeof(radTFieldKey) << " bytes" << std::endl;
    std::cout << "radTCompCriterium: " << sizeof(radTCompCriterium) << " bytes" << std::endl;

    std::cout << std::endl;
    std::cout << "radTField members:" << std::endl;
    std::cout << "  11 x TVector3d: " << 11 * sizeof(TVector3d) << " bytes" << std::endl;
    std::cout << "  2 x double: " << 2 * sizeof(double) << " bytes" << std::endl;
    std::cout << "  1 x int: " << sizeof(int) << " bytes" << std::endl;
    std::cout << "  1 x short: " << sizeof(short) << " bytes" << std::endl;
    std::cout << "  1 x radTFieldKey: " << sizeof(radTFieldKey) << " bytes" << std::endl;
    std::cout << "  1 x radTCompCriterium: " << sizeof(radTCompCriterium) << " bytes" << std::endl;

    return 0;
}
