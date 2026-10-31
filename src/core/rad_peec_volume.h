// Finite rectangular-volume Neumann kernel, linear nonmagnetic MQS.
#pragma once
#include <array>
#include <vector>
namespace radia {
struct VolumeBox { std::array<double,3> lo, hi; int axis; int sign; };
class VolumeGreenRule {
    std::vector<double> u_, weight_;
public:
    explicit VolumeGreenRule(int order);
    double Inductance(const VolumeBox& a, const VolumeBox& b) const;
};
}
