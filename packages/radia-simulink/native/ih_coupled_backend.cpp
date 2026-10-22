// Production link boundary. No canonical physical provider is qualified yet.
#include "ih_coupled_backend.hpp"

namespace radia_ih::detail {
std::unique_ptr<PhysicsBackend> make_coupled_physics_backend(const CoupledConfig& config) {
    // Direct callers receive the same structural gate as CoupledDomain callers.
    (void)validate_config(config);
    throw DomainError(ErrorCode::unsupported_physics,
        "IH backend independent_cpu is not implemented; no fallback. See packages/radia-simulink/docs/COUPLED_DOMAIN_INTERFACE.md");
}
} // namespace radia_ih::detail
