// Links the real production factory, never the synthetic scenario backend.
#include "ih_coupled_test_fixture.hpp"
#include "../../native/ih_coupled_backend.hpp"
#include <iostream>
#include <string>

namespace {
void require(bool value,const char* message) {
    if(!value)throw std::runtime_error(message);
}
void unavailable(const radia_ih::CoupledConfig& config) {
    using namespace radia_ih;
    (void)validate_config(config); // rejection must not be a malformed fixture
    bool factory_rejected=false;
    try { (void)detail::make_coupled_physics_backend(config); }
    catch(const DomainError& e) {
        const std::string message=e.what();
        factory_rejected=e.code()==ErrorCode::unsupported_physics &&
            message.find("independent_cpu")!=std::string::npos &&
            message.find("COUPLED_DOMAIN_INTERFACE.md")!=std::string::npos;
    }
    require(factory_rejected,"Production factory must raise unsupported_physics");
    bool constructor_rejected=false;
    try { CoupledDomain domain(config); }
    catch(const DomainError& e) {constructor_rejected=e.code()==ErrorCode::unsupported_physics;}
    require(constructor_rejected,"Unavailable domain must not publish an initialized owner");
}
}
int main() {
    try {
        auto config=ih_contract_fixture::fixture();
        unavailable(config);
        // A failed construction must not leave global state that changes the
        // next request. The main config accepts only explicit independent_cpu.
        unavailable(ih_contract_fixture::fixture());
        std::cout<<"IH production unavailable link/lifecycle contract PASS\n";
        return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
