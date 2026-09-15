#pragma once
#include <periodic.hpp>
#include <h1hofespace.hpp>
#include <map>
#include <vector>
#include <algorithm>

namespace radia {
// A single, bijective interface. Reconstruct face DOFs from all corner IDs,
// retaining Periodic's vertex-oriented finite elements and edge mappings.
class FaceConsistentPeriodicH1 : public ngcomp::PeriodicFESpace {
public:
    explicit FaceConsistentPeriodicH1(std::shared_ptr<ngcomp::FESpace> base)
        : PeriodicFESpace(base, base->GetFlags(),
                          std::make_shared<ngcore::Array<int>>()) {}

    void Update() override {
        if (!dynamic_cast<ngcomp::H1HighOrderFESpace*>(space.get()))
            throw ngcore::Exception("Face-consistent interface requires scalar H1");
        if (ma->GetDimension() != 3 || ma->GetNPeriodicIdentifications() != 1)
            throw ngcore::Exception("Exactly one 3D periodic interface is required");
        PeriodicFESpace::Update();
        std::map<std::vector<int>, size_t> faces;
        ngcore::Array<int> dofs, master;
        for (size_t face = 0; face < ma->GetNFaces(); ++face) {
            auto vertices = ma->GetFacePNums(face);
            std::vector<int> key;
            for (auto v : vertices) key.push_back(v);
            std::sort(key.begin(), key.end());
            if (!faces.emplace(key, face).second)
                throw ngcore::Exception("Duplicate face corner identity");
            space->GetDofNrs(ngcomp::NodeId(ngfem::NT_FACE, face), dofs);
            for (auto d : dofs) if (d >= 0) dofmap[d] = d;
        }
        size_t paired = 0;
        for (size_t face = 0; face < ma->GetNFaces(); ++face) {
            auto vertices = ma->GetFacePNums(face);
            std::vector<int> key;
            bool all_slave = true;
            for (auto v : vertices) {
                key.push_back(vertex_map[v]);
                all_slave = all_slave && (vertex_map[v] != v);
            }
            if (!all_slave) continue;
            std::sort(key.begin(), key.end());
            auto match = faces.find(key);
            if (match == faces.end() || match->second == face)
                throw ngcore::Exception("Interface face has no complete corner match");
            space->GetDofNrs(ngcomp::NodeId(ngfem::NT_FACE, face), dofs);
            space->GetDofNrs(ngcomp::NodeId(ngfem::NT_FACE, match->second), master);
            if (dofs.Size() != master.Size())
                throw ngcore::Exception("Interface face orders differ");
            for (size_t i = 0; i < dofs.Size(); ++i)
                if (dofs[i] >= 0) dofmap[dofs[i]] = master[i];
            ++paired;
        }
        if (!paired) throw ngcore::Exception("No paired interface faces");
        for (size_t d = 0; d < dofmap.Size(); ++d)
            ctofdof[d] = dofmap[d] == d ? space->GetDofCouplingType(d)
                                      : ngcomp::UNUSED_DOF;
    }
};
}
