#ifndef WRAPMESH_H
#define WRAPMESH_H

#include <OpenSim/Simulation/Wrap/WrapObject.h>
#include <OpenSim/Common/TriangularMesh.h>
#include <OpenSim/Common/STLFileAdapter.h>

namespace OpenSim {

class WrapMesh : public WrapObject {
    OpenSim_DECLARE_CONCRETE_OBJECT(WrapMesh, WrapObject);

public:
    // Constructors
    WrapMesh();
    WrapMesh(const std::string& stlFilePath);

    // Compute the wrapping path
    virtual void computePath(const SimTK::Vec3& pointA,
                             const SimTK::Vec3& pointB,
                             SimTK::Array_<SimTK::Vec3>& path) const override;

    // Load mesh from file
    void loadMesh(const std::string& stlFilePath);

protected:
    // Initialize properties
    void constructProperties();

private:
    TriangularMesh triangularMesh; // Custom mesh for wrapping computations.

    OpenSim_DECLARE_PROPERTY(mesh_file, std::string,
                             "Path to the .stl file representing the triangular mesh.");
};

} // namespace OpenSim

#endif // WRAPMESH_H