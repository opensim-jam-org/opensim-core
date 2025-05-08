#ifndef WRAPMESH_H
#define WRAPMESH_H

#include "WrapObject.h"
//#include <OpenSim/Simulation/Wrap/WrapObject.h>
#include <simmath/internal/ContactGeometry.h>
#include <SimTKcommon.h>
#include <vector>

using SimTK::Vec3;

namespace OpenSim {

    class PathWrap;
    class WrapResult;

    //=============================================================================
    //=============================================================================
    /**
     * A class implementing a triangular mesh for path wrapping.
     *
     * @author Colin Smith
     * 
     */
    class OSIMSIMULATION_API WrapMesh : public WrapObject {
        OpenSim_DECLARE_CONCRETE_OBJECT(WrapMesh, WrapObject);

    public:
        // Property declaration
        OpenSim_DECLARE_PROPERTY(mesh_file, std::string, "Path to the mesh file.");

        // Constructors
        //WrapMesh();
        WrapMesh(const std::string& meshFilePath);
        //virtual ~WrapMesh();

        SimTK::PolygonalMesh loadPolyMesh(const std::string& meshFilePath);

        const char* getWrapTypeName() const override
        {
            static const char* wrapTypeName = "mesh";
            return wrapTypeName;
        }

        std::string getDimensionsString() const override;

        // Override the lineWrap method
        int wrapLine(const SimTK::State& s, SimTK::Vec3& aPoint1, SimTK::Vec3& aPoint2,
            const PathWrap& aPathWrap, WrapResult& aWrapResult, bool& aFlag) const override;

        // Load the mesh
        void loadMesh(const std::string& meshFilePath);

        

    protected:
        void constructProperties();

        /// Implement generateDecorations to draw geometry in visualizer
        void generateDecorations(bool fixed, const ModelDisplayHints& hints, const SimTK::State& state,
            SimTK::Array_<SimTK::DecorativeGeometry>& appendToThis) const override;

    private:
        SimTK::ContactGeometry::TriangleMesh triangleMesh;

        // Helper methods
        void initializeKnots(const SimTK::Vec3& pointA, const SimTK::Vec3& pointB,
            SimTK::Vector_<Vec3>& knotPositions) const;
        double computePathLength(const SimTK::Vec3& aPoint1, const SimTK::Vec3& aPoint2, const SimTK::Vector_<Vec3>& knots) const;

        SimTK::Vector_<Vec3> solveBlockTridiagonal(
            SimTK::Vector_<SimTK::Mat33>& B_3d_mat, SimTK::Vector_<SimTK::Mat33>& A_3d_mat,
            SimTK::Vector_<SimTK::Mat33>& C_3d_mat, SimTK::Vector_<Vec3>& D_3d_mat) const;

        SimTK::Matrix computeStiffnessMatrix(const SimTK::Vector_<Vec3>& knots,
            double stiffness, double contactStiffness) const;
        SimTK::Vector_<Vec3> computeForces(const SimTK::Vector_<Vec3>& knots,
            double stiffness, double contactStiffness) const;
        SimTK::Vector_<Vec3> solveBlockTridiagonal(const SimTK::Matrix& K,
            const SimTK::Vector_<Vec3>& forces) const;
    };

} // namespace OpenSim

#endif // WRAPMESH_H