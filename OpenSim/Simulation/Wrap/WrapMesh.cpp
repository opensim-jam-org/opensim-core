#include "WrapMesh.h"
#include "PathWrap.h"
#include "WrapResult.h"
#include <OpenSim/Common/IO.h>
#include "OpenSim/Simulation/Model/Model.h"

//#include <SimTKcommon.h>
//#include <vector>

using SimTK::Vec3;

namespace OpenSim {

    // Default constructor
    /*WrapMesh::WrapMesh() : triangleMesh(SimTK::PolygonalMesh()) {
        constructProperties();
    }*/

    // Constructor with mesh file
    WrapMesh::WrapMesh(const std::string& meshFilePath) : triangleMesh(loadPolyMesh(meshFilePath)){
        constructProperties();
        set_mesh_file(meshFilePath);
        loadMesh(meshFilePath);
    }

    // Destructor.

   // WrapMesh::~WrapMesh()
    //{
    //}

    SimTK::PolygonalMesh WrapMesh::loadPolyMesh(const std::string& meshFilePath) {
        SimTK::PolygonalMesh mesh;
        mesh.loadFile(meshFilePath);
        return mesh;
    }

    // Define constructProperties
    void WrapMesh::constructProperties() {
        constructProperty_mesh_file("");
    }

    

    std::string WrapMesh::getDimensionsString() const
    {
        std::stringstream dimensions;
        dimensions << "meters";

        return dimensions.str();
    }


    // Load mesh from file
    void WrapMesh::loadMesh(const std::string& meshFilePath) {
        SimTK::PolygonalMesh mesh;
        mesh.loadFile(meshFilePath);
        triangleMesh = SimTK::ContactGeometry::TriangleMesh(mesh);
    }

    // Override lineWrap method
    int WrapMesh::wrapLine(const SimTK::State& s, SimTK::Vec3& aPoint1, SimTK::Vec3& aPoint2,
        const PathWrap& aPathWrap, WrapResult& aWrapResult, bool& aFlag) const {

            //noWrap,          // the path segment did not intersect the wrap object
            //insideRadius,    // one or both path points are inside the wrap object
            //wrapped,         // successful wrap, but may not be 'best' path
            //mandatoryWrap    // successful wrap that must be used (e.g., both tangent
                  // points are on the constrained side of the wrap object)

        aWrapResult.wrap_pts.setSize(0);
        aWrapResult.wrap_path_length = 0.0;

        // Check for noWrap
        SimTK::UnitVec3 dir_vec((aPoint2 - aPoint1).normalize());
        double distance;
        SimTK::UnitVec3 normal;
        if (!triangleMesh.intersectsRay(aPoint1, dir_vec, distance, normal)) {            
            return noWrap;
        }

        // Check for inside
        bool P1isInside = false;
        triangleMesh.findNearestPoint(aPoint1,P1isInside,normal);
        bool P2isInside = false;
        triangleMesh.findNearestPoint(aPoint2, P2isInside, normal);

       if (P1isInside || P2isInside) {
           return insideRadius;
       }

        //Find Wrap Path
        const double wrapStiffness = 1.0;
        const double contactStiffness = 10.0;
        const double tolerance = 1e-8;
        const int maxIterations = 20;
        const double wrapOffset = 1e-8;
        const int numKnots = 50; 
        const double damping = 0;// 10.0 / (numKnots * numKnots);

        // Build Wrap Stiffness Matrix
        SimTK::Mat33 K_wrap_mat(0.0);
        K_wrap_mat(0, 0) = 2 * wrapStiffness;
        K_wrap_mat(1, 1) = 2 * wrapStiffness;
        K_wrap_mat(2, 2) = 2 * wrapStiffness;

        SimTK::Vector_<Vec3> D_3d_mat(numKnots);

        SimTK::Vector_<SimTK::Mat33> B_3d_mat(numKnots);
        SimTK::Vector_<SimTK::Mat33> A_3d_mat(numKnots - 1);
        SimTK::Vector_<SimTK::Mat33> C_3d_mat(numKnots - 1);

        SimTK::Mat33 A_mat(1.0);
        SimTK::Mat33 C_mat(1.0);

        A_mat = A_mat * -wrapStiffness;
        C_mat = C_mat * -wrapStiffness;

        for (int k = 0; k < numKnots-1; ++k) {
            A_3d_mat(k) = A_mat;
            C_3d_mat(k) = C_mat;
        }


        SimTK::Vector_<Vec3> knotPositions(numKnots, Vec3(0));
        initializeKnots(aPoint1, aPoint2, knotPositions);

        SimTK::Vector_<Vec3> knotDisplacements(numKnots,Vec3(0));

        SimTK::Array_<bool> knotContact(numKnots,false);

        // Iteratively compute equilibrium knot positions
        double previousLength = computePathLength(aPoint1, aPoint2, knotPositions);
        double currentLength;
        
        for (int iter = 0; iter < maxIterations; ++iter) {
            //std::cout << "Iteration: " << iter << std::endl;
            for (int k = 0; k < numKnots; ++k) {

                // Wrap
                Vec3 forward_knot_pos;
                Vec3 backward_knot_pos;

                if (k == numKnots-1) { forward_knot_pos = aPoint2; }                    
                else{forward_knot_pos = knotPositions[k + 1];}                                   

                if (k == 0) {backward_knot_pos = aPoint1;}                    
                else {backward_knot_pos = knotPositions[k - 1];}

                Vec3 F_wrap = wrapStiffness *(forward_knot_pos - 2 * knotPositions[k] + backward_knot_pos);

                //std::cout << k << '\t' << (forward_knot_pos - 2 * knotPositions[k] + backward_knot_pos)<< std::endl;

                // Contact
                bool inside;
                SimTK::UnitVec3 surface_normal;
                double distance = 0.0;

                //std::cout << k << '\t' << knotPositions[k] << '\t' << knotPositions.size()  << std::endl;
                
                SimTK::Vec3 surface_point = triangleMesh.findNearestPoint(knotPositions[k], inside, surface_normal);
                
                SimTK::Vec3 knot_normal(0);

                if (inside) {
                    knot_normal = surface_point - knotPositions[k];
                    distance = -knot_normal.norm();
                    knot_normal = knot_normal.normalize();
                }
                  

                Vec3 F_contact(0);
                SimTK::Mat33 K_cnt_mat(0.0);

                if (inside) {
                    knotContact[k] = true;

                    F_contact = -contactStiffness * distance * knot_normal;
                   // K_cnt_mat = SimTK::Mat33(-10.0);
                   K_cnt_mat = (contactStiffness + damping) * SimTK::outer(knot_normal,knot_normal);                   
                   //std::cout << distance << "\t" << knot_normal << std::endl;
                }
                else {
                    knotContact[k] = false;
                }


                // Assemble Matrices
                Vec3 F_total = F_wrap + F_contact;

               // std::cout << k << '\t' << F_wrap << '\t' << F_contact << std::endl;

                D_3d_mat[k] = F_total;

                //std::cout << k << '\t' << K_wrap_mat << '\t' << K_cnt_mat << std::endl;

                SimTK::Mat33 B_mat = K_wrap_mat +K_cnt_mat;

                B_3d_mat[k] = B_mat;              


                //SimTK::Matrix blockTridiagonalK = computeStiffnessMatrix(knotPositions, stiffness, contactStiffness);
                //SimTK::Vector_<Vec3> forces = computeForces(knotPositions, stiffness, contactStiffness);
                //SimTK::Vector_<Vec3> displacements = solveBlockTridiagonal(blockTridiagonalK, forces);
            }

            //std::cout << B_3d_mat.size() << std::endl;
            //for (int k = 0; k < B_3d_mat.size(); ++k) {
            //    std::cout << k << '\t' << B_3d_mat[k] << std::endl;}

            knotDisplacements = solveBlockTridiagonal(B_3d_mat, A_3d_mat, C_3d_mat, D_3d_mat);

            

            for (int k = 0; k < numKnots; ++k) {
                //std::cout << k << '\t' << knotDisplacements[k] <<  std::endl;
                knotPositions[k] += knotDisplacements[k];
                //std::cout << k << '\t' << knotDisplacements[k] << std::endl;
            }


            // Compute Length and Convergence
            currentLength = computePathLength(aPoint1, aPoint2, knotPositions);
            //std::cout << "Length: " << currentLength << std::endl;
            double length_change = std::abs(currentLength - previousLength);
            if ( length_change < tolerance) {
                break;
            }
            
            previousLength = currentLength;
        }

        // Populate wrapResult with the computed wrapping path
        aWrapResult.wrap_path_length = currentLength;
        aWrapResult.wrap_pts.setSize(0);

        for (const auto& knot : knotPositions) {
            aWrapResult.wrap_pts.append(knot);
        }

        for (int k = 0; k < numKnots; ++k) {
            if (knotContact[k]) {
                aWrapResult.r1 = knotPositions[k];
                break;
            }
        }

        for (int k = numKnots-1; k >=0 ; --k) {
            if (knotContact[k]) {
                aWrapResult.r2 = knotPositions[k];
                break;
            }
        }
        //std::cout << "r1: " << aWrapResult.r1 << std::endl;
        //std::cout << "r2: " << aWrapResult.r2 << std::endl;

        //aWrapResult.wrap_path_length = computePathLength(knotPositions);

        return mandatoryWrap;
        //return wrapped;
    }

    // Helper methods
    void WrapMesh::initializeKnots(const SimTK::Vec3& pointA, const SimTK::Vec3& pointB,
        SimTK::Vector_<Vec3>& knotPositions) const {
        const int numKnots = knotPositions.size();
        SimTK::Vec3 direction = (pointB - pointA) / (numKnots+1);

        //std::cout << "A: " << pointA << std::endl;
        for (int i = 1; i <= numKnots; ++i) {
            knotPositions[i-1] = pointA + i * direction;
            //std::cout << i-1 << " " << knotPositions[i-1] << std::endl;
        }
        //std::cout << "B: " << pointB << std::endl;
    }

    double WrapMesh::computePathLength(const SimTK::Vec3& aPoint1, const SimTK::Vec3& aPoint2, const SimTK::Vector_<Vec3>& knots) const {
        
        double length = 0.0;
        for (int i = 0; i < knots.size() - 1; ++i) {
            length += (knots[i + 1] - knots[i]).norm();
        }

        length += (aPoint1 - knots[0]).norm();
        length += (knots[knots.size() - 1] - aPoint2).norm();
        return length;
    }


    SimTK::Vector_<SimTK::Vec3> WrapMesh::solveBlockTridiagonal(
        SimTK::Vector_<SimTK::Mat33>& A_3d_mat,
        SimTK::Vector_<SimTK::Mat33>& B_3d_mat,
        SimTK::Vector_<SimTK::Mat33>& C_3d_mat,
        SimTK::Vector_<SimTK::Vec3>& D_3d_mat) const
    {
        int n = A_3d_mat.size(); // Number of blocks
        //A : 3D array of diagonal blocks(m x m x n)
        // B : 3D array of upper diagonal blocks(m x m x n - 1)
        // C : 3D array of lower diagonal blocks(m x m x n - 1)
        // D : 2D array of right - hand side vectors(m x n)
        // 
        // Forward elimination
        SimTK::Mat33 L; // Multipliers
        SimTK::Vector_<SimTK::Vec3> x(n);   // Solution vector

        for (int i = 1; i < n-1; ++i) {
            L = C_3d_mat[i - 1] * A_3d_mat[i - 1].invert(); // Compute multiplier
            A_3d_mat[i] = A_3d_mat[i] - L * B_3d_mat[i - 1]; // Update diagonal block
            D_3d_mat[i] = D_3d_mat[i] - L * D_3d_mat[i - 1]; // Update RHS
        }

        // Back substitution
        x[n - 1] = A_3d_mat[n - 1].invert() * D_3d_mat[n - 1]; // Solve last block
        for (int i = n - 2; i >= 0; --i) {
            x[i] = A_3d_mat[i].invert() * (D_3d_mat[i] - B_3d_mat[i] * x[i + 1]);
        }

        return x;
    }

    SimTK::Matrix WrapMesh::computeStiffnessMatrix(const SimTK::Vector_<Vec3>& knots,
        double stiffness, double contactStiffness) const {
        int numKnots = knots.size();
        SimTK::Matrix K(3 * numKnots, 3 * numKnots, 0.0);

        for (size_t k = 1; k < numKnots - 1; ++k) {
            K[3 * k][3 * (k - 1)] = stiffness;
            K[3 * k][3 * k] = -2 * stiffness;
            K[3 * k][3 * (k + 1)] = stiffness;

            bool inside;
            SimTK::UnitVec3 surfaceNormal;
            SimTK::Vec3 surface_point = triangleMesh.findNearestPoint(knots[k], inside, surfaceNormal);
            if (inside) {                
                K[3 * k][3 * k] += contactStiffness;                
            }
        }
        return K;
    }

    SimTK::Vector_<Vec3> WrapMesh::computeForces(const SimTK::Vector_<Vec3>& knots,
        double stiffness, double contactStiffness) const {
        SimTK::Vector_<Vec3> forces(knots.size(), SimTK::Vec3(0));

        for (size_t k = 1; k < knots.size() - 1; ++k) {
            SimTK::Vec3 springForce = stiffness * (knots[k + 1] - 2 * knots[k] + knots[k - 1]);

            bool inside;
            SimTK::UnitVec3 surfaceNormal;
            SimTK::Vec3 contactForce(0);

            SimTK::Vec3 surface_point = triangleMesh.findNearestPoint(knots[k], inside, surfaceNormal);
            SimTK::Vec3 surf_vec = surface_point - knots[k];
            double distance = surf_vec.norm();

            if (inside) {
                contactForce = contactStiffness * distance * surfaceNormal;
            }
            

            forces[k] = springForce + contactForce;
        }
        return forces;
    }

    SimTK::Vector_<Vec3> WrapMesh::solveBlockTridiagonal(const SimTK::Matrix& K,
        const SimTK::Vector_<Vec3>& forces) const {
        SimTK::Vector_<Vec3> displacements(forces.size(), SimTK::Vec3(0));
        // Placeholder for block tridiagonal solve algorithm implementation
        return displacements;
    }

    void WrapMesh::generateDecorations(bool fixed, const ModelDisplayHints& hints, const SimTK::State& state,
        SimTK::Array_<SimTK::DecorativeGeometry>& appendToThis) const
    {

        Super::generateDecorations(fixed, hints, state, appendToThis);
        if (!fixed) return;

        if (hints.get_show_wrap_geometry()) {
            const Appearance& defaultAppearance = get_Appearance();
            if (!defaultAppearance.get_visible()) return;
            const Vec3 color = defaultAppearance.get_color();
            const auto X_BP = calcWrapGeometryTransformInBaseFrame();

            //std::cout << "FILE: " << get_mesh_file() << std::endl;
            /*
            appendToThis.push_back(
                SimTK::DecorativeMeshFile("lenhart2015-R-femur-bone_remesh.stl"));
            
            */
            appendToThis.push_back(
                SimTK::DecorativeMeshFile(get_mesh_file())
                .setTransform(X_BP)
                .setColor(color).setOpacity(defaultAppearance.get_opacity())
                .setRepresentation(defaultAppearance.get_representation())
                .setBodyId(getFrame().getMobilizedBodyIndex()));
                
        }


    }
} // namespace OpenSim