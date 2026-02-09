#ifndef MESH_SHAPE_DIAMETER_FUNCTION_H
#define MESH_SHAPE_DIAMETER_FUNCTION_H

#define M_PI 3.14159265358979323846

#include "Mesh_Types.h"

namespace MeshSDF {

    class ShapeDiameterFunction {
    public:
        ShapeDiameterFunction(const MeshTypes::Mesh3D& mesh) : mesh_(mesh) {}

        std::vector<double> ComputePerFaceSDF() {
            std::vector<double> sdf_values(mesh_.faces.size(), 0);

//#pragma omp parallel for
            for (long long face_idx = 0; face_idx < mesh_.faces.size(); ++face_idx) {
                const auto& face = mesh_.faces[face_idx];
                const MeshTypes::Vec3 v0 = mesh_.vertices[face[0]];
                const MeshTypes::Vec3 v1 = mesh_.vertices[face[1]];
                const MeshTypes::Vec3 v2 = mesh_.vertices[face[2]];
                const MeshTypes::Vec3 face_center = (v0 + v1 + v2) / 3.0;
                const MeshTypes::Vec3 normal = (v1 - v0).cross(v2 - v0).normalized();

                const auto dirs = GenerateConeDirections(normal, 10, 45.0);
                for (const auto& dir : dirs) {
                    const double t = IntersectRayWithMesh(face_center, dir, face_idx);
                    const double axis_align_t = t / 2.0 / dir.manhattan_norm();
                    sdf_values[face_idx] += axis_align_t;
                }
                sdf_values[face_idx] /= dirs.size();
            }
            return sdf_values;
        }

    private:
        const MeshTypes::Mesh3D& mesh_;

        std::vector<MeshTypes::Vec3> GenerateConeDirections(
            const MeshTypes::Vec3& cone_axis, size_t num_samples, double cone_angle_deg)
        {
            std::vector<MeshTypes::Vec3> directions;
            directions.reserve(num_samples);

            // deg to rad
            const double cone_angle_rad = cone_angle_deg * M_PI / 180.0;
            // cosine range: [cos(cone_angle), 1]
            const double cosThetaMax = std::cos(cone_angle_rad);
            // golden angle
            const double goldenAngle = M_PI * (3.0 - std::sqrt(5.0));

            // construct orthogonal basis
            MeshTypes::Vec3 u, v;
            if (std::fabs(cone_axis.z) < 0.999) {
                u = MeshTypes::Vec3(0, 0, 1).cross(cone_axis).normalized();
            }
            else {
                u = MeshTypes::Vec3(0, 1, 0).cross(cone_axis).normalized();
            }
            v = cone_axis.cross(u);

            // for each sampling
            for (size_t i = 0; i < num_samples; ++i) {
                // sample t uniformly in [0, 1]
                const double t = i / 1.0 / (num_samples - 1);
                // sample cosTheta uniformly in [cosThetaMax, 1]
                const double cosTheta = cosThetaMax + (1.0 - cosThetaMax) * t;
                const double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
                // azimuthal angle sampled with golden angle spiral
                const double phi = std::fmod(goldenAngle * i, 2.0 * M_PI);
                // orientation in the local coordinate system (centered on the z-axis)
                const MeshTypes::Vec3 localDir(
                    sinTheta * std::cos(phi),
                    sinTheta * std::sin(phi),
                    cosTheta
                );
                // convert local coordinates to global coordinates
                const MeshTypes::Vec3 dir =
                    u * localDir.x
                    + v * localDir.y
                    + cone_axis * localDir.z;
                directions.push_back(dir);
            }

            return directions;
        }

        double IntersectRayWithMesh(const MeshTypes::Vec3& origin,
            const MeshTypes::Vec3& direction,
            const size_t current_face_idx)
        {
            double min_t = std::numeric_limits<double>::max();
            for (size_t face_idx = 0; face_idx < mesh_.faces.size(); ++face_idx) {
                if (face_idx == current_face_idx) {
                    continue;
                }

                const auto& face = mesh_.faces[face_idx];
                const MeshTypes::Vec3 v0 = mesh_.vertices[face[0]];
                const MeshTypes::Vec3 v1 = mesh_.vertices[face[1]];
                const MeshTypes::Vec3 v2 = mesh_.vertices[face[2]];
                double t;

                if (RayTriangleIntersect(origin, direction, v0, v1, v2, t) && t < min_t) {
                    min_t = t;
                }
            }
            return min_t;
        }

        bool RayTriangleIntersect(const MeshTypes::Vec3& orig, const MeshTypes::Vec3& dir,
            const MeshTypes::Vec3& v0, const MeshTypes::Vec3& v1,
            const MeshTypes::Vec3& v2, double& t)
        {
            const MeshTypes::Vec3 edge1 = v1 - v0;
            const MeshTypes::Vec3 edge2 = v2 - v0;
            const MeshTypes::Vec3 pvec = dir.cross(edge2);
            const double det = edge1.dot(pvec);

            if (det == 0) {
                return false;
            }
            const double inv_det = 1.0 / det;

            const MeshTypes::Vec3 tvec = orig - v0;
            const double u = tvec.dot(pvec) * inv_det;
            if (u < 0.0 || u > 1.0) {
                return false;
            }

            const MeshTypes::Vec3 qvec = tvec.cross(edge1);
            const double v = dir.dot(qvec) * inv_det;
            if (v < 0.0 || u + v > 1.0) {
                return false;
            }

            t = std::fabs(edge2.dot(qvec) * inv_det);
            return true;
        }
    };

}

#endif