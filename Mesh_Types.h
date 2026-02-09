#ifndef MESH_TYPES_H
#define MESH_TYPES_H

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace MeshTypes {
#define HASH_EDGE(x, y) std::hash<std::string>{}(std::to_string(x) + " " + std::to_string(y))

#define HASH_SORTED_EDGE(x, y) (x > y ? HASH_EDGE(y, x) : HASH_EDGE(x, y))

#define HASH_VERTEX(x, y, z) std::hash<std::string>{}( \
    std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z))

#define HASH_FACE(face) ([&]() { \
    std::vector<size_t> smallest3(3); \
    std::partial_sort_copy((face).begin(), (face).end(), smallest3.begin(), smallest3.end()); \
    return HASH_VERTEX(smallest3[0], smallest3[1], smallest3[2]); \
})()

    struct Vec2 {
        double x, y;

        Vec2() : x(0.0), y(0.0) {
        }
        Vec2(const double x_, const double y_) : x(x_), y(y_) {
        }

        Vec2 operator-(const Vec2& other) const {
            return Vec2(x - other.x, y - other.y);
        }

        double cross(const Vec2& other) const {
            return x * other.y - y * other.x;
        }
    };

    struct Vec3 {
        double x, y, z;

        Vec3() : x(0.0), y(0.0), z(0.0) {
        }
        Vec3(const double x_, const double y_, const double z_) : x(x_), y(y_), z(z_) {
        }

        Vec3 operator+(const Vec3& other) const {
            return Vec3(x + other.x, y + other.y, z + other.z);
        }

        Vec3 operator-(const Vec3& other) const {
            return Vec3(x - other.x, y - other.y, z - other.z);
        }

        Vec3 operator*(const double scalar) const {
            return Vec3(x * scalar, y * scalar, z * scalar);
        }

        Vec3 operator/(const double scalar) const {
            if (scalar == 0) {
                return Vec3();
            }
            else {
                return Vec3(x / scalar, y / scalar, z / scalar);
            }
        }

        Vec3& operator+=(const Vec3& other) {
            x += other.x;
            y += other.y;
            z += other.z;
            return *this;
        }

        Vec3& operator-=(const Vec3& other) {
            x -= other.x;
            y -= other.y;
            z -= other.z;
            return *this;
        }

        Vec3& operator*=(const double scalar) {
            x *= scalar;
            y *= scalar;
            z *= scalar;
            return *this;
        }

        Vec3& operator/=(const double scalar) {
            if (scalar == 0) {
                x = 0;
                y = 0;
                z = 0;
            }
            else {
                x /= scalar;
                y /= scalar;
                z /= scalar;
            }
            return *this;
        }

        double dot(const Vec3& other) const {
            return x * other.x + y * other.y + z * other.z;
        }

        Vec3 cross(const Vec3& other) const {
            return Vec3(y * other.z - z * other.y,
                z * other.x - x * other.z,
                x * other.y - y * other.x);
        }

        double norm() const {
            return std::sqrt(x * x + y * y + z * z);
        }

        double manhattan_norm() const {
            return std::fabs(x) + std::fabs(y) + std::fabs(z);
        }

        Vec3 normalized() const {
            double norm = this->norm();
            return (norm == 0 ? Vec3() : *this / norm);
        }

        void normalize() {
            double norm = this->norm();
            if (norm == 0) {
                *this = Vec3();
            }
            else {
                *this /= norm;
            }
        }
    };

    inline Vec3 operator-(const Vec3& vec) {
        return Vec3(-vec.x, -vec.y, -vec.z);
    }

    inline Vec3 operator*(double scalar, const Vec3& vec) {
        return Vec3(vec.x * scalar, vec.y * scalar, vec.z * scalar);
    }

    struct VecNI {
        std::vector<size_t> indices;

        VecNI() = default;
        VecNI(std::initializer_list<size_t> list) : indices(list) {
        }
        // explicit constructor to prevent implicit conversions
        explicit VecNI(const std::vector<size_t>& vs) : indices(vs) {
        }

        size_t operator[](size_t index) const {
            return indices[index];
        }
        size_t& operator[](size_t index) {
            return indices[index];
        }

        size_t size() const {
            return indices.size();
        }

        void push_back(const size_t& value) {
            indices.push_back(value);
        }
        void push_back(size_t&& value) {
            indices.push_back(std::move(value));
        }

        auto begin() {
            return indices.begin();
        }
        auto end() {
            return indices.end();
        }
        auto begin() const {
            return indices.begin();
        }
        auto end() const {
            return indices.end();
        }
    };

    struct Mesh2D {
        std::vector<Vec2> vertices;
        std::vector<VecNI> faces;
    };

    struct Mesh3D {
        std::vector<Vec3> vertices;
        std::vector<VecNI> faces;
    };

    struct SharpFeature {
        std::vector<size_t> vertices = {};
        std::vector<std::pair<size_t, size_t>> edges = {};
    };
}

#endif