#pragma once

#include "Rational.hpp"
#include "Utils.hpp"

#include <Eigen/Core>

#include <array>
#include <vector>

namespace eccd
{

class TriPtF
{
public:
    TriPtF(const Vector3d &pts,
           const Vector3d &v1s, const Vector3d &v2s, const Vector3d &v3s,
           const Vector3d &pte,
           const Vector3d &v1e, const Vector3d &v2e, const Vector3d &v3e);

    Vector3r operator()(double u, double v, double t) const;

    std::array<Vector3r, 4> corners(int i) const;
    inline int n_patches() const { return 3; }

    std::vector<std::array<Vector3r, 3>> top_bottom_faces() const;

private:
    const Vector3d &pts_;
    const Vector3d &v1s_, &v2s_, &v3s_;

    const Vector3d &pte_;
    const Vector3d &v1e_, &v2e_, &v3e_;

    Vector3r ptsr_;
    Vector3r v1sr_, v2sr_, v3sr_;

    Vector3r pter_;
    Vector3r v1er_, v2er_, v3er_;
};

class EdgeEdgeF
{
public:
    //edge [a0, a1] and [b0, b1]
    EdgeEdgeF(const Vector3d &a0s, const Vector3d &a1s,
              const Vector3d &b0s, const Vector3d &b1s,
              const Vector3d &a0e, const Vector3d &a1e,
              const Vector3d &b0e, const Vector3d &b1e);

    Vector3r operator()(double ta, double tb, double t) const;

    std::array<Vector3r, 4> corners(int i) const;
    inline int n_patches() const { return 6; }

    std::vector<std::array<Vector3r, 3>> top_bottom_faces() const { return {}; }

private:
    const Vector3d &a0s_, &a1s_;
    const Vector3d &b0s_, &b1s_;

    const Vector3d &a0e_, &a1e_;
    const Vector3d &b0e_, &b1e_;

    Vector3r a0rs_, a1rs_;
    Vector3r b0rs_, b1rs_;

    Vector3r a0re_, a1re_;
    Vector3r b0re_, b1re_;
};

Rational phi(const Vector3r x, const std::array<Vector3r, 4> &corners);

bool is_origin_in_tet(const std::array<Vector3r, 4> &corners, const std::array<std::array<int, 3>, 4> &tet_faces);

int ray_flat_patch(const std::array<Vector3r, 4> &corners, const Vector3d &dir);

template <typename FuncF>
int ray_patch(const FuncF &func, int patch, const Vector3d &dir);

template <typename FuncF>
int ccd(const FuncF &func, const Vector3d &dir);

template <typename FuncF>
bool retrial_ccd(const FuncF &func);

bool vertexFaceCCD(const Vector3d &pts,
                   const Vector3d &v1s, const Vector3d &v2s, const Vector3d &v3s,
                   const Vector3d &pte,
                   const Vector3d &v1e, const Vector3d &v2e, const Vector3d &v3e);

bool edgeEdgeCCD(const Vector3d &a0s, const Vector3d &a1s,
                 const Vector3d &b0s, const Vector3d &b1s,
                 const Vector3d &a0e, const Vector3d &a1e,
                 const Vector3d &b0e, const Vector3d &b1e);

} // namespace eccd
