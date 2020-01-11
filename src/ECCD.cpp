#include "ECCD.hpp"

#include <cassert>

static const int max_trials = 8;

namespace eccd
{
TriPtF::TriPtF(const Vector3d &pts,
               const Vector3d &v1s, const Vector3d &v2s, const Vector3d &v3s,
               const Vector3d &pte,
               const Vector3d &v1e, const Vector3d &v2e, const Vector3d &v3e)
    : pts_(pts),
      v1s_(v1s), v2s_(v2s), v3s_(v3s),
      pte_(pte),
      v1e_(v1e), v2e_(v2e), v3e_(v3e)
{
    for (int i = 0; i < 3; ++i)
    {
        ptsr_[i] = Rational(pts_[i]);

        v1sr_[i] = Rational(v1s_[i]);
        v2sr_[i] = Rational(v2s_[i]);
        v3sr_[i] = Rational(v3s_[i]);

        pter_[i] = Rational(pte_[i]);

        v1er_[i] = Rational(v1e_[i]);
        v2er_[i] = Rational(v2e_[i]);
        v3er_[i] = Rational(v3e_[i]);
    }
}

Vector3r TriPtF::operator()(double u, double v, double t) const
{
    const Rational ur(u);
    const Rational vr(v);

    const Rational tr(t);
    return (1 - tr) * ptsr_ + tr * pter_ -
           ((1 - tr) * v1sr_ + t * v1er_) * (1 - ur - vr) -
           ((1 - tr) * v2sr_ + t * v2er_) * u -
           ((1 - tr) * v3sr_ + t * v3er_) * v;
}

std::array<Vector3r, 4> TriPtF::corners(int i) const
{
    //u=0
    if (i == 0)
    {
        return {{(*this)(0, 0, 0), (*this)(0, 1, 0), (*this)(0, 1, 1), (*this)(0, 0, 1)}};
    }
    //v = 0
    else if (i == 1)
    {
        return {{(*this)(0, 0, 0), (*this)(1, 0, 0), (*this)(1, 0, 1), (*this)(0, 0, 1)}};
    }
    //u + v= 1
    else if (i == 2)
    {
        return {{(*this)(1, 0, 0), (*this)(0, 1, 0), (*this)(0, 1, 1), (*this)(1, 0, 1)}};
    }
    else
    {
        assert(false);
    }
}

std::vector<std::array<Vector3r, 3>> TriPtF::top_bottom_faces() const
{
    return {
        {{(*this)(0, 0, 0), (*this)(0, 1, 0), (*this)(1, 1, 0)}},
        {{(*this)(0, 0, 1), (*this)(0, 1, 1), (*this)(1, 1, 1)}}};
}


///////////////////////////////////////////////////////
EdgeEdgeF::EdgeEdgeF(const Vector3d &a0s, const Vector3d &a1s,
                     const Vector3d &b0s, const Vector3d &b1s,
                     const Vector3d &a0e, const Vector3d &a1e,
                     const Vector3d &b0e, const Vector3d &b1e)
    : a0s_(a0s), a1s_(a1s),
      b0s_(b0s), b1s_(b1s),

      a0e_(a0e), a1e_(a1e),
      b0e_(b0e), b1e_(b1e)
{
    for (int i = 0; i < 3; ++i)
    {
        a0rs_[i] = Rational(a0s_[i]);
        a0re_[i] = Rational(a0e_[i]);

        a1rs_[i] = Rational(a1s_[i]);
        a1re_[i] = Rational(a1e_[i]);

        b0rs_[i] = Rational(b0s_[i]);
        b0re_[i] = Rational(b0e_[i]);

        b1rs_[i] = Rational(b1s_[i]);
        b1re_[i] = Rational(b1e_[i]);
    }
}

Vector3r EdgeEdgeF::operator()(double ta, double tb, double t) const
{
    const Rational tar(ta);
    const Rational tbr(tb);

    const Rational tr(t);
    return (1 - tar) * ((1 - tr) * a0rs_ + tr * a0re_) + tar * ((1 - tr) * a1rs_ + tr * a1re_) -
           ((1 - tbr) * ((1 - tr) * b0rs_ + tr * b0re_) + tbr * ((1 - tr) * b1rs_ + tr * b1re_));
}

std::array<Vector3r, 4> EdgeEdgeF::corners(int i) const
{
    //ta=0
    if (i == 0)
    {
        return {{(*this)(0, 0, 0), (*this)(0, 1, 0), (*this)(0, 1, 1), (*this)(0, 0, 1)}};
    }
    //tb = 0
    else if (i == 1)
    {
        return {{(*this)(0, 0, 0), (*this)(1, 0, 0), (*this)(1, 0, 1), (*this)(0, 0, 1)}};
    }
    //ta=1
    else if (i == 2)
    {
        return {{(*this)(1, 0, 0), (*this)(1, 1, 0), (*this)(1, 1, 1), (*this)(1, 0, 1)}};
    }
    //tb = 1
    else if (i == 3)
    {
        return {{(*this)(0, 1, 0), (*this)(1, 1, 0), (*this)(1, 1, 1), (*this)(0, 1, 1)}};
    }
    //t=0
    else if (i == 4)
    {
        return {{(*this)(0, 0, 0), (*this)(0, 1, 0), (*this)(1, 1, 0), (*this)(1, 0, 0)}};
    }
    //t=1
    else if (i == 5)
    {
        return {{(*this)(0, 0, 1), (*this)(0, 1, 1), (*this)(1, 1, 1), (*this)(1, 0, 1)}};
    }
    else
    {
        assert(false);
    }
}


namespace
{

Rational func_g(const Vector3r &x, const std::array<Vector3r, 4> &corners, const std::array<int, 3> &indices)
{
    const int p = indices[0];
    const int q = indices[1];
    const int r = indices[2];
    return (x - corners[p]).dot(cross(corners[q] - corners[p], corners[r] - corners[p]));
}

} // namespace

Rational phi(const Vector3r x, const std::array<Vector3r, 4> &corners)
{
    // const Rational ur(u);
    // const Rational vr(v);

    // const Rational tr(t);
    // (ur, vr, tr);

    const Rational g012 = func_g(x, corners, {{0, 1, 2}});
    const Rational g132 = func_g(x, corners, {{1, 3, 2}});
    const Rational g013 = func_g(x, corners, {{0, 1, 3}});
    const Rational g032 = func_g(x, corners, {{0, 3, 2}});

    const Rational h12 = g012 * g132;
    const Rational h03 = g013 * g032;

    const Rational phi = h12 - h03;

    return phi;
}

bool is_origin_in_tet(const std::array<Vector3r, 4> &corners, const std::array<std::array<int, 3>, 4> &tet_faces)
{
    Vector3d dir(0, 0, 1);
    int trials;

    for (trials = 0; trials < max_trials; ++trials)
    {
        int count = 0;
        for (int i = 0; i < 4; ++i)
        {
            const auto &f = tet_faces[i];
            const int res = origin_ray_triangle_inter(dir, corners[f[0]], corners[f[1]], corners[f[2]]);

            //bad luck
            if (res < 0)
            {
                count = -1;
                dir = Vector3d::Random();
                break;
            }

            if (res > 0)
                ++count;

            if (count > 1)
                break;
        }

        if (count == 1)
            return true;

        if (count >= 0)
            return false;
    }

    if (trials == max_trials)
    {
        std::cout << "All rays are on edges, increase trials" << std::endl;
        assert(false);
    }

    return false;
}


bool ray_flat_patch(const std::array<Vector3r, 4> &corners)
{
    int trials;
    assert(orient3d(corners[0], corners[1], corners[2], corners[3]) == 0);

    Vector3d dir(0, 0, 1);

    Vector3r inter;
    Vector3r n = cross(corners[0]-corners[1], corners[2] - corners[1]);
    bool inter_r;
    int dim;
    for(dim = 0; dim < 3; ++dim)
    {
        if (n[dim].get_sign() != 0){
            inter_r = segment_segment_inter(corners[0], corners[1], corners[2], corners[3], inter, dim);
            break;
        }
    }
    if(dim == 3)
    {
        std::cout << "n == 0" << std::endl;
        exit(0);
    }

    if (inter_r)
    {
        std::cout<<"butterfly 1"<<std::endl;
        exit(0);
    }

    for (dim = 0; dim < 3; ++dim)
    {
        if (n[dim].get_sign() != 0)
        {
            inter_r = segment_segment_inter(corners[1], corners[2], corners[3], corners[0], inter, dim);
            break;
        }
    }
    if (dim == 3)
    {
        std::cout << "n == 0" << std::endl;
        exit(0);
    }

    if (inter_r)
    {
        std::cout << "butterfly 2, cannot happend" << std::endl;
        exit(0);
    }

    for (trials = 0; trials < max_trials; ++trials)
    {
        int res0 = origin_ray_triangle_inter(dir, corners[0], corners[1], corners[2]);
        //hit
        if (res0 > 0)
            return true;

        int res1 = origin_ray_triangle_inter(dir, corners[0], corners[2], corners[3]);
        //hit
        if (res1 > 0)
            return true;

        //bad luck
        if (res1 < 0 || res1 < 0)
        {
            dir = Vector3d::Random();
            continue;
        }

        //both do not hit
        return false;
    }

    std::cout << "All rays are on edges, increase trials" << std::endl;
    return false;
}

template <typename FuncF>
bool ray_patch(const FuncF &func, int patch)
{
    static const std::array<std::array<int, 3>, 4> tet_faces = {{{{0, 1, 2}},
                                                                 {{0, 1, 3}},
                                                                 {{1, 2, 3}},
                                                                 {{2, 0, 3}}}};

    static const Vector3r zero(0, 0, 0);
    int trials;

    const std::array<Vector3r, 4> corners = func.corners(patch);

    if (orient3d(corners[0], corners[1], corners[2], corners[3]) == 0)
    {
        return ray_flat_patch(corners);
    }

    std::array<int, 2> Fp;
    std::array<int, 2> Fm;

    int ip = 0, im = 0;

    bool zero_inside = is_origin_in_tet(corners, tet_faces);

    for (int i = 0; i < 4; ++i)
    {
        const auto &f = tet_faces[i];

        const Vector3r bary = (corners[f[0]] + corners[f[1]] + corners[f[2]]) / 3;

        const auto phi_v = phi(bary, corners);
        assert(phi_v.get_sign() != 0);
        assert(ip <= 2);
        assert(im <= 2);

        if (phi_v.get_sign() > 0)
        {
            Fp[ip] = i;
            ip++;
        }
        else if (phi_v.get_sign() < 0)
        {
            Fm[im] = i;
            im++;
        }
        else
        {
            std::cout << "Barycenter has phi zero" << std::endl;
            assert(false);
        }
    }

    assert(ip == 2);
    assert(im == 2);


    if (zero_inside)
    {
        const auto phi_v = phi(zero, corners);

        Vector3d dir(0, 0, 1);
        const auto &F0 = phi_v.get_sign() > 0 ? tet_faces[Fm[0]] : tet_faces[Fp[0]];
        const auto &F1 = phi_v.get_sign() > 0 ? tet_faces[Fm[1]] : tet_faces[Fp[1]];

        for (trials = 0; trials < max_trials; ++trials)
        {
            int res0 = origin_ray_triangle_inter(dir, corners[F0[0]], corners[F0[1]], corners[F0[2]]);
            //hit
            if (res0 > 0)
                return true;

            int res1 = origin_ray_triangle_inter(dir, corners[F1[0]], corners[F1[1]], corners[F1[2]]);
            //hit
            if (res1 > 0)
                return true;

            //bad luck
            if (res1 < 0 || res1 < 0)
            {
                dir = Vector3d::Random();
                continue;
            }

            //both do not hit
            return false;
        }

        if (trials == max_trials)
        {
            std::cout << "All rays are on edges, increase trials" << std::endl;
            return false;
        }
    }
    else
    {
        Vector3d dir(0, 0, 1);
        const auto &Fp0 = tet_faces[Fp[0]];
        const auto &Fp1 = tet_faces[Fp[1]];

        for (trials = 0; trials < max_trials; ++trials)
        {
            int res0 = origin_ray_triangle_inter(dir, corners[Fp0[0]], corners[Fp0[1]], corners[Fp0[2]]);
            //bad luck
            if (res0 < 0)
            {
                dir = Vector3d::Random();
                continue;
            }

            int res1 = origin_ray_triangle_inter(dir, corners[Fp1[0]], corners[Fp1[1]], corners[Fp1[2]]);

            if (res1 < 0)
            {
                dir = Vector3d::Random();
                continue;
            }

            bool res0b = res0 == 1;
            bool res1b = res1 == 1;

            return res0b ^ res1b;
        }

        if (trials == max_trials)
        {
            std::cout << "All rays are on edges, increase trials" << std::endl;
            return false;
        }
    }

    //the other branches return
    assert(false);
    return false;
}

template <typename FuncF>
bool ccd(const FuncF &func)
{
    int S = 0;

    const int n_patches = func.n_patches();

    for (int patch = 0; patch < n_patches; ++patch)
    {
        bool is_ray_patch = ray_patch(func, patch);
        if (is_ray_patch)
            S++;
    }

    const auto caps = func.top_bottom_faces();
    int trials;

    for (const auto &tri : caps)
    {
        Vector3d dir(0, 0, 1);
        int res = -1;

        for (trials = 0; trials < max_trials; ++trials)
        {
            res = origin_ray_triangle_inter(dir, tri[0], tri[1], tri[2]);
            if (res >= 0)
                break;

            dir = Vector3d::Random();
        }

        if (trials == max_trials)
        {
            std::cout << "All rays are on edges, increase trials" << std::endl;
            return false;
        }

        if (res > 0)
            S++;
    }

    return S % 2 == 1;
}

bool vertexFaceCCD(const Vector3d &pts,
                   const Vector3d &v1s, const Vector3d &v2s, const Vector3d &v3s,
                   const Vector3d &pte,
                   const Vector3d &v1e, const Vector3d &v2e, const Vector3d &v3e)
{
    const auto tf = TriPtF(
        pts, v1s, v2s, v3s,
        pte, v1e, v2e, v3e);

    return ccd(tf);
}

bool edgeEdgeCCD(const Vector3d &a0s, const Vector3d &a1s,
                 const Vector3d &b0s, const Vector3d &b1s,
                 const Vector3d &a0e, const Vector3d &a1e,
                 const Vector3d &b0e, const Vector3d &b1e)
{
    const auto eef = EdgeEdgeF(
        a0s, a1s,
        b0s, b1s,
        a0e, a1e,
        b0e, b1e);

    return ccd(eef);
}

template bool ray_patch<TriPtF>(const TriPtF &, int);
template bool ccd<TriPtF>(const TriPtF &);

template bool ray_patch<EdgeEdgeF>(const EdgeEdgeF &, int);
template bool ccd<EdgeEdgeF>(const EdgeEdgeF &);

} // namespace eccd