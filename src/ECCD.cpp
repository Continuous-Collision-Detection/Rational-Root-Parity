#include "ECCD.hpp"

#include "Plots.hpp"

#include <fstream>
#include <cassert>

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
    init_rationals();
}

void TriPtF ::init_rationals()
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

void TriPtF::save(std::ostream &out) const
{
    write(pts_, out);

    write(v1s_, out);
    write(v2s_, out);
    write(v3s_, out);

    write(pte_, out);

    write(v1e_, out);
    write(v2e_, out);
    write(v3e_, out);
}

Vector3r TriPtF::operator()(double u, double v, double t) const
{
    const Rational ur(u);
    const Rational vr(v);

    const Rational tr(t);
    return (1 - tr) * ptsr_ + tr * pter_ -
           ((1 - tr) * v1sr_ + t * v1er_) * (1 - ur - vr) -
           ((1 - tr) * v2sr_ + t * v2er_) * ur -
           ((1 - tr) * v3sr_ + t * v3er_) * vr;
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
        throw "Prism invalid corner";
        assert(false);
    }
}

std::vector<std::array<Vector3r, 3>> TriPtF::top_bottom_faces() const
{
    return {
        {{(*this)(0, 0, 0), (*this)(0, 1, 0), (*this)(1, 0, 0)}},
        {{(*this)(0, 0, 1), (*this)(0, 1, 1), (*this)(1, 0, 1)}}};
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
    init_rationals();
}

void EdgeEdgeF::init_rationals()
{
    for (int i = 0; i < 3; ++i)
    {
        a0rs_[i] = Rational(a0s_[i]);
        a1rs_[i] = Rational(a1s_[i]);

        b0rs_[i] = Rational(b0s_[i]);
        b1rs_[i] = Rational(b1s_[i]);

        a0re_[i] = Rational(a0e_[i]);
        a1re_[i] = Rational(a1e_[i]);

        b0re_[i] = Rational(b0e_[i]);
        b1re_[i] = Rational(b1e_[i]);
    }
}

void EdgeEdgeF::save(std::ostream &out) const
{
    write(a0s_, out);
    write(a0e_, out);

    write(a1s_, out);
    write(a1e_, out);

    write(b0s_, out);
    write(b0e_, out);

    write(b1s_, out);
    write(b1e_, out);
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
        throw "Hex invalid corner";
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
    static const std::array<int, 4> vv = {{0, 1, 2, 3}};
    const Rational g012 = func_g(x, corners, {{vv[0], vv[1], vv[2]}});
    const Rational g132 = func_g(x, corners, {{vv[1], vv[3], vv[2]}});
    const Rational g013 = func_g(x, corners, {{vv[0], vv[1], vv[3]}});
    const Rational g032 = func_g(x, corners, {{vv[0], vv[3], vv[2]}});

    const Rational h12 = g012 * g032;
    const Rational h03 = g132 * g013;

    const Rational phi = h12 - h03;

    return phi;
}

bool is_origin_in_tet(const std::array<Vector3r, 4> &corners, const std::array<std::array<int, 3>, 4> &tet_faces)
{
    const int max_trials = 8;
    Vector3d dir(0, 0, 1);

    for (int trials = 0; trials < max_trials; ++trials)
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

    // std::cout << "All rays are on edges, increase trials" << std::endl;
    throw "All rays are on edges, increase trials, for point in tet";
    assert(false);

    return false;
}

int ray_flat_patch(const std::array<Vector3r, 4> &corners, const Vector3d &dir)
{
    assert(orient3d(corners[0], corners[1], corners[2], corners[3]) == 0);

    Vector3r inter;
    const Vector3r n012 = cross(corners[0] - corners[1], corners[2] - corners[1]);
    const Vector3r n023 = cross(corners[0] - corners[2], corners[3] - corners[2]);

    bool is_012_zero = n012[0].get_sign() == 0 && n012[1].get_sign() == 0 && n012[2].get_sign() == 0;
    bool is_023_zero = n023[0].get_sign() == 0 && n023[1].get_sign() == 0 && n023[2].get_sign() == 0;

    bool inter_r = false;
    bool ok = false;

    const Vector3r &n = is_012_zero ? n023 : n012;

    if (is_012_zero && is_023_zero)
    {
        return 0;
    }

    for (int dim = 0; dim < 3; ++dim)
    {
        if (n[dim].get_sign() != 0)
        {
            inter_r = segment_segment_inter(corners[0], corners[1], corners[2], corners[3], inter, dim);
            ok = true;
            break;
        }
    }
    if (!ok)
    {
        // std::cout << "n == 0" << std::endl;
        // throw "n == 0";
        // assert(false);
        return 0;
    }

    if (inter_r)
    {
        if (corners[0] == corners[2] || corners[3] == corners[1])
            throw "butterfly is two edges";
        int res0 = corners[0] == corners[3] ? 0 : origin_ray_triangle_inter(dir, corners[0], corners[3], inter);
        if (res0 < 0)
            return -1;

        int res1 = corners[1] == corners[2] ? 0 : origin_ray_triangle_inter(dir, corners[1], corners[2], inter);
        if (res1 < 0)
            return -1;

        if (res0 == 2 || res1 == 2)
            return 2;

        if (res0 > 0 || res1 > 0)
            return 1;

        //both no hit
        return 0;
    }

    ok = false;
    for (int dim = 0; dim < 3; ++dim)
    {
        if (n[dim].get_sign() != 0)
        {
            inter_r = segment_segment_inter(corners[1], corners[2], corners[3], corners[0], inter, dim);
            ok = true;
            break;
        }
    }
    if (!ok)
    {
        // std::cout << "n == 0, cannot happend" << std::endl;
        throw "n == 0, cannot happend";
        assert(false);
        return 0;
    }

    if (inter_r)
    {
        if (corners[1] == corners[3] || corners[2] == corners[0])
            throw "butterfly 2 is two edges";
        int res0 = corners[0] == corners[1] ? 0 : origin_ray_triangle_inter(dir, corners[0], corners[1], inter);
        if (res0 < 0)
            return -1;

        int res1 = corners[2] == corners[3] ? 0 : origin_ray_triangle_inter(dir, corners[2], corners[3], inter);
        if (res1 < 0)
            return -1;

        if (res0 == 2 || res1 == 2)
            return 2;

        if (res0 > 0 || res1 > 0)
            return 1;

        //both no hit
        return 0;
    }

    int res0 = is_012_zero ? 0 : origin_ray_triangle_inter(dir, corners[0], corners[1], corners[2]);
    if (res0 < 0)
        return -1;

    int res1 = is_023_zero ? 0 : origin_ray_triangle_inter(dir, corners[0], corners[2], corners[3]);
    if (res1 < 0)
        return -1;

    if (res0 == 2 || res1 == 2)
        return 2;

    if (res0 > 0 || res1 > 0)
        return 1;

    //both no hit
    return 0;
}

template <typename FuncF>
int ray_patch(const FuncF &func, int patch, const Vector3d &dir)
{
    //0 1 2
    //1 3 2
    //0 2 3
    //1 3 0
    static const std::array<std::array<int, 3>, 4> tet_faces = {{{{0, 1, 2}},
                                                                 {{1, 3, 2}},
                                                                 {{0, 2, 3}},
                                                                 {{1, 3, 0}}}};

    static const Vector3r zero(0, 0, 0);

    const std::array<Vector3r, 4> corners = func.corners(patch);

    if (orient3d(corners[0], corners[1], corners[2], corners[3]) == 0)
    {
        return ray_flat_patch(corners, dir);
    }

    std::array<int, 2> Fp;
    std::array<int, 2> Fm;

    int ip = 0, im = 0;

    bool zero_inside = is_origin_in_tet(corners, tet_faces);

    // const Vector3r e0 = (corners[0] + corners[1]) / Rational(2);
    // const Vector3r e1 = (corners[1] + corners[2]) / Rational(2);
    // const Vector3r e2 = (corners[2] + corners[3]) / Rational(2);
    // const Vector3r e3 = (corners[0] + corners[3]) / Rational(2);

    // std::cout << "phi "<<phi(e0, corners) << std::endl;
    // std::cout << "phi "<<phi(e1, corners) << std::endl;
    // std::cout << "phi "<<phi(e2, corners) << std::endl;
    // std::cout << "phi "<<phi(e3, corners) << std::endl;
    // exit(0);

    for (int i = 0; i < 4; ++i)
    {
        const auto &f = tet_faces[i];
        assert(f.size() == 3);

        const Vector3r bary = (corners[f[0]] + corners[f[1]] + corners[f[2]]) / Rational(3);

        // const Vector3r e0 = (corners[f[0]] + corners[f[2]]) / Rational(2);
        // const Vector3r e1 = (corners[f[1]] + corners[f[2]]) / Rational(2);
        // const Vector3r e2 = (corners[f[0]] + corners[f[2]]) / Rational(2);

        // std::cout << phi(e0, corners) << std::endl;
        // std::cout << phi(e1, corners) << std::endl;
        // std::cout << phi(e2, corners) << std::endl;
        // print(bary);

        const auto phi_v = phi(bary, corners);
        // std::cout << phi_v << std::endl;
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
            // std::cout << "Barycenter has phi zero" << std::endl;
            throw "Barycenter has phi zero";
            assert(false);
        }
    }

    assert(ip == 2);
    assert(im == 2);

    // {
    //     std::ofstream os("blaam.obj");
    //     os << "v " << corners[tet_faces[Fm[0]][0]][0] << " " << corners[tet_faces[Fm[0]][0]][1] << " " << corners[tet_faces[Fm[0]][0]][2] << "\n";
    //     os << "v " << corners[tet_faces[Fm[0]][1]][0] << " " << corners[tet_faces[Fm[0]][1]][1] << " " << corners[tet_faces[Fm[0]][1]][2] << "\n";
    //     os << "v " << corners[tet_faces[Fm[0]][2]][0] << " " << corners[tet_faces[Fm[0]][2]][1] << " " << corners[tet_faces[Fm[0]][2]][2] << "\n";
    //     os << "f 1 2 3\n";
    //     os << "v " << corners[tet_faces[Fm[1]][0]][0] << " " << corners[tet_faces[Fm[1]][0]][1] << " " << corners[tet_faces[Fm[1]][0]][2] << "\n";
    //     os << "v " << corners[tet_faces[Fm[1]][1]][0] << " " << corners[tet_faces[Fm[1]][1]][1] << " " << corners[tet_faces[Fm[1]][1]][2] << "\n";
    //     os << "v " << corners[tet_faces[Fm[1]][2]][0] << " " << corners[tet_faces[Fm[1]][2]][1] << " " << corners[tet_faces[Fm[1]][2]][2] << "\n";
    //     os << "f 4 5 6\n";
    //     os.close();
    // }
    // {
    //     std::ofstream os("blaap.obj");
    //     os << "v " << corners[tet_faces[Fp[0]][0]][0] << " " << corners[tet_faces[Fp[0]][0]][1] << " " << corners[tet_faces[Fp[0]][0]][2] << "\n";
    //     os << "v " << corners[tet_faces[Fp[0]][1]][0] << " " << corners[tet_faces[Fp[0]][1]][1] << " " << corners[tet_faces[Fp[0]][1]][2] << "\n";
    //     os << "v " << corners[tet_faces[Fp[0]][2]][0] << " " << corners[tet_faces[Fp[0]][2]][1] << " " << corners[tet_faces[Fp[0]][2]][2] << "\n";
    //     os << "f 1 2 3\n";
    //     os << "v " << corners[tet_faces[Fp[1]][0]][0] << " " << corners[tet_faces[Fp[1]][0]][1] << " " << corners[tet_faces[Fp[1]][0]][2] << "\n";
    //     os << "v " << corners[tet_faces[Fp[1]][1]][0] << " " << corners[tet_faces[Fp[1]][1]][1] << " " << corners[tet_faces[Fp[1]][1]][2] << "\n";
    //     os << "v " << corners[tet_faces[Fp[1]][2]][0] << " " << corners[tet_faces[Fp[1]][2]][1] << " " << corners[tet_faces[Fp[1]][2]][2] << "\n";
    //     os << "f 4 5 6\n";
    //     os.close();
    // }

    if (zero_inside)
    {
        const auto phi_v = phi(zero, corners);

        if (phi_v.get_sign() == 0)
            return 2;

        const auto &F0 = phi_v.get_sign() > 0 ? tet_faces[Fm[0]] : tet_faces[Fp[0]];
        const auto &F1 = phi_v.get_sign() > 0 ? tet_faces[Fm[1]] : tet_faces[Fp[1]];

        int res0 = origin_ray_triangle_inter(dir, corners[F0[0]], corners[F0[1]], corners[F0[2]]);
        if (res0 < 0)
            return -1;

        int res1 = origin_ray_triangle_inter(dir, corners[F1[0]], corners[F1[1]], corners[F1[2]]);
        if (res1 < 0)
            return -1;

        //hit
        if (res0 > 0 || res1 > 0)
            return 1;

        return 0;
    }
    else
    {
        const auto &Fp0 = tet_faces[Fp[0]];
        const auto &Fp1 = tet_faces[Fp[1]];

        int res0 = origin_ray_triangle_inter(dir, corners[Fp0[0]], corners[Fp0[1]], corners[Fp0[2]]);
        //bad luck
        if (res0 < 0)
            return -1;

        int res1 = origin_ray_triangle_inter(dir, corners[Fp1[0]], corners[Fp1[1]], corners[Fp1[2]]);
        if (res1 < 0)
            return -1;

        bool res0b = (res0 >= 1);
        bool res1b = (res1 >= 1);

        //res0b xor res1b
        return (!res0b != !res1b) ? 1 : 0;
    }

    //the other branches return
    assert(false);
    return 0;
}

template <typename FuncF>
int ccd(const FuncF &func, const Vector3d &dir)
{
    int S = 0;

    const int n_patches = func.n_patches();

    for (int patch = 0; patch < n_patches; ++patch)
    {
        int is_ray_patch = ray_patch(func, patch, dir);
        if (is_ray_patch == 2)
            return 1;

        if (is_ray_patch == -1)
            return -1;

        if (is_ray_patch == 1)
            S++;
    }

    const auto caps = func.top_bottom_faces();

    for (const auto &tri : caps)
    {
        int res = origin_ray_triangle_inter(dir, tri[0], tri[1], tri[2]);
        if (res == 2)
            return 1;
        if (res == -1)
            return -1;

        if (res > 0)
            S++;
    }

    return ((S % 2) == 1) ? 1 : 0;
}

template <typename FuncF>
bool retrial_ccd(const FuncF &func)
{
    static const int max_trials = 8;
    Vector3d dir(1, 0, 0);

    int res = -1;
    int trials;
    for (trials = 0; trials < max_trials; ++trials)
    {
        res = ccd(func, dir);

        if (res >= 0)
            break;

        dir = Vector3d::Random();
    }

    if (trials == max_trials)
    {

        // std::cout << "All rays are on edges, increase trials" << std::endl;
        throw "All rays are on edges, increase trials";
        return false;
    }

    return res >= 1;
}

bool vertexFaceCCD(const Vector3d &pts,
                   const Vector3d &v1s, const Vector3d &v2s, const Vector3d &v3s,
                   const Vector3d &pte,
                   const Vector3d &v1e, const Vector3d &v2e, const Vector3d &v3e)
{
    const auto tf = TriPtF(
        pts, v1s, v2s, v3s,
        pte, v1e, v2e, v3e);

    if (pts == pte && v1s == v1e && v2s == v2e && v3s == v3e)
    {
        return false;
    }

    if (
        (pts - pte) == (v1s - v1e) &&
        (pts - pte) == (v2s - v2e) &&
        (pts - pte) == (v3s - v3e) &&
        (v1s - v1e) == (v2s - v2e) &&
        (v3s - v3e) == (v2s - v2e))
    {
        return false;
    }

    try
    {
        bool ok = retrial_ccd(tf);
        return ok;
    }
    catch (const char *msg)
    {
        std::string name = "TF" + std::to_string(rand()) + ".txt";
        std::cout << "[ERROR] " << msg << "out written to: " << name << std::endl;
        std::ofstream ofs(name, std::ios::binary);
        tf.save(ofs);
        ofs.close();

        return false;
    }
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

    if (a0s == a0e && a1s == a1e && b0s == b0e && b1s == b1e)
    {
        return false;
    }

    if (
        (a0s - a0e) == (a1s - a1e) &&
        (a0s - a0e) == (b0s - b0e) &&
        (a0s - a0e) == (b1s - b1e) &&
        (a1s - a1e) == (b1s - b1e) &&
        (b1s - b1e) == (b0s - b0e))
    {
        return false;
    }

    try
    {
        bool ok = retrial_ccd(eef);
        return ok;
    }
    catch (const char *msg)
    {
        std::string name = "EE" + std::to_string(rand()) + ".txt";
        std::cout << "[ERROR] " << msg << "out written to: " << name << std::endl;
        std::ofstream ofs(name, std::ios::binary);
        eef.save(ofs);
        ofs.close();

        return false;
    }
}

template int ray_patch<TriPtF>(const TriPtF &, int, const Vector3d &);
template int ccd<TriPtF>(const TriPtF &, const Vector3d &);
template bool retrial_ccd<TriPtF>(const TriPtF &func);

template int ray_patch<EdgeEdgeF>(const EdgeEdgeF &, int, const Vector3d &);
template int ccd<EdgeEdgeF>(const EdgeEdgeF &, const Vector3d &);
template bool retrial_ccd<EdgeEdgeF>(const EdgeEdgeF &func);

} // namespace eccd