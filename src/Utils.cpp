#include "Utils.hpp"

#include <array>

namespace eccd
{

int orient3d(const Vector3r &a, const Vector3r &b, const Vector3r &c, const Vector3r &d)
{
    const Rational det = (a - d).dot(cross(b - d, c - d));
    return det.get_sign();
}

int origin_ray_triangle_inter(const Vector3d &dirf, const Vector3r &t1, const Vector3r &t2, const Vector3r &t3)
{
    const Vector3r dir(dirf[0], dirf[1], dirf[2]);
    Rational denom = dir[0] * t1[1] * t2[2] - dir[0] * t1[1] * t3[2] - dir[0] * t1[2] * t2[1] + dir[0] * t1[2] * t3[1] + dir[0] * t2[1] * t3[2] - dir[0] * t2[2] * t3[1] - dir[1] * t1[0] * t2[2] + dir[1] * t1[0] * t3[2] + dir[1] * t1[2] * t2[0] - dir[1] * t1[2] * t3[0] - dir[1] * t2[0] * t3[2] + dir[1] * t2[2] * t3[0] + dir[2] * t1[0] * t2[1] - dir[2] * t1[0] * t3[1] - dir[2] * t1[1] * t2[0] + dir[2] * t1[1] * t3[0] + dir[2] * t2[0] * t3[1] - dir[2] * t2[1] * t3[0];

    //infinite intersections
    if (denom.get_sign() == 0)
        return -1;

    if (denom.get_sign() < 0)
        denom = -denom;

    assert(denom.get_sign() > 0);

    const Rational u = dir[0] * t2[1] * t3[2] - dir[0] * t2[2] * t3[1] - dir[1] * t2[0] * t3[2] + dir[1] * t2[2] * t3[0] + dir[2] * t2[0] * t3[1] - dir[2] * t2[1] * t3[0];
    const Rational v = dir[0] * t1[2] * t3[1] - dir[0] * t1[1] * t3[2] + dir[1] * t1[0] * t3[2] - dir[1] * t1[2] * t3[0] - dir[2] * t1[0] * t3[1] + dir[2] * t1[1] * t3[0];
    if (u >= 0 && u <= denom && v >= 0 && v <= denom)
        return 1;

    return 0;
}

bool segment_segment_inter(const Vector3r &s0, const Vector3r &e0, const Vector3r &s1, const Vector3r &e1, Vector3r &res, int axis)
{
    const int i1 = (axis + 1) % 3;
    const int i2 = (axis + 2) % 3;

    Rational dd = e0[i1] * e1[i2] - e0[i1] * s1[i2] - e0[i2] * e1[i1] + e0[i2] * s1[i1] + e1[i1] * s0[i2] - e1[i2] * s0[i1] + s0[i1] * s1[i2] - s0[i2] * s1[i1];
    if (dd.get_sign() == 0)
    {
        std::array<int, 2> is = {{i1, i2}};
        int ik;
        for (ik = 0; ik < 2; ++ik)
        {
            const int ii = is[ik];

            const Rational ddd0 = e0[ii] - s0[ii];
            if (ddd0.get_sign() == 0)
                continue;
            const Rational tt0 = (s1[ii] - s0[ii]) / ddd0;
            if (tt0 < 0 || tt0 > 0)
                return false;

            const Rational ddd1 = e1[ii] - s1[ii];
            if (ddd1.get_sign() == 0)
                continue;
            const Rational tt1 = (-s1[ii] + s0[ii]) / ddd1;
            if (tt1 < 0 || tt1 > 0)
                return false;

            const Vector3r p0 = (1 - tt0) * s0 + tt0 * e0;
            const Vector3r p1 = (1 - tt1) * s1 + tt1 * e1;

            if (p0[0] == p1[0] && p0[1] == p1[1] && p0[2] == p1[2])
            {
                //TODO
                std::cout << "degenerate line line intersection" << std::endl;
            }
            else
                return false;
        }
        if(ik == 2)
            std::cout<<"shouldnt happend"<<std::endl;
        return false;
    }

    int sign = 1;
    if (dd.get_sign() < 0)
    {
        sign = -1;
        dd = -1 * dd;
    }

    assert(dd.get_sign() > 0);

    Rational t0 = e1[i1] * s0[i2] - e1[i1] * s1[i2] - e1[i2] * s0[i1] + e1[i2] * s1[i1] + s0[i1] * s1[i2] - s0[i2] * s1[i1];
    Rational t1 = e0[i1] * s0[i2] - e0[i1] * s1[i2] - e0[i2] * s0[i1] + e0[i2] * s1[i1] + s0[i1] * s1[i2] - s0[i2] * s1[i1];

    if (t0 < 0 || t0 > dd || t1 < 0 || t1 > dd)
    {
        return false;
    }

    t0 = sign * t0 / dd;

    res = (1 - t0) * s0 + t0 * e0;
#ifndef NDEBUG
    t1 = sign * t1 / dd;
    const Vector3r p1 = (1 - t1) * s1 + t1 * e1;

    assert(res[0] == p1[0] && res[1] == p1[1] && res[2] == p1[2]);
#endif
    return true;
}
} // namespace eccd
