#include "Utils.hpp"

#include <fstream>
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
    const Rational denom = dir[0] * t1[1] * t2[2] - dir[0] * t1[1] * t3[2] - dir[0] * t1[2] * t2[1] + dir[0] * t1[2] * t3[1] + dir[0] * t2[1] * t3[2] - dir[0] * t2[2] * t3[1] - dir[1] * t1[0] * t2[2] + dir[1] * t1[0] * t3[2] + dir[1] * t1[2] * t2[0] - dir[1] * t1[2] * t3[0] - dir[1] * t2[0] * t3[2] + dir[1] * t2[2] * t3[0] + dir[2] * t1[0] * t2[1] - dir[2] * t1[0] * t3[1] - dir[2] * t1[1] * t2[0] + dir[2] * t1[1] * t3[0] + dir[2] * t2[0] * t3[1] - dir[2] * t2[1] * t3[0];

    // const auto n = cross(t1 - t2, t3 - t2);
    // print(n);
    // std::cout<<denom<<std::endl;


    //infinite intersections
    if (denom.get_sign() == 0)
        return -1;


    // assert(denom.get_sign() > 0);



    const Rational u = (-1*dir[0] * t1[1] * t3[2] + dir[0] * t1[2] * t3[1] + dir[1] * t1[0] * t3[2] - dir[1] * t1[2] * t3[0] - dir[2] * t1[0] * t3[1] + dir[2] * t1[1] * t3[0])/denom;
    const Rational v = (dir[0] * t1[1] * t2[2] - dir[0] * t1[2] * t2[1] - dir[1] * t1[0] * t2[2] + dir[1] * t1[2] * t2[0] + dir[2] * t1[0] * t2[1] - dir[2] * t1[1] * t2[0])/denom;
    const Rational t = (t1[0] * t2[1] * t3[2] - t1[0] * t2[2] * t3[1] - t1[1] * t2[0] * t3[2] + t1[1] * t2[2] * t3[0] + t1[2] * t2[0] * t3[1] - t1[2] * t2[1] * t3[0])/denom;

    // std::ofstream os("blaa.obj");
    // os << "v " << t1[0] << " " << t1[1] << " " << t1[2] << "\n";
    // os << "v " << t2[0] << " " << t2[1] << " " << t2[2] << "\n";
    // os << "v " << t3[0] << " " << t3[1] << " " << t3[2] << "\n";
    // os << "f 1 2 3\n";
    // os.close();
    // std::ofstream os1("blaa1.obj");
    // os1 << "v " << dir[0]*t << " " << dir[1]*t << " " << dir[2]*t << "\n";
    // os1.close();

    // std::cout<<t<<std::endl;
    // std::cout << u << std::endl;
    // std::cout << v << std::endl;
    // std::cout << 1-u-v << std::endl;

    if (u >= 0 && u <= 1 && v >= 0 && v <= 1 && u+v<=1 && t>= 0){
        if (t.get_sign() == 0)
            return 2;
        //on a corner
        const auto w = 1-u-v;
        if (u.get_sign() == 0 || v.get_sign() == 0 || w.get_sign() == 0)
            return -1;

        return 1;
    }

    return 0;
}

bool segment_segment_inter(const Vector3r &s0, const Vector3r &e0, const Vector3r &s1, const Vector3r &e1, Vector3r &res, int axis)
{
    const int i1 = (axis + 1) % 3;
    const int i2 = (axis + 2) % 3;

    const Rational dd = e0[i1] * e1[i2] - e0[i1] * s1[i2] - e0[i2] * e1[i1] + e0[i2] * s1[i1] + e1[i1] * s0[i2] - e1[i2] * s0[i1] + s0[i1] * s1[i2] - s0[i2] * s1[i1];

    if (dd.get_sign() == 0)
    {
        return false;
    }

    const Rational t0 = (e1[i1] * s0[i2] - e1[i1] * s1[i2] - e1[i2] * s0[i1] + e1[i2] * s1[i1] + s0[i1] * s1[i2] - s0[i2] * s1[i1]) / dd;
    const Rational t1 = (e0[i1] * s0[i2] - e0[i1] * s1[i2] - e0[i2] * s0[i1] + e0[i2] * s1[i1] + s0[i1] * s1[i2] - s0[i2] * s1[i1]) / dd;

    //we exclude intersection on corners
    if (t0 <= 0 || t0 >= 1 || t1 <= 0 || t1 >= 1)
    {
        return false;
    }

    res = (1 - t0) * s0 + t0 * e0;
#ifndef NDEBUG
    const Vector3r p1 = (1 - t1) * s1 + t1 * e1;

    assert(res[0] == p1[0] && res[1] == p1[1] && res[2] == p1[2]);
#endif
    return true;
}

void write(const Vector3d &v, std::ostream &out)
{
    out.write(reinterpret_cast<const char *>(&v[0]), sizeof(v[0]));
    out.write(reinterpret_cast<const char *>(&v[1]), sizeof(v[1]));
    out.write(reinterpret_cast<const char *>(&v[2]), sizeof(v[2]));
}

Vector3d read(std::istream &in)
{
    Vector3d res;
    double tmp;
    in.read(reinterpret_cast<char *>(&tmp), sizeof(tmp));
    res[0] = tmp;

    in.read(reinterpret_cast<char *>(&tmp), sizeof(tmp));
    res[1] = tmp;

    in.read(reinterpret_cast<char *>(&tmp), sizeof(tmp));
    res[2] = tmp;

    return res;
}

int segment_triangle_inter(const Vector3d &ef0, const Vector3d &ef1, const Vector3d &tf1, const Vector3d &tf2, const Vector3d &tf3)
{
    Vector3r e0, e1, t1, t2, t3;

    for(int d = 0; d < 3; ++d)
    {
        e0[d] = ef0[d];
        e1[d] = ef1[d];

        t1[d] = tf1[d];
        t2[d] = tf2[d];
        t3[d] = tf3[d];
    }

    const Rational d = e0[0]*t1[1]*t2[2]-e0[0]*t1[1]*t3[2]-e0[0]*t1[2]*t2[1]+e0[0]*t1[2]*t3[1]+e0[0]*t2[1]*t3[2]-e0[0]*t2[2]*t3[1]-e0[1]*t1[0]*t2[2]+e0[1]*t1[0]*t3[2]+e0[1]*t1[2]*t2[0]-e0[1]*t1[2]*t3[0]-e0[1]*t2[0]*t3[2]+e0[1]*t2[2]*t3[0]+e0[2]*t1[0]*t2[1]-e0[2]*t1[0]*t3[1]-e0[2]*t1[1]*t2[0]+e0[2]*t1[1]*t3[0]+e0[2]*t2[0]*t3[1]-e0[2]*t2[1]*t3[0]-e1[0]*t1[1]*t2[2]+e1[0]*t1[1]*t3[2]+e1[0]*t1[2]*t2[1]-e1[0]*t1[2]*t3[1]-e1[0]*t2[1]*t3[2]+e1[0]*t2[2]*t3[1]+e1[1]*t1[0]*t2[2]-e1[1]*t1[0]*t3[2]-e1[1]*t1[2]*t2[0]+e1[1]*t1[2]*t3[0]+e1[1]*t2[0]*t3[2]-e1[1]*t2[2]*t3[0]-e1[2]*t1[0]*t2[1]+e1[2]*t1[0]*t3[1]+e1[2]*t1[1]*t2[0]-e1[2]*t1[1]*t3[0]-e1[2]*t2[0]*t3[1]+e1[2]*t2[1]*t3[0];
    if(d.get_sign() == 0)
        return -1;

    const Rational t = (e0[0] * t1[1] * t2[2] - e0[0] * t1[1] * t3[2] - e0[0] * t1[2] * t2[1] + e0[0] * t1[2] * t3[1] + e0[0] * t2[1] * t3[2] - e0[0] * t2[2] * t3[1] - e0[1] * t1[0] * t2[2] + e0[1] * t1[0] * t3[2] + e0[1] * t1[2] * t2[0] - e0[1] * t1[2] * t3[0] - e0[1] * t2[0] * t3[2] + e0[1] * t2[2] * t3[0] + e0[2] * t1[0] * t2[1] - e0[2] * t1[0] * t3[1] - e0[2] * t1[1] * t2[0] + e0[2] * t1[1] * t3[0] + e0[2] * t2[0] * t3[1] - e0[2] * t2[1] * t3[0] - t1[0] * t2[1] * t3[2] + t1[0] * t2[2] * t3[1] + t1[1] * t2[0] * t3[2] - t1[1] * t2[2] * t3[0] - t1[2] * t2[0] * t3[1] + t1[2] * t2[1] * t3[0]) / d;
    const Rational u = (-e0[0] * e1[1] * t1[2] + e0[0] * e1[1] * t3[2] + e0[0] * e1[2] * t1[1] - e0[0] * e1[2] * t3[1] - e0[0] * t1[1] * t3[2] + e0[0] * t1[2] * t3[1] + e0[1] * e1[0] * t1[2] - e0[1] * e1[0] * t3[2] - e0[1] * e1[2] * t1[0] + e0[1] * e1[2] * t3[0] + e0[1] * t1[0] * t3[2] - e0[1] * t1[2] * t3[0] - e0[2] * e1[0] * t1[1] + e0[2] * e1[0] * t3[1] + e0[2] * e1[1] * t1[0] - e0[2] * e1[1] * t3[0] - e0[2] * t1[0] * t3[1] + e0[2] * t1[1] * t3[0] + e1[0] * t1[1] * t3[2] - e1[0] * t1[2] * t3[1] - e1[1] * t1[0] * t3[2] + e1[1] * t1[2] * t3[0] + e1[2] * t1[0] * t3[1] - e1[2] * t1[1] * t3[0]) / d;
    const Rational v = (e0[0] * e1[1] * t1[2] - e0[0] * e1[1] * t2[2] - e0[0] * e1[2] * t1[1] + e0[0] * e1[2] * t2[1] + e0[0] * t1[1] * t2[2] - e0[0] * t1[2] * t2[1] - e0[1] * e1[0] * t1[2] + e0[1] * e1[0] * t2[2] + e0[1] * e1[2] * t1[0] - e0[1] * e1[2] * t2[0] - e0[1] * t1[0] * t2[2] + e0[1] * t1[2] * t2[0] + e0[2] * e1[0] * t1[1] - e0[2] * e1[0] * t2[1] - e0[2] * e1[1] * t1[0] + e0[2] * e1[1] * t2[0] + e0[2] * t1[0] * t2[1] - e0[2] * t1[1] * t2[0] - e1[0] * t1[1] * t2[2] + e1[0] * t1[2] * t2[1] + e1[1] * t1[0] * t2[2] - e1[1] * t1[2] * t2[0] - e1[2] * t1[0] * t2[1] + e1[2] * t1[1] * t2[0]) / d;

    // std::cout << t << std::endl;
    // std::cout << u << std::endl;
    // std::cout << v << std::endl;

    if(t < 0 || t > 1)
        return 0;

    if (u >= 0 && u <= 1 && v >= 0 && v <= 1 && u + v <= 1)
        return 1;

    return 0;

}
} // namespace eccd
