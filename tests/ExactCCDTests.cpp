#include <catch.hpp>

#include "ECCD.hpp"
#include "Plots.hpp"

#include <stdlib.h>

static const double EPSILON = std::numeric_limits<float>::epsilon();

TEST_CASE("Corner case0",
          "Corner case")
{
    const Eigen::Vector3d pts(0.0742893551979073, 0.0263742788498324, 0.0010474039494139426);
    const Eigen::Vector3d v1s(0.07290626476594385, 0.027311689801583078, 0.0013884370526333508);
    const Eigen::Vector3d v2s(0.07335896966730952, 0.02637427884983239, 0.001035981768919693);
    const Eigen::Vector3d v3s(0.07242914448368196, 0.026374278849832386, 0.001070228306583729);

    const Eigen::Vector3d pte(0.07562935519790728, 0.0263742788498324, 0.0010474039494139335);
    const Eigen::Vector3d v1e(0.07424626476594383, 0.027311689801583078, 0.0013884370526333339);
    const Eigen::Vector3d v2e(0.07469896966730949, 0.0263742788498324, 0.0010359817689196637);
    const Eigen::Vector3d v3e(0.07376914448368192, 0.0263742788498324, 0.001070228306583705);

    const eccd::TriPtF tf(
        pts,
        v1s, v2s, v3s,
        pte,
        v1e, v2e, v3e);
    const int n = 10;
    eccd::save_prism("prism.obj", tf, n);

    bool hit = eccd::vertexFaceCCD(pts, v1s, v2s, v3s, pte, v1e, v2e, v3e);
    // std::cout<<hit<<std::endl;
    CHECK(!hit);
    // exit(0);
}

TEST_CASE("Corner case1",
          "Corner case")
{
    const Eigen::Vector3d pts(1.0, 0.8660254037844386, 0.0);
    const Eigen::Vector3d v1s(0.0, 0.8660254037844386, 0.0);
    const Eigen::Vector3d v2s(0.0, 1.8660254037844388, 0.0);
    const Eigen::Vector3d v3s(1.0, 1.8660254037844388, 0.0);
    const Eigen::Vector3d pte(1.0, 0.7679589037844385, 0.0);
    const Eigen::Vector3d v1e(2.220446049250313E-16, 0.7679589037844385, 0.0);
    const Eigen::Vector3d v2e(2.220446049250313E-16, 1.7679589037844385, 2.220446049250313E-16);
    const Eigen::Vector3d v3e(1.0, 1.7679589037844385, 2.220446049250313E-16);

    // const eccd::TriPtF tf(
    //     pts,
    //     v1s, v2s, v3s,
    //     pte,
    //     v1e, v2e, v3e);
    // const int n = 10;
    // eccd::save_prism("prism.obj", tf, n);

    bool hit = eccd::vertexFaceCCD(pts, v1s, v2s, v3s, pte, v1e, v2e, v3e);
    CHECK(!hit);
    // exit(0);
}

TEST_CASE("Corner case",
          "Corner case")
{
    const double v1 = 0.5028730361513796;
    const double v11 = 0.5028730361513793;
    const double v2 = 2.220446049250313e-16;
    const double v3 = 0.27762640378443826;
    const double v31 = 0.2776264037844385;
    const double v32 = 1.0000000000000002;

    const Eigen::Vector3d pts(1, 0.5, 1);
    const Eigen::Vector3d v1s(1, v1, 1);
    const Eigen::Vector3d v2s(0, v1, v2);
    const Eigen::Vector3d v3s(1, v11, 0);
    const Eigen::Vector3d pte(1, 0.5, 1);
    const Eigen::Vector3d v1e(1, v31, 1);
    const Eigen::Vector3d v2e(0, v3, 0);
    const Eigen::Vector3d v3e(v32, v3, 0);

    // const eccd::TriPtF tf(
    //     pts,
    //     v1s, v2s, v3s,
    //     pte,
    //     v1e, v2e, v3e);
    // const int n = 10;
    // eccd::save_prism("prism.obj", tf, n);

    bool hit = eccd::vertexFaceCCD(pts, v1s, v2s, v3s, pte, v1e, v2e, v3e);
    CHECK(hit);
}

TEST_CASE("Test Edge-Edge Exact Continous Collision Detection",
          "[ccd][exact-ccd][edge-edge]")
{
    // EE Exact CCD unit test
    // e0 = (v0, v1)
    Eigen::Vector3d v0(-1, -1, 0);
    Eigen::Vector3d v1(1, -1, 0);
    // e2 = (v2, v3)
    double e1x = GENERATE(
        -1 - EPSILON, -1, -1 + EPSILON, -0.5, 0, 0.5, 1 - EPSILON, 1,
        1 + EPSILON);
    Eigen::Vector3d v2(e1x, 1, -1);
    Eigen::Vector3d v3(e1x, 1, 1);

    // displacements
    double y_displacement = GENERATE(
        -1.0, 0.0, 1 - EPSILON, 1.0, 1 + EPSILON, 2.0);
    Eigen::Vector3d u0(0, y_displacement, 0);
    Eigen::Vector3d u1(0, -y_displacement, 0);

    bool hit = eccd::edgeEdgeCCD(
        v0, v1, v2, v3,
        v0 + u0, v1 + u0, v2 + u1, v3 + u1);

    CAPTURE(y_displacement, e1x);
    CHECK(hit == (y_displacement >= 1.0 && e1x >= -1 && e1x <= 1));
}

TEST_CASE("Test Point-Triangle Exact Continous Collision Detection",
          "[ccd][exact-ccd][point-triangle]")
{
    srand(42);
    // PT Exact CCD unit test
    // point
    double v0z = GENERATE(0.0, -1.0);
    Eigen::Vector3d v0(0, 1, v0z);
    // triangle = (v1, v2, v3)
    Eigen::Vector3d v1(-1, 0, 1);
    Eigen::Vector3d v2(1, 0, 1);
    Eigen::Vector3d v3(0, 0, -1);

    // displacements
    double u0y = -GENERATE(-1.0, 0.0, 0.5 - EPSILON, 0.5, 0.5 + EPSILON, 1.0, 2.0);
    double u0z = GENERATE(-EPSILON, 0.0, EPSILON);
    Eigen::Vector3d u0(0, u0y, u0z);
    double u1y = GENERATE(-1.0, 0.0, 0.5 - EPSILON, 0.5, 0.5 + EPSILON, 1.0, 2.0);
    Eigen::Vector3d u1(0, u1y, 0);

    SECTION("Clockwise triangle")
    {
        // already in clockwise order
    }
    SECTION("Counter-clockwise triangle")
    {
        std::swap(v1, v2);
    }

    bool hit = eccd::vertexFaceCCD(
        v0, v1, v2, v3,
        v0 + u0, v1 + u1, v2 + u1, v3 + u1);

    CAPTURE(v0z, u0y, u1y, u0z, EPSILON);
    CHECK(hit == ((-u0y + u1y >= 1) && (v0z + u0z >= v3.z())));
}
