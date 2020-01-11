#include "ECCD.hpp"
#include "Plots.hpp"

using namespace eccd;

int main(int argc, char const *argv[])
{

    Eigen::Vector3d pt(-1, -1, 0);
    Eigen::Vector3d v1(1, -1, 0);
    // e2 = (v2, v3)
    double e1x = 0;
    // GENERATE(-1 - EPSILON, -1, -1 + EPSILON, -0.5, 0, 0.5, 1 - EPSILON, 1, 1 + EPSILON);
    Eigen::Vector3d v2(e1x, 1, -1);
    Eigen::Vector3d v3(e1x, 1, 1);

    // displacements
    double y_displacement = -1;
    //GENERATE(-1.0, 0.0, 1 - EPSILON, 1.0, 1 + EPSILON, 2.0);
    Eigen::Vector3d u0(0, y_displacement, 0);
    Eigen::Vector3d u1(0, -y_displacement, 0);
    Eigen::Vector3d up(1, -y_displacement, 0);

    const auto tf = TriPtF(
        pt, v1, v2, v3,
        pt+up, v1 + u0, v2 + u1, v3 + u1);

    const auto corners = tf.corners(0);

    const int n = 10;
    save_prism("prism.obj", tf, n);

    ccd(tf);

    return 0;
}
