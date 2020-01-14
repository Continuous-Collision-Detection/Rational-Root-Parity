#include "ECCD.hpp"
#include "Plots.hpp"

#include <stdlib.h>
#include <fstream>

using namespace eccd;

int main(int argc, char const *argv[])
{
    if (argc == 2)
    {
        const std::string path = argv[1];
        std::ifstream in(path, std::ios::binary);
        const bool is_tf = path.rfind("TF", 0) == 0;
        std::cout << "loading " << path << (is_tf ? " TF" : " EE") << std::endl;

        if (is_tf)
        {
            const Vector3d pts = read(in);

            const Vector3d v1s = read(in);
            const Vector3d v2s = read(in);
            const Vector3d v3s = read(in);

            const Vector3d pte = read(in);

            const Vector3d v1e = read(in);
            const Vector3d v2e = read(in);
            const Vector3d v3e = read(in);

            // std::cout<<"\n"<<pts.transpose()<<"\n"<<v1s.transpose()<<"\n"<<v2s.transpose()<<"\n"<<v3s.transpose()<<"\n"<<pte.transpose()<<"\n"<<v1e.transpose()<<"\n"<<v2e.transpose()<<"\n"<<v3e.transpose()<<std::endl;

            // std::cout<<(pts - pte)<<std::endl;
            // std::cout << (v1s - v1e) << std::endl;
            // std::cout << (v2s - v2e) << std::endl;
            // std::cout << (v3s - v3e) << std::endl;




            bool hit = eccd::vertexFaceCCD(
                pts,
                v1s, v2s, v3s,
                pte,
                v1e, v2e, v3e);

            const TriPtF tf(
                pts,
                v1s, v2s, v3s,
                pte,
                v1e, v2e, v3e);
            const int n = 1;
            save_prism("prism.obj", tf, n);

            std::cout << "hit? " << hit<<std::endl;
        }
        else
        {
            const Vector3d a0s = read(in);
            const Vector3d a0e = read(in);

            const Vector3d a1s = read(in);
            const Vector3d a1e = read(in);

            const Vector3d b0s = read(in);
            const Vector3d b0e = read(in);

            const Vector3d b1s = read(in);
            const Vector3d b1e = read(in);

            bool hit = eccd::edgeEdgeCCD(
                a0s, a1s,
                b0s, b1s,
                a0e, a1e,
                b0e, b1e);
        }

        in.close();
    }
    else
    {
        static const double EPSILON = std::numeric_limits<float>::epsilon();
        srand(42);
        // displacements
        double v0z = 4;
        double u0y = 0.5;
        double u1y = 2;
        double u0z = 6;

        Eigen::Vector3d v0(0, 1, v0z);
        // triangle = (v1, v2, v3)
        Eigen::Vector3d v1(-1, 0, 1);
        Eigen::Vector3d v2(1, 0, 1);
        Eigen::Vector3d v3(0, 0, -1);

        Eigen::Vector3d u0(0, u0y, u0z);

        Eigen::Vector3d u1(0.3, u1y, 1);

        bool hit = eccd::edgeEdgeCCD(
            v0, v1, v2, v3,
            v0 + u0, v1, v2, v3 + u1);

        const auto tf = EdgeEdgeF(
            v0, v1, v2, v3,
            v0 + u0, v1, v2, v3 + u1);
        const int n = 10;
        // save_prism("prism.obj", tf, n);
        save_hex("hex.obj", tf, n);

        bool check = ((-u0y + u1y >= 1) && (v0z + u0z >= v3.z()));
        std::cout << check << " code: " << hit << std::endl;
        assert(check == hit);

        // Eigen::Vector3d pt(-1, -1, 0);
        // Eigen::Vector3d v1(1, -1, 0);
        // // e2 = (v2, v3)
        // double e1x = 0;
        // // GENERATE(-1 - EPSILON, -1, -1 + EPSILON, -0.5, 0, 0.5, 1 - EPSILON, 1, 1 + EPSILON);
        // Eigen::Vector3d v2(e1x, 1, -1);
        // Eigen::Vector3d v3(e1x, 1, 1);

        // // displacements
        // double y_displacement = -1;
        // //GENERATE(-1.0, 0.0, 1 - EPSILON, 1.0, 1 + EPSILON, 2.0);
        // Eigen::Vector3d u0(0, y_displacement, 0);
        // Eigen::Vector3d u1(0, -y_displacement, 0);
        // Eigen::Vector3d up(1, -y_displacement, 0);

        // const auto tf = TriPtF(
        //     pt, v1, v2, v3,
        //     pt+up, v1 + u0, v2 + u1, v3 + u1);

        // const auto corners = tf.corners(0);

        // const int n = 10;
        // save_prism("prism.obj", tf, n);

        // ccd(tf);
    }
    return 0;
}
