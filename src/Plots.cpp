#include "Plots.hpp"

#include <vector>
#include <fstream>

namespace eccd
{
void save_prism(const std::string &name, const std::function<Vector3r(double, double, double)> &func, int n)
{
    std::ofstream out(name);

    if (!out.good())
    {
        out.close();
        return;
    }

    std::vector<std::vector<int>> indices(n + 1);
    int index = 1;

    /////////////////////////////////v=0
    for (int ui = 0; ui <= n; ++ui)
    {
        indices[ui].resize(n + 1);
        const double u = double(ui) / n;
        for (int ti = 0; ti <= n; ++ti)
        {
            const double t = double(ti) / n;
            const double v = 0;
            const auto res = func(u, v, t);

            out << "v "
                << res[0] << " " << res[1] << " " << res[2] << std::endl;
            indices[ui][ti] = index;
            ++index;
        }
    }

    for (int ui = 0; ui < n; ++ui)
    {
        for (int ti = 0; ti < n; ++ti)
        {
            out << "f " << indices[ui][ti] << " " << indices[ui + 1][ti] << " " << indices[ui + 1][ti + 1] << " " << indices[ui][ti + 1] << std::endl;
        }
    }


    /////////////////u=0
    for (int vi = 0; vi <= n; ++vi)
    {
        indices[vi].resize(n + 1);
        const double v = double(vi) / n;
        for (int ti = 0; ti <= n; ++ti)
        {
            const double t = double(ti) / n;
            const double u = 0;
            const auto res = func(u, v, t);

            out << "v "
                << res[0] << " " << res[1] << " " << res[2] << std::endl;
            indices[vi][ti] = index;
            ++index;
        }
    }

    for (int vi = 0; vi < n; ++vi)
    {
        for (int ti = 0; ti < n; ++ti)
        {
            out << "f " << indices[vi][ti] << " " << indices[vi + 1][ti] << " " << indices[vi + 1][ti + 1] << " " << indices[vi][ti + 1] << std::endl;
        }
    }

    /////////////////u+v=1
    for (int ki = 0; ki <= n; ++ki)
    {
        indices[ki].resize(n + 1);
        const double k = double(ki) / n;
        for (int ti = 0; ti <= n; ++ti)
        {
            const double t = double(ti) / n;
            const double u = k;
            const double v = 1-k;
            const auto res = func(u, v, t);

            out << "v "
                << res[0] << " " << res[1] << " " << res[2] << std::endl;
            indices[ki][ti] = index;
            ++index;
        }
    }

    for (int ki = 0; ki < n; ++ki)
    {
        for (int ti = 0; ti < n; ++ti)
        {
            out << "f " << indices[ki][ti] << " " << indices[ki + 1][ti] << " " << indices[ki + 1][ti + 1] << " " << indices[ki][ti + 1] << std::endl;
        }
    }
    out.close();
    {
    std::ofstream out("xx.obj");
    out << "v 0 0 0" << std::endl;
    out.close();
    }
}
} // namespace eccd
