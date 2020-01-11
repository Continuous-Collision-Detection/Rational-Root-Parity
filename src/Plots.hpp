#pragma once

#include "ECCD.hpp"

#include <functional>

namespace eccd
{
void save_prism(const std::string &name, const std::function<Vector3r(double, double, double)> &func, int n);
}
